//! CUDA executor: real GPU Smith-Waterman via [`cudarc`].
//!
//! This module is only compiled when the `gpu-cuda` feature is enabled
//! and `cudarc` is in the dependency tree.
//!
//! ## Lifecycle
//!
//! 1. [`CudaExecutor::new`] selects device 0, compiles the kernel source
//!    via NVRTC, loads the resulting PTX, and stashes [`CudaFunction`]
//!    handles for the kernels.
//! 2. [`CudaExecutor::score_batch`] sorts subjects by length (warp coherence),
//!    packs them with their offsets, allocates device buffers, copies inputs,
//!    launches the kernel, copies scores back, and un-sorts to caller order.
//! 3. The device handle (`Arc<CudaDevice>`) lives for the executor's
//!    lifetime; the PTX is loaded once.
//!
//! ## Memory model
//!
//! Per launch, the executor allocates:
//!   - `subjects_packed`: total subject residues (bytes)
//!   - `subj_offsets`: (num_subjects + 1) × u32
//!   - `query`: query_len bytes
//!   - `sub_matrix`: 24×24 i8 (protein only)
//!   - `work_H`, `work_E`: num_subjects × query_len × i16 each
//!   - `scores`: num_subjects × i32
//!
//! For a 500-aa query against 100k subjects, work_H + work_E is ~200 MB.
//! Callers should batch in chunks of ~50k subjects to stay under typical
//! consumer-GPU memory budgets.
//!
//! ## Known limits
//!
//! - Max query length: 65535 (int16 internal score range, fits typical proteins).
//! - Scores saturate at i16 max (32767); long highly-similar matches may clamp.
//!   For BLOSUM62 with default gap penalties this corresponds to ~150 aa of
//!   perfect identity, which is much longer than any single HSP.

use std::sync::Arc;

use anyhow::{anyhow, Context, Result};
use cudarc::driver::{CudaDevice, CudaFunction, CudaSlice, LaunchAsync, LaunchConfig};
use cudarc::nvrtc::compile_ptx;

use crate::algorithm::Scoring;

use super::encoding::{
    encode_nucleotide, flatten_matrix, pack_sequences, PROTEIN_ALPHABET_SIZE,
};
use super::{GpuExecutor, GpuScore};

const KERNEL_SOURCE: &str = include_str!("kernels/sw.cu");

const MODULE_NAME: &str = "saber_sw";
const KERNEL_PROTEIN: &str = "saber_sw_protein_score";
const KERNEL_NUCLEOTIDE: &str = "saber_sw_nucleotide_score";

/// Default block size. 128 is a reasonable starting point — Richard, tune
/// this on your hardware. Common values: 64 for older GPUs, 128 or 256 for
/// Ampere+. The kernel is register-bound, so larger blocks may spill.
const BLOCK_SIZE: u32 = 128;

/// Hard cap on query length. The kernel uses i16 for DP cell values; longer
/// queries risk score saturation for highly-similar regions. Realistically,
/// proteins are < 5000 aa so this is plenty.
const MAX_QUERY_LEN: usize = 16384;

pub struct CudaExecutor {
    device: Arc<CudaDevice>,
    protein_kernel: CudaFunction,
    nucleotide_kernel: CudaFunction,
    block_size: u32,
}

impl CudaExecutor {
    /// Initialize the executor on device 0.
    ///
    /// Compiles the CUDA kernel via NVRTC the first time. Subsequent calls
    /// (e.g. if the executor is dropped and recreated) will re-compile;
    /// for production use, hold the executor for the program's lifetime.
    pub fn new() -> Result<Self> {
        Self::new_on_device(0)
    }

    pub fn new_on_device(ordinal: usize) -> Result<Self> {
        log::info!("Initializing CUDA device {ordinal}");
        let device = CudaDevice::new(ordinal)
            .with_context(|| format!("could not open CUDA device {ordinal}"))?;

        log::info!("Compiling SABER GPU kernels via NVRTC...");
        let ptx = compile_ptx(KERNEL_SOURCE)
            .map_err(|e| {
                eprintln!("NVRTC ERROR:\n{:#?}", e);
                e
            })
            .context("NVRTC failed to compile sw.cu")?;

        device
            .load_ptx(ptx, MODULE_NAME, &[KERNEL_PROTEIN, KERNEL_NUCLEOTIDE])
            .context("failed to load SABER kernel PTX into CUDA context")?;

        let protein_kernel = device
            .get_func(MODULE_NAME, KERNEL_PROTEIN)
            .ok_or_else(|| anyhow!("kernel `{}` not found in loaded module", KERNEL_PROTEIN))?;
        let nucleotide_kernel = device
            .get_func(MODULE_NAME, KERNEL_NUCLEOTIDE)
            .ok_or_else(|| anyhow!("kernel `{}` not found in loaded module", KERNEL_NUCLEOTIDE))?;

        log::info!("CUDA executor ready on device {ordinal}");

        Ok(Self {
            device,
            protein_kernel,
            nucleotide_kernel,
            block_size: BLOCK_SIZE,
        })
    }

//    /// Free device memory budget in bytes. Useful for picking a batch size
//    /// before launching score_batch.
//    pub fn free_device_memory(&self) -> Result<usize> {
//        self.device
//            .free_memory()
//            .map_err(|e| anyhow!("could not query CUDA free memory: {e}"))
//    }
}

impl GpuExecutor for CudaExecutor {
    fn backend_name(&self) -> &str {
        "cuda"
    }

    fn score_batch(
        &self,
        query: &[u8],
        subjects: &[Vec<u8>],
        scoring: &Scoring,
        gap_open: i32,
        gap_extend: i32,
    ) -> Result<Vec<GpuScore>> {
        if query.is_empty() {
            return Ok(Vec::new());
        }
        if query.len() > MAX_QUERY_LEN {
            anyhow::bail!(
                "query length {} exceeds CUDA kernel limit {}; use CPU path for this query",
                query.len(),
                MAX_QUERY_LEN
            );
        }
        if subjects.is_empty() {
            return Ok(Vec::new());
        }

        // 1. Sort subjects by length (descending) for warp coherence.
        //    Each warp of 32 threads runs in lockstep; if subjects within a
        //    warp have wildly different lengths, the warp stalls on the
        //    longest. Sorting groups similar-length subjects together.
        let mut order: Vec<usize> = (0..subjects.len()).collect();
        order.sort_unstable_by_key(|&i| std::cmp::Reverse(subjects[i].len()));

        let sorted_subjects: Vec<&Vec<u8>> = order.iter().map(|&i| &subjects[i]).collect();

        // 2. Dispatch on scoring kind.
        let scores_sorted = match scoring {
            Scoring::Matrix(matrix) => self.run_protein(
                query,
                &sorted_subjects,
                matrix.as_ref(),
                gap_open,
                gap_extend,
            )?,
            Scoring::NucleotideSimple {
                match_score,
                mismatch,
            } => self.run_nucleotide(
                query,
                &sorted_subjects,
                *match_score,
                *mismatch,
                gap_open,
                gap_extend,
            )?,
        };

        // 3. Un-sort scores back into caller order.
        let mut out = vec![
            GpuScore {
                subject_idx: 0,
                score: 0,
            };
            subjects.len()
        ];
        for (sorted_idx, &orig_idx) in order.iter().enumerate() {
            out[orig_idx] = GpuScore {
                subject_idx: orig_idx,
                score: scores_sorted[sorted_idx],
            };
        }
        Ok(out)
    }
}

impl CudaExecutor {
    fn run_protein(
        &self,
        query: &[u8],
        subjects: &[&Vec<u8>],
        matrix: &crate::scoring::matrices::ScoreMatrix,
        gap_open: i32,
        gap_extend: i32,
    ) -> Result<Vec<i32>> {
        let query_enc = query.to_vec();
        let query_len = query_enc.len() as i32;
        let n_subj = subjects.len() as i32;

        let (subj_packed, subj_offsets) = pack_sequences(subjects, |s| s.to_vec());
        let sub_matrix = flatten_matrix(matrix);
        debug_assert_eq!(sub_matrix.len(), PROTEIN_ALPHABET_SIZE * PROTEIN_ALPHABET_SIZE);

        // Allocate + copy.
        let work_size = (subjects.len() * query_enc.len()) as usize;

        let d_query: CudaSlice<u8> = self.device.htod_copy(query_enc)?;
        let d_subjects: CudaSlice<u8> = self.device.htod_copy(subj_packed)?;
        let d_offsets: CudaSlice<u32> = self.device.htod_copy(subj_offsets)?;
        let d_sub_matrix: CudaSlice<i8> = self.device.htod_copy(sub_matrix)?;

        // alloc_zeros gets us NEG_INF-free init; the kernel zeros the rows itself,
        // so we don't actually need this to be zeroed — but alloc_zeros is the
        // simplest safe API in cudarc 0.9.
        let d_work_h: CudaSlice<i16> = self.device.alloc_zeros::<i16>(work_size)?;
        let d_work_e: CudaSlice<i16> = self.device.alloc_zeros::<i16>(work_size)?;
        let mut d_scores: CudaSlice<i32> = self.device.alloc_zeros::<i32>(subjects.len())?;

        // Launch config.
        let grid = (n_subj as u32 + self.block_size - 1) / self.block_size;
        let cfg = LaunchConfig {
            grid_dim: (grid, 1, 1),
            block_dim: (self.block_size, 1, 1),
            shared_mem_bytes: 0, // kernel declares its own __shared__ array
        };

        // SAFETY: kernel signature matches the argument tuple below.
        // If you change the .cu kernel signature, you MUST change this tuple too.
        unsafe {
            self.protein_kernel.clone().launch(
                cfg,
                (
                    &d_query,
                    query_len,
                    &d_subjects,
                    &d_offsets,
                    n_subj,
                    &d_sub_matrix,
                    gap_open,
                    gap_extend,
                    &d_work_h,
                    &d_work_e,
                    &mut d_scores,
                ),
            )?;
        }

        // Synchronize + copy back.
        self.device.synchronize()?;
        let scores: Vec<i32> = self.device.dtoh_sync_copy(&d_scores)?;

        Ok(scores)
    }

    fn run_nucleotide(
        &self,
        query: &[u8],
        subjects: &[&Vec<u8>],
        match_score: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Result<Vec<i32>> {
        let query_enc = encode_nucleotide(query);
        let query_len = query_enc.len() as i32;
        let n_subj = subjects.len() as i32;

        let (subj_packed, subj_offsets) = pack_sequences(subjects, |s| encode_nucleotide(s));

        let d_query: CudaSlice<u8> = self.device.htod_copy(query_enc.clone())?;
        let d_subjects: CudaSlice<u8> = self.device.htod_copy(subj_packed)?;
        let d_offsets: CudaSlice<u32> = self.device.htod_copy(subj_offsets)?;

        let work_size = subjects.len() * query_enc.len();
        let d_work_h: CudaSlice<i16> = self.device.alloc_zeros::<i16>(work_size)?;
        let d_work_e: CudaSlice<i16> = self.device.alloc_zeros::<i16>(work_size)?;
        let mut d_scores: CudaSlice<i32> = self.device.alloc_zeros::<i32>(subjects.len())?;

        let grid = (n_subj as u32 + self.block_size - 1) / self.block_size;
        let cfg = LaunchConfig {
            grid_dim: (grid, 1, 1),
            block_dim: (self.block_size, 1, 1),
            shared_mem_bytes: 0,
        };

        unsafe {
            self.nucleotide_kernel.clone().launch(
                cfg,
                (
                    &d_query,
                    query_len,
                    &d_subjects,
                    &d_offsets,
                    n_subj,
                    match_score,
                    mismatch,
                    gap_open,
                    gap_extend,
                    &d_work_h,
                    &d_work_e,
                    &mut d_scores,
                ),
            )?;
        }

        self.device.synchronize()?;
        let scores: Vec<i32> = self.device.dtoh_sync_copy(&d_scores)?;
        Ok(scores)
    }
}
