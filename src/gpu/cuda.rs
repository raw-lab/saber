//! CUDA implementation of [`ScoreBackend`] via `cudarc` (feature `gpu-cuda`).
//!
//! **Status: not executed on a GPU during v1.2 development** (no CUDA device
//! was available). The kernel logic is tested on the CPU by compiling
//! `kernels/sw.cu` as C++ against the CPU reference (`gpu/emulate/`); the
//! host code below type-checks against cudarc 0.9.15 but its runtime
//! behaviour on hardware is unverified. Treat GPU results as experimental
//! until validated with `saber --doctor --gpu` on a CUDA machine.
//!
//! Lifecycle: [`CudaBackend::new`] opens device 0 and compiles the kernel
//! with NVRTC once; each [`ScoreBackend::score_windows`] call sorts windows
//! by length (warp coherence), splits them into batches whose interleaved
//! DP state fits `work_budget_bytes`, launches one thread per window and
//! returns scores in input order. All DP arithmetic is int32 (exact).

use std::sync::{Arc, Mutex};

use anyhow::{anyhow, Context, Result};
use cudarc::driver::{CudaDevice, CudaFunction, CudaSlice, LaunchAsync, LaunchConfig};
use cudarc::nvrtc::compile_ptx;

use crate::scoring::stats::GapCosts;

use super::ScoreBackend;

const KERNEL_SOURCE: &str = include_str!("kernels/sw.cu");
const MODULE_NAME: &str = "saber_sw_v12";
const KERNEL: &str = "saber_sw_score_i32";
const BLOCK_SIZE: u32 = 128;

pub struct CudaBackend {
    device: Arc<CudaDevice>,
    kernel: CudaFunction,
    /// Upper bound for the interleaved H/E work arrays of one launch.
    work_budget_bytes: usize,
    /// Serialises launches from the rayon worker threads.
    lock: Mutex<()>,
}

impl CudaBackend {
    pub fn new() -> Result<Self> {
        Self::new_on_device(0)
    }

    pub fn new_on_device(ordinal: usize) -> Result<Self> {
        let device = CudaDevice::new(ordinal).with_context(|| format!("could not open CUDA device {ordinal}"))?;
        let ptx = compile_ptx(KERNEL_SOURCE).map_err(|e| anyhow!("NVRTC failed to compile sw.cu: {e:?}"))?;
        device.load_ptx(ptx, MODULE_NAME, &[KERNEL]).context("failed to load SABER kernel PTX")?;
        let kernel = device.get_func(MODULE_NAME, KERNEL).ok_or_else(|| anyhow!("kernel `{KERNEL}` not found"))?;
        Ok(Self { device, kernel, work_budget_bytes: 512 << 20, lock: Mutex::new(()) })
    }

    fn launch_batch(
        &self,
        query: &[u8],
        bias: &[i8],
        windows: &[&[u8]],
        matrix: &[i32],
        gaps: GapCosts,
    ) -> Result<Vec<i32>> {
        let n = windows.len();
        let total: usize = windows.iter().map(|w| w.len()).sum();
        anyhow::ensure!(total <= u32::MAX as usize, "batch too large for 32-bit offsets");
        let mut packed = Vec::with_capacity(total);
        let mut offsets = Vec::with_capacity(n + 1);
        offsets.push(0u32);
        for w in windows {
            packed.extend_from_slice(w);
            offsets.push(packed.len() as u32);
        }
        if packed.is_empty() {
            packed.push(22); // cudarc rejects zero-length copies
        }
        let work = n.checked_mul(query.len()).ok_or_else(|| anyhow!("work size overflow"))?;
        let d_query: CudaSlice<u8> = self.device.htod_copy(query.to_vec())?;
        let d_bias: CudaSlice<i8> = self.device.htod_copy(bias.to_vec())?;
        let d_subjects: CudaSlice<u8> = self.device.htod_copy(packed)?;
        let d_offsets: CudaSlice<u32> = self.device.htod_copy(offsets)?;
        let d_matrix: CudaSlice<i32> = self.device.htod_copy(matrix.to_vec())?;
        let d_h: CudaSlice<i32> = self.device.alloc_zeros::<i32>(work)?;
        let d_e: CudaSlice<i32> = self.device.alloc_zeros::<i32>(work)?;
        let mut d_scores: CudaSlice<i32> = self.device.alloc_zeros::<i32>(n)?;
        let grid = (n as u32 + BLOCK_SIZE - 1) / BLOCK_SIZE;
        let cfg = LaunchConfig { grid_dim: (grid, 1, 1), block_dim: (BLOCK_SIZE, 1, 1), shared_mem_bytes: 0 };
        // SAFETY: the tuple matches the kernel signature in kernels/sw.cu.
        unsafe {
            self.kernel.clone().launch(
                cfg,
                (
                    &d_query,
                    query.len() as i32,
                    &d_bias,
                    &d_subjects,
                    &d_offsets,
                    n as i32,
                    &d_matrix,
                    gaps.open,
                    gaps.extend,
                    &d_h,
                    &d_e,
                    &mut d_scores,
                ),
            )?;
        }
        self.device.synchronize()?;
        Ok(self.device.dtoh_sync_copy(&d_scores)?)
    }
}

impl ScoreBackend for CudaBackend {
    fn name(&self) -> &str {
        "cuda"
    }

    fn score_windows(
        &self,
        query: &[u8],
        bias: Option<&[i8]>,
        windows: &[&[u8]],
        matrix: &[[i32; 24]; 24],
        gaps: GapCosts,
    ) -> Result<Vec<i32>> {
        if windows.is_empty() {
            return Ok(Vec::new());
        }
        if query.is_empty() {
            return Ok(vec![0; windows.len()]);
        }
        let zeros;
        let bias = match bias {
            Some(b) => {
                anyhow::ensure!(b.len() == query.len(), "bias length mismatch");
                b
            }
            None => {
                zeros = vec![0i8; query.len()];
                &zeros
            }
        };
        let flat: Vec<i32> = matrix.iter().flat_map(|r| r.iter().copied()).collect();
        let mut order: Vec<usize> = (0..windows.len()).collect();
        order.sort_unstable_by_key(|&i| std::cmp::Reverse(windows[i].len()));
        let per_launch = (self.work_budget_bytes / (query.len() * 8)).clamp(1, 1 << 20);
        let _g = self.lock.lock().map_err(|_| anyhow!("CUDA backend lock poisoned"))?;
        let mut out = vec![0i32; windows.len()];
        for chunk in order.chunks(per_launch) {
            let ws: Vec<&[u8]> = chunk.iter().map(|&i| windows[i]).collect();
            let sc = self.launch_batch(query, bias, &ws, &flat, gaps)?;
            anyhow::ensure!(sc.len() == chunk.len(), "kernel returned {} scores for {} windows", sc.len(), chunk.len());
            for (k, &i) in chunk.iter().enumerate() {
                out[i] = sc[k];
            }
        }
        Ok(out)
    }
}
