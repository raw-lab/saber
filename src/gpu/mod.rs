//! GPU acceleration for Smith-Waterman alignment.
//!
//! The default build does NOT compile the GPU code. To enable it:
//!
//! ```text
//! cargo build --release --features gpu-cuda
//! ```
//!
//! At runtime, `try_init_gpu()` probes for a CUDA device. If none is
//! available, the caller can fall back to [`cpu_fallback::CpuExecutor`],
//! which runs the same Smith-Waterman through `rayon` on every core.
//!
//! ## Architecture
//!
//! 1. Host encodes query + subjects to compact byte arrays
//!    ([`encoding`]).
//! 2. (Optional) Host runs a k-mer prefilter on GPU to reduce the candidate
//!    subject set — closes the throughput gap to DIAMOND/MMseqs2 by avoiding
//!    full SW on obvious non-hits.
//! 3. Subjects are sorted by length on the host to minimize warp divergence.
//! 4. Host launches the [`saber_sw_protein_score`] or
//!    [`saber_sw_nucleotide_score`] kernel — one CUDA thread per
//!    (query, subject) pair, scoring only (no traceback on the GPU).
//! 5. Host runs CPU traceback only on the top-K hits using the existing
//!    [`crate::algorithm::SmithWaterman`] implementation.
//!
//! ## Why score-only on GPU?
//!
//! Smith-Waterman traceback requires the full O(m·n) DP matrix and a
//! divergent backtracking walk that's poorly suited to GPU. Standard
//! practice (CUDASW++, SW#) is to compute scores on GPU and re-align
//! on CPU for the small subset of subjects that pass the E-value filter.
//! For a 100k-sequence DB with ~50 hits above threshold, this means
//! 50 CPU re-alignments instead of 100k full traceback walks on GPU.

pub mod encoding;
pub mod cpu_fallback;

#[cfg(feature = "gpu-cuda")]
pub mod cuda;

use std::sync::Arc;

use crate::algorithm::Scoring;

/// One scored (query, subject) pair from the GPU.
#[derive(Debug, Clone)]
pub struct GpuScore {
    /// Index of the subject in the input array.
    pub subject_idx: usize,
    /// Raw Smith-Waterman score (matches what the CPU kernel would return).
    pub score: i32,
}

/// Trait implemented by GPU executors. Both [`cpu_fallback::CpuExecutor`]
/// and [`cuda::CudaExecutor`] (under feature flag) implement it.
pub trait GpuExecutor: Send + Sync {
    /// Human-readable backend name ("cuda", "cpu-fallback", ...).
    fn backend_name(&self) -> &str;

    /// Compute SW scores for one query against many subjects.
    ///
    /// Returns one [`GpuScore`] per input subject, in the same order.
    fn score_batch(
        &self,
        query: &[u8],
        subjects: &[Vec<u8>],
        scoring: &Scoring,
        gap_open: i32,
        gap_extend: i32,
    ) -> anyhow::Result<Vec<GpuScore>>;
}

/// Probe for an available GPU backend. Returns `None` and the caller
/// should use [`cpu_fallback::CpuExecutor::new`] instead.
#[cfg(feature = "gpu-cuda")]
pub fn try_init_gpu() -> Option<Arc<dyn GpuExecutor>> {
    match cuda::CudaExecutor::new() {
        Ok(exec) => Some(Arc::new(exec)),
        Err(e) => {
            log::warn!("CUDA GPU initialization failed: {e}. Falling back to CPU.");
            None
        }
    }
}

#[cfg(not(feature = "gpu-cuda"))]
pub fn try_init_gpu() -> Option<Arc<dyn GpuExecutor>> {
    log::debug!("Built without --features gpu-cuda; no GPU available.");
    None
}
