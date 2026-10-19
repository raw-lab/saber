//! CPU implementation of [`ScoreBackend`] (AVX2 inter-sequence kernel with
//! exact i32 fallback), plus the legacy [`CpuExecutor`].
//!
//! This is what `--gpu` uses when no CUDA device is available, and it is the
//! reference the CUDA backend must match bit for bit.

use std::cell::RefCell;
use std::sync::Arc;

use crate::algorithm::sw_simd::{score_many, KernelCounters, SimdScratch};
use crate::alphabet::ProteinCodes;
use crate::scoring::stats::GapCosts;

use super::{BackendExecutor, GpuExecutor, GpuScore, ScoreBackend};

thread_local! {
    static SCRATCH: RefCell<SimdScratch> = RefCell::new(SimdScratch::new());
}

/// Exact CPU score backend.
#[derive(Debug, Default, Clone, Copy)]
pub struct CpuBackend;

impl ScoreBackend for CpuBackend {
    fn name(&self) -> &str {
        "cpu-simd"
    }

    fn score_windows(
        &self,
        query: &[u8],
        bias: Option<&[i8]>,
        windows: &[&[u8]],
        matrix: &[[i32; 24]; 24],
        gaps: GapCosts,
    ) -> anyhow::Result<Vec<i32>> {
        if let Some(b) = bias {
            anyhow::ensure!(b.len() == query.len(), "bias length {} != query length {}", b.len(), query.len());
        }
        let views: Vec<ProteinCodes<'_>> = windows.iter().map(|w| ProteinCodes::new_unchecked(w)).collect();
        let mut kc = KernelCounters::default();
        Ok(SCRATCH.with(|s| score_many(ProteinCodes::new_unchecked(query), bias, &views, matrix, gaps, &mut s.borrow_mut(), &mut kc)))
    }
}

/// Legacy CPU executor (ASCII input; see [`GpuExecutor`]).
pub struct CpuExecutor {
    inner: BackendExecutor,
}

impl CpuExecutor {
    pub fn new() -> Self {
        Self { inner: BackendExecutor::new(Arc::new(CpuBackend)) }
    }
}

impl Default for CpuExecutor {
    fn default() -> Self {
        Self::new()
    }
}

impl GpuExecutor for CpuExecutor {
    fn backend_name(&self) -> &str {
        "cpu-simd"
    }

    fn score_batch(
        &self,
        query: &[u8],
        subjects: &[Vec<u8>],
        scoring: &crate::algorithm::Scoring,
        gap_open: i32,
        gap_extend: i32,
    ) -> anyhow::Result<Vec<GpuScore>> {
        self.inner.score_batch(query, subjects, scoring, gap_open, gap_extend)
    }
}
