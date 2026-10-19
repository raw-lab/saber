//! Optional accelerator backends for the region-scoring stage.
//!
//! The default build does NOT compile the CUDA code. To enable it:
//!
//! ```text
//! cargo build --release --features gpu-cuda
//! ```
//!
//! ## Contract (v1.2)
//!
//! [`ScoreBackend::score_windows`] scores one query against many subject
//! windows (local Smith-Waterman, affine gaps, score only):
//!
//! * sequences are **canonical alignment codes** ([`crate::alphabet`]):
//!   protein `ARNDCQEGHILKMFPSTWYVBZX*` → 0..23, nucleotide `ACGTN` → 0..4.
//!   Backends never re-encode; codes ≥ 24 are scored as X (22);
//! * gap convention is SABER's: a gap of length k costs `open + (k-1)·extend`;
//! * `bias[i]` (optional) is added to every aligned column that involves
//!   query position `i` (composition-bias correction);
//! * scores are **exact** `i32` values identical to the CPU kernel
//!   ([`crate::algorithm::sw_simd::score_many`]); a backend must return an
//!   error rather than a saturated or approximate score.
//!
//! The search engine ([`crate::search::engine::Engine::with_backend`]) calls
//! the backend once per query with the candidate windows; on any error the
//! batch is re-scored on the CPU and the fallback is counted.
//!
//! The legacy [`GpuExecutor`] API (ASCII in, scores out) is kept for library
//! users and is implemented on top of the same backends, so the CPU and CUDA
//! executors now share one input contract (v1.1 encoded protein input twice
//! on one path and not at all on another).

pub mod cpu_fallback;
pub mod encoding;

#[cfg(feature = "gpu-cuda")]
pub mod cuda;

use std::sync::Arc;

use crate::algorithm::Scoring;
use crate::scoring::stats::GapCosts;

/// Exact local-alignment scoring of one query against many windows.
pub trait ScoreBackend: Send + Sync {
    /// Backend name ("cuda", "cpu-simd", ...).
    fn name(&self) -> &str;

    /// One score per window, in input order. See the module docs for the
    /// encoding, gap and bias conventions.
    fn score_windows(
        &self,
        query: &[u8],
        bias: Option<&[i8]>,
        windows: &[&[u8]],
        matrix: &[[i32; 24]; 24],
        gaps: GapCosts,
    ) -> anyhow::Result<Vec<i32>>;
}

/// Probe for a GPU score backend (CUDA device 0). `None` when the binary was
/// built without `gpu-cuda` or no device could be initialised.
#[cfg(feature = "gpu-cuda")]
pub fn try_init_backend() -> Option<Arc<dyn ScoreBackend>> {
    match cuda::CudaBackend::new() {
        Ok(b) => Some(Arc::new(b)),
        Err(e) => {
            log::warn!("CUDA initialisation failed: {e}");
            None
        }
    }
}

#[cfg(not(feature = "gpu-cuda"))]
pub fn try_init_backend() -> Option<Arc<dyn ScoreBackend>> {
    None
}

// ---------------------------------------------------------------------------
// Legacy API (v1.0/v1.1 library compatibility)
// ---------------------------------------------------------------------------

/// One scored (query, subject) pair.
#[derive(Debug, Clone)]
pub struct GpuScore {
    /// Index of the subject in the input array.
    pub subject_idx: usize,
    /// Raw Smith-Waterman score (identical to the CPU kernel).
    pub score: i32,
}

/// Legacy executor interface. **Inputs are ASCII sequences** for every
/// implementation; they are encoded once (protein or nucleotide according to
/// `scoring`) and passed to a [`ScoreBackend`].
pub trait GpuExecutor: Send + Sync {
    fn backend_name(&self) -> &str;

    fn score_batch(
        &self,
        query: &[u8],
        subjects: &[Vec<u8>],
        scoring: &Scoring,
        gap_open: i32,
        gap_extend: i32,
    ) -> anyhow::Result<Vec<GpuScore>>;
}

/// Adapter: legacy ASCII executor over any [`ScoreBackend`].
pub struct BackendExecutor {
    backend: Arc<dyn ScoreBackend>,
}

impl BackendExecutor {
    pub fn new(backend: Arc<dyn ScoreBackend>) -> Self {
        Self { backend }
    }
}

impl GpuExecutor for BackendExecutor {
    fn backend_name(&self) -> &str {
        self.backend.name()
    }

    fn score_batch(
        &self,
        query: &[u8],
        subjects: &[Vec<u8>],
        scoring: &Scoring,
        gap_open: i32,
        gap_extend: i32,
    ) -> anyhow::Result<Vec<GpuScore>> {
        let (matrix, q, subj): ([[i32; 24]; 24], Vec<u8>, Vec<Vec<u8>>) = match scoring {
            Scoring::Matrix(m) => (
                m.matrix,
                crate::alphabet::EncodedProtein::from_ascii(query).as_slice().to_vec(),
                subjects.iter().map(|s| crate::alphabet::EncodedProtein::from_ascii(s).as_slice().to_vec()).collect(),
            ),
            Scoring::NucleotideSimple { match_score, mismatch } => (
                crate::scoring::matrices::ScoreMatrix::nucleotide(*match_score, *mismatch).matrix,
                crate::alphabet::encode_nucleotide(query),
                subjects.iter().map(|s| crate::alphabet::encode_nucleotide(s)).collect(),
            ),
        };
        if gap_extend < 1 || gap_open < gap_extend {
            anyhow::bail!("invalid gap penalties open={gap_open} extend={gap_extend}");
        }
        let gaps = GapCosts { open: gap_open, extend: gap_extend };
        let views: Vec<&[u8]> = subj.iter().map(|s| s.as_slice()).collect();
        let scores = self.backend.score_windows(&q, None, &views, &matrix, gaps)?;
        Ok(scores.into_iter().enumerate().map(|(i, score)| GpuScore { subject_idx: i, score }).collect())
    }
}

/// Legacy probe: a CUDA-backed [`GpuExecutor`], or `None`.
pub fn try_init_gpu() -> Option<Arc<dyn GpuExecutor>> {
    try_init_backend().map(|b| Arc::new(BackendExecutor::new(b)) as Arc<dyn GpuExecutor>)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::reference::sw_reference;
    use crate::scoring::{ScoreMatrix, ScoringSystem};

    #[test]
    fn legacy_cpu_executor_takes_ascii_for_both_kinds() {
        let ex = cpu_fallback::CpuExecutor::new();
        let m = std::sync::Arc::new(ScoreMatrix::new(ScoringSystem::Blosum62));
        let q = b"MKTAYIAKQRQISFVKSHFSRQ".to_vec();
        let subs = vec![b"MKTAYIAKQRQISFVKSHFSRQ".to_vec(), b"PPPPGGGG".to_vec(), b"mktayiakqr".to_vec()];
        let got = ex.score_batch(&q, &subs, &Scoring::Matrix(m.clone()), 11, 1).unwrap();
        let qc = crate::alphabet::EncodedProtein::from_ascii(&q);
        for (i, s) in subs.iter().enumerate() {
            let sc = crate::alphabet::EncodedProtein::from_ascii(s);
            let (r, _, _) = sw_reference(qc.as_slice(), sc.as_slice(), &m.matrix, None, 11, 1);
            assert_eq!(got[i].score as i64, r);
        }
        let n = ex
            .score_batch(b"ACGTACGTAC", &[b"ACGTACGTAC".to_vec(), b"acgtnnacgt".to_vec()], &Scoring::NucleotideSimple { match_score: 2, mismatch: -3 }, 7, 2)
            .unwrap();
        assert_eq!(n[0].score, 20);
        assert_eq!(n[1].score, 8);
    }
}
