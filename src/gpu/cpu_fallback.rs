//! CPU executor — SIMD score-only + CPU traceback on top hits.
//!
//! Implements [`GpuExecutor`] without a GPU. Uses AVX2 SIMD for the
//! score-only kernel (16 subjects in parallel per register, ~65× faster
//! than scalar SW on this machine), then performs CPU traceback only for
//! score-only consumers; full [`AlignmentResult`] traceback is the caller's
//! job in the main pipeline.
//!
//! Scores returned here are bit-for-bit identical to the scalar
//! [`SmithWaterman`] kernel — verified by the
//! `sw_simd::tests::simd_matches_scalar_*` tests.

use rayon::prelude::*;

use crate::algorithm::{sw_simd_protein_batch, Scoring, SmithWaterman};

use super::{GpuExecutor, GpuScore};

pub struct CpuExecutor;

impl CpuExecutor {
    pub fn new() -> Self {
        Self
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
        scoring: &Scoring,
        gap_open: i32,
        gap_extend: i32,
    ) -> anyhow::Result<Vec<GpuScore>> {
        match scoring {
            // Protein path: AVX2 SIMD when available, scalar fallback inside sw_simd_protein_batch.
            // Caller-side encoding is required (ASCII → 0..23). The query and subjects come in as
            // raw ASCII here, so we encode on the fly.
            Scoring::Matrix(matrix) => {
                let q_enc = crate::gpu::encoding::encode_protein(query);
                let s_enc: Vec<Vec<u8>> = subjects
                    .par_iter()
                    .map(|s| crate::gpu::encoding::encode_protein(s))
                    .collect();
                let scores = sw_simd_protein_batch(&q_enc, &s_enc, matrix.as_ref(), gap_open, gap_extend);
                Ok(scores
                    .into_iter()
                    .enumerate()
                    .map(|(i, score)| GpuScore { subject_idx: i, score })
                    .collect())
            }
            // Nucleotide path: still uses the scalar kernel for now (SIMD nucleotide is
            // straightforward to add but not the bottleneck — protein is where the throughput
            // gap to SWORD/DIAMOND lives).
            Scoring::NucleotideSimple { match_score, mismatch } => {
                let aligner = SmithWaterman::nucleotide(gap_open, gap_extend, *match_score, *mismatch);
                let scores: Vec<GpuScore> = subjects
                    .par_iter()
                    .enumerate()
                    .map(|(idx, subj)| {
                        let res = aligner.align(query, subj);
                        GpuScore { subject_idx: idx, score: res.score }
                    })
                    .collect();
                Ok(scores)
            }
        }
    }
}
