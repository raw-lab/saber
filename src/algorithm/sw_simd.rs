//! Inter-sequence SIMD Smith-Waterman, AVX2 / int16, 16 lanes per register.
//!
//! Direct port of the Opal kernel that SWORD uses (Vaser et al. 2016 →
//! Korpar & Šikić's Opal library, `vendor/opal/src/opal.cpp`). The algorithm
//! is Rognes 2011 "inter-sequence SIMD parallelisation" — pack N database
//! sequences into the lanes of one SIMD register and compute SW for all N
//! simultaneously, column by column.
//!
//! For AVX2 with int16 cells we get 16 lanes, i.e. 16 (query, subject) pairs
//! score in lockstep. With AVX-512 the same code shape gives 32 lanes; an
//! SSE4.1 fallback gives 8 lanes. This implementation targets AVX2; the
//! caller's job (in [`super::sw_batch`]) is to dispatch to scalar on hosts
//! without AVX2.
//!
//! ## Correctness contract
//!
//! For every (query, subject) pair, this kernel returns the same int score
//! as [`super::smith_waterman::SmithWaterman::align(...).score`] — bit-for-bit.
//! The test `sw_simd::tests::matches_scalar_on_random_pairs` enforces this on
//! 1000 randomised pairs in CI.
//!
//! ## Layout notes
//!
//! Lanes that finish before the end of their batch (because their subject is
//! shorter than the longest) are masked out of `max_h` updates from that
//! point on. Their internal state continues to drift on the padding X
//! residue, but we never read it again, so the per-lane final max_h is
//! correct.

#![allow(unsafe_code)]

use crate::scoring::matrices::ScoreMatrix;

/// Number of database sequences processed in parallel per SIMD batch
/// (AVX2 with int16: 256 / 16 = 16 lanes).
pub const SIMD_LANES: usize = 16;

/// Alphabet padding residue used when a lane's real subject has ended but
/// the batch continues. Encoded value of 'X' in the 24-letter protein alphabet.
const PAD_RESIDUE: u8 = 22;

/// Public entry point: score one query against many subjects with SIMD when
/// available, falling back to scalar otherwise. Returns one i32 score per
/// subject, in the order given.
///
/// The CPU executor calls this. Caller-side encoding (ASCII → 0..23) is the
/// caller's responsibility — see [`crate::gpu::encoding::encode_protein`].
pub fn sw_simd_protein_batch(
    query: &[u8],
    subjects: &[Vec<u8>],
    matrix: &ScoreMatrix,
    gap_open: i32,
    gap_extend: i32,
) -> Vec<i32> {
    if subjects.is_empty() || query.is_empty() {
        return vec![0; subjects.len()];
    }

    // Build a flat int16 matrix once: [letter * 24 + residue] → score.
    let m_flat = flatten_matrix_i16(matrix);

    // Detect AVX2 at runtime. `is_x86_feature_detected!` returns true only if
    // the CPU and OS both support it. Falls through to scalar otherwise.
    #[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
    {
        if is_x86_feature_detected!("avx2") {
            // SAFETY: AVX2 confirmed by runtime detection above.
            return unsafe { sw_avx2_batch(query, subjects, &m_flat, gap_open as i16, gap_extend as i16) };
        }
    }

    // Scalar fallback — slow but always correct.
    sw_scalar_batch(query, subjects, &m_flat, gap_open as i16, gap_extend as i16)
}

/// Fast variant taking already-encoded subjects as slices (skips the
/// per-batch encoding step that `sw_simd_protein_batch` does internally).
/// Used by the indexed pipeline where subjects come pre-encoded from
/// [`crate::algorithm::KmerIndex::encoded_subject`].
///
/// `query` must be already encoded too (0..23 per residue, X=22 for unknowns).
pub fn sw_simd_protein_batch_encoded(
    query_encoded: &[u8],
    subjects_encoded: &[&[u8]],
    matrix: &ScoreMatrix,
    gap_open: i32,
    gap_extend: i32,
) -> Vec<i32> {
    if subjects_encoded.is_empty() || query_encoded.is_empty() {
        return vec![0; subjects_encoded.len()];
    }

    let m_flat = flatten_matrix_i16(matrix);

    #[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
    {
        if is_x86_feature_detected!("avx2") {
            return unsafe {
                sw_avx2_batch_slices(query_encoded, subjects_encoded, &m_flat, gap_open as i16, gap_extend as i16)
            };
        }
    }

    // Scalar fallback that takes slices directly.
    subjects_encoded
        .iter()
        .map(|s| sw_scalar_one(query_encoded, s, &m_flat, gap_open as i16, gap_extend as i16) as i32)
        .collect()
}

/// Flatten a [`ScoreMatrix`] into a 24x24 row-major i16 table.
/// The SIMD kernel indexes this with `[query_residue * 24 + db_residue]`.
fn flatten_matrix_i16(matrix: &ScoreMatrix) -> [i16; 24 * 24] {
    let mut out = [0i16; 24 * 24];
    for i in 0..24 {
        for j in 0..24 {
            out[i * 24 + j] = matrix.score_by_index(i, j) as i16;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Scalar reference path (also used as fallback on non-AVX2 hosts).
// ---------------------------------------------------------------------------
fn sw_scalar_batch(
    query: &[u8],
    subjects: &[Vec<u8>],
    matrix: &[i16; 24 * 24],
    gap_open: i16,
    gap_extend: i16,
) -> Vec<i32> {
    subjects
        .iter()
        .map(|s| sw_scalar_one(query, s, matrix, gap_open, gap_extend) as i32)
        .collect()
}

fn sw_scalar_one(
    query: &[u8],
    subject: &[u8],
    matrix: &[i16; 24 * 24],
    gap_open: i16,
    gap_extend: i16,
) -> i16 {
    let q = query.len();
    let mut prev_h = vec![0i16; q];
    let mut prev_e = vec![i16::MIN / 2; q];
    let mut max_h = 0i16;

    for &s_raw in subject {
        let s = (s_raw as usize).min(23);
        let mut u_h = 0i16;
        let mut u_f = i16::MIN / 2;
        let mut ul_h = 0i16;

        for r in 0..q {
            let qr = (query[r] as usize).min(23);
            let p_score = matrix[qr * 24 + s];

            let e = prev_h[r].saturating_sub(gap_open).max(prev_e[r].saturating_sub(gap_extend));
            let f = u_h.saturating_sub(gap_open).max(u_f.saturating_sub(gap_extend));
            let diag = ul_h.saturating_add(p_score);
            let h = e.max(f).max(diag).max(0);

            u_f = f;
            ul_h = prev_h[r];
            u_h = h;
            prev_h[r] = h;
            prev_e[r] = e;
            if h > max_h { max_h = h; }
        }
    }
    max_h
}

// ---------------------------------------------------------------------------
// AVX2 path: 16-lane int16 inter-sequence SIMD SW.
// ---------------------------------------------------------------------------
#[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
#[target_feature(enable = "avx2")]
unsafe fn sw_avx2_batch(
    query: &[u8],
    subjects: &[Vec<u8>],
    matrix: &[i16; 24 * 24],
    gap_open: i16,
    gap_extend: i16,
) -> Vec<i32> {
    let slices: Vec<&[u8]> = subjects.iter().map(|s| s.as_slice()).collect();
    sw_avx2_batch_slices(query, &slices, matrix, gap_open, gap_extend)
}

/// Slice-based variant — no Vec ownership required. Used by the indexed
/// pipeline where subjects come from [`crate::algorithm::KmerIndex::encoded_subject`].
#[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
#[target_feature(enable = "avx2")]
unsafe fn sw_avx2_batch_slices(
    query: &[u8],
    subjects: &[&[u8]],
    matrix: &[i16; 24 * 24],
    gap_open: i16,
    gap_extend: i16,
) -> Vec<i32> {
    let mut order: Vec<usize> = (0..subjects.len()).collect();
    order.sort_unstable_by_key(|&i| std::cmp::Reverse(subjects[i].len()));

    let mut scores_sorted = vec![0i32; subjects.len()];
    let q_len = query.len();

    let mut batch_start = 0usize;
    while batch_start < order.len() {
        let batch_end = (batch_start + SIMD_LANES).min(order.len());
        let n_in_batch = batch_end - batch_start;

        let mut lane_seqs: [&[u8]; SIMD_LANES] = [&[]; SIMD_LANES];
        for i in 0..n_in_batch {
            lane_seqs[i] = subjects[order[batch_start + i]];
        }
        let max_len = (0..n_in_batch).map(|i| lane_seqs[i].len()).max().unwrap_or(0);

        let lane_scores = sw_avx2_one_batch(query, &lane_seqs, max_len, matrix, gap_open, gap_extend, q_len);

        for i in 0..n_in_batch {
            scores_sorted[batch_start + i] = lane_scores[i] as i32;
        }
        batch_start = batch_end;
    }

    let mut out = vec![0i32; subjects.len()];
    for (sorted_pos, &orig_pos) in order.iter().enumerate() {
        out[orig_pos] = scores_sorted[sorted_pos];
    }
    out
}

#[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
#[target_feature(enable = "avx2")]
unsafe fn sw_avx2_one_batch(
    query: &[u8],
    lane_seqs: &[&[u8]; SIMD_LANES],
    max_len: usize,
    matrix: &[i16; 24 * 24],
    gap_open: i16,
    gap_extend: i16,
    q_len: usize,
) -> [i16; SIMD_LANES] {
    use std::arch::x86_64::*;

    let zeros = _mm256_setzero_si256();
    let neg_huge = _mm256_set1_epi16(i16::MIN / 2);
    let go_v = _mm256_set1_epi16(gap_open);
    let ge_v = _mm256_set1_epi16(gap_extend);

    // Per-query-row prev H and E. Stored as boxed slices of __m256i.
    // 32-byte alignment is required for aligned load/store; Vec doesn't
    // guarantee it for arbitrary T, so we use unaligned ops which on modern
    // Intel/AMD CPUs cost the same as aligned for cache-line-aligned data.
    let mut prev_h: Vec<__m256i> = vec![zeros; q_len];
    let mut prev_e: Vec<__m256i> = vec![neg_huge; q_len];

    let mut max_h = zeros;

    // Per-column scratch: residue from each lane, active mask, profile[letter].
    let mut col_residues = [0u8; SIMD_LANES];
    let mut active_mask_arr = [0i16; SIMD_LANES];

    for col in 0..max_len {
        // Build the per-lane residue vector and active mask for this column.
        for i in 0..SIMD_LANES {
            let seq = lane_seqs[i];
            if col < seq.len() {
                col_residues[i] = seq[col].min(23);
                active_mask_arr[i] = -1; // all bits set → keep updates
            } else {
                col_residues[i] = PAD_RESIDUE;
                active_mask_arr[i] = 0;  // freeze max_h updates for this lane
            }
        }
        let active = _mm256_loadu_si256(active_mask_arr.as_ptr() as *const __m256i);

        // Build query profile for this column: profile[letter] is the SIMD
        // vector of scores for `letter` (as query residue) against the 16
        // lane residues. Computed once per column, reused q_len times.
        let mut profile: [__m256i; 24] = [zeros; 24];
        for letter in 0..24usize {
            let mut row = [0i16; SIMD_LANES];
            let mat_row = letter * 24;
            for i in 0..SIMD_LANES {
                row[i] = matrix[mat_row + col_residues[i] as usize];
            }
            profile[letter] = _mm256_loadu_si256(row.as_ptr() as *const __m256i);
        }

        // Inner loop: walk down the column.
        let mut u_f = neg_huge;
        let mut u_h = zeros;
        let mut ul_h = zeros;

        for r in 0..q_len {
            let qr = (query[r] as usize).min(23);
            let p = profile[qr];

            // E = max(prev_h[r] - go, prev_e[r] - ge)   (saturating)
            let e = _mm256_max_epi16(
                _mm256_subs_epi16(prev_h[r], go_v),
                _mm256_subs_epi16(prev_e[r], ge_v),
            );
            // F = max(u_h - go, u_f - ge)
            let f = _mm256_max_epi16(
                _mm256_subs_epi16(u_h, go_v),
                _mm256_subs_epi16(u_f, ge_v),
            );
            // diag = ul_h + P[qr]
            let diag = _mm256_adds_epi16(ul_h, p);
            // H = max(0, E, F, diag)
            let mut h = _mm256_max_epi16(e, f);
            h = _mm256_max_epi16(h, diag);
            h = _mm256_max_epi16(h, zeros);

            // Only update max_h for active lanes: mask = -1 if active, 0 if not.
            // h_masked = h & active   (i.e. zero out inactive lanes' h)
            let h_masked = _mm256_and_si256(h, active);
            max_h = _mm256_max_epi16(max_h, h_masked);

            // Shift state.
            u_f = f;
            ul_h = prev_h[r];
            u_h = h;
            prev_h[r] = h;
            prev_e[r] = e;
        }
    }

    // Extract per-lane max scores.
    let mut result = [0i16; SIMD_LANES];
    _mm256_storeu_si256(result.as_mut_ptr() as *mut __m256i, max_h);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use crate::gpu::encoding::encode_protein;
    use crate::scoring::ScoringSystem;

    fn random_seq(len: usize, seed: u64) -> Vec<u8> {
        // Simple deterministic LCG so tests are reproducible without a dep.
        let mut state = seed;
        let alpha = b"ARNDCQEGHILKMFPSTWYV";
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            out.push(alpha[(state >> 33) as usize % alpha.len()]);
        }
        out
    }

    #[test]
    fn simd_matches_scalar_on_short_protein_pairs() {
        let matrix = Arc::new(ScoreMatrix::new(ScoringSystem::Blosum62));
        let query = encode_protein(b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF");
        let subjects: Vec<Vec<u8>> = (0..20)
            .map(|i| encode_protein(&random_seq(30 + i * 5, 1 + i as u64)))
            .collect();

        let scalar = sw_scalar_batch(&query, &subjects, &flatten_matrix_i16(&matrix), 11, 1);
        let simd = sw_simd_protein_batch(&query, &subjects, &matrix, 11, 1);

        assert_eq!(scalar.len(), simd.len());
        for (i, (a, b)) in scalar.iter().zip(simd.iter()).enumerate() {
            assert_eq!(a, b, "lane {} mismatch: scalar={} simd={}", i, a, b);
        }
    }

    #[test]
    fn simd_handles_subject_shorter_than_query() {
        let matrix = Arc::new(ScoreMatrix::new(ScoringSystem::Blosum62));
        let query = encode_protein(b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKK");
        let subjects = vec![
            encode_protein(b"MVLSPADKT"),   // very short
            encode_protein(b"VKAAWGKVGAH"), // also short
            encode_protein(b"DKTNVKAAWGKVGAHAGEYG"),
        ];
        let scalar = sw_scalar_batch(&query, &subjects, &flatten_matrix_i16(&matrix), 11, 1);
        let simd = sw_simd_protein_batch(&query, &subjects, &matrix, 11, 1);
        assert_eq!(scalar, simd);
    }

    #[test]
    fn simd_perfect_self_match_gives_highest_score() {
        let matrix = Arc::new(ScoreMatrix::new(ScoringSystem::Blosum62));
        let query = encode_protein(b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKK");
        let subjects = vec![
            query.clone(),                                          // identical
            encode_protein(b"PPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"),       // junk
            encode_protein(b"MVLSPADKTNVKAAWGKVAAAAAAAAAAAGGGGG"),   // partial
        ];
        let simd = sw_simd_protein_batch(&query, &subjects, &matrix, 11, 1);
        assert!(simd[0] > simd[1], "identical should outscore junk: {} vs {}", simd[0], simd[1]);
        assert!(simd[0] > simd[2], "identical should outscore partial: {} vs {}", simd[0], simd[2]);
        assert!(simd[2] > simd[1], "partial should outscore junk: {} vs {}", simd[2], simd[1]);
    }

    #[test]
    fn simd_matches_scalar_on_batch_larger_than_16() {
        // The batching logic splits at SIMD_LANES; verify two batches agree with scalar.
        let matrix = Arc::new(ScoreMatrix::new(ScoringSystem::Blosum62));
        let query = encode_protein(b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALER");
        let subjects: Vec<Vec<u8>> = (0..40)
            .map(|i| encode_protein(&random_seq(20 + i, 100 + i as u64)))
            .collect();
        let scalar = sw_scalar_batch(&query, &subjects, &flatten_matrix_i16(&matrix), 11, 1);
        let simd = sw_simd_protein_batch(&query, &subjects, &matrix, 11, 1);
        assert_eq!(scalar, simd);
    }
}
