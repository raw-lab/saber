//! Ungapped X-drop extension — BLAST's HSP filter.
//!
//! Between seeding and gapped Smith-Waterman, BLAST extends each seed
//! ungapped (along a single diagonal, no gaps allowed) using the X-drop
//! algorithm to find the best high-scoring segment pair (HSP). Only
//! candidates whose HSP score exceeds a threshold get the expensive
//! gapped SW. This is what lets BLAST dismiss 99% of seed-finding
//! candidates without running gapped DP.
//!
//! ## Why it works
//!
//! True homologs have stretches of high local similarity even before
//! you allow gaps — that's what ungapped extension captures. Random
//! sequences that happen to share a few k-mers don't extend well:
//! the score drops fast once you leave the immediate seed neighborhood.
//! So ungapped HSP score is a strong predictor of whether gapped SW
//! will score above threshold, at a tiny fraction of the cost.
//!
//! ## Algorithm (per seed)
//!
//! Start at the seed position (q0, s0). Walk forward along the diagonal
//! `d = s0 - q0` (i.e., increment both q and s in lock-step), accumulating
//! score from the substitution matrix. Track `best_right` = max score
//! reached so far. Stop when `current < best_right - X` (the score has
//! dropped too far below its peak — extending further is unlikely to
//! recover). Symmetrically extend left.
//!
//! Final HSP score = `best_right + best_left` (the seed position itself
//! is counted in one of them; we use the convention that the seed
//! residue is included in `best_right`).
//!
//! ## X-drop value
//!
//! BLAST blastp uses X = 7 bits ≈ 16 in raw matrix units for BLOSUM62
//! during the ungapped stage (`-xdrop_ungap`). We default to X = 20
//! to be slightly more permissive — false negatives at this stage are
//! permanent.
//!
//! ## HSP threshold
//!
//! BLAST derives this from K, λ, query length and database size
//! (`S_min` such that E ≤ 10 for the most permissive setting). For
//! BLOSUM62 with typical queries that lands near raw score 30-45. We
//! default to 25, configurable via `--hsp-threshold`. Bench data on
//! random-sequence DBs shows this preserves the full hit set found by
//! `--no-index`.

use crate::scoring::matrices::ScoreMatrix;

/// Default X-drop value (raw score units, BLOSUM62-tuned).
pub const DEFAULT_X_DROP: i32 = 20;

/// Default minimum ungapped HSP score for a candidate to survive to
/// gapped SW. Conservative default — tighten with `--hsp-threshold` for
/// faster, slightly less sensitive searches.
pub const DEFAULT_HSP_THRESHOLD: i32 = 25;

/// Run ungapped X-drop extension from a seed on diagonal `d = s_anchor - q_anchor`.
///
/// `q` and `s` must be encoded into 0..23 (the substitution-matrix index
/// space; X = 22 for unknown residues). The caller is responsible for that.
///
/// Returns the best HSP score along the diagonal passing through the seed.
/// Returns 0 if the seed position is out of bounds.
pub fn extend_ungapped(
    q: &[u8],
    s: &[u8],
    q_anchor: usize,
    s_anchor: usize,
    matrix: &ScoreMatrix,
    x_drop: i32,
) -> i32 {
    extend_ungapped_with_range(q, s, q_anchor, s_anchor, matrix, x_drop).0
}

/// Like [`extend_ungapped`] but also returns the HSP range as `(score,
/// q_start, q_end)`. The range is in query coordinates (subject coordinates
/// are `q_start + diagonal` and `q_end + diagonal`). Range is empty when
/// score = 0 (no HSP).
pub fn extend_ungapped_with_range(
    q: &[u8],
    s: &[u8],
    q_anchor: usize,
    s_anchor: usize,
    matrix: &ScoreMatrix,
    x_drop: i32,
) -> (i32, usize, usize) {
    if q_anchor >= q.len() || s_anchor >= s.len() {
        return (0, 0, 0);
    }

    let qi = q[q_anchor] as usize;
    let si = s[s_anchor] as usize;
    if qi >= 24 || si >= 24 {
        return (0, 0, 0);
    }
    let seed_score = matrix.matrix[qi][si];

    // ---- Extend right ----
    let mut best_right: i32 = seed_score;
    let mut best_right_qp: usize = q_anchor; // q-coord of best-so-far on right
    let mut cur: i32 = seed_score;
    let mut qp = q_anchor + 1;
    let mut sp = s_anchor + 1;
    while qp < q.len() && sp < s.len() {
        let qi = q[qp] as usize;
        let si = s[sp] as usize;
        if qi >= 24 || si >= 24 {
            break;
        }
        cur += matrix.matrix[qi][si];
        if cur > best_right {
            best_right = cur;
            best_right_qp = qp;
        } else if cur < best_right - x_drop {
            break;
        }
        qp += 1;
        sp += 1;
    }

    // ---- Extend left ----
    let mut best_left: i32 = 0;
    let mut best_left_qp: usize = q_anchor; // q-coord of best-so-far on left
    let mut cur: i32 = 0;
    let mut qp = q_anchor as isize - 1;
    let mut sp = s_anchor as isize - 1;
    while qp >= 0 && sp >= 0 {
        let qi = q[qp as usize] as usize;
        let si = s[sp as usize] as usize;
        if qi >= 24 || si >= 24 {
            break;
        }
        cur += matrix.matrix[qi][si];
        if cur > best_left {
            best_left = cur;
            best_left_qp = qp as usize;
        } else if cur < best_left - x_drop {
            break;
        }
        qp -= 1;
        sp -= 1;
    }

    (best_right + best_left, best_left_qp, best_right_qp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::ScoringSystem;
    use crate::gpu::encoding::encode_protein;

    #[test]
    fn identical_sequences_extend_far() {
        let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
        let q = encode_protein(b"MVLSPADKTNVKAAWGKVGAH");
        let s = encode_protein(b"MVLSPADKTNVKAAWGKVGAH");
        // Anchor near the middle; should extend the full length both ways.
        let score = extend_ungapped(&q, &s, 10, 10, &matrix, DEFAULT_X_DROP);
        // BLOSUM62 self-scores sum to a large positive number for this peptide.
        assert!(score > 50, "expected high score for identical match, got {}", score);
    }

    #[test]
    fn unrelated_sequences_dont_extend() {
        let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
        // Two totally random sequences with a single matching residue.
        let q = encode_protein(b"AAAAAAMAAAAAA");
        let s = encode_protein(b"DDDDDDMDDDDDD");
        // The M at position 6 in each matches with M-M = 5 in BLOSUM62.
        // Extension hits A-D and M-D scores which are bad and drop fast.
        let score = extend_ungapped(&q, &s, 6, 6, &matrix, DEFAULT_X_DROP);
        // Should be roughly the M-M self-score (5) plus a small amount of noise.
        assert!(score < 20, "expected low score, got {}", score);
        assert!(score >= 5, "expected at least the seed match score, got {}", score);
    }

    #[test]
    fn partial_homology_scores_in_between() {
        let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
        // Conserved core with random tails.
        let q = encode_protein(b"AAAAAMVLSPADKAAAAA");
        let s = encode_protein(b"DDDDDMVLSPADKDDDDD");
        // Anchor in the middle of the conserved region.
        let score = extend_ungapped(&q, &s, 9, 9, &matrix, DEFAULT_X_DROP);
        // Should pick up the MVLSPADK core (positive scores).
        assert!(score > 25, "expected medium-high score, got {}", score);
    }

    #[test]
    fn handles_out_of_bounds_anchor() {
        let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
        let q = encode_protein(b"MVLSP");
        let s = encode_protein(b"MVLSP");
        let score = extend_ungapped(&q, &s, 100, 100, &matrix, DEFAULT_X_DROP);
        assert_eq!(score, 0);
    }

    #[test]
    fn handles_ambiguous_residues() {
        // Unknown bytes encoded as X (idx 22 in the BLOSUM matrix).
        // Should still produce some score without crashing.
        let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
        let q = encode_protein(b"MVL*SP");
        let s = encode_protein(b"MVL*SP");
        let _score = extend_ungapped(&q, &s, 1, 1, &matrix, DEFAULT_X_DROP);
        // Just shouldn't panic.
    }

    #[test]
    fn x_drop_stops_extension_at_bad_region() {
        let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
        // Good seed, then trash on both sides.
        let q = encode_protein(b"DDDDDDDDDMVLSPADKDDDDDDDDD");
        let s = encode_protein(b"NNNNNNNNNMVLSPADKNNNNNNNNN");
        let anchor = 12; // middle of MVLSPADK
        let score_loose = extend_ungapped(&q, &s, anchor, anchor, &matrix, 100); // huge X-drop
        let score_tight = extend_ungapped(&q, &s, anchor, anchor, &matrix, 5);    // small X-drop
        // Tight X-drop should stop earlier → lower or equal score.
        assert!(score_tight <= score_loose);
        // Both should at least capture the conserved core.
        assert!(score_tight > 15, "got {}", score_tight);
    }
}
