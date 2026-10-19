//! Independent, deliberately naive reference implementations used to test the
//! optimised kernels (SIMD, banded, traceback). Nothing here is fast; every
//! function keeps the full O(m·n) matrices in `i64` so it cannot overflow and
//! shares no code with the production kernels.

use crate::alphabet::code_index;

/// Full-matrix Gotoh Smith–Waterman score in i64.
///
/// Conventions (identical to SABER's): the first residue of a gap costs
/// `gap_open`, every further residue `gap_extend`; `bias[i]` (if given) is
/// added to every aligned column involving query position `i`.
/// Returns `(score, query_end_exclusive, subject_end_exclusive)` of the first
/// (row-major) maximal cell.
pub fn sw_reference(
    q: &[u8],
    s: &[u8],
    matrix: &[[i32; 24]; 24],
    bias: Option<&[i8]>,
    gap_open: i64,
    gap_extend: i64,
) -> (i64, usize, usize) {
    let m = q.len();
    let n = s.len();
    if m == 0 || n == 0 {
        return (0, 0, 0);
    }
    let neg = i64::MIN / 4;
    let w = n + 1;
    let mut h = vec![0i64; (m + 1) * w];
    let mut e = vec![neg; (m + 1) * w];
    let mut f = vec![neg; (m + 1) * w];
    let mut best = (0i64, 0usize, 0usize);
    for i in 1..=m {
        let qi = code_index(q[i - 1]);
        let b = bias.map_or(0, |b| b[i - 1] as i64);
        for j in 1..=n {
            let sj = code_index(s[j - 1]);
            let ee = (h[i * w + j - 1] - gap_open).max(e[i * w + j - 1] - gap_extend);
            let ff = (h[(i - 1) * w + j] - gap_open).max(f[(i - 1) * w + j] - gap_extend);
            let dd = h[(i - 1) * w + j - 1] + matrix[qi][sj] as i64 + b;
            let v = 0.max(ee).max(ff).max(dd);
            e[i * w + j] = ee;
            f[i * w + j] = ff;
            h[i * w + j] = v;
            if v > best.0 {
                best = (v, i, j);
            }
        }
    }
    best
}

/// Re-score an alignment described by start coordinates and SAM-convention
/// operations (`M` aligned pair, `I` query residue against a gap, `D` subject
/// residue against a gap). Used to prove that a reported score, its
/// coordinates and its CIGAR are mutually consistent.
pub fn score_alignment(
    q: &[u8],
    s: &[u8],
    q_start: usize,
    s_start: usize,
    ops: &[(u8, u32)],
    matrix: &[[i32; 24]; 24],
    bias: Option<&[i8]>,
    gap_open: i64,
    gap_extend: i64,
) -> Option<(i64, usize, usize)> {
    let (mut i, mut j) = (q_start, s_start);
    let mut score = 0i64;
    for &(op, len) in ops {
        let len = len as usize;
        match op {
            b'M' => {
                for _ in 0..len {
                    if i >= q.len() || j >= s.len() {
                        return None;
                    }
                    let b = bias.map_or(0, |b| b[i] as i64);
                    score += matrix[code_index(q[i])][code_index(s[j])] as i64 + b;
                    i += 1;
                    j += 1;
                }
            }
            b'I' => {
                if len == 0 || i + len > q.len() {
                    return None;
                }
                score -= gap_open + (len as i64 - 1) * gap_extend;
                i += len;
            }
            b'D' => {
                if len == 0 || j + len > s.len() {
                    return None;
                }
                score -= gap_open + (len as i64 - 1) * gap_extend;
                j += len;
            }
            _ => return None,
        }
    }
    Some((score, i, j))
}

/// Exhaustive ungapped best segment on diagonal `d = j − i` (i64, no X-drop).
pub fn best_ungapped_on_diagonal(q: &[u8], s: &[u8], d: i64, matrix: &[[i32; 24]; 24], bias: Option<&[i8]>) -> i64 {
    let mut best = 0i64;
    let mut cur = 0i64;
    for i in 0..q.len() as i64 {
        let j = i + d;
        if j < 0 || j >= s.len() as i64 {
            continue;
        }
        let b = bias.map_or(0, |b| b[i as usize] as i64);
        cur = (cur + matrix[code_index(q[i as usize])][code_index(s[j as usize])] as i64 + b).max(0);
        best = best.max(cur);
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::EncodedProtein;
    use crate::scoring::{ScoreMatrix, ScoringSystem};

    #[test]
    fn reference_self_score() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let q = EncodedProtein::from_ascii(b"WAGNER");
        let (sc, i, j) = sw_reference(q.as_slice(), q.as_slice(), &m.matrix, None, 11, 1);
        assert_eq!(sc, 37);
        assert_eq!((i, j), (6, 6));
    }

    #[test]
    fn rescoring_a_gapped_alignment() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let q = EncodedProtein::from_ascii(b"WWWWAAAWWWW");
        let s = EncodedProtein::from_ascii(b"WWWWWWWW");
        // 4 W pairs, 3 query residues against gaps, 4 W pairs.
        let ops = [(b'M', 4), (b'I', 3), (b'M', 4)];
        let (sc, ie, je) = score_alignment(q.as_slice(), s.as_slice(), 0, 0, &ops, &m.matrix, None, 11, 1).unwrap();
        assert_eq!(sc, 8 * 11 - 13);
        assert_eq!((ie, je), (11, 8));
    }
}
