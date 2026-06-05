//! Smith-Waterman local sequence alignment with affine gap penalties.
//!
//! Three scoring backends are supported, selected when constructing the aligner:
//!
//! * **Match/mismatch** — for nucleotide-nucleotide alignment
//! * **Scoring matrix** — for protein-protein alignment (BLOSUM/PAM)
//! * Both can be combined with affine gap penalties via `gap_open` + `gap_extend`
//!
//! The DP is the Gotoh affine-gap formulation:
//!
//! ```text
//! H(i,j) = max( 0,
//!               H(i-1,j-1) + s(qᵢ, sⱼ),
//!               E(i,j),
//!               F(i,j) )
//! E(i,j) = max( H(i,j-1) - gap_open,   E(i,j-1) - gap_extend )
//! F(i,j) = max( H(i-1,j) - gap_open,   F(i-1,j) - gap_extend )
//! ```
//!
//! Traceback consults all three matrices to recover the optimal local path —
//! the original implementation used H only and produced incorrect gap structures.

use serde::{Deserialize, Serialize};

use crate::scoring::ScoreMatrix;

/// One local alignment between a query and a subject.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentResult {
    pub query_id: String,
    pub subject_id: String,
    /// Raw Smith-Waterman score.
    pub score: i32,
    /// E-value (Karlin-Altschul). Filled in by the caller after scoring.
    pub evalue: f64,
    /// Bit score. Filled in by the caller after scoring.
    pub bitscore: f64,
    /// Percent identity, 0-100.
    pub identity: f64,
    pub query_coverage: f64,
    pub subject_coverage: f64,
    /// 1-based inclusive coordinates.
    pub query_start: usize,
    pub query_end: usize,
    pub subject_start: usize,
    pub subject_end: usize,
    pub align_length: usize,
    pub matches: usize,
    pub mismatches: usize,
    /// Number of gap columns (insertions + deletions).
    pub gaps: usize,
    /// CIGAR string of the alignment (M/I/D).
    pub cigar: String,
    pub query_seq: String,
    pub subject_seq: String,
    pub match_string: String,
}

impl Default for AlignmentResult {
    fn default() -> Self {
        Self {
            query_id: String::new(),
            subject_id: String::new(),
            score: 0,
            evalue: f64::INFINITY,
            bitscore: 0.0,
            identity: 0.0,
            query_coverage: 0.0,
            subject_coverage: 0.0,
            query_start: 0,
            query_end: 0,
            subject_start: 0,
            subject_end: 0,
            align_length: 0,
            matches: 0,
            mismatches: 0,
            gaps: 0,
            cigar: String::new(),
            query_seq: String::new(),
            subject_seq: String::new(),
            match_string: String::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentSet {
    pub query_id: String,
    pub query_length: usize,
    pub alignments: Vec<AlignmentResult>,
}

/// What kind of pairwise score the aligner uses.
#[derive(Clone)]
pub enum Scoring {
    /// Fixed match/mismatch — used for nucleotide alignment.
    NucleotideSimple { match_score: i32, mismatch: i32 },
    /// Substitution matrix — used for protein alignment.
    Matrix(std::sync::Arc<ScoreMatrix>),
}

/// Smith-Waterman aligner with affine gaps.
///
/// `gap_open` is the penalty paid the first time a gap is opened. `gap_extend`
/// is paid for every additional gap character. Both are positive numbers (the
/// DP subtracts them).
pub struct SmithWaterman {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub scoring: Scoring,
}

impl SmithWaterman {
    /// New nucleotide aligner.
    pub fn nucleotide(gap_open: i32, gap_extend: i32, match_score: i32, mismatch: i32) -> Self {
        Self {
            gap_open,
            gap_extend,
            scoring: Scoring::NucleotideSimple { match_score, mismatch },
        }
    }

    /// New protein aligner with a substitution matrix.
    pub fn protein(gap_open: i32, gap_extend: i32, matrix: std::sync::Arc<ScoreMatrix>) -> Self {
        Self { gap_open, gap_extend, scoring: Scoring::Matrix(matrix) }
    }

    /// Score one residue pair using whatever backend is configured.
    #[inline]
    fn pair_score(&self, a: u8, b: u8) -> i32 {
        match &self.scoring {
            Scoring::NucleotideSimple { match_score, mismatch } => {
                // Treat as nucleotide; A/C/G/T compared case-insensitively, anything else mismatches.
                let au = a.to_ascii_uppercase();
                let bu = b.to_ascii_uppercase();
                if au == bu && matches!(au, b'A' | b'C' | b'G' | b'T' | b'U') {
                    *match_score
                } else {
                    *mismatch
                }
            }
            Scoring::Matrix(m) => m.score(a, b),
        }
    }

    /// Run Smith-Waterman on the two byte slices. Returns the best local alignment.
    pub fn align(&self, query: &[u8], subject: &[u8]) -> AlignmentResult {
        let m = query.len();
        let n = subject.len();
        if m == 0 || n == 0 {
            return AlignmentResult::default();
        }

        let gop = self.gap_open;
        let gex = self.gap_extend;
        let neg = i32::MIN / 2;

        // LINEAR-SPACE forward pass.
        //
        // The classic Gotoh SW keeps three full (m+1)×(n+1) i32 matrices
        // (H, E, F) = 12 bytes/cell, which dominates peak memory for long
        // sequences (a 2949-residue self-alignment needs ~100 MB). Here we
        // keep only two score *rows* (H/F of the previous row, H/E/F of the
        // current row) and a single 1-byte direction matrix used for
        // traceback. That's ~1 byte/cell instead of 12 — a 12× reduction —
        // and the score recurrence is unchanged, so results are identical.
        //
        // Direction byte layout (per cell):
        //   bits 0-1  H origin: 0=stop(v==0) 1=diagonal 2=from-E 3=from-F
        //   bit  2    E came from E-extend (1) vs H-open (0)
        //   bit  3    F came from F-extend (1) vs H-open (0)
        let mut prev_h = vec![0i32; n + 1];
        let mut prev_f = vec![neg; n + 1];
        let mut cur_h = vec![0i32; n + 1];
        let mut cur_e = vec![neg; n + 1];
        let mut cur_f = vec![neg; n + 1];

        let mut dir = vec![0u8; (m + 1) * (n + 1)];
        let row = n + 1;

        let mut max_score = 0i32;
        let mut max_i = 0usize;
        let mut max_j = 0usize;

        for i in 1..=m {
            let qi = query[i - 1];
            cur_h[0] = 0;
            cur_e[0] = neg;
            cur_f[0] = neg;
            let dbase = i * row;
            for j in 1..=n {
                let sj = subject[j - 1];
                let s = self.pair_score(qi, sj);

                // E: gap consuming a subject column (query gets '-').
                let e_open = cur_h[j - 1] - gop;
                let e_ext = cur_e[j - 1] - gex;
                let (e_ij, e_from_ext) = if e_ext > e_open { (e_ext, true) } else { (e_open, false) };

                // F: gap consuming a query row (subject gets '-').
                let f_open = prev_h[j] - gop;
                let f_ext = prev_f[j] - gex;
                let (f_ij, f_from_ext) = if f_ext > f_open { (f_ext, true) } else { (f_open, false) };

                let diag = prev_h[j - 1] + s;

                let mut v = 0i32;
                let mut horigin: u8 = 0;
                if diag > v { v = diag; horigin = 1; }
                if e_ij > v { v = e_ij; horigin = 2; }
                if f_ij > v { v = f_ij; horigin = 3; }

                cur_h[j] = v;
                cur_e[j] = e_ij;
                cur_f[j] = f_ij;

                let mut d = horigin;
                if e_from_ext { d |= 0b0100; }
                if f_from_ext { d |= 0b1000; }
                dir[dbase + j] = d;

                if v > max_score {
                    max_score = v;
                    max_i = i;
                    max_j = j;
                }
            }
            std::mem::swap(&mut prev_h, &mut cur_h);
            std::mem::swap(&mut prev_f, &mut cur_f);
        }

        if max_score == 0 {
            return AlignmentResult::default();
        }

        self.traceback(&dir, query, subject, max_i, max_j, max_score, n)
    }

    /// Walk back from the maximum-scoring cell using the direction matrix
    /// written during the forward pass. `state` tracks which Gotoh plane we
    /// are in: H (diagonal/start), E (gap consuming subject columns), or F
    /// (gap consuming query rows). This reproduces the exact alignment the
    /// full-matrix traceback would, but reads 1 byte/cell instead of three
    /// i32 matrices.
    fn traceback(
        &self,
        dir: &[u8],
        query: &[u8],
        subject: &[u8],
        mut i: usize,
        mut j: usize,
        score: i32,
        n: usize,
    ) -> AlignmentResult {
        let row = n + 1;

        let mut q_aln: Vec<u8> = Vec::new();
        let mut s_aln: Vec<u8> = Vec::new();
        let mut m_aln: Vec<u8> = Vec::new();
        let mut cigar_rev: Vec<u8> = Vec::new();

        let mut matches = 0usize;
        let mut mismatches = 0usize;
        let mut gaps = 0usize;

        let query_end = i;
        let subject_end = j;

        // 0 = H plane, 1 = E plane (consume column j), 2 = F plane (consume row i).
        let mut state: u8 = 0;

        loop {
            let d = dir[i * row + j];
            match state {
                0 => {
                    let horigin = d & 0b11;
                    if horigin == 0 {
                        break; // SW stop cell (v == 0)
                    } else if horigin == 1 {
                        // Diagonal: match or mismatch.
                        let qi = query[i - 1];
                        let sj = subject[j - 1];
                        let sc = self.pair_score(qi, sj);
                        q_aln.push(qi);
                        s_aln.push(sj);
                        if qi.to_ascii_uppercase() == sj.to_ascii_uppercase() && qi != b'-' {
                            m_aln.push(b'|');
                            matches += 1;
                        } else if sc > 0 {
                            m_aln.push(b'+');
                            mismatches += 1;
                        } else {
                            m_aln.push(b' ');
                            mismatches += 1;
                        }
                        cigar_rev.push(b'M');
                        i -= 1;
                        j -= 1;
                    } else if horigin == 2 {
                        state = 1; // move into E plane at (i, j)
                    } else {
                        state = 2; // move into F plane at (i, j)
                    }
                }
                1 => {
                    // E plane: gap consuming subject column j; query gets '-'.
                    let sj = subject[j - 1];
                    q_aln.push(b'-');
                    s_aln.push(sj);
                    m_aln.push(b' ');
                    cigar_rev.push(b'I');
                    gaps += 1;
                    let extend = (d & 0b0100) != 0;
                    j -= 1;
                    if !extend {
                        state = 0; // gap opened from H at the new (i, j)
                    }
                }
                _ => {
                    // F plane: gap consuming query row i; subject gets '-'.
                    let qi = query[i - 1];
                    q_aln.push(qi);
                    s_aln.push(b'-');
                    m_aln.push(b' ');
                    cigar_rev.push(b'D');
                    gaps += 1;
                    let extend = (d & 0b1000) != 0;
                    i -= 1;
                    if !extend {
                        state = 0;
                    }
                }
            }
            if i == 0 || j == 0 {
                break;
            }
        }

        q_aln.reverse();
        s_aln.reverse();
        m_aln.reverse();
        cigar_rev.reverse();

        let query_start = i + 1;
        let subject_start = j + 1;
        let align_length = q_aln.len();
        let identity = if align_length > 0 {
            (matches as f64 / align_length as f64) * 100.0
        } else {
            0.0
        };

        AlignmentResult {
            query_id: String::new(),
            subject_id: String::new(),
            score,
            evalue: f64::INFINITY,
            bitscore: 0.0,
            identity,
            query_coverage: 0.0,
            subject_coverage: 0.0,
            query_start,
            query_end,
            subject_start,
            subject_end,
            align_length,
            matches,
            mismatches,
            gaps,
            cigar: compact_cigar(&cigar_rev),
            query_seq: String::from_utf8_lossy(&q_aln).into_owned(),
            subject_seq: String::from_utf8_lossy(&s_aln).into_owned(),
            match_string: String::from_utf8_lossy(&m_aln).into_owned(),
        }
    }
}

/// Compact a per-column CIGAR like "MMMIDD" into "3M1I2D".
fn compact_cigar(ops: &[u8]) -> String {
    let mut out = String::new();
    if ops.is_empty() {
        return out;
    }
    let mut cur = ops[0];
    let mut run = 1usize;
    for &c in &ops[1..] {
        if c == cur {
            run += 1;
        } else {
            out.push_str(&format!("{}{}", run, cur as char));
            cur = c;
            run = 1;
        }
    }
    out.push_str(&format!("{}{}", run, cur as char));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{ScoreMatrix, ScoringSystem};
    use std::sync::Arc;

    #[test]
    fn identical_nucleotide_sequences_score_nonzero() {
        let sw = SmithWaterman::nucleotide(10, 1, 2, -3);
        let r = sw.align(b"ACGTACGT", b"ACGTACGT");
        assert_eq!(r.matches, 8);
        assert_eq!(r.mismatches, 0);
        assert_eq!(r.gaps, 0);
        assert_eq!(r.score, 16); // 8 * 2
        assert!(r.identity > 99.9);
    }

    #[test]
    fn distant_nucleotide_sequences_score_lower() {
        let sw = SmithWaterman::nucleotide(10, 1, 2, -3);
        let same = sw.align(b"ACGTACGTACGT", b"ACGTACGTACGT");
        let diff = sw.align(b"ACGTACGTACGT", b"TGCATGCATGCA");
        assert!(same.score > diff.score, "identical should beat unrelated");
    }

    #[test]
    fn protein_identical_uses_blosum_diagonal() {
        let m = Arc::new(ScoreMatrix::new(ScoringSystem::Blosum62));
        let sw = SmithWaterman::protein(11, 1, m);
        let r = sw.align(b"WAGNER", b"WAGNER");
        // BLOSUM62 diagonal: W=11, A=4, G=6, N=6, E=5, R=5 → 37
        assert_eq!(r.score, 37);
        assert!(r.identity > 99.9);
    }

    #[test]
    fn protein_distant_scores_lower_than_identical() {
        let m = Arc::new(ScoreMatrix::new(ScoringSystem::Blosum62));
        let sw = SmithWaterman::protein(11, 1, m);
        let identical = sw.align(b"WAGNERSEQ", b"WAGNERSEQ");
        let distant = sw.align(b"WAGNERSEQ", b"PPPPPPPPP");
        assert!(identical.score > distant.score);
    }

    #[test]
    fn empty_sequence_returns_zero() {
        let sw = SmithWaterman::nucleotide(10, 1, 2, -3);
        let r = sw.align(b"", b"ACGT");
        assert_eq!(r.score, 0);
    }

    #[test]
    fn affine_gap_traceback_yields_compact_cigar() {
        let sw = SmithWaterman::nucleotide(2, 1, 2, -3);
        let r = sw.align(b"ACGTACGT", b"ACGTACGT");
        assert_eq!(r.cigar, "8M");
    }

    #[test]
    fn single_insertion_in_subject_is_one_gap() {
        // Subject has an extra residue. Optimal: 8 matches, one length-1 gap
        // in the query. The exact column of the gap is co-optimal (it can sit
        // on either side of the duplicated residue), so we assert the
        // invariants — score, matches, gap count, and a single "1I" in the
        // CIGAR — rather than a fixed gap position.
        let sw = SmithWaterman::nucleotide(2, 1, 2, -3);
        let r = sw.align(b"ACGTACGT", b"ACGTTACGT");
        assert_eq!(r.matches, 8, "all 8 query residues align as matches");
        assert_eq!(r.gaps, 1, "exactly one gap column");
        assert!(r.cigar.contains("1I"), "one insertion op, got {}", r.cigar);
        assert!(!r.cigar.contains('D'), "no deletions, got {}", r.cigar);
        assert_eq!(r.score, 14, "8*2 - gap_open(2)");
    }

    #[test]
    fn single_deletion_in_query_is_one_gap() {
        // Query has the extra residue → one length-1 gap in the subject (D).
        let sw = SmithWaterman::nucleotide(2, 1, 2, -3);
        let r = sw.align(b"ACGTTACGT", b"ACGTACGT");
        assert_eq!(r.matches, 8);
        assert_eq!(r.gaps, 1);
        assert!(r.cigar.contains("1D"), "one deletion op, got {}", r.cigar);
        assert!(!r.cigar.contains('I'), "no insertions, got {}", r.cigar);
        assert_eq!(r.score, 14);
    }

    #[test]
    fn longer_gap_uses_affine_extend() {
        // A 3-residue insertion must cost gap_open + 2*gap_extend (affine),
        // NOT three separate opens. With gap_open=4, gap_extend=1 that's
        // 4 + 2 = 6, so score = 8*2 - 6 = 10. If the traceback wrongly
        // decomposed the gap into opens it would still score the same
        // (scoring is from the forward DP), but the gap MUST be a single
        // contiguous 3-wide run, which we check via "3I".
        let sw = SmithWaterman::nucleotide(4, 1, 2, -3);
        let r = sw.align(b"AAAACCCC", b"AAAAGGGCCCC");
        assert_eq!(r.matches, 8);
        assert_eq!(r.gaps, 3, "three gap columns");
        assert_eq!(r.cigar, "4M3I4M", "single contiguous 3-wide gap");
        assert_eq!(r.score, 10, "16 - (gap_open 4 + 2*gap_extend 1)");
    }
}

