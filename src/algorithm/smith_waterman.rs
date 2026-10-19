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
    /// CIGAR string of the alignment, SAM convention with the subject as the
    /// reference: `M` aligned pair, `I` query residue against a gap, `D`
    /// subject residue against a gap. (Before v1.2 `I`/`D` were swapped.)
    pub cigar: String,
    pub query_seq: String,
    pub subject_seq: String,
    pub match_string: String,
    /// Number of gap openings (BLAST tabular column 6).
    #[serde(default)]
    pub gap_opens: usize,
    /// Aligned pairs with a positive substitution score.
    #[serde(default)]
    pub positives: usize,
    #[serde(default)]
    pub query_length: usize,
    #[serde(default)]
    pub subject_length: usize,
    /// Translated searches: frame of the query / subject (±1..3), 0 otherwise.
    #[serde(default)]
    pub query_frame: i8,
    #[serde(default)]
    pub subject_frame: i8,
    /// False when the E-value was computed with UNCALIBRATED parameters.
    #[serde(default = "default_true")]
    pub evalue_calibrated: bool,
    /// False when the traceback band hit its memory cap before reproducing
    /// the score of the searched window (explicit fallback).
    #[serde(default = "default_true")]
    pub exact: bool,
}

fn default_true() -> bool {
    true
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
            gap_opens: 0,
            positives: 0,
            query_length: 0,
            subject_length: 0,
            query_frame: 0,
            subject_frame: 0,
            evalue_calibrated: true,
            exact: true,
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

/// Reusable buffers for [`SmithWaterman::align_with`].
///
/// Holds the Gotoh score rows, the packed direction matrix, and the
/// residue→matrix-index translations. All grow to a high-water mark and are
/// then reused, which is what removes the per-alignment allocation and
/// zero-page traffic from the hot loop.
#[derive(Default)]
pub struct AlignScratch {
    prev_h: Vec<i32>,
    prev_f: Vec<i32>,
    cur_h: Vec<i32>,
    cur_e: Vec<i32>,
    cur_f: Vec<i32>,
    dir: Vec<u8>,
    qidx: Vec<u8>,
    sidx: Vec<u8>,
}

impl AlignScratch {
    pub fn new() -> Self {
        Self::default()
    }

    /// Release the direction matrix if an unusually large window grew it.
    ///
    /// Reuse is what makes the arena fast, but retention has a cost: the
    /// traceback window is capped at 2048x2048, so one outlier alignment can
    /// leave 4 MB resident *per thread*. On a 64-core run that is 256 MB held
    /// for buffers that are typically twenty times smaller. Anything above the
    /// retention threshold is handed back; ordinary alignments stay well under
    /// it and never allocate again after the first few calls.
    #[inline]
    fn trim(&mut self) {
        const RETAIN_CELLS: usize = 1 << 20; // 1 MB
        if self.dir.len() > RETAIN_CELLS {
            self.dir = Vec::new();
        }
    }
}

thread_local! {
    static SW_SCRATCH: std::cell::RefCell<AlignScratch> =
        std::cell::RefCell::new(AlignScratch::default());
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
    ///
    /// Uses a thread-local scratch arena so repeated calls (the normal case —
    /// one per surviving candidate) do not re-allocate the DP rows or the
    /// direction matrix. Use [`Self::align_with`] to supply your own arena.
    pub fn align(&self, query: &[u8], subject: &[u8]) -> AlignmentResult {
        SW_SCRATCH.with(|sc| self.align_with(query, subject, &mut sc.borrow_mut()))
    }

    /// Run Smith-Waterman reusing a caller-owned scratch arena.
    ///
    /// The forward pass is identical to the original in arithmetic and in the
    /// direction bits it writes; the differences are all about *where the
    /// bytes live* and *how the score is fetched*:
    ///
    /// * **No per-call allocation.** The original allocated five score rows
    ///   plus a full `(m+1)·(n+1)` direction matrix on every call. That
    ///   matrix is ~100 KB for a typical 325×313 alignment, and `vec![0u8; N]`
    ///   routes to `alloc_zeroed`, so each call took fresh zero pages from the
    ///   kernel and faulted them in one by one. Here the buffers grow to a
    ///   high-water mark and are then reused untouched.
    ///
    /// * **No zeroing of the direction matrix.** Only cells with `i>=1, j>=1`
    ///   are ever written, and traceback stops at `i==0 || j==0` or on a stop
    ///   bit, so row 0 and column 0 are never read. Growing the buffer is
    ///   therefore the only time we pay for initialisation.
    ///
    /// * **No per-cell branch on the scoring backend.** The original called
    ///   `pair_score`, which matched on the `Scoring` enum and then did two
    ///   lookups into a `[usize; 256]` table (2 KB, competing with the DP rows
    ///   for L1) before indexing the matrix. Residues are now translated to
    ///   matrix indices once per call, and the substitution row for the
    ///   current query residue is hoisted out of the inner loop, so a cell
    ///   costs one `i32` load from a 24-entry slice.
    pub fn align_with(
        &self,
        query: &[u8],
        subject: &[u8],
        sc: &mut AlignScratch,
    ) -> AlignmentResult {
        let m = query.len();
        let n = subject.len();
        if m == 0 || n == 0 {
            return AlignmentResult::default();
        }

        let gop = self.gap_open;
        let gex = self.gap_extend;
        let neg = i32::MIN / 2;
        let row = n + 1;
        let need = row.checked_mul(m + 1).expect("alignment window too large");

        // Grow-only buffers: after the first few alignments these are already
        // large enough and the resize is a no-op.
        if sc.dir.len() < need {
            sc.dir.resize(need, 0);
        }
        for b in [&mut sc.prev_h, &mut sc.cur_h, &mut sc.cur_e, &mut sc.cur_f, &mut sc.prev_f] {
            if b.len() < row {
                b.resize(row, 0);
            }
        }
        sc.prev_h[..row].fill(0);
        sc.prev_f[..row].fill(neg);

        // Translate residues to matrix indices once, instead of twice per cell.
        sc.qidx.clear();
        sc.sidx.clear();
        let use_matrix = matches!(self.scoring, Scoring::Matrix(_));
        if use_matrix {
            if let Scoring::Matrix(mx) = &self.scoring {
                sc.qidx.extend(query.iter().map(|&c| mx.char_to_idx[c as usize] as u8));
                sc.sidx.extend(subject.iter().map(|&c| mx.char_to_idx[c as usize] as u8));
            }
        }

        let mut max_score = 0i32;
        let mut max_i = 0usize;
        let mut max_j = 0usize;

        for i in 1..=m {
            sc.cur_h[0] = 0;
            sc.cur_e[0] = neg;
            sc.cur_f[0] = neg;
            let dbase = i * row;

            // Hoist the substitution row for this query residue. For the
            // matrix backend that is a 24-entry slice; for nucleotides we
            // fall back to the generic per-cell scorer.
            let mrow: Option<&[i32; 24]> = if use_matrix {
                match &self.scoring {
                    Scoring::Matrix(mx) => Some(&mx.matrix[sc.qidx[i - 1] as usize]),
                    _ => None,
                }
            } else {
                None
            };
            let qi_raw = query[i - 1];

            for j in 1..=n {
                let s = match mrow {
                    Some(r) => r[sc.sidx[j - 1] as usize],
                    None => self.pair_score(qi_raw, subject[j - 1]),
                };

                // E: gap consuming a subject column (query gets '-').
                let e_open = sc.cur_h[j - 1] - gop;
                let e_ext = sc.cur_e[j - 1] - gex;
                let (e_ij, e_from_ext) =
                    if e_ext > e_open { (e_ext, true) } else { (e_open, false) };

                // F: gap consuming a query row (subject gets '-').
                let f_open = sc.prev_h[j] - gop;
                let f_ext = sc.prev_f[j] - gex;
                let (f_ij, f_from_ext) =
                    if f_ext > f_open { (f_ext, true) } else { (f_open, false) };

                let diag = sc.prev_h[j - 1] + s;

                let mut v = 0i32;
                let mut horigin: u8 = 0;
                if diag > v { v = diag; horigin = 1; }
                if e_ij > v { v = e_ij; horigin = 2; }
                if f_ij > v { v = f_ij; horigin = 3; }

                sc.cur_h[j] = v;
                sc.cur_e[j] = e_ij;
                sc.cur_f[j] = f_ij;

                let mut d = horigin;
                if e_from_ext { d |= 0b0100; }
                if f_from_ext { d |= 0b1000; }
                sc.dir[dbase + j] = d;

                if v > max_score {
                    max_score = v;
                    max_i = i;
                    max_j = j;
                }
            }
            std::mem::swap(&mut sc.prev_h, &mut sc.cur_h);
            std::mem::swap(&mut sc.prev_f, &mut sc.cur_f);
        }

        if max_score == 0 {
            sc.trim();
            return AlignmentResult::default();
        }

        let res = self.traceback(&sc.dir, query, subject, max_i, max_j, max_score, n);
        sc.trim();
        res
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
                    // query gap: subject residue present -> SAM 'D'
                    cigar_rev.push(b'D');
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
                    // subject gap: query residue present -> SAM 'I'
                    cigar_rev.push(b'I');
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

        let gap_opens = cigar_rev.windows(2).filter(|w| w[1] != b'M' && w[1] != w[0]).count()
            + usize::from(cigar_rev.first().map_or(false, |&c| c != b'M'));
        let positives = m_aln.iter().filter(|&&c| c == b'|' || c == b'+').count();
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
            gap_opens,
            positives,
            query_length: query.len(),
            subject_length: subject.len(),
            ..AlignmentResult::default()
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
        // Subject has the extra residue: SAM 'D' (reference residue vs query gap).
        assert!(r.cigar.contains("1D"), "one deletion op, got {}", r.cigar);
        assert!(!r.cigar.contains('I'), "no insertions, got {}", r.cigar);
        assert_eq!(r.score, 14, "8*2 - gap_open(2)");
        assert_eq!(r.gap_opens, 1);
    }

    #[test]
    fn single_deletion_in_query_is_one_gap() {
        // Query has the extra residue → one length-1 gap in the subject (D).
        let sw = SmithWaterman::nucleotide(2, 1, 2, -3);
        let r = sw.align(b"ACGTTACGT", b"ACGTACGT");
        assert_eq!(r.matches, 8);
        assert_eq!(r.gaps, 1);
        // Query has the extra residue: SAM 'I' (query residue vs reference gap).
        assert!(r.cigar.contains("1I"), "one insertion op, got {}", r.cigar);
        assert!(!r.cigar.contains('D'), "no deletions, got {}", r.cigar);
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
        assert_eq!(r.cigar, "4M3D4M", "single contiguous 3-wide gap (subject residues vs query gap)");
        assert_eq!(r.score, 10, "16 - (gap_open 4 + 2*gap_extend 1)");
    }
}

