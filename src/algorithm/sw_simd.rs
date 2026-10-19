//! Inter-sequence SIMD Smith–Waterman score kernel (AVX2, 16 × i16 lanes).
//!
//! Rognes (2011) "inter-sequence" layout, as in SWIPE/Opal: sixteen database
//! subjects occupy the sixteen 16-bit lanes of one AVX2 register and are
//! aligned against the same query column by column.
//!
//! ## Correctness contract (v1.2)
//!
//! * Inputs are typed [`ProteinCodes`] — encoded residues 0..=23. ASCII can
//!   not reach this kernel.
//! * The per-query-position composition-bias correction (see
//!   [`crate::scoring::bias`]) is added to every substitution score, exactly
//!   as in the scalar, banded and traceback kernels.
//! * 16-bit saturation is detected and repaired: every H, E and F value of a
//!   local alignment is bounded above by the optimal score, so the i16
//!   computation is exact whenever the optimum is below `i16::MAX`, and a lane
//!   whose maximum equals `i16::MAX` is recomputed exactly in i32 by the
//!   scalar kernel. The number of such retries is reported in
//!   [`KernelCounters::saturation_retries`].
//! * Hosts without AVX2 use the scalar i32 kernel for every subject (same
//!   results, slower).
//!
//! `tests::simd_matches_independent_reference` compares this kernel with the
//! independent i64 reference in [`super::reference`] — not with another
//! saturating implementation.

#![allow(unsafe_code)]

use super::banded::{banded_score, Band, BandScratch, Scheme};
use crate::alphabet::ProteinCodes;
use crate::scoring::matrices::ScoreMatrix;
use crate::scoring::stats::GapCosts;

/// Lanes per AVX2 register (i16).
pub const SIMD_LANES: usize = 16;

/// Code used to pad finished lanes (X). Their updates are masked out.
const PAD_RESIDUE: u8 = 22;

/// Kernel work counters (summed into the search statistics).
#[derive(Debug, Default, Clone, Copy)]
pub struct KernelCounters {
    /// DP cells evaluated by the SIMD kernel (query length × subject length).
    pub simd_cells: u64,
    /// DP cells evaluated by scalar kernels (fallback / saturation retries).
    pub scalar_cells: u64,
    /// Lanes whose 16-bit score saturated and were recomputed in 32 bits.
    pub saturation_retries: u64,
    /// Subjects whose 8-bit first-pass score reached the 8-bit limit and were
    /// re-scored by the 16-bit kernel.
    pub u8_rescored: u64,
    /// SIMD batches executed.
    pub batches: u64,
}

impl KernelCounters {
    pub fn add(&mut self, o: &KernelCounters) {
        self.simd_cells += o.simd_cells;
        self.scalar_cells += o.scalar_cells;
        self.saturation_retries += o.saturation_retries;
        self.u8_rescored += o.u8_rescored;
        self.batches += o.batches;
    }
}

#[repr(C, align(32))]
#[derive(Clone, Copy, Default)]
struct Lane16([i16; 16]);

#[repr(C, align(32))]
#[derive(Clone, Copy, Default)]
struct Lane8([u8; 32]);

/// Re-usable buffers for [`score_many`].
#[derive(Default)]
pub struct SimdScratch {
    h: Vec<Lane16>,
    e: Vec<Lane16>,
    bias: Vec<Lane16>,
    h8: Vec<Lane8>,
    e8: Vec<Lane8>,
    bias8: Vec<Lane8>,
    order: Vec<usize>,
    band: BandScratch,
}

impl SimdScratch {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Exact SW scores of `query` against every subject (i32, one per subject,
/// caller order). `bias` (if any) must have the query's length.
pub fn score_many(
    query: ProteinCodes<'_>,
    bias: Option<&[i8]>,
    subjects: &[ProteinCodes<'_>],
    matrix: &[[i32; 24]; 24],
    gaps: GapCosts,
    sc: &mut SimdScratch,
    ctr: &mut KernelCounters,
) -> Vec<i32> {
    let mut out = vec![0i32; subjects.len()];
    if subjects.is_empty() || query.is_empty() {
        return out;
    }
    if let Some(b) = bias {
        assert_eq!(b.len(), query.len(), "bias length must equal query length");
    }
    let scheme = Scheme { matrix, gap_open: gaps.open, gap_extend: gaps.extend };

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            // Pass 1: 8-bit lanes (32 subjects per register); subjects whose
            // score may have saturated are re-scored by the 16-bit kernel.
            // SAFETY: AVX2 availability checked at runtime.
            if let Some(redo) = unsafe { score_many_u8_avx2(query.as_slice(), bias, subjects, matrix, gaps, sc, ctr, &mut out) } {
                if !redo.is_empty() {
                    ctr.u8_rescored += redo.len() as u64;
                    let sub: Vec<ProteinCodes<'_>> = redo.iter().map(|&i| subjects[i]).collect();
                    let mut o16 = vec![0i32; sub.len()];
                    // SAFETY: as above.
                    unsafe { score_many_avx2(query.as_slice(), bias, &sub, matrix, gaps, sc, ctr, &mut o16) };
                    for (k, &i) in redo.iter().enumerate() {
                        out[i] = o16[k];
                    }
                }
            } else {
                // SAFETY: AVX2 availability checked at runtime.
                unsafe { score_many_avx2(query.as_slice(), bias, subjects, matrix, gaps, sc, ctr, &mut out) };
            }
            // Repair saturated lanes exactly in i32.
            for (i, s) in subjects.iter().enumerate() {
                if out[i] >= i16::MAX as i32 {
                    let r = banded_score(query.as_slice(), s.as_slice(), bias, scheme,
                                         Band::full(query.len(), s.len()), &mut sc.band);
                    ctr.saturation_retries += 1;
                    ctr.scalar_cells += r.cells;
                    out[i] = r.score;
                }
            }
            return out;
        }
    }
    for (i, s) in subjects.iter().enumerate() {
        let r = banded_score(query.as_slice(), s.as_slice(), bias, scheme, Band::full(query.len(), s.len()), &mut sc.band);
        ctr.scalar_cells += r.cells;
        out[i] = r.score;
    }
    out
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn score_many_avx2(
    query: &[u8],
    bias: Option<&[i8]>,
    subjects: &[ProteinCodes<'_>],
    matrix: &[[i32; 24]; 24],
    gaps: GapCosts,
    sc: &mut SimdScratch,
    ctr: &mut KernelCounters,
    out: &mut [i32],
) {
    use std::arch::x86_64::*;
    let m = query.len();
    sc.order.clear();
    sc.order.extend(0..subjects.len());
    // Longest first; ties by index so the batch composition is deterministic.
    sc.order.sort_unstable_by(|&a, &b| subjects[b].len().cmp(&subjects[a].len()).then(a.cmp(&b)));

    sc.h.clear();
    sc.h.resize(m, Lane16::default());
    sc.e.clear();
    sc.e.resize(m, Lane16::default());
    sc.bias.clear();
    match bias {
        Some(b) => sc.bias.extend(b.iter().map(|&v| Lane16([v as i16; 16]))),
        None => sc.bias.resize(m, Lane16::default()),
    }

    // Per query letter: 16-byte shuffle tables for subject letters 0..15 and
    // 16..23 (scores fit in i8 for every supported matrix).
    let mut tlo = [_mm_setzero_si128(); 24];
    let mut thi = [_mm_setzero_si128(); 24];
    for a in 0..24 {
        let mut lo = [0i8; 16];
        let mut hi = [0i8; 16];
        for x in 0..16 {
            lo[x] = matrix[a][x].clamp(-128, 127) as i8;
        }
        for x in 16..24 {
            hi[x - 16] = matrix[a][x].clamp(-128, 127) as i8;
        }
        tlo[a] = _mm_loadu_si128(lo.as_ptr() as *const __m128i);
        thi[a] = _mm_loadu_si128(hi.as_ptr() as *const __m128i);
    }
    let go_v = _mm256_set1_epi16(gaps.open as i16);
    let ge_v = _mm256_set1_epi16(gaps.extend as i16);
    let zero = _mm256_setzero_si256();
    let neg = _mm256_set1_epi16(i16::MIN / 2);
    let fifteen = _mm_set1_epi8(15);
    let low_nib = _mm_set1_epi8(0x0f);

    let hp = sc.h.as_mut_ptr() as *mut __m256i;
    let ep = sc.e.as_mut_ptr() as *mut __m256i;
    let bp = sc.bias.as_ptr() as *const __m256i;

    let mut start = 0usize;
    while start < sc.order.len() {
        let end = (start + SIMD_LANES).min(sc.order.len());
        let nlanes = end - start;
        let mut seqs: [&[u8]; SIMD_LANES] = [&[]; SIMD_LANES];
        let mut max_len = 0usize;
        for l in 0..nlanes {
            seqs[l] = subjects[sc.order[start + l]].as_slice();
            max_len = max_len.max(seqs[l].len());
        }
        for r in 0..m {
            _mm256_store_si256(hp.add(r), zero);
            _mm256_store_si256(ep.add(r), neg);
        }
        let mut max_h = zero;
        let mut resid = [PAD_RESIDUE; 16];
        let mut act = [0i16; 16];
        let mut profile = [zero; 24];
        for col in 0..max_len {
            for l in 0..SIMD_LANES {
                let s = seqs[l];
                if col < s.len() {
                    resid[l] = crate::alphabet::code_index(s[col]) as u8;
                    act[l] = -1;
                } else {
                    resid[l] = PAD_RESIDUE;
                    act[l] = 0;
                }
            }
            let active = _mm256_loadu_si256(act.as_ptr() as *const __m256i);
            let rv = _mm_loadu_si128(resid.as_ptr() as *const __m128i);
            let idx = _mm_and_si128(rv, low_nib);
            let is_hi = _mm_cmpgt_epi8(rv, fifteen);
            for a in 0..24 {
                let lo = _mm_shuffle_epi8(tlo[a], idx);
                let hi = _mm_shuffle_epi8(thi[a], idx);
                profile[a] = _mm256_cvtepi8_epi16(_mm_blendv_epi8(lo, hi, is_hi));
            }
            let mut u_f = neg;
            let mut u_h = zero;
            let mut ul_h = zero;
            for r in 0..m {
                let qr = crate::alphabet::code_index(*query.get_unchecked(r));
                let p = _mm256_adds_epi16(*profile.get_unchecked(qr), _mm256_load_si256(bp.add(r)));
                let ph = _mm256_load_si256(hp.add(r));
                let e = _mm256_max_epi16(_mm256_subs_epi16(ph, go_v), _mm256_subs_epi16(_mm256_load_si256(ep.add(r)), ge_v));
                let f = _mm256_max_epi16(_mm256_subs_epi16(u_h, go_v), _mm256_subs_epi16(u_f, ge_v));
                let diag = _mm256_adds_epi16(ul_h, p);
                let h = _mm256_max_epi16(_mm256_max_epi16(e, f), _mm256_max_epi16(diag, zero));
                max_h = _mm256_max_epi16(max_h, _mm256_and_si256(h, active));
                u_f = f;
                ul_h = ph;
                u_h = h;
                _mm256_store_si256(hp.add(r), h);
                _mm256_store_si256(ep.add(r), e);
            }
        }
        let mut res = [0i16; 16];
        _mm256_storeu_si256(res.as_mut_ptr() as *mut __m256i, max_h);
        for l in 0..nlanes {
            out[sc.order[start + l]] = res[l] as i32;
            ctr.simd_cells += (seqs[l].len() * m) as u64;
        }
        ctr.batches += 1;
        start = end;
    }
}

/// 8-bit first pass of [`score_many`]: unsigned saturating lanes with a
/// substitution offset. Writes the score of every subject whose lane maximum
/// stayed below the 8-bit exactness limit and returns the indices of the
/// others (to be re-scored exactly). `None` when the scoring scheme does not
/// fit 8-bit lanes (the caller then uses the 16-bit kernel for everything).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn score_many_u8_avx2(
    query: &[u8],
    bias: Option<&[i8]>,
    subjects: &[ProteinCodes<'_>],
    matrix: &[[i32; 24]; 24],
    gaps: GapCosts,
    sc: &mut SimdScratch,
    ctr: &mut KernelCounters,
    out: &mut [i32],
) -> Option<Vec<usize>> {
    use std::arch::x86_64::*;
    let m = query.len();
    let (mut smin, mut smax) = (i32::MAX, i32::MIN);
    for row in matrix.iter() {
        for &v in row.iter() {
            smin = smin.min(v);
            smax = smax.max(v);
        }
    }
    let (bmin, bmax) = match bias {
        Some(b) => (b.iter().map(|&v| v as i32).min().unwrap_or(0), b.iter().map(|&v| v as i32).max().unwrap_or(0)),
        None => (0, 0),
    };
    let off = -(smin + bmin.min(0)).min(0);
    let pmax = smax + bmax.max(0) + off;
    if pmax >= 128 || off > 127 || !(0..=255).contains(&gaps.open) || !(0..=255).contains(&gaps.extend) {
        return None;
    }
    // A lane maximum below this cannot have saturated (every value in the
    // lane is bounded by its maximum, and an addition saturates only when
    // the operand H is at least 255 - pmax).
    let limit = 255 - pmax;
    const L: usize = 32;
    sc.order.clear();
    sc.order.extend(0..subjects.len());
    sc.order.sort_unstable_by(|&a, &b| subjects[b].len().cmp(&subjects[a].len()).then(a.cmp(&b)));
    sc.h8.clear();
    sc.h8.resize(m, Lane8::default());
    sc.e8.clear();
    sc.e8.resize(m, Lane8::default());
    sc.bias8.clear();
    match bias {
        Some(b) => sc.bias8.extend(b.iter().map(|&v| Lane8([v as u8; 32]))),
        None => sc.bias8.resize(m, Lane8::default()),
    }
    let mut tlo = [_mm256_setzero_si256(); 24];
    let mut thi = [_mm256_setzero_si256(); 24];
    for a in 0..24 {
        let mut t0 = [0u8; 16];
        let mut t1 = [0u8; 16];
        for x in 0..16 {
            t0[x] = (matrix[a][x] + off) as u8;
        }
        for x in 16..24 {
            t1[x - 16] = (matrix[a][x] + off) as u8;
        }
        tlo[a] = _mm256_broadcastsi128_si256(_mm_loadu_si128(t0.as_ptr() as *const __m128i));
        thi[a] = _mm256_broadcastsi128_si256(_mm_loadu_si128(t1.as_ptr() as *const __m128i));
    }
    let go = _mm256_set1_epi8(gaps.open as u8 as i8);
    let ge = _mm256_set1_epi8(gaps.extend as u8 as i8);
    let offv = _mm256_set1_epi8(off as u8 as i8);
    let zero = _mm256_setzero_si256();
    let nib = _mm256_set1_epi8(0x0f);
    let fifteen = _mm256_set1_epi8(15);
    let hp = sc.h8.as_mut_ptr() as *mut __m256i;
    let ep = sc.e8.as_mut_ptr() as *mut __m256i;
    let bp = sc.bias8.as_ptr() as *const __m256i;
    let qcodes: Vec<u8> = query.iter().map(|&c| crate::alphabet::code_index(c) as u8).collect();
    let mut redo = Vec::new();
    let mut start = 0usize;
    while start < sc.order.len() {
        let end = (start + L).min(sc.order.len());
        let nlanes = end - start;
        let mut seqs: [&[u8]; L] = [&[]; L];
        let mut max_len = 0usize;
        for l in 0..nlanes {
            seqs[l] = subjects[sc.order[start + l]].as_slice();
            max_len = max_len.max(seqs[l].len());
        }
        for r in 0..m {
            _mm256_store_si256(hp.add(r), zero);
            _mm256_store_si256(ep.add(r), zero);
        }
        let mut max_h = zero;
        let mut resid = [PAD_RESIDUE; L];
        let mut act = [0u8; L];
        let mut profile = [zero; 24];
        for col in 0..max_len {
            for l in 0..L {
                let s = seqs[l];
                if col < s.len() {
                    resid[l] = crate::alphabet::code_index(*s.get_unchecked(col)) as u8;
                    act[l] = 0xFF;
                } else {
                    resid[l] = PAD_RESIDUE;
                    act[l] = 0;
                }
            }
            let active = _mm256_loadu_si256(act.as_ptr() as *const __m256i);
            let rv = _mm256_loadu_si256(resid.as_ptr() as *const __m256i);
            let idx = _mm256_and_si256(rv, nib);
            let ishi = _mm256_cmpgt_epi8(rv, fifteen);
            for a in 0..24 {
                profile[a] = _mm256_blendv_epi8(_mm256_shuffle_epi8(tlo[a], idx), _mm256_shuffle_epi8(thi[a], idx), ishi);
            }
            let mut u_f = zero;
            let mut u_h = zero;
            let mut ul_h = zero;
            for r in 0..m {
                let p = _mm256_add_epi8(*profile.get_unchecked(*qcodes.get_unchecked(r) as usize), _mm256_load_si256(bp.add(r)));
                let ph = _mm256_load_si256(hp.add(r));
                let e = _mm256_max_epu8(_mm256_subs_epu8(ph, go), _mm256_subs_epu8(_mm256_load_si256(ep.add(r)), ge));
                let f = _mm256_max_epu8(_mm256_subs_epu8(u_h, go), _mm256_subs_epu8(u_f, ge));
                let diag = _mm256_subs_epu8(_mm256_adds_epu8(ul_h, p), offv);
                let h = _mm256_max_epu8(_mm256_max_epu8(e, f), diag);
                max_h = _mm256_max_epu8(max_h, _mm256_and_si256(h, active));
                u_f = f;
                ul_h = ph;
                u_h = h;
                _mm256_store_si256(hp.add(r), h);
                _mm256_store_si256(ep.add(r), e);
            }
        }
        let mut res = [0u8; L];
        _mm256_storeu_si256(res.as_mut_ptr() as *mut __m256i, max_h);
        for l in 0..nlanes {
            let k = sc.order[start + l];
            if res[l] as i32 >= limit {
                redo.push(k);
            } else {
                out[k] = res[l] as i32;
            }
            ctr.simd_cells += (seqs[l].len() * m) as u64;
        }
        ctr.batches += 1;
        start = end;
    }
    redo.sort_unstable();
    Some(redo)
}

// ---------------------------------------------------------------------------
// Backwards-compatible entry points (pre-v1.2 signatures)
// ---------------------------------------------------------------------------

/// Score one encoded query against many encoded subjects (no bias), using the
/// matrix's own default affine costs passed in. Kept for API compatibility;
/// panics if a byte is not a protein code.
pub fn sw_simd_protein_batch(
    query: &[u8],
    subjects: &[Vec<u8>],
    matrix: &ScoreMatrix,
    gap_open: i32,
    gap_extend: i32,
) -> Vec<i32> {
    let q = ProteinCodes::new(query).expect("sw_simd_protein_batch requires encoded protein codes");
    let subs: Vec<ProteinCodes<'_>> = subjects
        .iter()
        .map(|s| ProteinCodes::new(s).expect("sw_simd_protein_batch requires encoded protein codes"))
        .collect();
    let gaps = GapCosts { open: gap_open, extend: gap_extend };
    let mut sc = SimdScratch::new();
    let mut ctr = KernelCounters::default();
    score_many(q, None, &subs, &matrix.matrix, gaps, &mut sc, &mut ctr)
}

/// Slice-based variant of [`sw_simd_protein_batch`].
pub fn sw_simd_protein_batch_encoded(
    query_encoded: &[u8],
    subjects_encoded: &[&[u8]],
    matrix: &ScoreMatrix,
    gap_open: i32,
    gap_extend: i32,
) -> Vec<i32> {
    let q = ProteinCodes::new(query_encoded).expect("encoded protein codes required");
    let subs: Vec<ProteinCodes<'_>> = subjects_encoded
        .iter()
        .map(|s| ProteinCodes::new(s).expect("encoded protein codes required"))
        .collect();
    let gaps = GapCosts { open: gap_open, extend: gap_extend };
    let mut sc = SimdScratch::new();
    let mut ctr = KernelCounters::default();
    score_many(q, None, &subs, &matrix.matrix, gaps, &mut sc, &mut ctr)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::reference::sw_reference;
    use crate::alphabet::EncodedProtein;
    use crate::scoring::{ScoreMatrix, ScoringSystem};

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self, m: u64) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 33) % m
        }
    }

    #[test]
    fn simd_matches_independent_reference() {
        for sys in [ScoringSystem::Blosum62, ScoringSystem::Pam30, ScoringSystem::Blosum45] {
            let m = ScoreMatrix::new(sys);
            let mut r = Lcg(99);
            for trial in 0..20 {
                let ql = 1 + r.next(200) as usize;
                let q: Vec<u8> = (0..ql).map(|_| r.next(24) as u8).collect();
                let bias: Vec<i8> = (0..ql).map(|_| r.next(9) as i8 - 4).collect();
                let n = 1 + r.next(40) as usize;
                let subs: Vec<Vec<u8>> = (0..n)
                    .map(|_| {
                        let l = 1 + r.next(300) as usize;
                        (0..l).map(|_| r.next(24) as u8).collect()
                    })
                    .collect();
                let views: Vec<ProteinCodes> = subs.iter().map(|s| ProteinCodes::new(s).unwrap()).collect();
                for (go, ge) in [(11, 1), (7, 2), (13, 3)] {
                    let gaps = GapCosts { open: go, extend: ge };
                    let use_bias = trial % 2 == 0;
                    let b = if use_bias { Some(&bias[..]) } else { None };
                    let mut sc = SimdScratch::new();
                    let mut ctr = KernelCounters::default();
                    let got = score_many(ProteinCodes::new(&q).unwrap(), b, &views, &m.matrix, gaps, &mut sc, &mut ctr);
                    for (k, s) in subs.iter().enumerate() {
                        let (want, _, _) = sw_reference(&q, s, &m.matrix, b, go as i64, ge as i64);
                        assert_eq!(got[k] as i64, want, "{sys} trial {trial} subj {k} gaps {go}/{ge}");
                    }
                }
            }
        }
    }

    #[test]
    fn saturation_is_detected_and_repaired() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let w = EncodedProtein::from_ascii(&vec![b'W'; 4000]);
        let short = EncodedProtein::from_ascii(b"MKVLAAGIVG");
        let subs = vec![w.codes(), short.codes(), w.codes()];
        let mut sc = SimdScratch::new();
        let mut ctr = KernelCounters::default();
        let got = score_many(w.codes(), None, &subs, &m.matrix, GapCosts { open: 11, extend: 1 }, &mut sc, &mut ctr);
        assert_eq!(got[0], 44000);
        assert_eq!(got[2], 44000);
        assert_eq!(ctr.saturation_retries, 2);
    }

    #[test]
    fn ambiguous_residues_score_as_x_not_stop() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let q = EncodedProtein::from_ascii(b"MKVLXXXXAAGIVG");
        let s = EncodedProtein::from_ascii(b"MKVLXXXXAAGIVG");
        let mut sc = SimdScratch::new();
        let mut ctr = KernelCounters::default();
        let got = score_many(q.codes(), None, &[s.codes()], &m.matrix, GapCosts { open: 11, extend: 1 }, &mut sc, &mut ctr);
        let (want, _, _) = sw_reference(q.as_slice(), s.as_slice(), &m.matrix, None, 11, 1);
        assert_eq!(got[0] as i64, want);
        // X-X is -1 in BLOSUM62; a stop-symbol encoding would give +1 each
        // ('*' vs '*') or -4 against residues — the score must reflect X.
        let as_stop = EncodedProtein::from_ascii(b"MKVL****AAGIVG");
        let (stop_score, _, _) = sw_reference(as_stop.as_slice(), as_stop.as_slice(), &m.matrix, None, 11, 1);
        assert_ne!(want, stop_score);
    }
}
