//! Inter-sequence banded Smith–Waterman score kernel (AVX2, 32 × u8 lanes).
//!
//! Scores many *banded* jobs of one query at once: thirty-two subject windows
//! occupy the 8-bit lanes of one AVX2 register (Rognes' inter-sequence
//! layout) and every lane is restricted to its own diagonal band. Each job is
//! placed in "lane coordinates" `c = j - lo`, so all bands start on the same
//! lane diagonal and a column only visits the rows that can be inside some
//! lane's band; cells outside a lane's band or outside its subject are forced
//! to zero (no path), which for local alignment is exactly "not allowed".
//!
//! ## Exactness
//!
//! Unsigned saturating arithmetic with a substitution offset `B` (so that
//! `s + bias + B >= 0`). Every H/E/F value of a local alignment is bounded by
//! the lane's optimum, so a lane is exact whenever its maximum stays below
//! `255 - max(s + bias + B)`; a lane that reaches that value, a job whose band
//! is wider than the kernel supports, and every job on a host without AVX2
//! are scored by [`banded_score`] (exact, i16/i32). The result for each job
//! equals `banded_score(query, subject, bias, scheme, Band { lo, hi })`;
//! `tests::banded_u8_matches_banded_score` checks this on random inputs.

#![allow(unsafe_code)]

use super::banded::{banded_score, Band, BandScratch, Scheme};
use super::sw_simd::KernelCounters;

/// One banded scoring job: a subject window and the band (diagonals
/// `j - i`, window coordinates) its alignment must stay in.
#[derive(Clone, Copy, Debug)]
pub struct BandJob<'a> {
    pub subject: &'a [u8],
    pub band: Band,
}

/// Lanes per register.
pub const U8_LANES: usize = 32;
/// Widest band (hi - lo) handled by the 8-bit kernel (signed byte compare).
pub const U8_MAX_WIDTH: i64 = 126;
/// Residue code marking a cell column outside a lane's subject (the lane's
/// cells there are forced to zero).
const INVALID: u8 = 31;

#[repr(C, align(32))]
#[derive(Clone, Copy)]
struct V32([u8; 32]);

/// Re-usable buffers for [`score_banded_many`].
#[derive(Default)]
pub struct Band8Scratch {
    h: Vec<V32>,
    e: Vec<V32>,
    bias: Vec<V32>,
    cols: Vec<V32>,
    /// Lane-major residues before the transpose (32 rows of padded length).
    rows: Vec<u8>,
    order: Vec<usize>,
    band: BandScratch,
}

impl Band8Scratch {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Counters specific to the 8-bit banded kernel.
#[derive(Debug, Default, Clone, Copy)]
pub struct Band8Counters {
    /// Lanes re-scored by the exact kernel after reaching the 8-bit limit.
    pub u8_overflows: u64,
    /// Jobs scored by the exact kernel directly (band too wide / no AVX2).
    pub wide_jobs: u64,
}

/// Exact banded SW scores of `query` against every job (one per job, caller
/// order). `bias` (if any) must have the query's length.
#[allow(clippy::too_many_arguments)]
pub fn score_banded_many(
    query: &[u8],
    bias: Option<&[i8]>,
    jobs: &[BandJob<'_>],
    scheme: Scheme<'_>,
    sc: &mut Band8Scratch,
    kc: &mut KernelCounters,
    bc: &mut Band8Counters,
) -> Vec<i32> {
    let m = query.len();
    let mut out = vec![0i32; jobs.len()];
    if jobs.is_empty() || m == 0 {
        return out;
    }
    if let Some(b) = bias {
        assert_eq!(b.len(), m, "bias length must equal query length");
    }
    // Clip every band to its rectangle; empty bands score 0.
    let clipped: Vec<Option<Band>> = jobs
        .iter()
        .map(|j| {
            let n = j.subject.len();
            if n == 0 || j.band.hi < j.band.lo {
                return None;
            }
            let b = j.band.clip(m, n);
            if b.hi < b.lo {
                None
            } else {
                Some(b)
            }
        })
        .collect();
    let mut exact: Vec<usize> = Vec::new();
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            // Substitution offset and the 8-bit exactness limit for this query.
            let (mut smin, mut smax) = (i32::MAX, i32::MIN);
            for a in 0..24 {
                for x in 0..24 {
                    smin = smin.min(scheme.matrix[a][x]);
                    smax = smax.max(scheme.matrix[a][x]);
                }
            }
            let (bmin, bmax) = match bias {
                Some(b) => (b.iter().map(|&v| v as i32).min().unwrap_or(0), b.iter().map(|&v| v as i32).max().unwrap_or(0)),
                None => (0, 0),
            };
            let off = -(smin + bmin.min(0)).min(0);
            let pmax = smax + bmax.max(0) + off;
            let gaps_ok = (0..=255).contains(&scheme.gap_open) && (0..=255).contains(&scheme.gap_extend);
            if pmax < 128 && gaps_ok && off <= 127 {
                let limit = 255 - pmax; // a lane maximum >= limit may have saturated
                sc.order.clear();
                for (k, b) in clipped.iter().enumerate() {
                    match b {
                        Some(b) if b.hi - b.lo <= U8_MAX_WIDTH => sc.order.push(k),
                        Some(_) => exact.push(k),
                        None => {}
                    }
                }
                // Narrow bands together, then by subject length (deterministic).
                let order = std::mem::take(&mut sc.order);
                let mut order = order;
                order.sort_unstable_by_key(|&k| {
                    let b = clipped[k].unwrap();
                    (b.hi - b.lo, jobs[k].subject.len(), k)
                });
                let mut over: Vec<usize> = Vec::new();
                for chunk in order.chunks(U8_LANES) {
                    // SAFETY: AVX2 detected above.
                    let maxes = unsafe { batch_u8(query, bias, jobs, &clipped, chunk, scheme, off, sc, kc) };
                    for (l, &k) in chunk.iter().enumerate() {
                        if maxes[l] as i32 >= limit {
                            over.push(k);
                        } else {
                            out[k] = maxes[l] as i32;
                        }
                    }
                }
                sc.order = order;
                bc.u8_overflows += over.len() as u64;
                exact.extend(over);
            } else {
                exact.extend(clipped.iter().enumerate().filter(|(_, b)| b.is_some()).map(|(k, _)| k));
                bc.wide_jobs += exact.len() as u64;
            }
            for &k in &exact {
                let b = clipped[k].unwrap();
                let r = banded_score(query, jobs[k].subject, bias, scheme, b, &mut sc.band);
                kc.scalar_cells += r.cells;
                out[k] = r.score;
            }
            return out;
        }
    }
    for (k, b) in clipped.iter().enumerate() {
        if let Some(b) = b {
            let r = banded_score(query, jobs[k].subject, bias, scheme, *b, &mut sc.band);
            kc.scalar_cells += r.cells;
            out[k] = r.score;
            bc.wide_jobs += 1;
        }
    }
    let _ = &mut exact;
    out
}

/// Transpose a 16 x 16 byte block: row `i` of the source (stride `sstride`)
/// becomes column `i` of the destination (rows of stride `dstride`).
#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn transpose16x16(src: *const u8, sstride: usize, dst: *mut u8, dstride: usize) {
    use std::arch::x86_64::*;
    let mut r = [_mm_setzero_si128(); 16];
    for (i, x) in r.iter_mut().enumerate() {
        *x = _mm_loadu_si128(src.add(i * sstride) as *const __m128i);
    }
    // bytes of row pairs, 16-bit groups of row quads, 32-bit groups of row
    // octets, then the two octets: each step interleaves twice as many rows
    let mut a = [_mm_setzero_si128(); 16];
    for k in 0..8 {
        a[2 * k] = _mm_unpacklo_epi8(r[2 * k], r[2 * k + 1]);
        a[2 * k + 1] = _mm_unpackhi_epi8(r[2 * k], r[2 * k + 1]);
    }
    let mut b = [_mm_setzero_si128(); 16];
    for q in 0..4 {
        for h in 0..2 {
            b[4 * q + 2 * h] = _mm_unpacklo_epi16(a[4 * q + h], a[4 * q + 2 + h]);
            b[4 * q + 2 * h + 1] = _mm_unpackhi_epi16(a[4 * q + h], a[4 * q + 2 + h]);
        }
    }
    let mut c = [_mm_setzero_si128(); 16];
    for p in 0..2 {
        for h in 0..2 {
            for g in 0..2 {
                c[8 * p + 4 * h + 2 * g] = _mm_unpacklo_epi32(b[8 * p + 2 * h + g], b[8 * p + 4 + 2 * h + g]);
                c[8 * p + 4 * h + 2 * g + 1] = _mm_unpackhi_epi32(b[8 * p + 2 * h + g], b[8 * p + 4 + 2 * h + g]);
            }
        }
    }
    for h in 0..2 {
        for g in 0..2 {
            for e in 0..2 {
                let col = 8 * h + 4 * g + 2 * e;
                let (x, y) = (c[4 * h + 2 * g + e], c[8 + 4 * h + 2 * g + e]);
                _mm_storeu_si128(dst.add(col * dstride) as *mut __m128i, _mm_unpacklo_epi64(x, y));
                _mm_storeu_si128(dst.add((col + 1) * dstride) as *mut __m128i, _mm_unpackhi_epi64(x, y));
            }
        }
    }
}

/// One batch of up to 32 jobs; returns the per-lane maxima (u8).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn batch_u8(
    query: &[u8],
    bias: Option<&[i8]>,
    jobs: &[BandJob<'_>],
    clipped: &[Option<Band>],
    chunk: &[usize],
    scheme: Scheme<'_>,
    off: i32,
    sc: &mut Band8Scratch,
    kc: &mut KernelCounters,
) -> [u8; 32] {
    use std::arch::x86_64::*;
    let m = query.len();
    let nl = chunk.len();
    // Lane geometry: lane column c <-> window column j = c + lo_k; the band
    // of lane k is 0 <= c - i <= w_k.
    let mut lo = [0i64; 32];
    let mut w = [0i64; 32];
    let mut wmax = 0i64;
    for (l, &k) in chunk.iter().enumerate() {
        let b = clipped[k].unwrap();
        lo[l] = b.lo;
        w[l] = b.hi - b.lo;
        wmax = wmax.max(w[l]);
    }
    let ncols = m + wmax as usize;
    // Lane-major residues in lane coordinates (INVALID outside the lane's
    // subject), then one 32 x ncols byte transpose.
    let ncp = (ncols + 15) & !15;
    sc.rows.clear();
    sc.rows.resize(32 * ncp, INVALID);
    let mut cells = 0u64;
    for (l, &k) in chunk.iter().enumerate() {
        let s = jobs[k].subject;
        let n = s.len() as i64;
        // c in [max(0, -lo), min(ncols, n - lo))
        let c0 = (-lo[l]).max(0);
        let c1 = (n - lo[l]).min(ncols as i64);
        if c1 > c0 {
            let row = &mut sc.rows[l * ncp + c0 as usize..l * ncp + c1 as usize];
            for (d, &x) in row.iter_mut().zip(&s[(c0 + lo[l]) as usize..(c1 + lo[l]) as usize]) {
                // codes outside the alphabet score as X (code_index)
                *d = if x < 24 { x } else { crate::alphabet::AA_X };
            }
        }
        // cells inside the clipped band (for the counters): per diagonal
        let b = clipped[k].unwrap();
        for d in b.lo..=b.hi {
            let (ilo, ihi) = ((-d).max(0), (m as i64).min(n - d));
            if ihi > ilo {
                cells += (ihi - ilo) as u64;
            }
        }
    }
    sc.cols.clear();
    sc.cols.resize(ncp, V32([INVALID; 32]));
    {
        let src = sc.rows.as_ptr();
        let dst = sc.cols.as_mut_ptr() as *mut u8;
        let mut c = 0;
        while c < ncp {
            transpose16x16(src.add(c), ncp, dst.add(c * 32), 32);
            transpose16x16(src.add(16 * ncp + c), ncp, dst.add(c * 32 + 16), 32);
            c += 16;
        }
    }
    kc.simd_cells += cells;
    kc.batches += 1;
    let _ = nl;
    // Per-lane band width + 1 (signed bytes; lanes without a job get 0 = no cell).
    let mut wp1 = [0i8; 32];
    for l in 0..nl {
        wp1[l] = (w[l] + 1) as i8;
    }
    let wp1v = _mm256_loadu_si256(wp1.as_ptr() as *const __m256i);
    // Bias rows (two's complement bytes; profile entries carry +off).
    sc.bias.clear();
    match bias {
        Some(b) => sc.bias.extend(b.iter().map(|&v| V32([v as u8; 32]))),
        None => sc.bias.resize(m, V32([0u8; 32])),
    }
    sc.h.clear();
    sc.h.resize(m, V32([0u8; 32]));
    sc.e.clear();
    sc.e.resize(m, V32([0u8; 32]));
    // Shuffle tables: s(a, x) + off as u8 for x in 0..16 and 16..24.
    let mut tlo = [_mm256_setzero_si256(); 24];
    let mut thi = [_mm256_setzero_si256(); 24];
    for a in 0..24 {
        let mut t0 = [0u8; 16];
        let mut t1 = [0u8; 16];
        for x in 0..16 {
            t0[x] = (scheme.matrix[a][x] + off) as u8;
        }
        for x in 16..24 {
            t1[x - 16] = (scheme.matrix[a][x] + off) as u8;
        }
        tlo[a] = _mm256_broadcastsi128_si256(_mm_loadu_si128(t0.as_ptr() as *const __m128i));
        thi[a] = _mm256_broadcastsi128_si256(_mm_loadu_si128(t1.as_ptr() as *const __m128i));
    }
    let go = _mm256_set1_epi8(scheme.gap_open as u8 as i8);
    let ge = _mm256_set1_epi8(scheme.gap_extend as u8 as i8);
    let offv = _mm256_set1_epi8(off as u8 as i8);
    let zero = _mm256_setzero_si256();
    let nib = _mm256_set1_epi8(0x0f);
    let fifteen = _mm256_set1_epi8(15);
    let hp = sc.h.as_mut_ptr() as *mut __m256i;
    let ep = sc.e.as_mut_ptr() as *mut __m256i;
    let bp = sc.bias.as_ptr() as *const __m256i;
    let qcodes: Vec<u8> = query.iter().map(|&q| crate::alphabet::code_index(q) as u8).collect();
    let invalid = _mm256_set1_epi8(INVALID as i8);
    let mut max_h = zero;
    let mut profile = [zero; 24];
    for c in 0..ncols {
        let r0 = (c as i64 - wmax).max(0) as usize;
        let r1 = c.min(m - 1);
        if r0 > r1 {
            continue;
        }
        let rv = _mm256_loadu_si256(sc.cols.as_ptr().add(c) as *const __m256i);
        let colv = _mm256_cmpgt_epi8(invalid, rv);
        let idx = _mm256_and_si256(rv, nib);
        let ishi = _mm256_cmpgt_epi8(rv, fifteen);
        for a in 0..24 {
            profile[a] = _mm256_blendv_epi8(_mm256_shuffle_epi8(tlo[a], idx), _mm256_shuffle_epi8(thi[a], idx), ishi);
        }
        let mut ul_h = if r0 > 0 { _mm256_load_si256(hp.add(r0 - 1)) } else { zero };
        let mut u_h = zero;
        let mut u_f = zero;
        // d = c - i, decreasing by one per row
        let mut dv = _mm256_set1_epi8((c - r0) as i8);
        let one = _mm256_set1_epi8(1);
        for i in r0..=r1 {
            let valid = _mm256_and_si256(colv, _mm256_cmpgt_epi8(wp1v, dv));
            let p = _mm256_add_epi8(*profile.get_unchecked(*qcodes.get_unchecked(i) as usize), _mm256_load_si256(bp.add(i)));
            let ph = _mm256_load_si256(hp.add(i));
            let e = _mm256_and_si256(_mm256_max_epu8(_mm256_subs_epu8(ph, go), _mm256_subs_epu8(_mm256_load_si256(ep.add(i)), ge)), valid);
            let f = _mm256_and_si256(_mm256_max_epu8(_mm256_subs_epu8(u_h, go), _mm256_subs_epu8(u_f, ge)), valid);
            let diag = _mm256_subs_epu8(_mm256_adds_epu8(ul_h, p), offv);
            let h = _mm256_and_si256(_mm256_max_epu8(_mm256_max_epu8(e, f), diag), valid);
            max_h = _mm256_max_epu8(max_h, h);
            ul_h = ph;
            u_h = h;
            u_f = f;
            _mm256_store_si256(hp.add(i), h);
            _mm256_store_si256(ep.add(i), e);
            dv = _mm256_sub_epi8(dv, one);
        }
    }
    let mut res = [0u8; 32];
    _mm256_storeu_si256(res.as_mut_ptr() as *mut __m256i, max_h);
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{ScoreMatrix, ScoringSystem};

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self, m: u64) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 33) % m
        }
    }

    #[test]
    fn banded_u8_matches_banded_score() {
        for sys in [ScoringSystem::Blosum62, ScoringSystem::Pam30, ScoringSystem::Blosum45] {
            let mx = ScoreMatrix::new(sys);
            let mut r = Lcg(7);
            for trial in 0..40 {
                let ql = 1 + r.next(220) as usize;
                let q: Vec<u8> = (0..ql).map(|_| r.next(24) as u8).collect();
                let bias: Vec<i8> = (0..ql).map(|_| r.next(9) as i8 - 4).collect();
                let n = 1 + r.next(70) as usize;
                let subs: Vec<Vec<u8>> = (0..n)
                    .map(|k| {
                        let l = 1 + r.next(320) as usize;
                        let mut s: Vec<u8> = (0..l).map(|_| r.next(24) as u8).collect();
                        if k % 3 == 0 {
                            // a related stretch so some lanes score high
                            let o = r.next(l as u64) as usize;
                            for t in 0..l.saturating_sub(o).min(ql) {
                                if r.next(4) > 0 {
                                    s[o + t] = q[t];
                                }
                            }
                        }
                        s
                    })
                    .collect();
                let jobs: Vec<BandJob> = subs
                    .iter()
                    .map(|s| {
                        let c = r.next((s.len() + ql) as u64) as i64 - ql as i64;
                        let w = r.next(if trial % 5 == 0 { 200 } else { 90 }) as i64;
                        BandJob { subject: s, band: Band { lo: c - w / 2, hi: c + (w + 1) / 2 } }
                    })
                    .collect();
                for (go, ge) in [(11, 1), (7, 2)] {
                    let scheme = Scheme { matrix: &mx.matrix, gap_open: go, gap_extend: ge };
                    let b = if trial % 2 == 0 { Some(&bias[..]) } else { None };
                    let mut sc = Band8Scratch::new();
                    let (mut kc, mut bc) = (KernelCounters::default(), Band8Counters::default());
                    let got = score_banded_many(&q, b, &jobs, scheme, &mut sc, &mut kc, &mut bc);
                    let mut bs = BandScratch::new();
                    for (k, j) in jobs.iter().enumerate() {
                        let want = banded_score(&q, j.subject, b, scheme, j.band, &mut bs).score;
                        assert_eq!(got[k], want, "{sys} trial {trial} job {k} band {:?} gaps {go}/{ge}", j.band);
                    }
                }
            }
        }
    }

    #[test]
    fn overflowing_lanes_are_rescored_exactly() {
        let mx = ScoreMatrix::new(ScoringSystem::Blosum62);
        let q = crate::alphabet::EncodedProtein::from_ascii(&vec![b'W'; 300]);
        let other = crate::alphabet::EncodedProtein::from_ascii(b"MKVLAAGIVGKLLA");
        let jobs = vec![
            BandJob { subject: q.as_slice(), band: Band { lo: -5, hi: 5 } },
            BandJob { subject: other.as_slice(), band: Band { lo: -3, hi: 3 } },
        ];
        let scheme = Scheme { matrix: &mx.matrix, gap_open: 11, gap_extend: 1 };
        let mut sc = Band8Scratch::new();
        let (mut kc, mut bc) = (KernelCounters::default(), Band8Counters::default());
        let got = score_banded_many(q.as_slice(), None, &jobs, scheme, &mut sc, &mut kc, &mut bc);
        assert_eq!(got[0], 300 * 11);
        assert_eq!(bc.u8_overflows, 1);
    }
}
