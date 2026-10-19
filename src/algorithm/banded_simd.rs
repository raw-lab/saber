//! AVX2 implementation of the banded affine-gap DP (score-only and with
//! direction bits for traceback). Eight band cells of one query row are
//! computed per vector (i32 lanes, so no overflow concerns). The horizontal
//! (query-gap) recurrence `E[k] = max(H[k−1] − go, E[k−1] − ge)` is resolved
//! inside each vector with a max-plus prefix scan:
//!
//! ```text
//! E[k] = max_{k'<k} (H'[k'] + k'·ge) − go − (k−1)·ge,   H' = max(D, F, 0)
//! ```
//! (valid because go ≥ ge, so re-opening a gap is never better than
//! extending it). The scan of a chunk is computed without the value entering
//! from the previous chunk, which is then applied as `max(E, carry − t·ge)`;
//! the chunk-to-chunk dependency is two vector operations, so successive
//! chunks overlap in the pipeline. Direction bits are packed to bytes with
//! vector packs and the row maximum is resolved once per row. Results —
//! scores, end cells, direction bits and therefore
//! alignments — are identical to the scalar kernel in [`super::banded`];
//! `banded::tests::simd_and_scalar_kernels_agree` checks this on random and
//! related pairs with and without bias.

#![allow(unsafe_code)]

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use super::banded::{Band, Scheme};

pub(crate) const NEG: i32 = i32::MIN / 4;

/// Per-call outputs of the SIMD kernel.
pub(crate) struct SimdBandOut {
    pub score: i32,
    pub best_i: usize,
    pub best_k: usize,
    pub cells: u64,
}

/// `[prev[7], v0, v1, ..., v6]` (lane-crossing shift with a vector carry).
#[cfg(target_arch = "x86_64")]
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn shift_in_vec(v: __m256i, prev: __m256i) -> __m256i {
    // t = [prev.hi128, v.lo128]; per 128-bit lane (v : t) >> 12 bytes.
    let t = _mm256_permute2x128_si256(prev, v, 0x21);
    _mm256_alignr_epi8(v, t, 12)
}

/// `[v1, ..., v7, next[0]]` (i32 lanes): the previous row's values one band
/// cell to the right, built from two loads aligned with the previous row's
/// stores (an unaligned reload of just-stored data would defeat store
/// forwarding).
#[cfg(target_arch = "x86_64")]
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn shift_out(v: __m256i, next: __m256i) -> __m256i {
    let t = _mm256_permute2x128_si256(v, next, 0x21); // [v.hi, next.lo]
    _mm256_alignr_epi8(t, v, 4)
}

#[cfg(target_arch = "x86_64")]
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn prefix_max(v: __m256i) -> __m256i {
    let neg = _mm256_set1_epi32(NEG);
    let s1 = _mm256_blend_epi32(_mm256_permutevar8x32_epi32(v, _mm256_setr_epi32(0, 0, 1, 2, 3, 4, 5, 6)), neg, 0b0000_0001);
    let v = _mm256_max_epi32(v, s1);
    let s2 = _mm256_blend_epi32(_mm256_permutevar8x32_epi32(v, _mm256_setr_epi32(0, 0, 0, 1, 2, 3, 4, 5)), neg, 0b0000_0011);
    let v = _mm256_max_epi32(v, s2);
    let s4 = _mm256_blend_epi32(_mm256_permutevar8x32_epi32(v, _mm256_setr_epi32(0, 0, 0, 0, 0, 1, 2, 3)), neg, 0b0000_1111);
    _mm256_max_epi32(v, s4)
}

#[cfg(target_arch = "x86_64")]
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn hmax(v: __m256i) -> i32 {
    let a = _mm_max_epi32(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1));
    let b = _mm_max_epi32(a, _mm_shuffle_epi32(a, 0b01_00_11_10));
    let c = _mm_max_epi32(b, _mm_shuffle_epi32(b, 0b10_11_00_01));
    _mm_cvtsi128_si32(c)
}

/// First `k` in `kmin..=kmax` with `h[k] == v` (vectorised; `h` holds at
/// least `kmax + 8` entries). Returns `kmin` if none (cannot happen for the
/// row maximum).
#[cfg(target_arch = "x86_64")]
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn first_eq32(h: &[i32], kmin: usize, kmax: usize, v: i32) -> usize {
    let vv = _mm256_set1_epi32(v);
    let mut k = kmin;
    while k <= kmax {
        let m = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(_mm256_loadu_si256(h.as_ptr().add(k) as *const __m256i), vv))) as u32;
        if m != 0 {
            let p = k + m.trailing_zeros() as usize;
            return if p <= kmax { p } else { kmin };
        }
        k += 8;
    }
    kmin
}

/// 16-bit version of [`first_eq32`] (`h` holds at least `kmax + 16` entries).
#[cfg(target_arch = "x86_64")]
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn first_eq16(h: &[i16], kmin: usize, kmax: usize, v: i16) -> usize {
    let vv = _mm256_set1_epi16(v);
    let mut k = kmin;
    while k <= kmax {
        let m = _mm256_movemask_epi8(_mm256_cmpeq_epi16(_mm256_loadu_si256(h.as_ptr().add(k) as *const __m256i), vv)) as u32;
        if m != 0 {
            let p = k + (m.trailing_zeros() / 2) as usize;
            return if p <= kmax { p } else { kmin };
        }
        k += 16;
    }
    kmin
}

/// Banded DP over rows `0..m`. `hp`, `fp`, `hc`, `fc` must hold at least
/// `w + 16` entries; `dir` (when `TB`) at least `m*w + 16` bytes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn banded_avx2<const TB: bool>(
    q: &[u8],
    s: &[u8],
    bias: Option<&[i8]>,
    scheme: Scheme<'_>,
    band: Band,
    w: usize,
    hp: &mut Vec<i32>,
    fp: &mut Vec<i32>,
    hc: &mut Vec<i32>,
    fc: &mut Vec<i32>,
    dir: &mut [u8],
    stop_at: i32,
) -> SimdBandOut {
    let m = q.len();
    let n = s.len();
    let go = scheme.gap_open;
    let ge = scheme.gap_extend;
    let wp = w + 16;
    for v in [&mut *hp, &mut *hc] {
        v.clear();
        v.resize(wp, 0);
    }
    for v in [&mut *fp, &mut *fc] {
        v.clear();
        v.resize(wp, NEG);
    }
    // Row -1 is the boundary: H = 0, F = -inf.
    // Per query letter: shuffle tables (i8) for subject letters 0..15 / 16..23.
    let mut tlo = [_mm_setzero_si128(); 24];
    let mut thi = [_mm_setzero_si128(); 24];
    for a in 0..24 {
        let mut lo = [0i8; 16];
        let mut hi = [0i8; 16];
        for x in 0..16 {
            lo[x] = scheme.matrix[a][x].clamp(-128, 127) as i8;
        }
        for x in 16..24 {
            hi[x - 16] = scheme.matrix[a][x].clamp(-128, 127) as i8;
        }
        tlo[a] = _mm_loadu_si128(lo.as_ptr() as *const __m128i);
        thi[a] = _mm_loadu_si128(hi.as_ptr() as *const __m128i);
    }
    let low_nib = _mm_set1_epi8(0x0f);
    let fifteen = _mm_set1_epi8(15);
    let zero = _mm256_setzero_si256();
    let negv = _mm256_set1_epi32(NEG);
    let gov = _mm256_set1_epi32(go);
    let gev = _mm256_set1_epi32(ge);
    let lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    let ramp = _mm256_mullo_epi32(lane, gev); // t·ge
    let one8 = _mm256_set1_epi32(1);
    let two8 = _mm256_set1_epi32(2);
    let three8 = _mm256_set1_epi32(3);
    let four8 = _mm256_set1_epi32(4);
    let eight8 = _mm256_set1_epi32(8);
    let ge8 = _mm256_set1_epi32(8 * ge);
    let idx7 = _mm256_set1_epi32(7);

    let mut best = 0i32;
    let mut best_i = 0usize;
    let mut best_k = 0usize;
    let mut cells = 0u64;
    let mut sbuf = [23u8; 16];

    for i in 0..m {
        let ii = i as i64;
        let kmin = (-(ii + band.lo)).max(0);
        let kmax = (n as i64 - 1 - ii - band.lo).min(w as i64 - 1);
        // Reset current row: boundary 0 left of the band start, -inf elsewhere.
        if kmin > kmax {
            for k in 0..wp {
                hc[k] = if (k as i64) < kmin { 0 } else { NEG };
                fc[k] = NEG;
            }
            std::mem::swap(hp, hc);
            std::mem::swap(fp, fc);
            continue;
        }
        let (kmin, kmax) = (kmin as usize, kmax as usize);
        for k in 0..kmin {
            hc[k] = 0;
            fc[k] = NEG;
        }
        let qi = crate::alphabet::code_index(q[i]);
        let tl = tlo[qi];
        let th = thi[qi];
        let bv = _mm256_set1_epi32(bias.map_or(0, |b| b[i] as i32));
        // E entering the next chunk (broadcast). Within a chunk
        //   E = max(E_local, carry - t·ge)
        // where E_local is the in-chunk scan without the carry, so the only
        // loop-carried dependency between chunks is
        //   carry' = max(max(H'[7] - go, E_local[7] - ge), carry - 8·ge)
        // (two vector ops), and the scans of successive chunks overlap.
        let mut carry = negv;
        let mut h_prev = negv; // previous chunk's H (lane 7: H of the previous cell)
        let mut e_prev = negv; // previous chunk's E
        let mut rowmax = zero;
        let base = ii + band.lo;
        let dbase = i * w;
        let mut k0 = kmin;
        let mut hcur = _mm256_loadu_si256(hp.as_ptr().add(k0) as *const __m256i);
        let mut fcur = _mm256_loadu_si256(fp.as_ptr().add(k0) as *const __m256i);
        while k0 <= kmax {
            let nvalid = (kmax + 1 - k0).min(8);
            let full = nvalid == 8;
            // Subject residues j = base + k0 .. +8
            let j0 = (base + k0 as i64) as usize;
            let sv = if j0 + 16 <= n {
                _mm_loadu_si128(s.as_ptr().add(j0) as *const __m128i)
            } else {
                for t in 0..16 {
                    sbuf[t] = if j0 + t < n { crate::alphabet::code_index(s[j0 + t]) as u8 } else { 23 };
                }
                _mm_loadu_si128(sbuf.as_ptr() as *const __m128i)
            };
            let idx = _mm_and_si128(sv, low_nib);
            let is_hi = _mm_cmpgt_epi8(sv, fifteen);
            let sub8 = _mm_blendv_epi8(_mm_shuffle_epi8(tl, idx), _mm_shuffle_epi8(th, idx), is_hi);
            let sub = _mm256_cvtepi8_epi32(sub8);
            let hpk = hcur;
            let hnext = _mm256_loadu_si256(hp.as_ptr().add(k0 + 8) as *const __m256i);
            let fnext = _mm256_loadu_si256(fp.as_ptr().add(k0 + 8) as *const __m256i);
            let hp1 = shift_out(hpk, hnext);
            let fp1 = shift_out(fcur, fnext);
            hcur = hnext;
            fcur = fnext;
            // The cell (i-1, j) for the last band diagonal is outside the
            // band: hp/fp hold -inf there by construction (index w).
            let d = _mm256_add_epi32(_mm256_add_epi32(hpk, sub), bv);
            let fo = _mm256_sub_epi32(hp1, gov);
            let fe = _mm256_sub_epi32(fp1, gev);
            let f = _mm256_max_epi32(fo, fe);
            let hprime = _mm256_max_epi32(_mm256_max_epi32(d, f), zero);
            let valid = if full { _mm256_set1_epi32(-1) } else { _mm256_cmpgt_epi32(_mm256_set1_epi32(nvalid as i32), lane) };
            let hprime_v = if full { hprime } else { _mm256_blendv_epi8(negv, hprime, valid) };
            // E scan: V = [carry, H'0-go, ..., H'6-go]; E = prefmax(V + t·ge) - t·ge
            let open = _mm256_sub_epi32(hprime_v, gov);
            let e_loc = _mm256_sub_epi32(prefix_max(_mm256_add_epi32(shift_in_vec(open, negv), ramp)), ramp);
            let e = _mm256_max_epi32(e_loc, _mm256_sub_epi32(carry, ramp));
            let h = _mm256_max_epi32(hprime, e);
            let (h_v, f_v) = if full { (h, f) } else { (_mm256_blendv_epi8(negv, h, valid), _mm256_blendv_epi8(negv, f, valid)) };
            _mm256_storeu_si256(hc.as_mut_ptr().add(k0) as *mut __m256i, h_v);
            _mm256_storeu_si256(fc.as_mut_ptr().add(k0) as *mut __m256i, f_v);
            if TB {
                // origin bits: diag(1) > 0, then E(2) strictly better, then F(3)
                let c1 = _mm256_cmpgt_epi32(d, zero);
                let h1 = _mm256_max_epi32(d, zero);
                let c2 = _mm256_cmpgt_epi32(e, h1);
                let h2 = _mm256_max_epi32(h1, e);
                let c3 = _mm256_cmpgt_epi32(f, h2);
                let mut org = _mm256_and_si256(c1, one8);
                org = _mm256_blendv_epi8(org, two8, c2);
                org = _mm256_blendv_epi8(org, three8, c3);
                // e_ext: E[k-1]-ge > H[k-1]-go (previous cell of the row)
                let hs = shift_in_vec(h_v, h_prev);
                let es = shift_in_vec(e, e_prev);
                let eext = _mm256_cmpgt_epi32(_mm256_sub_epi32(es, gev), _mm256_sub_epi32(hs, gov));
                let fext = _mm256_cmpgt_epi32(fe, fo);
                let bits = _mm256_or_si256(org, _mm256_or_si256(_mm256_and_si256(eext, four8), _mm256_and_si256(fext, eight8)));
                let drow = dir.as_mut_ptr().add(dbase + k0);
                if full {
                    // 8 x i32 (values 0..15) -> 8 bytes
                    let p16 = _mm256_packus_epi32(bits, bits);
                    let p8 = _mm256_packus_epi16(p16, p16);
                    let packed = _mm_unpacklo_epi32(_mm256_castsi256_si128(p8), _mm256_extracti128_si256(p8, 1));
                    _mm_storel_epi64(drow as *mut __m128i, packed);
                } else {
                    let mut ba = [0i32; 8];
                    _mm256_storeu_si256(ba.as_mut_ptr() as *mut __m256i, bits);
                    for t in 0..nvalid {
                        *drow.add(t) = ba[t] as u8;
                    }
                }
            }
            rowmax = _mm256_max_epi32(rowmax, h_v);
            if full {
                let a = _mm256_max_epi32(_mm256_sub_epi32(hprime_v, gov), _mm256_sub_epi32(e_loc, gev));
                carry = _mm256_max_epi32(_mm256_permutevar8x32_epi32(a, idx7), _mm256_sub_epi32(carry, ge8));
                h_prev = h_v;
                e_prev = e;
            }
            k0 += 8;
        }
        cells += (kmax - kmin + 1) as u64;
        // Best cell (first occurrence in row-major order).
        let rm = hmax(rowmax);
        if rm > best {
            best = rm;
            best_i = i;
            best_k = first_eq32(hc, kmin, kmax, rm);
        }
        for k in kmax + 1..wp {
            hc[k] = NEG;
            fc[k] = NEG;
        }
        std::mem::swap(hp, hc);
        std::mem::swap(fp, fc);
        if best >= stop_at {
            // No later cell can score higher and ties keep the first cell.
            break;
        }
    }
    SimdBandOut { score: best, best_i, best_k, cells }
}

/// Sentinel for unreachable cells in the 16-bit kernel. Far enough from
/// `i16::MIN` that subtracting any gap cost (≤ 1000) does not saturate, and
/// far enough below every reachable value (≥ −1000) that adding the scan
/// ramp (≤ 15·1000) keeps it below them.
pub(crate) const NEG16: i16 = -31_000;

#[cfg(target_arch = "x86_64")]
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn shift_in16(v: __m256i, prev: __m256i) -> __m256i {
    // [prev[15], v0, ..., v14]
    let t = _mm256_permute2x128_si256(prev, v, 0x21);
    _mm256_alignr_epi8(v, t, 14)
}

/// `[v1, ..., v15, next[0]]` (i16 lanes); see [`shift_out`].
#[cfg(target_arch = "x86_64")]
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn shift_out16(v: __m256i, next: __m256i) -> __m256i {
    let t = _mm256_permute2x128_si256(v, next, 0x21);
    _mm256_alignr_epi8(t, v, 2)
}

#[cfg(target_arch = "x86_64")]
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn prefix_max16(v: __m256i, negv: __m256i) -> __m256i {
    // within 128-bit halves: shifts by 1, 2, 4 lanes (NEG fill)
    let v = _mm256_max_epi16(v, _mm256_alignr_epi8(v, negv, 14));
    let v = _mm256_max_epi16(v, _mm256_alignr_epi8(v, negv, 12));
    let v = _mm256_max_epi16(v, _mm256_alignr_epi8(v, negv, 8));
    // high half: max with the low half's last lane (word 7 broadcast)
    let t = _mm256_permute2x128_si256(v, v, 0x08); // [0, v.lo]
    let b = _mm256_shuffle_epi32(_mm256_shufflehi_epi16(t, 0xFF), 0xFF);
    _mm256_max_epi16(v, _mm256_blend_epi32(negv, b, 0xF0))
}

#[cfg(target_arch = "x86_64")]
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn hmax16(v: __m256i) -> i16 {
    let a = _mm_max_epi16(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1));
    let b = _mm_max_epi16(a, _mm_shuffle_epi32(a, 0b01_00_11_10));
    let c = _mm_max_epi16(b, _mm_shuffle_epi32(b, 0b10_11_00_01));
    let d = _mm_max_epi16(c, _mm_shufflelo_epi16(c, 0b10_11_00_01));
    _mm_extract_epi16(d, 0) as u16 as i16
}

/// 16-bit variant of [`banded_avx2`]: sixteen band cells per vector with
/// saturating arithmetic. Returns `None` as soon as a row's best score
/// exceeds the range in which every intermediate value provably fits
/// (`32767 − 16·ge − 192`); the caller then reruns the 32-bit kernel. When it
/// returns `Some`, scores, end cells and direction bits are identical to the
/// 32-bit and scalar kernels: every comparison that decides a result
/// involves a reachable value (H ≥ 0, E/F ≥ −go), which is exact, and
/// unreachable values only ever compare below them. `hp`, `fp`, `hc`, `fc`
/// must hold `w + 32` entries; `dir` (when `TB`) at least `m*w + 16` bytes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn banded_avx2_16<const TB: bool>(
    q: &[u8],
    s: &[u8],
    bias: Option<&[i8]>,
    scheme: Scheme<'_>,
    band: Band,
    w: usize,
    hp: &mut Vec<i16>,
    fp: &mut Vec<i16>,
    hc: &mut Vec<i16>,
    fc: &mut Vec<i16>,
    prof: &mut Vec<i8>,
    dir: &mut [u8],
    stop_at: i32,
) -> Option<SimdBandOut> {
    let m = q.len();
    let n = s.len();
    let go = scheme.gap_open;
    let ge = scheme.gap_extend;
    if !(0..=1000).contains(&go) || !(0..=1000).contains(&ge) || ge > go {
        return None;
    }
    let limit = 32_767 - 16 * ge - 192;
    let wp = w + 32;
    for v in [&mut *hp, &mut *hc] {
        v.clear();
        v.resize(wp, 0);
    }
    for v in [&mut *fp, &mut *fc] {
        v.clear();
        v.resize(wp, NEG16);
    }
    // Substitution profile of the subject for the query letters present:
    // prof[a·stride + j] = S(a, s_j), padded with S(a, '*') like the scalar
    // kernel's out-of-range handling. One 16-byte load per chunk then
    // replaces the per-chunk table lookups (and frees their registers).
    let stride = (n + 32 + 15) & !15;
    let mut present = 0u32;
    for &c in q {
        present |= 1 << crate::alphabet::code_index(c);
    }
    prof.clear();
    prof.resize(25 * stride, 0);
    {
        // row 24 holds the subject codes (clamped like `code_index`, padded with '*')
        let (rows, codes) = prof.split_at_mut(24 * stride);
        for (j, c) in codes.iter_mut().enumerate() {
            *c = if j < n { crate::alphabet::code_index(s[j]) as i8 } else { 23 };
        }
        let low_nib = _mm_set1_epi8(0x0f);
        let fifteen = _mm_set1_epi8(15);
        for a in 0..24 {
            if present & (1 << a) == 0 {
                continue;
            }
            let mut lo = [0i8; 16];
            let mut hi = [0i8; 16];
            for x in 0..16 {
                lo[x] = scheme.matrix[a][x].clamp(-128, 127) as i8;
            }
            for x in 16..24 {
                hi[x - 16] = scheme.matrix[a][x].clamp(-128, 127) as i8;
            }
            let tl = _mm_loadu_si128(lo.as_ptr() as *const __m128i);
            let th = _mm_loadu_si128(hi.as_ptr() as *const __m128i);
            let dst = rows.as_mut_ptr().add(a * stride);
            let mut j = 0;
            while j < stride {
                let sv = _mm_loadu_si128(codes.as_ptr().add(j) as *const __m128i);
                let idx = _mm_and_si128(sv, low_nib);
                let is_hi = _mm_cmpgt_epi8(sv, fifteen);
                let v = _mm_blendv_epi8(_mm_shuffle_epi8(tl, idx), _mm_shuffle_epi8(th, idx), is_hi);
                _mm_storeu_si128(dst.add(j) as *mut __m128i, v);
                j += 16;
            }
        }
    }
    let zero = _mm256_setzero_si256();
    let negv = _mm256_set1_epi16(NEG16);
    let gov = _mm256_set1_epi16(go as i16);
    let gev = _mm256_set1_epi16(ge as i16);
    let lane = _mm256_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
    let ramp = _mm256_mullo_epi16(lane, gev); // t·ge ≤ 15 000
    let ge16 = _mm256_set1_epi16((16 * ge) as i16);

    let mut best = 0i32;
    let mut best_i = 0usize;
    let mut best_k = 0usize;
    let mut cells = 0u64;

    for i in 0..m {
        let ii = i as i64;
        let kmin = (-(ii + band.lo)).max(0);
        let kmax = (n as i64 - 1 - ii - band.lo).min(w as i64 - 1);
        if kmin > kmax {
            for k in 0..wp {
                hc[k] = if (k as i64) < kmin { 0 } else { NEG16 };
                fc[k] = NEG16;
            }
            std::mem::swap(hp, hc);
            std::mem::swap(fp, fc);
            continue;
        }
        let (kmin, kmax) = (kmin as usize, kmax as usize);
        for k in 0..kmin {
            hc[k] = 0;
            fc[k] = NEG16;
        }
        let prow = prof.as_ptr().add(crate::alphabet::code_index(q[i]) * stride);
        let bv = _mm256_set1_epi16(bias.map_or(0, |b| b[i] as i16));
        let mut carry = negv;
        let mut h_prev = negv;
        let mut e_prev = negv;
        let mut rowmax = zero;
        let base = ii + band.lo;
        let dbase = i * w;
        let mut k0 = kmin;
        let mut hcur = _mm256_loadu_si256(hp.as_ptr().add(k0) as *const __m256i);
        let mut fcur = _mm256_loadu_si256(fp.as_ptr().add(k0) as *const __m256i);
        while k0 <= kmax {
            let nvalid = (kmax + 1 - k0).min(16);
            let full = nvalid == 16;
            let j0 = (base + k0 as i64) as usize;
            let sub = _mm256_cvtepi8_epi16(_mm_loadu_si128(prow.add(j0) as *const __m128i));
            let hpk = hcur;
            let hnext = _mm256_loadu_si256(hp.as_ptr().add(k0 + 16) as *const __m256i);
            let fnext = _mm256_loadu_si256(fp.as_ptr().add(k0 + 16) as *const __m256i);
            let hp1 = shift_out16(hpk, hnext);
            let fp1 = shift_out16(fcur, fnext);
            hcur = hnext;
            fcur = fnext;
            let d = _mm256_adds_epi16(_mm256_adds_epi16(hpk, sub), bv);
            let fo = _mm256_subs_epi16(hp1, gov);
            let fe = _mm256_subs_epi16(fp1, gev);
            let f = _mm256_max_epi16(fo, fe);
            let hprime = _mm256_max_epi16(_mm256_max_epi16(d, f), zero);
            let valid = if full { _mm256_set1_epi16(-1) } else { _mm256_cmpgt_epi16(_mm256_set1_epi16(nvalid as i16), lane) };
            let hprime_v = if full { hprime } else { _mm256_blendv_epi8(negv, hprime, valid) };
            let open = _mm256_subs_epi16(hprime_v, gov);
            let e_loc = _mm256_subs_epi16(prefix_max16(_mm256_adds_epi16(shift_in16(open, negv), ramp), negv), ramp);
            let e = _mm256_max_epi16(e_loc, _mm256_subs_epi16(carry, ramp));
            let h = _mm256_max_epi16(hprime, e);
            let (h_v, f_v) = if full { (h, f) } else { (_mm256_blendv_epi8(negv, h, valid), _mm256_blendv_epi8(negv, f, valid)) };
            _mm256_storeu_si256(hc.as_mut_ptr().add(k0) as *mut __m256i, h_v);
            _mm256_storeu_si256(fc.as_mut_ptr().add(k0) as *mut __m256i, f_v);
            if TB {
                // origin: 1 diagonal (> 0), 2 E strictly better, 3 F strictly
                // better — priority equals value, so a max of the masked codes
                // (mask >> 15 = 1, (mask >> 15) << 1 = 2, mask >> 14 = 3).
                let c1 = _mm256_cmpgt_epi16(d, zero);
                let h1 = _mm256_max_epi16(d, zero);
                let c2 = _mm256_cmpgt_epi16(e, h1);
                let h2 = _mm256_max_epi16(h1, e);
                let c3 = _mm256_cmpgt_epi16(f, h2);
                let org = _mm256_max_epi16(
                    _mm256_max_epi16(_mm256_srli_epi16(c1, 15), _mm256_slli_epi16(_mm256_srli_epi16(c2, 15), 1)),
                    _mm256_srli_epi16(c3, 14),
                );
                let hs = shift_in16(h_v, h_prev);
                let es = shift_in16(e, e_prev);
                let eext = _mm256_cmpgt_epi16(_mm256_subs_epi16(es, gev), _mm256_subs_epi16(hs, gov));
                let fext = _mm256_cmpgt_epi16(fe, fo);
                let bits = _mm256_or_si256(
                    org,
                    _mm256_or_si256(_mm256_slli_epi16(_mm256_srli_epi16(eext, 15), 2), _mm256_slli_epi16(_mm256_srli_epi16(fext, 15), 3)),
                );
                let drow = dir.as_mut_ptr().add(dbase + k0);
                let p8 = _mm256_permute4x64_epi64(_mm256_packus_epi16(bits, bits), 0b00_00_10_00);
                if full {
                    _mm_storeu_si128(drow as *mut __m128i, _mm256_castsi256_si128(p8));
                } else {
                    let mut ba = [0u8; 16];
                    _mm_storeu_si128(ba.as_mut_ptr() as *mut __m128i, _mm256_castsi256_si128(p8));
                    std::ptr::copy_nonoverlapping(ba.as_ptr(), drow, nvalid);
                }
            }
            rowmax = _mm256_max_epi16(rowmax, h_v);
            if full {
                let a = _mm256_max_epi16(_mm256_subs_epi16(hprime_v, gov), _mm256_subs_epi16(e_loc, gev));
                let a15 = _mm256_shufflehi_epi16(_mm256_shufflelo_epi16(_mm256_permute4x64_epi64(a, 0xFF), 0xFF), 0xFF);
                carry = _mm256_max_epi16(a15, _mm256_subs_epi16(carry, ge16));
                h_prev = h_v;
                e_prev = e;
            }
            k0 += 16;
        }
        cells += (kmax - kmin + 1) as u64;
        let rm = hmax16(rowmax) as i32;
        if rm > limit {
            return None;
        }
        if rm > best {
            best = rm;
            best_i = i;
            best_k = first_eq16(hc, kmin, kmax, rm as i16);
        }
        // cells right of the band read by the next row (≤ kmax + 16)
        _mm256_storeu_si256(hc.as_mut_ptr().add(kmax + 1) as *mut __m256i, negv);
        _mm256_storeu_si256(fc.as_mut_ptr().add(kmax + 1) as *mut __m256i, negv);
        std::mem::swap(hp, hc);
        std::mem::swap(fp, fc);
        if best >= stop_at {
            break;
        }
    }
    Some(SimdBandOut { score: best, best_i, best_k, cells })
}
