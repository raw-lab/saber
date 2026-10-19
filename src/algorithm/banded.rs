//! Genuine banded affine-gap Smith–Waterman (score-only and traceback).
//!
//! The band is expressed in diagonal coordinates `d = j − i` (subject position
//! minus query position, 0-based): only cells with `band.lo ≤ j − i ≤ band.hi`
//! are evaluated. Cells outside the band are treated as unreachable (−∞), so
//! the result is the optimal local alignment *constrained to the band*; when
//! the band covers every diagonal of the rectangle the result is identical to
//! full Smith–Waterman (same recurrences, same tie-breaking as
//! [`super::reference::sw_reference`]).
//!
//! Memory: score-only uses O(band width); traceback stores one direction byte
//! per evaluated cell, i.e. O(rows × band width) instead of O(rows × columns).
//!
//! Recurrences (SABER convention, first gap residue costs `gap_open`):
//! ```text
//! E(i,j) = max(H(i,j−1) − go, E(i,j−1) − ge)        query gap   (SAM 'D')
//! F(i,j) = max(H(i−1,j) − go, F(i−1,j) − ge)        subject gap (SAM 'I')
//! H(i,j) = max(0, H(i−1,j−1) + S(q_i,s_j) + bias_i, E(i,j), F(i,j))
//! ```

/// Substitution matrix + affine gap costs (SABER convention).
#[derive(Clone, Copy)]
pub struct Scheme<'a> {
    pub matrix: &'a [[i32; 24]; 24],
    pub gap_open: i32,
    pub gap_extend: i32,
}

/// Inclusive diagonal band `lo ≤ j − i ≤ hi`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Band {
    pub lo: i64,
    pub hi: i64,
}

impl Band {
    /// Band covering every diagonal of an `m × n` rectangle.
    pub fn full(m: usize, n: usize) -> Self {
        Self { lo: -(m as i64), hi: n as i64 }
    }
    pub fn width(&self) -> usize {
        (self.hi - self.lo + 1).max(0) as usize
    }
    /// Clip to the diagonals that exist in an `m × n` rectangle.
    pub fn clip(&self, m: usize, n: usize) -> Self {
        Self { lo: self.lo.max(-(m as i64 - 1)), hi: self.hi.min(n as i64 - 1) }
    }
    pub fn covers(&self, m: usize, n: usize) -> bool {
        self.lo <= -(m as i64 - 1) && self.hi >= n as i64 - 1
    }
    pub fn widen(&self, by: i64) -> Self {
        Self { lo: self.lo - by, hi: self.hi + by }
    }
}

/// Re-usable buffers (one per worker thread).
#[derive(Default)]
pub struct BandScratch {
    h_prev: Vec<i32>,
    f_prev: Vec<i32>,
    h_cur: Vec<i32>,
    f_cur: Vec<i32>,
    h16_prev: Vec<i16>,
    f16_prev: Vec<i16>,
    h16_cur: Vec<i16>,
    f16_cur: Vec<i16>,
    prof: Vec<i8>,
    dir: Vec<u8>,
    ops: Vec<u8>,
}

impl BandScratch {
    pub fn new() -> Self {
        Self::default()
    }
    /// Release buffers that grew beyond `keep` bytes (after an outlier).
    pub fn trim(&mut self, keep: usize) {
        if self.dir.capacity() > keep {
            self.dir = Vec::new();
        }
        if self.ops.capacity() > keep {
            self.ops = Vec::new();
        }
        if self.prof.capacity() > keep {
            self.prof = Vec::new();
        }
    }
}

/// Result of a score-only banded pass.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BandScore {
    pub score: i32,
    /// Exclusive end coordinates of the (first, row-major) optimal cell.
    pub q_end: usize,
    pub s_end: usize,
    /// Cells evaluated.
    pub cells: u64,
}

/// A local alignment in SAM-convention operations relative to the subject:
/// `M` aligned pair, `I` query residue vs gap, `D` subject residue vs gap.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct LocalAlignment {
    pub score: i32,
    /// 0-based, half-open.
    pub q_start: usize,
    pub q_end: usize,
    pub s_start: usize,
    pub s_end: usize,
    pub ops: Vec<(u8, u32)>,
    pub cells: u64,
}

impl LocalAlignment {
    pub fn columns(&self) -> usize {
        self.ops.iter().map(|&(_, l)| l as usize).sum()
    }
}

const NEG: i32 = i32::MIN / 4;

#[inline]
fn row_range(i: usize, band: Band, n: usize) -> Option<(usize, usize)> {
    // k ranges over [0, width); j = i + lo + k
    let w = band.width() as i64;
    let i = i as i64;
    let kmin = (-(i + band.lo)).max(0);
    let kmax = (n as i64 - 1 - i - band.lo).min(w - 1);
    if kmin > kmax { None } else { Some((kmin as usize, kmax as usize)) }
}

/// Score-only banded SW over `q × s` (codes 0..=23), portable scalar kernel.
pub fn banded_score_scalar(
    q: &[u8],
    s: &[u8],
    bias: Option<&[i8]>,
    scheme: Scheme<'_>,
    band: Band,
    sc: &mut BandScratch,
) -> BandScore {
    let m = q.len();
    let n = s.len();
    let mut out = BandScore { score: 0, q_end: 0, s_end: 0, cells: 0 };
    if m == 0 || n == 0 || band.hi < band.lo {
        return out;
    }
    let band = band.clip(m, n);
    let w = band.width();
    if w == 0 {
        return out;
    }
    for v in [&mut sc.h_prev, &mut sc.f_prev, &mut sc.h_cur, &mut sc.f_cur] {
        v.clear();
        v.resize(w + 1, NEG);
    }
    let go = scheme.gap_open;
    let ge = scheme.gap_extend;
    let mut cells = 0u64;
    for i in 0..m {
        let Some((kmin, kmax)) = row_range(i, band, n) else {
            for k in 0..w {
                sc.h_cur[k] = NEG;
                sc.f_cur[k] = NEG;
            }
            std::mem::swap(&mut sc.h_prev, &mut sc.h_cur);
            std::mem::swap(&mut sc.f_prev, &mut sc.f_cur);
            continue;
        };
        let row = &scheme.matrix[crate::alphabet::code_index(q[i])];
        let b = bias.map_or(0, |b| b[i] as i32);
        let mut e = NEG;
        let mut h_left = NEG;
        let base = i as i64 + band.lo;
        // Local slices: lets the compiler keep buffer pointers in registers
        // (field accesses through `&mut sc` would be reloaded after every
        // store because the heap writes might alias the struct).
        {
            let hp = &sc.h_prev[..w + 1];
            let fp = &sc.f_prev[..w + 1];
            let hc = &mut sc.h_cur[..w + 1];
            let fc = &mut sc.f_cur[..w + 1];
            let mut best = out.score;
            let mut bq = out.q_end;
            let mut bs = out.s_end;
            for k in kmin..=kmax {
                let j = (base + k as i64) as usize;
                let ecand = (h_left - go).max(e - ge);
                e = if k > kmin { ecand } else { NEG };
                // hp/fp have w+1 entries; index w holds NEG (never written).
                let f = if i >= 1 { (hp[k + 1] - go).max(fp[k + 1] - ge) } else { NEG };
                let hd = if i >= 1 && j >= 1 { hp[k] } else { 0 };
                let diag = hd + row[crate::alphabet::code_index(s[j])] + b;
                let h = diag.max(e).max(f).max(0);
                hc[k] = h;
                fc[k] = f;
                h_left = h;
                if h > best {
                    best = h;
                    bq = i + 1;
                    bs = j + 1;
                }
            }
            out.score = best;
            out.q_end = bq;
            out.s_end = bs;
        }
        cells += (kmax - kmin + 1) as u64;
        // Cells of this row outside [kmin, kmax] must read as unreachable in
        // the next row.
        for k in 0..kmin {
            sc.h_cur[k] = NEG;
            sc.f_cur[k] = NEG;
        }
        for k in kmax + 1..w {
            sc.h_cur[k] = NEG;
            sc.f_cur[k] = NEG;
        }
        std::mem::swap(&mut sc.h_prev, &mut sc.h_cur);
        std::mem::swap(&mut sc.f_prev, &mut sc.f_cur);
    }
    out.cells = cells;
    out
}

/// Portable scalar banded SW with traceback (reference for the SIMD path).
pub fn banded_align_scalar(
    q: &[u8],
    s: &[u8],
    bias: Option<&[i8]>,
    scheme: Scheme<'_>,
    band: Band,
    sc: &mut BandScratch,
    max_dir_bytes: usize,
) -> Result<Option<LocalAlignment>, usize> {
    let m = q.len();
    let n = s.len();
    if m == 0 || n == 0 || band.hi < band.lo {
        return Ok(None);
    }
    let band = band.clip(m, n);
    let w = band.width();
    if w == 0 {
        return Ok(None);
    }
    let need = m.checked_mul(w).unwrap_or(usize::MAX);
    if need > max_dir_bytes {
        return Err(need);
    }
    for v in [&mut sc.h_prev, &mut sc.f_prev, &mut sc.h_cur, &mut sc.f_cur] {
        v.clear();
        v.resize(w + 1, NEG);
    }
    if sc.dir.len() < need {
        sc.dir.resize(need, 0);
    }
    let go = scheme.gap_open;
    let ge = scheme.gap_extend;
    let mut best = (0i32, 0usize, 0usize); // score, i, k
    let mut cells = 0u64;
    for i in 0..m {
        let dbase = i * w;
        let Some((kmin, kmax)) = row_range(i, band, n) else {
            for k in 0..w {
                sc.h_cur[k] = NEG;
                sc.f_cur[k] = NEG;
            }
            std::mem::swap(&mut sc.h_prev, &mut sc.h_cur);
            std::mem::swap(&mut sc.f_prev, &mut sc.f_cur);
            continue;
        };
        let row = &scheme.matrix[crate::alphabet::code_index(q[i])];
        let b = bias.map_or(0, |b| b[i] as i32);
        let mut e = NEG;
        let mut h_left = NEG;
        let base = i as i64 + band.lo;
        // Branch-free cell update: the comparisons become flag bits instead
        // of unpredictable branches (identical results and tie-breaking:
        // diagonal > E > F, extension only when strictly better). Local
        // slices keep buffer pointers in registers.
        {
            let hp = &sc.h_prev[..w + 1];
            let fp = &sc.f_prev[..w + 1];
            let hc = &mut sc.h_cur[..w + 1];
            let fc = &mut sc.f_cur[..w + 1];
            let dirrow = &mut sc.dir[dbase..dbase + w];
            for k in kmin..=kmax {
                let j = (base + k as i64) as usize;
                let eo = h_left - go;
                let ee = e - ge;
                let e_ext = (ee > eo) as u8;
                e = if k > kmin { ee.max(eo) } else { NEG };
                let fo = hp[k + 1] - go;
                let fe = fp[k + 1] - ge;
                let f_ext = (fe > fo) as u8;
                let f = if i >= 1 { fe.max(fo) } else { NEG };
                let hd = if i >= 1 && j >= 1 { hp[k] } else { 0 };
                let diag = hd + row[crate::alphabet::code_index(s[j])] + b;
                let c1 = (diag > 0) as u8;
                let h1 = diag.max(0);
                let c2 = (e > h1) as u8;
                let h2 = h1.max(e);
                let c3 = (f > h2) as u8;
                let h = h2.max(f);
                let org = if c3 == 1 { 3 } else if c2 == 1 { 2 } else { c1 };
                hc[k] = h;
                fc[k] = f;
                dirrow[k] = org | (e_ext << 2) | (f_ext << 3);
                h_left = h;
                if h > best.0 {
                    best = (h, i, k);
                }
            }
        }
        cells += (kmax - kmin + 1) as u64;
        for k in 0..kmin {
            sc.h_cur[k] = NEG;
            sc.f_cur[k] = NEG;
        }
        for k in kmax + 1..w {
            sc.h_cur[k] = NEG;
            sc.f_cur[k] = NEG;
        }
        std::mem::swap(&mut sc.h_prev, &mut sc.h_cur);
        std::mem::swap(&mut sc.f_prev, &mut sc.f_cur);
    }
    if best.0 <= 0 {
        return Ok(None);
    }

    Ok(Some(walk(&sc.dir, &mut sc.ops, w, band, best, cells)))
}


/// Follow the direction bits from the best cell back to the alignment start.
fn walk(dir: &[u8], opbuf: &mut Vec<u8>, w: usize, band: Band, best: (i32, usize, usize), cells: u64) -> LocalAlignment {
    let (score, mut i, mut k) = best;
    let q_end = i + 1;
    let s_end = (i as i64 + band.lo + k as i64) as usize + 1;
    opbuf.clear();
    let mut state = 0u8; // 0 = H, 1 = E (query gap, 'D'), 2 = F (subject gap, 'I')
    loop {
        let d = dir[i * w + k];
        match state {
            0 => match d & 0b11 {
                0 => break,
                1 => {
                    opbuf.push(b'M');
                    let j = i as i64 + band.lo + k as i64;
                    if i == 0 || j == 0 {
                        // aligned pair at the matrix border starts the alignment
                        i = usize::MAX;
                        break;
                    }
                    i -= 1;
                }
                2 => state = 1,
                _ => state = 2,
            },
            1 => {
                opbuf.push(b'D');
                let ext = d & 0b0100 != 0;
                // move left: (i, j-1) -> k-1
                k -= 1;
                if !ext {
                    state = 0;
                }
            }
            _ => {
                opbuf.push(b'I');
                let ext = d & 0b1000 != 0;
                // move up: (i-1, j) -> k+1
                i -= 1;
                k += 1;
                if !ext {
                    state = 0;
                }
            }
        }
    }
    // Starting coordinates: the last emitted op consumed the first residue(s).
    opbuf.reverse();
    let (mut qs, mut ss) = (q_end, s_end);
    for &op in opbuf.iter() {
        match op {
            b'M' => {
                qs -= 1;
                ss -= 1;
            }
            b'I' => qs -= 1,
            _ => ss -= 1,
        }
    }
    let _ = i;
    let mut ops: Vec<(u8, u32)> = Vec::new();
    for &op in opbuf.iter() {
        match ops.last_mut() {
            Some((o, l)) if *o == op => *l += 1,
            _ => ops.push((op, 1)),
        }
    }
    LocalAlignment { score, q_start: qs, q_end, s_start: ss, s_end, ops, cells }
}

/// Score-only banded SW over `q × s` (codes 0..=23). Uses the AVX2 kernels
/// when available (16-bit lanes, or 32-bit lanes when scores could exceed the
/// 16-bit range; identical results), otherwise the scalar kernel.
pub fn banded_score(
    q: &[u8],
    s: &[u8],
    bias: Option<&[i8]>,
    scheme: Scheme<'_>,
    band: Band,
    sc: &mut BandScratch,
) -> BandScore {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            let m = q.len();
            let n = s.len();
            if m == 0 || n == 0 || band.hi < band.lo {
                return BandScore { score: 0, q_end: 0, s_end: 0, cells: 0 };
            }
            let band = band.clip(m, n);
            let w = band.width();
            if w == 0 {
                return BandScore { score: 0, q_end: 0, s_end: 0, cells: 0 };
            }
            let mut nodir: [u8; 0] = [];
            // SAFETY: AVX2 checked at runtime; buffers sized inside.
            let o = unsafe {
                super::banded_simd::banded_avx2_16::<false>(
                    q, s, bias, scheme, band, w, &mut sc.h16_prev, &mut sc.f16_prev, &mut sc.h16_cur, &mut sc.f16_cur, &mut sc.prof, &mut nodir, i32::MAX,
                )
                .unwrap_or_else(|| {
                    super::banded_simd::banded_avx2::<false>(
                        q, s, bias, scheme, band, w, &mut sc.h_prev, &mut sc.f_prev, &mut sc.h_cur, &mut sc.f_cur, &mut nodir, i32::MAX,
                    )
                })
            };
            if o.score <= 0 {
                return BandScore { score: 0, q_end: 0, s_end: 0, cells: o.cells };
            }
            let j = (o.best_i as i64 + band.lo + o.best_k as i64) as usize;
            return BandScore { score: o.score, q_end: o.best_i + 1, s_end: j + 1, cells: o.cells };
        }
    }
    banded_score_scalar(q, s, bias, scheme, band, sc)
}

/// Banded SW with traceback. Returns `Ok(None)` when the best score is 0 and
/// `Err(bytes)` when the direction matrix would exceed `max_dir_bytes` (the
/// caller must then report the fallback explicitly). AVX2 when available.
pub fn banded_align(
    q: &[u8],
    s: &[u8],
    bias: Option<&[i8]>,
    scheme: Scheme<'_>,
    band: Band,
    sc: &mut BandScratch,
    max_dir_bytes: usize,
) -> Result<Option<LocalAlignment>, usize> {
    banded_align_until(q, s, bias, scheme, band, sc, max_dir_bytes, i32::MAX)
}

/// [`banded_align`] for a caller that knows an upper bound `stop_at` on the
/// score (e.g. the exact score of the enclosing window): rows after the first
/// one whose best cell reaches `stop_at` are not computed — no later cell can
/// score higher and ties keep the first cell, so the alignment is identical;
/// only the reported `cells` count is smaller.
#[allow(clippy::too_many_arguments)]
pub fn banded_align_until(
    q: &[u8],
    s: &[u8],
    bias: Option<&[i8]>,
    scheme: Scheme<'_>,
    band: Band,
    sc: &mut BandScratch,
    max_dir_bytes: usize,
    stop_at: i32,
) -> Result<Option<LocalAlignment>, usize> {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            let m = q.len();
            let n = s.len();
            if m == 0 || n == 0 || band.hi < band.lo {
                return Ok(None);
            }
            let band = band.clip(m, n);
            let w = band.width();
            if w == 0 {
                return Ok(None);
            }
            let need = m.checked_mul(w).unwrap_or(usize::MAX);
            if need > max_dir_bytes {
                return Err(need);
            }
            if sc.dir.len() < need + 16 {
                sc.dir.resize(need + 16, 0);
            }
            // SAFETY: AVX2 checked at runtime; dir sized above.
            let o = unsafe {
                match super::banded_simd::banded_avx2_16::<true>(
                    q, s, bias, scheme, band, w, &mut sc.h16_prev, &mut sc.f16_prev, &mut sc.h16_cur, &mut sc.f16_cur, &mut sc.prof, &mut sc.dir, stop_at,
                ) {
                    Some(o) => o,
                    None => super::banded_simd::banded_avx2::<true>(
                        q, s, bias, scheme, band, w, &mut sc.h_prev, &mut sc.f_prev, &mut sc.h_cur, &mut sc.f_cur, &mut sc.dir, stop_at,
                    ),
                }
            };
            if o.score <= 0 {
                return Ok(None);
            }
            return Ok(Some(walk(&sc.dir, &mut sc.ops, w, band, (o.score, o.best_i, o.best_k), o.cells)));
        }
    }
    let _ = stop_at;
    banded_align_scalar(q, s, bias, scheme, band, sc, max_dir_bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::reference::{score_alignment, sw_reference};
    use crate::alphabet::EncodedProtein;
    use crate::scoring::{ScoreMatrix, ScoringSystem};

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self, m: u64) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 33) % m
        }
    }

    fn rand_seq(r: &mut Lcg, len: usize) -> Vec<u8> {
        (0..len).map(|_| r.next(24) as u8).collect()
    }

    fn mutate(r: &mut Lcg, s: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        let mut i = 0;
        while i < s.len() {
            let x = r.next(100);
            if x < 6 {
                // insertion of 1-12 residues
                for _ in 0..1 + r.next(12) {
                    out.push(r.next(20) as u8);
                }
            } else if x < 12 {
                i += 1 + r.next(10) as usize; // deletion
                continue;
            } else if x < 45 {
                out.push(r.next(20) as u8);
            } else {
                out.push(s[i]);
            }
            i += 1;
        }
        out
    }

    #[test]
    fn full_band_equals_reference_on_random_and_related_pairs() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let scheme = Scheme { matrix: &m.matrix, gap_open: 11, gap_extend: 1 };
        let mut r = Lcg(42);
        let mut sc = BandScratch::new();
        for t in 0..300 {
            let la = 1 + r.next(160) as usize;
            let a = rand_seq(&mut r, la);
            let lb = 1 + r.next(160) as usize;
            let b = if t % 2 == 0 { mutate(&mut r, &a) } else { rand_seq(&mut r, lb) };
            let bias: Vec<i8> = (0..a.len()).map(|_| r.next(7) as i8 - 3).collect();
            let biasopt = if t % 3 == 0 { Some(&bias[..]) } else { None };
            let (rs, ri, rj) = sw_reference(&a, &b, &m.matrix, biasopt, 11, 1);
            let full = Band::full(a.len(), b.len());
            let bs = banded_score(&a, &b, biasopt, scheme, full, &mut sc);
            assert_eq!(bs.score as i64, rs, "score t={t}");
            if rs > 0 {
                assert_eq!((bs.q_end, bs.s_end), (ri, rj), "end cell t={t}");
            }
            let al = banded_align(&a, &b, biasopt, scheme, full, &mut sc, usize::MAX).unwrap();
            if rs == 0 {
                assert!(al.is_none());
                continue;
            }
            let al = al.unwrap();
            assert_eq!(al.score as i64, rs);
            let (rescored, qe, se) = score_alignment(&a, &b, al.q_start, al.s_start, &al.ops, &m.matrix, biasopt, 11, 1).unwrap();
            assert_eq!(rescored, rs, "rescored CIGAR t={t} ops={:?}", al.ops);
            assert_eq!((qe, se), (al.q_end, al.s_end));
        }
    }

    #[test]
    fn band_never_exceeds_full_and_is_exact_when_path_inside() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let scheme = Scheme { matrix: &m.matrix, gap_open: 11, gap_extend: 1 };
        let mut r = Lcg(7);
        let mut sc = BandScratch::new();
        for _ in 0..200 {
            let la = 50 + r.next(200) as usize;
            let a = rand_seq(&mut r, la);
            let b = mutate(&mut r, &a);
            let (rs, _, _) = sw_reference(&a, &b, &m.matrix, None, 11, 1);
            for half in [2i64, 8, 32] {
                let bs = banded_score(&a, &b, None, scheme, Band { lo: -half, hi: half }, &mut sc);
                assert!(bs.score as i64 <= rs);
                let al = banded_align(&a, &b, None, scheme, Band { lo: -half, hi: half }, &mut sc, usize::MAX).unwrap();
                if let Some(al) = al {
                    assert_eq!(al.score, bs.score);
                    let (rescored, _, _) = score_alignment(&a, &b, al.q_start, al.s_start, &al.ops, &m.matrix, None, 11, 1).unwrap();
                    assert_eq!(rescored, al.score as i64);
                    // every column stays inside the band
                    let (mut i, mut j) = (al.q_start as i64, al.s_start as i64);
                    for &(op, l) in &al.ops {
                        for _ in 0..l {
                            assert!(j - i >= -half && j - i <= half);
                            match op { b'M' => { i += 1; j += 1 } b'I' => i += 1, _ => j += 1 }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn shifted_diagonal_and_long_gap_need_band_on_that_diagonal() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let scheme = Scheme { matrix: &m.matrix, gap_open: 11, gap_extend: 1 };
        let core = EncodedProtein::from_ascii(b"MKVLAAGIVGHHWWPQRSTVYCCMKDEFGHIKLMNPQRSTVWY");
        let mut subj = vec![0u8; 300];
        subj.extend_from_slice(core.as_slice());
        let (rs, _, _) = sw_reference(core.as_slice(), &subj, &m.matrix, None, 11, 1);
        let mut sc = BandScratch::new();
        let near = banded_score(core.as_slice(), &subj, None, scheme, Band { lo: -10, hi: 10 }, &mut sc);
        assert!(near.score < rs as i32 / 2);
        let on = banded_score(core.as_slice(), &subj, None, scheme, Band { lo: 290, hi: 310 }, &mut sc);
        assert_eq!(on.score as i64, rs);
    }

    #[test]
    fn simd_and_scalar_kernels_agree() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let mut r = Lcg(2024);
        let mut sc1 = BandScratch::new();
        let mut sc2 = BandScratch::new();
        for t in 0..600 {
            let (go, ge) = [(11, 1), (10, 2), (6, 3), (15, 1), (8, 8)][t % 5];
            let scheme = Scheme { matrix: &m.matrix, gap_open: go, gap_extend: ge };
            let la = 1 + r.next(if t % 7 == 0 { 700 } else { 220 }) as usize;
            let a = rand_seq(&mut r, la);
            let lb = 1 + r.next(220) as usize;
            let b = if t % 2 == 0 { mutate(&mut r, &a) } else { rand_seq(&mut r, lb) };
            let bias: Vec<i8> = (0..a.len()).map(|_| r.next(9) as i8 - 4).collect();
            let bo = if t % 3 == 0 { Some(&bias[..]) } else { None };
            let lo = -(r.next(60) as i64);
            let hi = r.next(60) as i64 + if t % 5 == 0 { 300 } else { 0 };
            for band in [Band { lo, hi }, Band::full(a.len(), b.len()), Band { lo: 5, hi: 40 }] {
                let x = banded_score(&a, &b, bo, scheme, band, &mut sc1);
                let y = banded_score_scalar(&a, &b, bo, scheme, band, &mut sc2);
                assert_eq!((x.score, x.q_end, x.s_end, x.cells), (y.score, y.q_end, y.s_end, y.cells), "score t={t} band={band:?}");
                let p = banded_align(&a, &b, bo, scheme, band, &mut sc1, usize::MAX).unwrap();
                let q = banded_align_scalar(&a, &b, bo, scheme, band, &mut sc2, usize::MAX).unwrap();
                assert_eq!(p, q, "align t={t} band={band:?}");
            }
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn i16_and_i32_simd_kernels_agree_and_overflow_falls_back() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        use crate::algorithm::banded_simd::{banded_avx2, banded_avx2_16};
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let mut r = Lcg(99);
        let (mut h1, mut f1, mut h2, mut f2) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let mut pr = Vec::new();
        let (mut g1, mut e1, mut g2, mut e2) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for t in 0..400 {
            let (go, ge) = [(11, 1), (10, 2), (6, 3), (15, 1), (8, 8)][t % 5];
            let scheme = Scheme { matrix: &m.matrix, gap_open: go, gap_extend: ge };
            let la = 1 + r.next(400) as usize;
            let a = rand_seq(&mut r, la);
            let lb = 1 + r.next(400) as usize;
            let b = if t % 2 == 0 { mutate(&mut r, &a) } else { rand_seq(&mut r, lb) };
            let bias: Vec<i8> = (0..a.len()).map(|_| r.next(9) as i8 - 4).collect();
            let bo = if t % 3 == 0 { Some(&bias[..]) } else { None };
            let lo = -(r.next(80) as i64) - if t % 4 == 0 { 200 } else { 0 };
            let hi = r.next(80) as i64;
            let band = Band { lo, hi }.clip(a.len(), b.len());
            if band.hi < band.lo {
                continue;
            }
            let w = band.width();
            let mut d16 = vec![0u8; a.len() * w + 16];
            let mut d32 = vec![0u8; a.len() * w + 16];
            let x = unsafe { banded_avx2_16::<true>(&a, &b, bo, scheme, band, w, &mut h1, &mut f1, &mut h2, &mut f2, &mut pr, &mut d16, i32::MAX) }
                .expect("scores far below the 16-bit limit");
            let y = unsafe { banded_avx2::<true>(&a, &b, bo, scheme, band, w, &mut g1, &mut e1, &mut g2, &mut e2, &mut d32, i32::MAX) };
            assert_eq!((x.score, x.best_i, x.best_k, x.cells), (y.score, y.best_i, y.best_k, y.cells), "t={t}");
            for i in 0..a.len() {
                if let Some((kmin, kmax)) = row_range(i, band, b.len()) {
                    assert_eq!(&d16[i * w + kmin..=i * w + kmax], &d32[i * w + kmin..=i * w + kmax], "dir row {i} t={t}");
                }
            }
            let mut nod: [u8; 0] = [];
            let xs = unsafe { banded_avx2_16::<false>(&a, &b, bo, scheme, band, w, &mut h1, &mut f1, &mut h2, &mut f2, &mut pr, &mut nod, i32::MAX) }.unwrap();
            assert_eq!((xs.score, xs.best_i, xs.best_k, xs.cells), (y.score, y.best_i, y.best_k, y.cells));
        }
        // Scores above the 16-bit range: the 16-bit kernel declines and the
        // dispatcher's 32-bit result equals the scalar kernel.
        let scheme = Scheme { matrix: &m.matrix, gap_open: 11, gap_extend: 1 };
        let a = vec![17u8; 4000]; // W x 4000: score 44 000
        let band = Band { lo: -20, hi: 20 };
        let w = band.width();
        let mut d = vec![0u8; a.len() * w + 16];
        assert!(unsafe { banded_avx2_16::<true>(&a, &a, None, scheme, band, w, &mut h1, &mut f1, &mut h2, &mut f2, &mut pr, &mut d, i32::MAX) }.is_none());
        let mut sc1 = BandScratch::new();
        let mut sc2 = BandScratch::new();
        let p = banded_align(&a, &a, None, scheme, band, &mut sc1, usize::MAX).unwrap().unwrap();
        let q = banded_align_scalar(&a, &a, None, scheme, band, &mut sc2, usize::MAX).unwrap().unwrap();
        assert_eq!(p.score, 44_000);
        assert_eq!(p, q);
        assert_eq!(banded_score(&a, &a, None, scheme, band, &mut sc1), banded_score_scalar(&a, &a, None, scheme, band, &mut sc2));
    }

    #[test]
    fn stopping_at_the_known_score_gives_the_same_alignment() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let scheme = Scheme { matrix: &m.matrix, gap_open: 11, gap_extend: 1 };
        let mut r = Lcg(5);
        let mut sc = BandScratch::new();
        let mut saved = 0u64;
        for _ in 0..200 {
            let la = 20 + r.next(300) as usize;
            let a = rand_seq(&mut r, la);
            let mut b = rand_seq(&mut r, 40);
            b.extend(mutate(&mut r, &a[..la / 2]));
            b.extend(rand_seq(&mut r, 40));
            let band = Band { lo: 0, hi: 90 };
            let full = banded_align(&a, &b, None, scheme, band, &mut sc, usize::MAX).unwrap();
            let Some(full) = full else { continue };
            let early = banded_align_until(&a, &b, None, scheme, band, &mut sc, usize::MAX, full.score).unwrap().unwrap();
            assert_eq!((early.score, early.q_start, early.q_end, early.s_start, early.s_end, &early.ops), (full.score, full.q_start, full.q_end, full.s_start, full.s_end, &full.ops));
            assert!(early.cells <= full.cells);
            saved += full.cells - early.cells;
        }
        assert!(saved > 0, "alignments ending before the last row stop early");
    }

    #[test]
    fn traceback_respects_memory_cap() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let scheme = Scheme { matrix: &m.matrix, gap_open: 11, gap_extend: 1 };
        let a = vec![17u8; 1000];
        let mut sc = BandScratch::new();
        let r = banded_align(&a, &a, None, scheme, Band { lo: -500, hi: 500 }, &mut sc, 1000);
        assert!(r.is_err());
    }
}
