//! Low-complexity / tandem-repeat masking for *seeding*.
//!
//! SABER implements the tantan model of Frith (2011), "A new repeat-masking
//! method enables specific detection of homologous sequences", Nucleic Acids
//! Research 39:e23 — the masker DIAMOND applies by default (`--masking
//! tantan`). Masked residues are excluded from seed (k-mer) formation on both
//! the query and the database side; the original residues are always kept for
//! ungapped extension, gapped alignment, traceback and output, so masking can
//! only remove *seeds*, never change a reported score.
//!
//! ## Model
//!
//! A hidden Markov model with one background state `B` and repeat states
//! `F_1 … F_W`, where `F_i` emits a letter related to the letter `i`
//! positions upstream:
//!
//! ```text
//! B → B   1 − r            B → F_i  r · f_i,  f_i = (1−d)/(1−d^W) · d^(i−1)
//! F_i → F_i  1 − e         F_i → B  e
//! emission ratio of F_i at position p:  exp(λ · S(x_p, x_{p−i}))
//! ```
//!
//! The posterior probability of being in any `F_i` is computed with the
//! scaled forward–backward algorithm and positions with probability ≥ 0.5 are
//! masked. Defaults are tantan's protein defaults (`tantan -p`): r = 0.005,
//! e = 0.05, W = 50, d = 0.9, BLOSUM62, threshold 0.5, no gaps.
//! `λ` is the scale implied by the matrix itself (the root of
//! Σ (e^{λS})⁻¹ = 1, Yu, Wootton & Altschul 2003), which is how tantan/LAST
//! derive it. `bench/scripts/validate_masking.sh` compares SABER's posteriors
//! with the tantan 49 binary (`bench/results/masking_vs_tantan.txt`).

use crate::scoring::matrices::ScoreMatrix;

/// Parameters of the tantan model.
#[derive(Debug, Clone, Copy)]
pub struct TantanParams {
    pub repeat_start: f64,
    pub repeat_end: f64,
    pub max_period: usize,
    pub decay: f64,
    pub min_prob: f64,
}

impl Default for TantanParams {
    fn default() -> Self {
        Self { repeat_start: 0.005, repeat_end: 0.05, max_period: 50, decay: 0.9, min_prob: 0.5 }
    }
}

/// Repeat states of the default model (tantan's protein default), for which
/// the AVX2 path is written.
const W50: usize = 50;

/// Which masking algorithm to apply to seeds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaskingMode {
    None,
    Tantan,
}

impl std::str::FromStr for MaskingMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "none" | "0" | "off" | "no" => Ok(Self::None),
            "tantan" | "1" | "on" | "yes" => Ok(Self::Tantan),
            o => Err(format!("unknown masking mode '{o}' (expected tantan|none)")),
        }
    }
}

impl std::fmt::Display for MaskingMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::None => f.write_str("none"),
            Self::Tantan => f.write_str("tantan"),
        }
    }
}

/// Reusable tantan masker for protein sequences in alignment encoding.
///
/// The forward–backward recursions run in single precision (eight lanes
/// per AVX2 vector, half the memory traffic of the per-position tables).
/// Every forward row is renormalised, so values stay in [0, 1] scale and
/// rounding moves posteriors by ~1e-6; agreement with the reference tantan
/// binary is reported in `bench/results/v1.2.0/masking_vs_tantan.txt`.
/// With AVX2 and the default 50 repeat states, [`Tantan::compute_avx2_w50`]
/// computes the same posteriors bit for bit with explicit vectors.
pub struct Tantan {
    p: TantanParams,
    /// Likelihood ratios exp(λ S(a,b)) over the 24-letter alphabet; columns
    /// 24..32 are 0 (the ratio of a position without a partner residue,
    /// looked up by the AVX2 path).
    lr: [[f32; 32]; 24],
    lambda: f64,
    /// b2f_i for i = 1..=W.
    b2f: Vec<f32>,
}

/// Scratch buffers for one masking call (re-used across sequences).
#[derive(Default)]
pub struct TantanScratch {
    fwd_scale: Vec<f32>,
    fwd: Vec<f32>, // (n) × (W+1) row-major: [B, F_1..F_W]
    post: Vec<f32>,
    codes: Vec<u8>,
    emit: Vec<f32>, // (n) × W emission ratios, shared by both passes
    /// AVX2 path: the codes in reverse order followed by code 24 (no
    /// partner), forward B states, forward F-state rows of 56 floats (states
    /// 50..56 zero) and the padded b2f / zero / beta vectors.
    rev: Vec<u8>,
    fwd_b: Vec<f32>,
    fwd_f: Vec<f32>,
    work: Vec<f32>,
}

/// Sum with eight independent accumulators (one AVX2 vector; the result
/// differs from a sequential sum only by rounding).
#[inline(always)]
fn sum8(v: &[f32]) -> f32 {
    let mut a = [0.0f32; 8];
    let mut it = v.chunks_exact(8);
    for c in &mut it {
        for k in 0..8 {
            a[k] += c[k];
        }
    }
    let mut t = ((a[0] + a[1]) + (a[2] + a[3])) + ((a[4] + a[5]) + (a[6] + a[7]));
    for &x in it.remainder() {
        t += x;
    }
    t
}

#[inline(always)]
fn dot8(a: &[f32], b: &[f32]) -> f32 {
    let mut acc = [0.0f32; 8];
    let n = a.len().min(b.len());
    let (a, b) = (&a[..n], &b[..n]);
    let mut ia = a.chunks_exact(8);
    let mut ib = b.chunks_exact(8);
    for (x, y) in (&mut ia).zip(&mut ib) {
        for k in 0..8 {
            acc[k] += x[k] * y[k];
        }
    }
    let mut t = ((acc[0] + acc[1]) + (acc[2] + acc[3])) + ((acc[4] + acc[5]) + (acc[6] + acc[7]));
    for (x, y) in ia.remainder().iter().zip(ib.remainder()) {
        t += x * y;
    }
    t
}

impl Tantan {
    /// Masker with tantan's protein defaults for the given matrix.
    pub fn protein(m: &ScoreMatrix) -> Self {
        Self::with_params(m, TantanParams::default(), None)
    }

    /// Fully specified masker; `lambda = None` derives λ from the matrix.
    pub fn with_params(m: &ScoreMatrix, p: TantanParams, lambda: Option<f64>) -> Self {
        let lambda = lambda.unwrap_or_else(|| implied_lambda(m).unwrap_or(0.3176));
        let mut lr = [[0.0f32; 32]; 24];
        for a in 0..24 {
            for b in 0..24 {
                lr[a][b] = (lambda * m.matrix[a][b] as f64).exp() as f32;
            }
        }
        let w = p.max_period.max(1);
        let first = if (p.decay - 1.0).abs() < 1e-12 {
            1.0 / w as f64
        } else {
            (1.0 - p.decay) / (1.0 - p.decay.powi(w as i32))
        };
        let b2f: Vec<f32> = (0..w).map(|i| (p.repeat_start * first * p.decay.powi(i as i32)) as f32).collect();
        Self { p, lr, lambda, b2f }
    }

    pub fn lambda(&self) -> f64 {
        self.lambda
    }

    /// Posterior repeat probability for every position of `seq` (codes 0..=23).
    pub fn repeat_probs(&self, seq: &[u8], sc: &mut TantanScratch) -> Vec<f32> {
        self.compute(seq, sc);
        sc.post.clone()
    }

    /// Set `mask[i] = true` for masked positions (mask is resized to seq.len()).
    pub fn mask(&self, seq: &[u8], sc: &mut TantanScratch, mask: &mut Vec<bool>) -> usize {
        self.compute(seq, sc);
        mask.clear();
        mask.extend(sc.post.iter().map(|&pr| pr as f64 >= self.p.min_prob));
        mask.iter().filter(|&&m| m).count()
    }

    fn compute(&self, seq: &[u8], sc: &mut TantanScratch) {
        #[cfg(target_arch = "x86_64")]
        {
            if self.b2f.len() == W50 && is_x86_feature_detected!("avx2") {
                // SAFETY: AVX2 detected; 50 repeat states.
                unsafe { self.compute_avx2_w50(seq, sc) };
                return;
            }
            if is_x86_feature_detected!("avx2") {
                // SAFETY: AVX2 detected. Same arithmetic (no FMA), wider registers.
                unsafe { self.compute_avx2(seq, sc) };
                return;
            }
        }
        self.compute_body(seq, sc);
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn compute_avx2(&self, seq: &[u8], sc: &mut TantanScratch) {
        self.compute_body(seq, sc);
    }

    /// [`Tantan::compute_body`] for the default 50 repeat states with
    /// explicit AVX2 vectors (about 2.5× faster). Every value is computed by
    /// the same operations in the same order — the sums keep `sum8`/`dot8`'s
    /// eight accumulators, their combination tree and the sequential
    /// remainder (states 48, 49); no FMA — so the posteriors are
    /// bit-identical (`tests::avx2_masker_is_bit_identical_to_the_scalar_body`).
    /// Only the layout differs: the F states of a position are a row of seven
    /// vectors (states 50..56 stay zero), the B state is stored apart, and
    /// the emission ratios lr[x_p][x_{p-1-i}] are looked up in the forward
    /// step (in-register permutes over the 32-entry row of x_p, indexed by
    /// the preceding codes in reverse order) and stored for the backward pass.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn compute_avx2_w50(&self, seq: &[u8], sc: &mut TantanScratch) {
        use std::arch::x86_64::*;
        const WS: usize = 56;
        debug_assert_eq!(self.b2f.len(), W50);
        let n = seq.len();
        sc.post.clear();
        sc.post.resize(n, 0.0);
        if n == 0 {
            return;
        }
        let b2b = (1.0 - self.p.repeat_start) as f32;
        let f2f = (1.0 - self.p.repeat_end) as f32;
        let f2b = self.p.repeat_end as f32;
        sc.codes.clear();
        sc.codes.extend(seq.iter().map(|&c| (c as usize).min(23) as u8));
        // rev[n - p + i] = x_{p-1-i} for i < p; code 24 (ratio 0) for i >= p
        sc.rev.clear();
        sc.rev.extend(sc.codes.iter().rev());
        sc.rev.resize(n + WS + 8, 24);
        // work: b2f | zeros (initial forward row) | beta | next beta
        sc.work.clear();
        sc.work.resize(4 * WS, 0.0);
        sc.work[..W50].copy_from_slice(&self.b2f);
        for j in 0..W50 {
            sc.work[2 * WS + j] = f2b;
        }
        if sc.emit.len() < n * WS {
            sc.emit.resize(n * WS, 0.0);
        }
        if sc.fwd_f.len() < n * WS {
            sc.fwd_f.resize(n * WS, 0.0);
        }
        if sc.fwd_b.len() < n {
            sc.fwd_b.resize(n, 0.0);
        }
        if sc.fwd_scale.len() < n {
            sc.fwd_scale.resize(n, 1.0);
        }
        let codes = sc.codes.as_ptr();
        let rev = sc.rev.as_ptr();
        let emit = sc.emit.as_mut_ptr();
        let fwd_f = sc.fwd_f.as_mut_ptr();
        let fwd_b = sc.fwd_b.as_mut_ptr();
        let scale = sc.fwd_scale.as_mut_ptr();
        let post = sc.post.as_mut_ptr();
        let b2f = sc.work.as_ptr();
        let lrp = self.lr.as_ptr() as *const f32;
        let f2f_v = _mm256_set1_ps(f2f);
        let keep2 = _mm256_castsi256_ps(_mm256_setr_epi32(-1, -1, 0, 0, 0, 0, 0, 0));

        // ((a0+a1)+(a2+a3))+((a4+a5)+(a6+a7)), as sum8/dot8 combine their accumulators
        macro_rules! tree {
            ($acc:expr) => {{
                let acc = $acc;
                let q = _mm_hadd_ps(_mm256_castps256_ps128(acc), _mm256_extractf128_ps::<1>(acc));
                let r = _mm_hadd_ps(q, q);
                _mm_cvtss_f32(_mm_add_ss(r, _mm_movehdup_ps(r)))
            }};
        }

        // ---- forward ----
        let mut prev_b = 1.0f32;
        let mut prev_f: *const f32 = b2f.add(WS); // zeros
        for p in 0..n {
            // emission ratios of position p
            let e = emit.add(p * WS);
            let lrow = lrp.add(32 * (*codes.add(p) as usize));
            let (r0, r1, r2, r3) = (_mm256_loadu_ps(lrow), _mm256_loadu_ps(lrow.add(8)), _mm256_loadu_ps(lrow.add(16)), _mm256_loadu_ps(lrow.add(24)));
            let src = rev.add(n - p);
            macro_rules! em {
                ($c:expr) => {{
                    let idx = _mm256_cvtepu8_epi32(_mm_loadl_epi64(src.add(8 * $c) as *const __m128i));
                    // bit 3 of the code selects within 0..16 / 16..32, bit 4 between them
                    let s3 = _mm256_castsi256_ps(_mm256_slli_epi32::<28>(idx));
                    let s4 = _mm256_castsi256_ps(_mm256_slli_epi32::<27>(idx));
                    let ab = _mm256_blendv_ps(_mm256_permutevar8x32_ps(r0, idx), _mm256_permutevar8x32_ps(r1, idx), s3);
                    let cd = _mm256_blendv_ps(_mm256_permutevar8x32_ps(r2, idx), _mm256_permutevar8x32_ps(r3, idx), s3);
                    _mm256_blendv_ps(ab, cd, s4)
                }};
            }
            let (m0, m1, m2, m3, m4, m5) = (em!(0), em!(1), em!(2), em!(3), em!(4), em!(5));
            let m6 = _mm256_and_ps(em!(6), keep2);
            _mm256_storeu_ps(e, m0);
            _mm256_storeu_ps(e.add(8), m1);
            _mm256_storeu_ps(e.add(16), m2);
            _mm256_storeu_ps(e.add(24), m3);
            _mm256_storeu_ps(e.add(32), m4);
            _mm256_storeu_ps(e.add(40), m5);
            _mm256_storeu_ps(e.add(48), m6);
            // sum_f = sum8(prev_f)
            let mut acc = _mm256_add_ps(_mm256_setzero_ps(), _mm256_loadu_ps(prev_f));
            acc = _mm256_add_ps(acc, _mm256_loadu_ps(prev_f.add(8)));
            acc = _mm256_add_ps(acc, _mm256_loadu_ps(prev_f.add(16)));
            acc = _mm256_add_ps(acc, _mm256_loadu_ps(prev_f.add(24)));
            acc = _mm256_add_ps(acc, _mm256_loadu_ps(prev_f.add(32)));
            acc = _mm256_add_ps(acc, _mm256_loadu_ps(prev_f.add(40)));
            let mut sum_f = tree!(acc);
            sum_f += *prev_f.add(48);
            sum_f += *prev_f.add(49);
            let cur_b = prev_b * b2b + sum_f * f2b;
            // cur_f[i] = (prev_b * b2f[i] + prev_f[i] * f2f) * emit[i]
            let pb = _mm256_set1_ps(prev_b);
            macro_rules! cur {
                ($c:expr, $m:expr) => {
                    _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(pb, _mm256_loadu_ps(b2f.add(8 * $c))), _mm256_mul_ps(_mm256_loadu_ps(prev_f.add(8 * $c)), f2f_v)), $m)
                };
            }
            let (c0, c1, c2, c3, c4, c5, c6) = (cur!(0, m0), cur!(1, m1), cur!(2, m2), cur!(3, m3), cur!(4, m4), cur!(5, m5), cur!(6, m6));
            // tot = cur_b + sum8(cur_f)
            let mut acc = _mm256_add_ps(_mm256_setzero_ps(), c0);
            acc = _mm256_add_ps(acc, c1);
            acc = _mm256_add_ps(acc, c2);
            acc = _mm256_add_ps(acc, c3);
            acc = _mm256_add_ps(acc, c4);
            acc = _mm256_add_ps(acc, c5);
            let mut tsum = tree!(acc);
            let last = _mm256_castps256_ps128(c6);
            tsum += _mm_cvtss_f32(last);
            tsum += _mm_cvtss_f32(_mm_movehdup_ps(last));
            let tot = cur_b + tsum;
            let s = if tot > 0.0 { 1.0 / tot } else { 1.0 };
            *scale.add(p) = s;
            let rb = cur_b * s;
            *fwd_b.add(p) = rb;
            let sv = _mm256_set1_ps(s);
            let row = fwd_f.add(p * WS);
            _mm256_storeu_ps(row, _mm256_mul_ps(c0, sv));
            _mm256_storeu_ps(row.add(8), _mm256_mul_ps(c1, sv));
            _mm256_storeu_ps(row.add(16), _mm256_mul_ps(c2, sv));
            _mm256_storeu_ps(row.add(24), _mm256_mul_ps(c3, sv));
            _mm256_storeu_ps(row.add(32), _mm256_mul_ps(c4, sv));
            _mm256_storeu_ps(row.add(40), _mm256_mul_ps(c5, sv));
            _mm256_storeu_ps(row.add(48), _mm256_mul_ps(c6, sv));
            prev_b = rb;
            prev_f = row;
        }

        // ---- backward (same scale factors) ----
        let mut nb = 1.0f32;
        let mut nf = sc.work.as_mut_ptr().add(2 * WS);
        let mut nnf = sc.work.as_mut_ptr().add(3 * WS);
        for p in (0..n).rev() {
            let row = fwd_f.add(p * WS) as *const f32;
            let pb = *fwd_b.add(p) * nb;
            // pf = dot8(row, nf)
            let mut acc = _mm256_add_ps(_mm256_setzero_ps(), _mm256_mul_ps(_mm256_loadu_ps(row), _mm256_loadu_ps(nf)));
            acc = _mm256_add_ps(acc, _mm256_mul_ps(_mm256_loadu_ps(row.add(8)), _mm256_loadu_ps(nf.add(8))));
            acc = _mm256_add_ps(acc, _mm256_mul_ps(_mm256_loadu_ps(row.add(16)), _mm256_loadu_ps(nf.add(16))));
            acc = _mm256_add_ps(acc, _mm256_mul_ps(_mm256_loadu_ps(row.add(24)), _mm256_loadu_ps(nf.add(24))));
            acc = _mm256_add_ps(acc, _mm256_mul_ps(_mm256_loadu_ps(row.add(32)), _mm256_loadu_ps(nf.add(32))));
            acc = _mm256_add_ps(acc, _mm256_mul_ps(_mm256_loadu_ps(row.add(40)), _mm256_loadu_ps(nf.add(40))));
            let mut pf = tree!(acc);
            pf += *row.add(48) * *nf.add(48);
            pf += *row.add(49) * *nf.add(49);
            let z = pb + pf;
            *post.add(p) = if z > 0.0 { pf / z } else { 0.0 };
            if p == 0 {
                break;
            }
            let e = emit.add(p * WS) as *const f32;
            let s = *scale.add(p);
            // be[i] = emit[i] * nf[i]
            macro_rules! be {
                ($c:expr) => {
                    _mm256_mul_ps(_mm256_loadu_ps(e.add(8 * $c)), _mm256_loadu_ps(nf.add(8 * $c)))
                };
            }
            let (e0, e1, e2, e3, e4, e5, e6) = (be!(0), be!(1), be!(2), be!(3), be!(4), be!(5), be!(6));
            // to_fg = dot8(b2f, be)
            let mut acc = _mm256_add_ps(_mm256_setzero_ps(), _mm256_mul_ps(_mm256_loadu_ps(b2f), e0));
            acc = _mm256_add_ps(acc, _mm256_mul_ps(_mm256_loadu_ps(b2f.add(8)), e1));
            acc = _mm256_add_ps(acc, _mm256_mul_ps(_mm256_loadu_ps(b2f.add(16)), e2));
            acc = _mm256_add_ps(acc, _mm256_mul_ps(_mm256_loadu_ps(b2f.add(24)), e3));
            acc = _mm256_add_ps(acc, _mm256_mul_ps(_mm256_loadu_ps(b2f.add(32)), e4));
            acc = _mm256_add_ps(acc, _mm256_mul_ps(_mm256_loadu_ps(b2f.add(40)), e5));
            let mut to_fg = tree!(acc);
            let last = _mm256_castps256_ps128(e6);
            to_fg += *b2f.add(48) * _mm_cvtss_f32(last);
            to_fg += *b2f.add(49) * _mm_cvtss_f32(_mm_movehdup_ps(last));
            // nnf[i] = (f2f * be[i] + f2b * nb) * s
            let fbv = _mm256_set1_ps(f2b * nb);
            let sv = _mm256_set1_ps(s);
            _mm256_storeu_ps(nnf, _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(f2f_v, e0), fbv), sv));
            _mm256_storeu_ps(nnf.add(8), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(f2f_v, e1), fbv), sv));
            _mm256_storeu_ps(nnf.add(16), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(f2f_v, e2), fbv), sv));
            _mm256_storeu_ps(nnf.add(24), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(f2f_v, e3), fbv), sv));
            _mm256_storeu_ps(nnf.add(32), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(f2f_v, e4), fbv), sv));
            _mm256_storeu_ps(nnf.add(40), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(f2f_v, e5), fbv), sv));
            _mm256_storeu_ps(nnf.add(48), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(f2f_v, e6), fbv), sv));
            nb = (b2b * nb + to_fg) * s;
            std::mem::swap(&mut nf, &mut nnf);
        }
    }

    #[inline(always)]
    fn compute_body(&self, seq: &[u8], sc: &mut TantanScratch) {
        let n = seq.len();
        let w = self.b2f.len();
        let ns = w + 1;
        sc.post.clear();
        sc.post.resize(n, 0.0);
        if n == 0 {
            return;
        }
        let b2b = (1.0 - self.p.repeat_start) as f32;
        let f2f = (1.0 - self.p.repeat_end) as f32;
        let f2b = self.p.repeat_end as f32;
        sc.codes.clear();
        sc.codes.extend(seq.iter().map(|&c| (c as usize).min(23) as u8));
        let codes = &sc.codes;

        // Every element of `fwd`, `fwd_scale` and `emit` below index n·ns /
        // n / n·w is written before it is read, so the buffers are only
        // grown, never cleared (zero-filling them per sequence was a
        // measurable share of index construction).
        if sc.fwd.len() < n * ns {
            sc.fwd.resize(n * ns, 0.0);
        }
        if sc.fwd_scale.len() < n {
            sc.fwd_scale.resize(n, 1.0);
        }
        // Emission ratio of repeat state i (period i+1) at position p:
        // lr[x_p][x_{p-i-1}], 0 where no partner residue exists.
        if sc.emit.len() < n * w {
            sc.emit.resize(n * w, 0.0);
        }
        for p in 0..n {
            let lr = &self.lr[codes[p] as usize];
            let maxi = p.min(w);
            let row = &mut sc.emit[p * w..(p + 1) * w];
            for (dst, &c) in row[..maxi].iter_mut().zip(codes[p - maxi..p].iter().rev()) {
                // codes are clamped to 0..=23 above
                *dst = lr[(c as usize).min(23)];
            }
            for dst in &mut row[maxi..] {
                *dst = 0.0;
            }
        }

        // ---- forward ----
        let mut prev_b = 1.0f32;
        let mut prev_f = vec![0.0f32; w];
        let mut cur_f = vec![0.0f32; w];
        for p in 0..n {
            let emit = &sc.emit[p * w..(p + 1) * w];
            let sum_f = sum8(&prev_f);
            let cur_b = prev_b * b2b + sum_f * f2b;
            for i in 0..w {
                cur_f[i] = (prev_b * self.b2f[i] + prev_f[i] * f2f) * emit[i];
            }
            let tot = cur_b + sum8(&cur_f);
            let s = if tot > 0.0 { 1.0 / tot } else { 1.0 };
            sc.fwd_scale[p] = s;
            let row = &mut sc.fwd[p * ns..(p + 1) * ns];
            row[0] = cur_b * s;
            for (r, &f) in row[1..].iter_mut().zip(cur_f.iter()) {
                *r = f * s;
            }
            prev_b = row[0];
            prev_f.copy_from_slice(&row[1..]);
        }

        // ---- backward (same scale factors) ----
        // The sequence must end in the background state (as in tantan): a
        // repeat state reaches the end only through its F->B transition.
        let mut nb = 1.0f32; // beta_B(p+1)
        let mut nf = vec![f2b; w]; // beta_Fi(p+1)
        let mut nnf = vec![0.0f32; w];
        let mut be = vec![0.0f32; w];
        for p in (0..n).rev() {
            let row = &sc.fwd[p * ns..(p + 1) * ns];
            let pb = row[0] * nb;
            let pf = dot8(&row[1..], &nf);
            let z = pb + pf;
            sc.post[p] = if z > 0.0 { pf / z } else { 0.0 };
            if p == 0 {
                break;
            }
            let emit = &sc.emit[p * w..(p + 1) * w];
            let s = sc.fwd_scale[p];
            for i in 0..w {
                be[i] = emit[i] * nf[i];
            }
            let to_fg = dot8(&self.b2f, &be);
            let fb = f2b * nb;
            for i in 0..w {
                nnf[i] = (f2f * be[i] + fb) * s;
            }
            nb = (b2b * nb + to_fg) * s;
            std::mem::swap(&mut nf, &mut nnf);
        }
    }
}

/// λ implied by a score matrix over the 20 standard letters: the root of
/// Σ_ij [(e^{λS})⁻¹]_ij = 1 for which the implied letter frequencies are all
/// positive (Yu, Wootton & Altschul 2003; the method of LAST/tantan).
pub fn implied_lambda(m: &ScoreMatrix) -> Option<f64> {
    let f = |l: f64| -> Option<(f64, bool)> {
        let mut a = [[0.0f64; 20]; 20];
        for i in 0..20 {
            for j in 0..20 {
                a[i][j] = (l * m.matrix[i][j] as f64).exp();
            }
        }
        let inv = invert20(a)?;
        let mut total = 0.0;
        let mut positive = true;
        for i in 0..20 {
            let rs: f64 = inv[i].iter().sum();
            let cs: f64 = (0..20).map(|k| inv[k][i]).sum();
            if rs <= 0.0 || cs <= 0.0 {
                positive = false;
            }
            total += rs;
        }
        Some((total - 1.0, positive))
    };
    let mut prev: Option<(f64, f64)> = None;
    let mut l = 0.02;
    while l < 3.0 {
        if let Some((v, _)) = f(l) {
            if let Some((pl, pv)) = prev {
                if (pv > 0.0) != (v > 0.0) {
                    // bisection
                    let (mut lo, mut hi, mut vlo) = (pl, l, pv);
                    for _ in 0..80 {
                        let mid = 0.5 * (lo + hi);
                        let vm = f(mid)?.0;
                        if (vm > 0.0) == (vlo > 0.0) {
                            lo = mid;
                            vlo = vm;
                        } else {
                            hi = mid;
                        }
                    }
                    let root = 0.5 * (lo + hi);
                    if let Some((_, pos)) = f(root) {
                        if pos {
                            return Some(root);
                        }
                    }
                }
            }
            prev = Some((l, v));
        }
        l += 0.01;
    }
    None
}

fn invert20(mut a: [[f64; 20]; 20]) -> Option<[[f64; 20]; 20]> {
    let mut inv = [[0.0f64; 20]; 20];
    for (i, row) in inv.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for col in 0..20 {
        let mut piv = col;
        for r in col + 1..20 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-300 {
            return None;
        }
        a.swap(col, piv);
        inv.swap(col, piv);
        let d = a[col][col];
        for j in 0..20 {
            a[col][j] /= d;
            inv[col][j] /= d;
        }
        for r in 0..20 {
            if r != col {
                let f = a[r][col];
                if f != 0.0 {
                    for j in 0..20 {
                        a[r][j] -= f * a[col][j];
                        inv[r][j] -= f * inv[col][j];
                    }
                }
            }
        }
    }
    Some(inv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::EncodedProtein;
    use crate::scoring::ScoringSystem;

    #[test]
    fn implied_lambda_of_blosum62_is_half_bit_scale() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let l = implied_lambda(&m).unwrap();
        // BLOSUM62 is in half-bit units: λ ≈ ln2/2 ≈ 0.347 (rounding shifts it).
        assert!(l > 0.30 && l < 0.36, "lambda {l}");
    }

    #[test]
    fn homopolymer_is_masked_diverse_sequence_is_not() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let t = Tantan::protein(&m);
        let mut s = b"MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPTDPRRRSRNLGKVIDTLTCGFADLMGYIPLVGAPLGGAARALAHGVRVLEDGVNYATGNLPGCSFSIFLLALLSCLTVPASA".to_vec();
        let start = s.len();
        s.extend(std::iter::repeat(b'Q').take(30));
        s.extend_from_slice(b"MKVLAAGIVALLLAAGCSSSKEETTEEKPAEEAPKAEEKPAEEAPK");
        let e = EncodedProtein::from_ascii(&s);
        let mut sc = TantanScratch::default();
        let mut mask = Vec::new();
        t.mask(e.as_slice(), &mut sc, &mut mask);
        let polyq_masked = mask[start..start + 30].iter().filter(|&&b| b).count();
        assert!(polyq_masked >= 25, "poly-Q masked {polyq_masked}/30");
        let head_masked = mask[..60].iter().filter(|&&b| b).count();
        assert!(head_masked <= 6, "diverse head masked {head_masked}/60");
    }

    /// The AVX2 path (default 50 repeat states) gives the scalar body's
    /// posteriors bit for bit, on random, low-complexity and periodic
    /// sequences of many lengths (including empty and shorter than a vector).
    #[test]
    fn avx2_masker_is_bit_identical_to_the_scalar_body() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let t = Tantan::protein(&m);
        let mut x: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut rnd = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let mut seqs: Vec<Vec<u8>> = Vec::new();
        for &n in &[0usize, 1, 2, 3, 7, 8, 9, 48, 49, 50, 51, 56, 57, 300, 2000] {
            // random residues, including X / B / Z / * codes (20..23)
            seqs.push((0..n).map(|_| (rnd() % 24) as u8).collect());
        }
        for period in 1..=60usize {
            // periodic repeats with point mutations, flanked by random sequence
            let unit: Vec<u8> = (0..period).map(|_| (rnd() % 20) as u8).collect();
            let mut s: Vec<u8> = (0..40).map(|_| (rnd() % 20) as u8).collect();
            for k in 0..(400 / period).max(3) * period {
                s.push(if rnd() % 10 == 0 { (rnd() % 20) as u8 } else { unit[k % period] });
            }
            s.extend((0..40).map(|_| (rnd() % 20) as u8));
            seqs.push(s);
        }
        seqs.push(vec![5u8; 1000]); // homopolymer
        seqs.push((0..1500).map(|i| [0u8, 7, 15][i % 3]).collect());
        let (mut sa, mut sb) = (TantanScratch::default(), TantanScratch::default());
        let mut checked = 0usize;
        for s in &seqs {
            t.compute_body(s, &mut sa);
            t.compute(s, &mut sb);
            assert_eq!(sa.post.len(), s.len());
            assert_eq!(sb.post.len(), s.len());
            for (i, (a, b)) in sa.post.iter().zip(&sb.post).enumerate() {
                assert_eq!(a.to_bits(), b.to_bits(), "length {} position {i}: {a} vs {b}", s.len());
            }
            checked += s.len();
        }
        assert!(checked > 20_000);
        // a masker with another number of repeat states takes the generic path
        let t2 = Tantan::with_params(&m, TantanParams { max_period: 20, ..TantanParams::default() }, None);
        let s = &seqs[seqs.len() - 1];
        t2.compute_body(s, &mut sa);
        t2.compute(s, &mut sb);
        assert!(sa.post.iter().zip(&sb.post).all(|(a, b)| a.to_bits() == b.to_bits()));
    }

    #[test]
    fn probabilities_are_in_unit_interval() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let t = Tantan::protein(&m);
        let e = EncodedProtein::from_ascii(b"AAAAAAAAAAPPPPPPPPGSGSGSGSGSGSMKVLHHHHHHHHHH");
        let mut sc = TantanScratch::default();
        let p = t.repeat_probs(e.as_slice(), &mut sc);
        assert!(p.iter().all(|&x| (0.0..=1.0001).contains(&x)));
    }
}
