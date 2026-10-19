//! Local composition-bias correction of substitution scores.
//!
//! Compositionally biased regions (low complexity, coiled coils, transmembrane
//! stretches, His/Gln/Pro-rich linkers, …) produce high Smith–Waterman scores
//! against unrelated sequences of similar composition, which is a major source
//! of false positives. SABER uses the *documented* per-position correction of
//! Hauser, Steinegger & Söding (2016, Bioinformatics 32:1323) as implemented in
//! MMseqs2 (`SubstitutionMatrix::calcLocalAaBiasCorrection`) and in DIAMOND's
//! default `--comp-based-stats 1` (`HauserCorrection`, src/stats/
//! hauser_correction.cpp):
//!
//! ```text
//! c_i = Σ_a p_a · S(q_i, a)  −  (Σ_{j∈W(i), j≠i} S(q_i, q_j)) / (|W(i)| − 1)
//! ```
//!
//! where `W(i)` is the window of up to 41 query positions centred on `i`
//! (`i−20 … i+20`, clipped at the sequence ends) and `p_a` are the background
//! frequencies (Robinson & Robinson 1991, NCBI `Robinson_prob`). The second
//! term is the score `q_i` would receive against its own local neighbourhood;
//! when the neighbourhood is enriched in residues similar to `q_i` the term is
//! large and `c_i` becomes negative, damping every substitution score in that
//! query row. `c_i` is rounded to the nearest integer (as DIAMOND does) and
//! added to `S(q_i, s_j)` for every aligned (non-gap) column in *all* scoring
//! stages — ungapped extension, SIMD score-only DP, banded DP and traceback —
//! so the final score, the E-value pre-filter and the reported alignment all
//! agree. Only standard residues receive a correction.

use super::matrices::ScoreMatrix;

/// Robinson & Robinson (1991) background amino-acid frequencies, per 1000,
/// in NCBI order ARNDCQEGHILKMFPSTWYV (from NCBI blast_stat.c `Robinson_prob`).
pub const ROBINSON_PER_1000: [f64; 20] = [
    78.05, // A
    51.29, // R
    44.87, // N
    53.64, // D
    19.25, // C
    42.64, // Q
    62.95, // E
    73.77, // G
    21.99, // H
    51.42, // I
    90.19, // L
    57.44, // K
    22.43, // M
    38.56, // F
    52.03, // P
    71.20, // S
    58.41, // T
    13.30, // W
    32.16, // Y
    64.41, // V
];

/// Background frequencies normalised to sum to 1.
pub fn background_freqs() -> [f64; 20] {
    let s: f64 = ROBINSON_PER_1000.iter().sum();
    let mut out = [0.0; 20];
    for i in 0..20 {
        out[i] = ROBINSON_PER_1000[i] / s;
    }
    out
}

/// Half-width of the correction window (DIAMOND `--cbs-window 40`).
pub const BIAS_WINDOW_HALF: usize = 20;

/// Expected score of residue `r` against the background: Σ_a p_a S(r, a).
pub fn background_scores(m: &ScoreMatrix) -> [f64; 20] {
    let p = background_freqs();
    let mut out = [0.0; 20];
    for r in 0..20 {
        out[r] = (0..20).map(|a| p[a] * m.matrix[r][a] as f64).sum();
    }
    out
}

/// Compute the per-position integer correction for an encoded query
/// (codes 0..=23). Positions holding non-standard residues get 0.
pub fn composition_bias(query: &[u8], m: &ScoreMatrix, out: &mut Vec<i8>) {
    out.clear();
    out.resize(query.len(), 0);
    let n = query.len();
    if n < 2 {
        return;
    }
    let bg = background_scores(m);
    // Rolling sum of S(r, q_j) for all 20 r over the window.
    let mut acc = [0i32; 20];
    let add = |acc: &mut [i32; 20], c: u8, sign: i32| {
        let c = (c as usize).min(23);
        for r in 0..20 {
            acc[r] += sign * m.matrix[r][c];
        }
    };
    let half = BIAS_WINDOW_HALF;
    let mut lo = 0usize; // window [lo, hi)
    let mut hi = 0usize;
    for i in 0..n {
        let want_hi = (i + half + 1).min(n);
        let want_lo = i.saturating_sub(half);
        while hi < want_hi {
            add(&mut acc, query[hi], 1);
            hi += 1;
        }
        while lo < want_lo {
            add(&mut acc, query[lo], -1);
            lo += 1;
        }
        let r = query[i] as usize;
        if r < 20 {
            let w = (hi - lo) as i32;
            if w >= 2 {
                let own = m.matrix[r][r];
                let neigh = (acc[r] - own) as f64 / (w - 1) as f64;
                let f = bg[r] - neigh;
                let v = if f < 0.0 { f - 0.5 } else { f + 0.5 } as i32;
                out[i] = v.clamp(-127, 127) as i8;
            }
        }
    }
}

/// Convenience: allocate and return the correction vector.
pub fn composition_bias_vec(query: &[u8], m: &ScoreMatrix) -> Vec<i8> {
    let mut v = Vec::new();
    composition_bias(query, m, &mut v);
    v
}

/// Ungapped λ for a score matrix and residue frequencies `p` (rows) and `q`
/// (columns): the positive root of Σ p_i q_j e^{λ s_ij} = 1. Returns `None`
/// when the expected score is not negative (no positive root exists).
pub fn ungapped_lambda(m: &ScoreMatrix, p: &[f64; 20], q: &[f64; 20]) -> Option<f64> {
    let mut expected = 0.0;
    let mut maxs = i32::MIN;
    for i in 0..20 {
        for j in 0..20 {
            expected += p[i] * q[j] * m.matrix[i][j] as f64;
            if p[i] > 0.0 && q[j] > 0.0 {
                maxs = maxs.max(m.matrix[i][j]);
            }
        }
    }
    if expected >= 0.0 || maxs <= 0 {
        return None;
    }
    let f = |l: f64| -> f64 {
        let mut s = 0.0;
        for i in 0..20 {
            for j in 0..20 {
                s += p[i] * q[j] * (l * m.matrix[i][j] as f64).exp();
            }
        }
        s - 1.0
    };
    let (mut lo, mut hi) = (1e-6f64, 0.5f64);
    while f(hi) < 0.0 {
        hi *= 2.0;
        if hi > 50.0 {
            return None;
        }
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 { hi = mid } else { lo = mid }
    }
    Some(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::EncodedProtein;
    use crate::scoring::ScoringSystem;

    #[test]
    fn robinson_lambda_reproduces_ncbi_ungapped_blosum62() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let p = background_freqs();
        let l = ungapped_lambda(&m, &p, &p).unwrap();
        // NCBI reports 0.3176 for BLOSUM62 with Robinson frequencies.
        assert!((l - 0.3176).abs() < 0.002, "lambda {l}");
    }

    #[test]
    fn low_complexity_region_is_penalised() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let mut s = b"MKTAYIAKQRQISFVKSHFSRQ".to_vec();
        s.extend(std::iter::repeat(b'Q').take(40));
        s.extend_from_slice(b"LEDNWETLNDNLKVIEKADNAAQVKDALTKMRAAALDAQKATPPKLEDKSPDSPEMKDFRHGFDILVGQIDDALKLANEGKVKEAQAAAEQLKTTRNAYHQKYR");
        let e = EncodedProtein::from_ascii(&s);
        let c = composition_bias_vec(e.as_slice(), &m);
        let polyq_mid = 22 + 20;
        assert!(c[polyq_mid] <= -3, "poly-Q correction {}", c[polyq_mid]);
        // Diverse region: corrections are small.
        let tail_mean: f64 = c[80..120].iter().map(|&v| v as f64).sum::<f64>() / 40.0;
        assert!(tail_mean.abs() < 1.5, "diverse region mean correction {tail_mean}");
    }

    #[test]
    fn nonstandard_residues_get_zero() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let e = EncodedProtein::from_ascii(b"XXXXBZ*XXXX");
        assert!(composition_bias_vec(e.as_slice(), &m).iter().all(|&v| v == 0));
    }
}
