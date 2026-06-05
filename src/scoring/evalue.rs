//! Karlin-Altschul E-value and bit-score calculation.
//!
//! Given a raw alignment score `S`, search space size `m * n` (query length
//! times database length, both in residues/bases), and the matrix-specific
//! parameters λ and K:
//!
//! ```text
//! E = K * m * n * exp(-λ * S)
//! S' = (λ * S - ln K) / ln 2     // bit score
//! E = m * n * 2^(-S')
//! ```
//!
//! The old implementation hard-coded `evalue = 0` and `bitscore = 0` in the
//! traceback and never plumbed real parameters through. This module fixes that
//! end-to-end so `-m 8` (BLAST tab) output carries real, comparable numbers.

use std::f64::consts::LN_2;

#[derive(Debug, Clone, Copy)]
pub struct EValue {
    /// Karlin-Altschul λ
    pub lambda: f64,
    /// Karlin-Altschul K
    pub kappa: f64,
    /// Database size in residues (sum of all subject lengths).
    pub db_size: usize,
    /// Effective query length (after edge correction).
    pub query_size: usize,
}

impl Default for EValue {
    fn default() -> Self {
        // BLOSUM62 ungapped defaults.
        Self { lambda: 0.3176, kappa: 0.134, db_size: 1_000_000, query_size: 1_000 }
    }
}

impl EValue {
    /// Build with explicit matrix parameters.
    pub fn new(lambda: f64, kappa: f64, db_size: usize, query_size: usize) -> Self {
        Self { lambda, kappa, db_size, query_size }
    }

    /// Effective search space `m * n`.
    #[inline]
    pub fn search_space(&self) -> f64 {
        (self.query_size as f64) * (self.db_size as f64)
    }

    /// E = K * m * n * exp(-λ S).
    pub fn score_to_evalue(&self, raw_score: i32) -> f64 {
        if raw_score <= 0 {
            return f64::INFINITY;
        }
        let s = raw_score as f64;
        let e = self.kappa * self.search_space() * (-self.lambda * s).exp();
        // Clamp absurd underflow to a tiny positive so M8 formatters don't print 0.
        if e == 0.0 { f64::MIN_POSITIVE } else { e }
    }

    /// S' = (λ S - ln K) / ln 2.
    pub fn score_to_bitscore(&self, raw_score: i32) -> f64 {
        if raw_score <= 0 {
            return 0.0;
        }
        let s = raw_score as f64;
        (self.lambda * s - self.kappa.ln()) / LN_2
    }

    /// E = m * n * 2^(-S').
    pub fn bitscore_to_evalue(&self, bitscore: f64) -> f64 {
        self.search_space() * (-bitscore * LN_2).exp()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct BitScore {
    pub raw_score: i32,
    pub value: f64,
    pub evalue: f64,
}

impl BitScore {
    pub fn calculate(raw_score: i32, ev: &EValue) -> Self {
        Self {
            raw_score,
            value: ev.score_to_bitscore(raw_score),
            evalue: ev.score_to_evalue(raw_score),
        }
    }

    pub fn passes_evalue_threshold(&self, threshold: f64) -> bool {
        self.evalue <= threshold
    }

    pub fn meets_bitscore_minimum(&self, minimum: f64) -> bool {
        self.value >= minimum
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn higher_score_means_lower_evalue() {
        let ev = EValue::new(0.3176, 0.134, 1_000_000, 1_000);
        let low = ev.score_to_evalue(20);
        let high = ev.score_to_evalue(100);
        assert!(low > high, "score 100 should beat score 20: {} vs {}", low, high);
    }

    #[test]
    fn zero_score_yields_infinite_evalue() {
        let ev = EValue::default();
        assert!(ev.score_to_evalue(0).is_infinite());
        assert!(ev.score_to_evalue(-5).is_infinite());
    }

    #[test]
    fn bitscore_monotonic_in_score() {
        let ev = EValue::default();
        let a = ev.score_to_bitscore(50);
        let b = ev.score_to_bitscore(100);
        assert!(b > a);
    }

    #[test]
    fn roundtrip_bitscore_to_evalue_close_to_direct() {
        let ev = EValue::new(0.3176, 0.134, 1_000_000, 1_000);
        let bs = ev.score_to_bitscore(80);
        let e_via_bits = ev.bitscore_to_evalue(bs);
        let e_direct = ev.score_to_evalue(80);
        let ratio = e_via_bits / e_direct;
        assert!((ratio - 1.0).abs() < 1e-6, "ratio {}", ratio);
    }
}
