//! Karlin–Altschul statistics for *gapped* local alignments.
//!
//! ## Provenance
//!
//! Every parameter set below is copied from the NCBI C++ Toolkit,
//! `src/algo/blast/core/blast_stat.c` (arrays `blosum45_values`,
//! `blosum50_values`, `blosum62_values`, `blosum80_values`, `blosum90_values`,
//! `pam30_values`, `pam70_values`, `pam250_values`, `blastn_values_*`), fetched
//! from github.com/ncbi/ncbi-cxx-toolkit-public on 2026-10-18
//! (sha256 `d3f28396…a3`). The length-adjustment routine is a line-by-line port
//! of `BLAST_ComputeLengthAdjustment` from the same file.
//!
//! ## Gap-penalty conventions
//!
//! NCBI tables are keyed by (gap *existence* `a`, gap *extension* `b`) and a
//! gap of length `k` costs `a + b·k`. SABER's `--gap-open G --gap-extend E`
//! has always meant: the first gap residue costs `G`, every further residue
//! costs `E`, i.e. `G + (k−1)·E`. The two agree when `a = G − E`, `b = E`, so
//! SABER's historical default `-g 11 -e 1` is BLAST's `-gapopen 10
//! -gapextend 1` (and BLAST/DIAMOND's default 11/1 is SABER `-g 12 -e 1`).
//! [`GapConvention::Blast`] lets users type BLAST-style numbers directly.
//!
//! ## What happens with unsupported combinations
//!
//! Only (matrix, gap) combinations present in the NCBI tables have published
//! gapped parameters. Anything else is rejected by default with a list of the
//! supported values. With `allow_uncalibrated` the search runs using the
//! *ungapped* parameters and the result is flagged `calibrated = false`, which
//! the CLI reports on stderr and in structured outputs. Parameters are never
//! substituted silently.

use std::f64::consts::LN_2;

use super::matrices::ScoringSystem;

/// How the user-facing gap penalties should be interpreted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapConvention {
    /// SABER: cost(k) = open + (k−1)·extend   (historical SABER semantics)
    Saber,
    /// BLAST/DIAMOND: cost(k) = open + k·extend
    Blast,
}

impl std::str::FromStr for GapConvention {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "saber" => Ok(Self::Saber),
            "blast" | "diamond" | "ncbi" => Ok(Self::Blast),
            o => Err(format!("unknown gap convention '{o}' (expected saber|blast)")),
        }
    }
}

/// Affine gap costs in SABER's internal convention:
/// first gap residue costs `open`, each further residue `extend`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GapCosts {
    pub open: i32,
    pub extend: i32,
}

impl GapCosts {
    /// Convert user input to SABER's internal convention, validating range.
    pub fn from_user(open: i32, extend: i32, conv: GapConvention) -> Result<Self, String> {
        if extend < 1 {
            return Err(format!("gap extension penalty must be >= 1 (got {extend})"));
        }
        if open < 0 {
            return Err(format!("gap open penalty must be >= 0 (got {open})"));
        }
        let open_internal = match conv {
            GapConvention::Saber => open,
            GapConvention::Blast => open
                .checked_add(extend)
                .ok_or_else(|| "gap penalties overflow".to_string())?,
        };
        if open_internal < extend {
            return Err(format!(
                "gap open ({open_internal} in SABER convention) must be >= gap extend ({extend})"
            ));
        }
        // Keep everything comfortably inside i16 so the SIMD kernels cannot
        // wrap on the penalty itself.
        if open_internal > 1000 || extend > 1000 {
            return Err("gap penalties above 1000 are not supported".into());
        }
        Ok(Self { open: open_internal, extend })
    }

    /// NCBI (existence, extension) pair for table lookup.
    #[inline]
    pub fn ncbi(&self) -> (i32, i32) {
        (self.open - self.extend, self.extend)
    }

    /// Cost of a gap of length `k >= 1`.
    #[inline]
    pub fn cost(&self, k: u32) -> i64 {
        self.open as i64 + (k.saturating_sub(1)) as i64 * self.extend as i64
    }
}

/// One row of an NCBI table: (existence, extension, λ, K, H, α, β).
#[derive(Debug, Clone, Copy)]
struct Row(i32, i32, f64, f64, f64, f64, f64);

// --- Protein tables (blast_stat.c). Row 0 of each NCBI array is the ungapped
//     set, stored separately here as `UNGAPPED_*`. ---------------------------
const UNGAPPED_BLOSUM45: Row = Row(0, 0, 0.2291, 0.0924, 0.2514, 0.9113, -5.7);
const BLOSUM45: &[Row] = &[
    Row(13, 3, 0.207, 0.049, 0.14, 1.5, -22.0),
    Row(12, 3, 0.199, 0.039, 0.11, 1.8, -34.0),
    Row(11, 3, 0.190, 0.031, 0.095, 2.0, -38.0),
    Row(10, 3, 0.179, 0.023, 0.075, 2.4, -51.0),
    Row(16, 2, 0.210, 0.051, 0.14, 1.5, -24.0),
    Row(15, 2, 0.203, 0.041, 0.12, 1.7, -31.0),
    Row(14, 2, 0.195, 0.032, 0.10, 1.9, -36.0),
    Row(13, 2, 0.185, 0.024, 0.084, 2.2, -45.0),
    Row(12, 2, 0.171, 0.016, 0.061, 2.8, -65.0),
    Row(19, 1, 0.205, 0.040, 0.11, 1.9, -43.0),
    Row(18, 1, 0.198, 0.032, 0.10, 2.0, -43.0),
    Row(17, 1, 0.189, 0.024, 0.079, 2.4, -57.0),
    Row(16, 1, 0.176, 0.016, 0.063, 2.8, -67.0),
];
const UNGAPPED_BLOSUM50: Row = Row(0, 0, 0.2318, 0.112, 0.3362, 0.6895, -4.0);
const BLOSUM50: &[Row] = &[
    Row(13, 3, 0.212, 0.063, 0.19, 1.1, -16.0),
    Row(12, 3, 0.206, 0.055, 0.17, 1.2, -18.0),
    Row(11, 3, 0.197, 0.042, 0.14, 1.4, -25.0),
    Row(10, 3, 0.186, 0.031, 0.11, 1.7, -34.0),
    Row(9, 3, 0.172, 0.022, 0.082, 2.1, -48.0),
    Row(16, 2, 0.215, 0.066, 0.20, 1.05, -15.0),
    Row(15, 2, 0.210, 0.058, 0.17, 1.2, -20.0),
    Row(14, 2, 0.202, 0.045, 0.14, 1.4, -27.0),
    Row(13, 2, 0.193, 0.035, 0.12, 1.6, -32.0),
    Row(12, 2, 0.181, 0.025, 0.095, 1.9, -41.0),
    Row(19, 1, 0.212, 0.057, 0.18, 1.2, -21.0),
    Row(18, 1, 0.207, 0.050, 0.15, 1.4, -28.0),
    Row(17, 1, 0.198, 0.037, 0.12, 1.6, -33.0),
    Row(16, 1, 0.186, 0.025, 0.10, 1.9, -42.0),
    Row(15, 1, 0.171, 0.015, 0.063, 2.7, -76.0),
];
const UNGAPPED_BLOSUM62: Row = Row(0, 0, 0.3176, 0.134, 0.4012, 0.7916, -3.2);
const BLOSUM62: &[Row] = &[
    Row(11, 2, 0.297, 0.082, 0.27, 1.1, -10.0),
    Row(10, 2, 0.291, 0.075, 0.23, 1.3, -15.0),
    Row(9, 2, 0.279, 0.058, 0.19, 1.5, -19.0),
    Row(8, 2, 0.264, 0.045, 0.15, 1.8, -26.0),
    Row(7, 2, 0.239, 0.027, 0.10, 2.5, -46.0),
    Row(6, 2, 0.201, 0.012, 0.061, 3.3, -58.0),
    Row(13, 1, 0.292, 0.071, 0.23, 1.2, -11.0),
    Row(12, 1, 0.283, 0.059, 0.19, 1.5, -19.0),
    Row(11, 1, 0.267, 0.041, 0.14, 1.9, -30.0),
    Row(10, 1, 0.243, 0.024, 0.10, 2.5, -44.0),
    Row(9, 1, 0.206, 0.010, 0.052, 4.0, -87.0),
];
const UNGAPPED_BLOSUM80: Row = Row(0, 0, 0.3430, 0.177, 0.6568, 0.5222, -1.6);
const BLOSUM80: &[Row] = &[
    Row(25, 2, 0.342, 0.17, 0.66, 0.52, -1.6),
    Row(13, 2, 0.336, 0.15, 0.57, 0.59, -3.0),
    Row(9, 2, 0.319, 0.11, 0.42, 0.76, -6.0),
    Row(8, 2, 0.308, 0.090, 0.35, 0.89, -9.0),
    Row(7, 2, 0.293, 0.070, 0.27, 1.1, -14.0),
    Row(6, 2, 0.268, 0.045, 0.19, 1.4, -19.0),
    Row(11, 1, 0.314, 0.095, 0.35, 0.90, -9.0),
    Row(10, 1, 0.299, 0.071, 0.27, 1.1, -14.0),
    Row(9, 1, 0.279, 0.048, 0.20, 1.4, -19.0),
];
const UNGAPPED_BLOSUM90: Row = Row(0, 0, 0.3346, 0.190, 0.7547, 0.4434, -1.4);
const BLOSUM90: &[Row] = &[
    Row(9, 2, 0.310, 0.12, 0.46, 0.67, -6.0),
    Row(8, 2, 0.300, 0.099, 0.39, 0.76, -7.0),
    Row(7, 2, 0.283, 0.072, 0.30, 0.93, -11.0),
    Row(6, 2, 0.259, 0.048, 0.22, 1.2, -16.0),
    Row(11, 1, 0.302, 0.093, 0.39, 0.78, -8.0),
    Row(10, 1, 0.290, 0.075, 0.28, 1.04, -15.0),
    Row(9, 1, 0.265, 0.044, 0.20, 1.3, -19.0),
];
const UNGAPPED_PAM30: Row = Row(0, 0, 0.3400, 0.283, 1.754, 0.1938, -0.3);
const PAM30: &[Row] = &[
    Row(7, 2, 0.305, 0.15, 0.87, 0.35, -3.0),
    Row(6, 2, 0.287, 0.11, 0.68, 0.42, -4.0),
    Row(5, 2, 0.264, 0.079, 0.45, 0.59, -7.0),
    Row(10, 1, 0.309, 0.15, 0.88, 0.35, -3.0),
    Row(9, 1, 0.294, 0.11, 0.61, 0.48, -6.0),
    Row(8, 1, 0.270, 0.072, 0.40, 0.68, -10.0),
    Row(15, 3, 0.339, 0.28, 1.70, 0.20, -0.5),
    Row(14, 2, 0.337, 0.27, 1.62, 0.21, -0.8),
    Row(14, 1, 0.333, 0.27, 1.43, 0.23, -1.4),
    Row(13, 3, 0.338, 0.27, 1.69, 0.20, -0.5),
];
const UNGAPPED_PAM70: Row = Row(0, 0, 0.3345, 0.229, 1.029, 0.3250, -0.7);
const PAM70: &[Row] = &[
    Row(8, 2, 0.301, 0.12, 0.54, 0.56, -5.0),
    Row(7, 2, 0.286, 0.093, 0.43, 0.67, -7.0),
    Row(6, 2, 0.264, 0.064, 0.29, 0.90, -12.0),
    Row(11, 1, 0.305, 0.12, 0.52, 0.59, -6.0),
    Row(10, 1, 0.291, 0.091, 0.41, 0.71, -9.0),
    Row(9, 1, 0.270, 0.060, 0.28, 0.97, -14.0),
    Row(11, 2, 0.323, 0.186, 0.80, 1.32, -27.0),
    Row(12, 3, 0.330, 0.219, 0.93, 0.82, -16.0),
];
const UNGAPPED_PAM250: Row = Row(0, 0, 0.2252, 0.0868, 0.2223, 0.98, -5.0);
const PAM250: &[Row] = &[
    Row(15, 3, 0.205, 0.049, 0.13, 1.6, -23.0),
    Row(14, 3, 0.200, 0.043, 0.12, 1.7, -26.0),
    Row(13, 3, 0.194, 0.036, 0.10, 1.9, -31.0),
    Row(12, 3, 0.186, 0.029, 0.085, 2.2, -41.0),
    Row(11, 3, 0.174, 0.020, 0.070, 2.5, -48.0),
    Row(17, 2, 0.204, 0.047, 0.12, 1.7, -28.0),
    Row(16, 2, 0.198, 0.038, 0.11, 1.8, -29.0),
    Row(15, 2, 0.191, 0.031, 0.087, 2.2, -44.0),
    Row(14, 2, 0.182, 0.024, 0.073, 2.5, -53.0),
    Row(13, 2, 0.171, 0.017, 0.059, 2.9, -64.0),
    Row(21, 1, 0.205, 0.045, 0.11, 1.8, -34.0),
    Row(20, 1, 0.199, 0.037, 0.10, 1.9, -35.0),
    Row(19, 1, 0.192, 0.029, 0.083, 2.3, -52.0),
    Row(18, 1, 0.183, 0.021, 0.070, 2.6, -60.0),
    Row(17, 1, 0.171, 0.014, 0.052, 3.3, -86.0),
];

// --- Nucleotide tables (blast_stat.c `blastn_values_R_P`). The NCBI "0, 0"
//     rows denote non-affine (megablast-style) gapping, which SABER does not
//     implement, so they are omitted. ------------------------------------------
struct NtTable {
    reward: i32,
    penalty: i32,
    /// 2/-3, 2/-5, 2/-7: odd scores must be rounded down before E-values.
    even_only: bool,
    rows: &'static [Row],
}
const NT_TABLES: &[NtTable] = &[
    NtTable { reward: 1, penalty: -5, even_only: false, rows: &[Row(3, 3, 1.39, 0.747, 1.38, 1.00, 0.0)] },
    NtTable { reward: 1, penalty: -4, even_only: false, rows: &[
        Row(1, 2, 1.36, 0.67, 1.2, 1.1, 0.0), Row(0, 2, 1.26, 0.43, 0.90, 1.4, -1.0),
        Row(2, 1, 1.35, 0.61, 1.1, 1.2, -1.0), Row(1, 1, 1.22, 0.35, 0.72, 1.7, -3.0)] },
    NtTable { reward: 2, penalty: -7, even_only: true, rows: &[
        Row(2, 4, 0.68, 0.67, 1.2, 0.55, 0.0), Row(0, 4, 0.63, 0.43, 0.90, 0.7, -1.0),
        Row(4, 2, 0.675, 0.62, 1.1, 0.6, -1.0), Row(2, 2, 0.61, 0.35, 0.72, 1.7, -3.0)] },
    NtTable { reward: 1, penalty: -3, even_only: false, rows: &[
        Row(2, 2, 1.37, 0.70, 1.2, 1.1, 0.0), Row(1, 2, 1.35, 0.64, 1.1, 1.2, -1.0),
        Row(0, 2, 1.25, 0.42, 0.83, 1.5, -2.0), Row(2, 1, 1.34, 0.60, 1.1, 1.2, -1.0),
        Row(1, 1, 1.21, 0.34, 0.71, 1.7, -2.0)] },
    NtTable { reward: 2, penalty: -5, even_only: true, rows: &[
        Row(2, 4, 0.67, 0.59, 1.1, 0.6, -1.0), Row(0, 4, 0.62, 0.39, 0.78, 0.8, -2.0),
        Row(4, 2, 0.67, 0.61, 1.0, 0.65, -2.0), Row(2, 2, 0.56, 0.32, 0.59, 0.95, -4.0)] },
    NtTable { reward: 1, penalty: -2, even_only: false, rows: &[
        Row(2, 2, 1.33, 0.62, 1.1, 1.2, 0.0), Row(1, 2, 1.30, 0.52, 0.93, 1.4, -2.0),
        Row(0, 2, 1.19, 0.34, 0.66, 1.8, -3.0), Row(3, 1, 1.32, 0.57, 1.0, 1.3, -1.0),
        Row(2, 1, 1.29, 0.49, 0.92, 1.4, -1.0), Row(1, 1, 1.14, 0.26, 0.52, 2.2, -5.0)] },
    NtTable { reward: 2, penalty: -3, even_only: true, rows: &[
        Row(4, 4, 0.63, 0.42, 0.84, 0.75, -2.0), Row(2, 4, 0.615, 0.37, 0.72, 0.85, -3.0),
        Row(0, 4, 0.55, 0.21, 0.46, 1.2, -5.0), Row(3, 3, 0.615, 0.37, 0.68, 0.9, -3.0),
        Row(6, 2, 0.63, 0.42, 0.84, 0.75, -2.0), Row(5, 2, 0.625, 0.41, 0.78, 0.8, -2.0),
        Row(4, 2, 0.61, 0.35, 0.68, 0.9, -3.0), Row(2, 2, 0.515, 0.14, 0.33, 1.55, -9.0)] },
    NtTable { reward: 3, penalty: -4, even_only: false, rows: &[
        Row(6, 3, 0.389, 0.25, 0.56, 0.7, -5.0), Row(5, 3, 0.375, 0.21, 0.47, 0.8, -6.0),
        Row(4, 3, 0.351, 0.14, 0.35, 1.0, -9.0), Row(6, 2, 0.362, 0.16, 0.45, 0.8, -4.0),
        Row(5, 2, 0.330, 0.092, 0.28, 1.2, -13.0), Row(4, 2, 0.281, 0.046, 0.16, 1.8, -23.0)] },
    NtTable { reward: 4, penalty: -5, even_only: false, rows: &[
        Row(6, 5, 0.28, 0.21, 0.47, 0.6, -7.0), Row(5, 5, 0.27, 0.17, 0.39, 0.7, -9.0),
        Row(4, 5, 0.25, 0.10, 0.31, 0.8, -10.0), Row(3, 5, 0.23, 0.065, 0.25, 0.9, -11.0)] },
    NtTable { reward: 1, penalty: -1, even_only: false, rows: &[
        Row(3, 2, 1.09, 0.31, 0.55, 2.0, -2.0), Row(2, 2, 1.07, 0.27, 0.49, 2.2, -3.0),
        Row(1, 2, 1.02, 0.21, 0.36, 2.8, -6.0), Row(0, 2, 0.80, 0.064, 0.17, 4.8, -16.0),
        Row(4, 1, 1.08, 0.28, 0.54, 2.0, -2.0), Row(3, 1, 1.06, 0.25, 0.46, 2.3, -4.0),
        Row(2, 1, 0.99, 0.17, 0.30, 3.3, -10.0)] },
    NtTable { reward: 3, penalty: -2, even_only: false, rows: &[Row(5, 5, 0.208, 0.030, 0.072, 2.9, -47.0)] },
    NtTable { reward: 5, penalty: -4, even_only: false, rows: &[
        Row(10, 6, 0.163, 0.068, 0.16, 1.0, -19.0), Row(8, 6, 0.146, 0.039, 0.11, 1.3, -29.0)] },
];

/// Karlin–Altschul parameters actually used for a search.
#[derive(Debug, Clone, PartialEq)]
pub struct KarlinParams {
    pub lambda: f64,
    pub k: f64,
    pub log_k: f64,
    pub h: f64,
    /// α/λ for the length adjustment (BLAST uses 1/H for ungapped).
    pub alpha_d_lambda: f64,
    pub beta: f64,
    /// True when (matrix, gap costs) has published gapped parameters.
    pub calibrated: bool,
    /// Round odd raw scores down before computing E-values (blastn 2/-x).
    pub even_scores_only: bool,
    /// Human-readable provenance, e.g. "NCBI blast_stat.c BLOSUM62 10/1".
    pub source: String,
}

impl KarlinParams {
    fn from_row(r: Row, calibrated: bool, even: bool, source: String) -> Self {
        let Row(_, _, lambda, k, h, alpha, beta) = r;
        Self {
            lambda,
            k,
            log_k: k.ln(),
            h,
            alpha_d_lambda: alpha / lambda,
            beta,
            calibrated,
            even_scores_only: even,
            source,
        }
    }
}

/// Error for (matrix, gap) combinations without published gapped parameters.
#[derive(Debug, Clone)]
pub struct UnsupportedStats {
    pub message: String,
}
impl std::fmt::Display for UnsupportedStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}
impl std::error::Error for UnsupportedStats {}

fn protein_table(m: ScoringSystem) -> (Row, &'static [Row]) {
    match m {
        ScoringSystem::Blosum45 => (UNGAPPED_BLOSUM45, BLOSUM45),
        ScoringSystem::Blosum50 => (UNGAPPED_BLOSUM50, BLOSUM50),
        ScoringSystem::Blosum62 => (UNGAPPED_BLOSUM62, BLOSUM62),
        ScoringSystem::Blosum80 => (UNGAPPED_BLOSUM80, BLOSUM80),
        ScoringSystem::Blosum90 => (UNGAPPED_BLOSUM90, BLOSUM90),
        ScoringSystem::Pam30 => (UNGAPPED_PAM30, PAM30),
        ScoringSystem::Pam70 => (UNGAPPED_PAM70, PAM70),
        ScoringSystem::Pam250 => (UNGAPPED_PAM250, PAM250),
    }
}

/// Ungapped (λ, K, H) for a protein matrix (NCBI row 0).
pub fn ungapped_protein(m: ScoringSystem) -> KarlinParams {
    let (u, _) = protein_table(m);
    let mut p = KarlinParams::from_row(u, true, false, format!("NCBI blast_stat.c {m} ungapped"));
    // BLAST uses 1/H (and beta = 0) for the ungapped length adjustment.
    p.alpha_d_lambda = 1.0 / p.h;
    p.beta = 0.0;
    p
}

/// List supported gap pairs for a matrix, in SABER and BLAST conventions.
pub fn supported_protein_gaps(m: ScoringSystem) -> Vec<(i32, i32)> {
    protein_table(m).1.iter().map(|r| (r.0, r.1)).collect()
}

/// Gapped parameters for a protein search.
pub fn gapped_protein(
    m: ScoringSystem,
    gaps: GapCosts,
    allow_uncalibrated: bool,
) -> Result<KarlinParams, UnsupportedStats> {
    let (ungapped, rows) = protein_table(m);
    let (a, b) = gaps.ncbi();
    if let Some(r) = rows.iter().find(|r| r.0 == a && r.1 == b) {
        return Ok(KarlinParams::from_row(
            *r,
            true,
            false,
            format!(
                "NCBI blast_stat.c {m} existence={a} extension={b} (SABER -g {} -e {})",
                gaps.open, gaps.extend
            ),
        ));
    }
    let supported: Vec<String> = rows
        .iter()
        .map(|r| format!("-g {} -e {}", r.0 + r.1, r.1))
        .collect();
    let msg = format!(
        "no published gapped Karlin-Altschul parameters for {m} with SABER gap costs \
         -g {} -e {} (= BLAST existence {a}, extension {b}). Supported for {m} in SABER \
         convention: {}. Use one of these, or pass --allow-uncalibrated-evalues to run \
         with E-values explicitly marked UNCALIBRATED.",
        gaps.open,
        gaps.extend,
        supported.join(", ")
    );
    if allow_uncalibrated {
        let mut p = KarlinParams::from_row(
            ungapped,
            false,
            false,
            format!("UNCALIBRATED: {m} ungapped parameters used for gapped scores -g {} -e {}", gaps.open, gaps.extend),
        );
        p.alpha_d_lambda = 1.0 / p.h;
        p.beta = 0.0;
        Ok(p)
    } else {
        Err(UnsupportedStats { message: msg })
    }
}

/// Gapped parameters for a nucleotide (match/mismatch) search.
pub fn gapped_nucleotide(
    reward: i32,
    penalty: i32,
    gaps: GapCosts,
    allow_uncalibrated: bool,
) -> Result<KarlinParams, UnsupportedStats> {
    let (a, b) = gaps.ncbi();
    let table = NT_TABLES.iter().find(|t| t.reward == reward && t.penalty == penalty);
    if let Some(t) = table {
        if let Some(r) = t.rows.iter().find(|r| r.0 == a && r.1 == b) {
            return Ok(KarlinParams::from_row(
                *r,
                true,
                t.even_only,
                format!("NCBI blast_stat.c blastn {reward}/{penalty} existence={a} extension={b}"),
            ));
        }
    }
    let supported = table
        .map(|t| t.rows.iter().map(|r| format!("-g {} -e {}", r.0 + r.1, r.1)).collect::<Vec<_>>().join(", "))
        .unwrap_or_else(|| "none (unsupported match/mismatch pair)".into());
    let msg = format!(
        "no published gapped parameters for nucleotide scoring {reward}/{penalty} with SABER gap \
         costs -g {} -e {} (= BLAST existence {a}, extension {b}). Supported: {}. Pass \
         --allow-uncalibrated-evalues to run anyway with E-values marked UNCALIBRATED.",
        gaps.open, gaps.extend, supported
    );
    if allow_uncalibrated {
        // Approximate ungapped λ for uniform base composition by solving
        // 1/4·e^{λr} + 3/4·e^{λp} = 1; K is set to 0.1 and flagged.
        let lambda = solve_lambda_match_mismatch(reward as f64, penalty as f64);
        Ok(KarlinParams {
            lambda,
            k: 0.1,
            log_k: 0.1f64.ln(),
            h: 1.0,
            alpha_d_lambda: 1.0,
            beta: 0.0,
            calibrated: false,
            even_scores_only: false,
            source: format!("UNCALIBRATED: nucleotide {reward}/{penalty} approximate ungapped lambda"),
        })
    } else {
        Err(UnsupportedStats { message: msg })
    }
}

fn solve_lambda_match_mismatch(r: f64, p: f64) -> f64 {
    let f = |l: f64| 0.25 * (l * r).exp() + 0.75 * (l * p).exp() - 1.0;
    let (mut lo, mut hi) = (1e-6, 10.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 { hi = mid } else { lo = mid }
    }
    0.5 * (lo + hi)
}

/// Port of NCBI `BLAST_ComputeLengthAdjustment` (blast_stat.c).
///
/// Returns `(length_adjustment, converged)`.
pub fn compute_length_adjustment(
    k: f64,
    log_k: f64,
    alpha_d_lambda: f64,
    beta: f64,
    query_length: i64,
    db_length: i64,
    db_num_seqs: i64,
) -> (i64, bool) {
    const MAX_ITER: i32 = 20;
    let m = query_length as f64;
    let n = db_length as f64;
    let nn = db_num_seqs as f64;
    let mut ell_min = 0.0f64;
    let mut ell_max;
    let mut converged = false;
    let mut ell_next = 0.0f64;
    {
        let a = nn;
        let mb = m * nn + n;
        let c = n * m - m.max(n) / k;
        if c < 0.0 {
            return (0, false);
        }
        ell_max = 2.0 * c / (mb + (mb * mb - 4.0 * a * c).sqrt());
    }
    for i in 1..=MAX_ITER {
        let ell = ell_next;
        let ss = (m - ell) * (n - nn * ell);
        let ell_bar = alpha_d_lambda * (log_k + ss.ln()) + beta;
        if ell_bar >= ell {
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            if ell_min == ell_max {
                break;
            }
        } else {
            ell_max = ell;
        }
        if ell_min <= ell_bar && ell_bar <= ell_max {
            ell_next = ell_bar;
        } else {
            ell_next = if i == 1 { ell_max } else { (ell_min + ell_max) / 2.0 };
        }
    }
    if converged {
        let mut adj = ell_min as i64;
        let ell = ell_min.ceil();
        if ell <= ell_max {
            let ss = (m - ell) * (n - nn * ell);
            if alpha_d_lambda * (log_k + ss.ln()) + beta >= ell {
                adj = ell as i64;
            }
        }
        (adj, true)
    } else {
        (ell_min as i64, false)
    }
}

/// How an E-value accounts for the database.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvalueModel {
    /// One search space for the whole database (BLAST):
    /// `E = K·(m−ℓ)·(n_db − N·ℓ)·e^{−λS}`, ℓ the length adjustment of the
    /// query against the database. Every hit of a query with the same score
    /// gets the same E-value, whatever the length of its subject.
    Database,
    /// Per subject (FASTA/SSEARCH, HMMER-style per-sequence E-values): the
    /// pairwise search space of the query and that subject, with BLAST's
    /// length adjustment for the pair, times the number of database
    /// sequences N: `E = N·K·(m−ℓ)·(n−ℓ)·e^{−λS}`. A score is less likely by
    /// chance in a short subject than in a long one; on a database of equal
    /// lengths both models agree.
    PerSubject,
}

impl EvalueModel {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "per-subject" | "subject" => Some(Self::PerSubject),
            "database" | "db" => Some(Self::Database),
            _ => None,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::PerSubject => "per-subject",
            Self::Database => "database",
        }
    }
}

/// Database-level statistical context: parameters plus database size.
#[derive(Debug, Clone)]
pub struct SearchStatistics {
    pub params: KarlinParams,
    /// Effective database length in residues of the *alignment* alphabet
    /// (for a translated database: total nucleotides / 3).
    pub db_length: u64,
    /// Number of database sequences (BLAST's N in the length adjustment).
    pub db_num_seqs: u64,
    pub model: EvalueModel,
}

/// Per-query statistics: effective search space for one query.
#[derive(Debug, Clone, Copy)]
pub struct QueryStatistics {
    pub lambda: f64,
    pub k: f64,
    pub log_k: f64,
    pub length_adjustment: i64,
    pub eff_query_len: f64,
    pub eff_db_len: f64,
    /// Search space of the database model (also the reference for the
    /// per-subject model's region heuristics).
    pub search_space: f64,
    pub even_scores_only: bool,
    pub model: EvalueModel,
    query_len: i64,
    num_seqs: f64,
    alpha_d_lambda: f64,
    beta: f64,
}

impl SearchStatistics {
    /// Database model (BLAST E-values).
    pub fn new(params: KarlinParams, db_length: u64, db_num_seqs: u64) -> Self {
        Self::with_model(params, db_length, db_num_seqs, EvalueModel::Database)
    }

    pub fn with_model(params: KarlinParams, db_length: u64, db_num_seqs: u64, model: EvalueModel) -> Self {
        Self { params, db_length: db_length.max(1), db_num_seqs: db_num_seqs.max(1), model }
    }

    /// Effective search space for a query of `query_len` residues, using the
    /// BLAST length adjustment: (m − ℓ)·(n − N·ℓ).
    pub fn for_query(&self, query_len: usize) -> QueryStatistics {
        let p = &self.params;
        let m = query_len.max(1) as i64;
        let (mut adj, _conv) = compute_length_adjustment(
            p.k,
            p.log_k,
            p.alpha_d_lambda,
            p.beta,
            m,
            self.db_length as i64,
            self.db_num_seqs as i64,
        );
        // Mirror BLAST: never let the effective query length drop below 1/K.
        let min_q = (1.0 / p.k).ceil() as i64;
        if m - adj < min_q {
            adj = (m - min_q).max(0);
        }
        let eff_q = ((m - adj) as f64).max(1.0);
        let eff_db = (self.db_length as f64 - self.db_num_seqs as f64 * adj as f64).max(1.0);
        QueryStatistics {
            lambda: p.lambda,
            k: p.k,
            log_k: p.log_k,
            length_adjustment: adj,
            eff_query_len: eff_q,
            eff_db_len: eff_db,
            search_space: eff_q * eff_db,
            even_scores_only: p.even_scores_only,
            model: self.model,
            query_len: m,
            num_seqs: self.db_num_seqs as f64,
            alpha_d_lambda: p.alpha_d_lambda,
            beta: p.beta,
        }
    }
}

impl QueryStatistics {
    #[inline]
    fn adj_score(&self, raw: i64) -> f64 {
        let s = if self.even_scores_only && raw % 2 != 0 { raw - 1 } else { raw };
        s as f64
    }

    /// Search space of a hit to a subject of `subject_len` residues under
    /// the statistics' model (the database search space for
    /// [`EvalueModel::Database`]).
    pub fn space_for(&self, subject_len: usize) -> f64 {
        match self.model {
            EvalueModel::Database => self.search_space,
            EvalueModel::PerSubject => {
                let m = self.query_len.max(1);
                let n = subject_len.max(1) as i64;
                let (adj, _) = compute_length_adjustment(self.k, self.log_k, self.alpha_d_lambda, self.beta, m, n, 1);
                // As BLAST does for the query: effective lengths not below 1/K.
                let min_l = (1.0 / self.k).ceil() as i64;
                let adj = adj.min((m - min_l).max(0)).min((n - min_l).max(0)).max(0);
                self.num_seqs * ((m - adj).max(1) as f64) * ((n - adj).max(1) as f64)
            }
        }
    }

    /// A lower bound on [`Self::space_for`] that costs one logarithm: the
    /// length adjustment of a pair never exceeds α/λ·(ln K + ln(m·n)) + β.
    pub fn space_lower_bound(&self, subject_len: usize) -> f64 {
        match self.model {
            EvalueModel::Database => self.search_space,
            EvalueModel::PerSubject => {
                let m = self.query_len.max(1) as f64;
                let n = subject_len.max(1) as f64;
                let up = (self.alpha_d_lambda * (self.log_k + (m * n).ln()) + self.beta).max(0.0).ceil();
                self.num_seqs * (m - up).max(1.0) * (n - up).max(1.0)
            }
        }
    }

    /// A lower bound on the smallest raw score that reaches the E-value
    /// cut-off and the bit-score floor for a subject of `subject_len`
    /// residues (`db_floor`: the database model's exact floor).
    pub fn floor_lower_bound(&self, evalue_cut: f64, min_bits: Option<f64>, subject_len: usize, db_floor: i64) -> f64 {
        match self.model {
            EvalueModel::Database => db_floor as f64,
            EvalueModel::PerSubject => {
                let mut f = 1.0f64;
                if evalue_cut > 0.0 && evalue_cut.is_finite() {
                    f = f.max((self.k * self.space_lower_bound(subject_len) / evalue_cut).ln() / self.lambda);
                }
                if let Some(b) = min_bits {
                    f = f.max((b * LN_2 + self.log_k) / self.lambda);
                }
                if self.even_scores_only {
                    f -= 1.0;
                }
                f
            }
        }
    }

    /// E-value of a hit to a subject of `subject_len` residues.
    pub fn evalue_for(&self, raw: i64, subject_len: usize) -> f64 {
        self.evalue_in(raw, self.space_for(subject_len))
    }

    /// E-value for a given search space.
    pub fn evalue_in(&self, raw: i64, space: f64) -> f64 {
        if raw <= 0 {
            return f64::INFINITY;
        }
        let s = self.adj_score(raw);
        let e = self.k * space * (-self.lambda * s).exp();
        if e == 0.0 { f64::MIN_POSITIVE } else { e }
    }

    /// E = K·m'·n'·exp(−λS) with the database search space.
    pub fn evalue(&self, raw: i64) -> f64 {
        if raw <= 0 {
            return f64::INFINITY;
        }
        let s = self.adj_score(raw);
        let e = self.k * self.search_space * (-self.lambda * s).exp();
        if e == 0.0 { f64::MIN_POSITIVE } else { e }
    }

    /// Bit score S' = (λS − ln K)/ln 2.
    pub fn bitscore(&self, raw: i64) -> f64 {
        if raw <= 0 {
            return 0.0;
        }
        (self.lambda * self.adj_score(raw) - self.log_k) / LN_2
    }

    /// Smallest raw score that passes both an E-value cutoff and an optional
    /// bit-score floor. Exact (monotone inversion), so filtering on this
    /// floor before traceback never changes the reported set.
    pub fn min_score(&self, evalue_cut: f64, min_bits: Option<f64>) -> i64 {
        let mut floor = 1.0f64;
        if evalue_cut > 0.0 && evalue_cut.is_finite() {
            floor = floor.max(((self.k * self.search_space) / evalue_cut).ln() / self.lambda);
        }
        if let Some(b) = min_bits {
            floor = floor.max((b * LN_2 + self.log_k) / self.lambda);
        }
        if !floor.is_finite() {
            return 1;
        }
        let mut s = floor.ceil().max(1.0) as i64;
        // Guard against floating rounding at the boundary.
        while s > 1 && self.passes(s - 1, evalue_cut, min_bits) {
            s -= 1;
        }
        while !self.passes(s, evalue_cut, min_bits) && s < i64::MAX / 2 {
            s += 1;
        }
        s
    }

    #[inline]
    fn passes(&self, s: i64, evalue_cut: f64, min_bits: Option<f64>) -> bool {
        self.evalue(s) <= evalue_cut && min_bits.map_or(true, |b| self.bitscore(s) >= b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b62_11_1() -> KarlinParams {
        gapped_protein(ScoringSystem::Blosum62, GapCosts::from_user(11, 1, GapConvention::Blast).unwrap(), false).unwrap()
    }

    #[test]
    fn database_model_is_unchanged_by_the_subject_length() {
        let st = SearchStatistics::with_model(b62_11_1(), 38_000_000, 91_000, EvalueModel::Database);
        let q = st.for_query(170);
        for n in [30usize, 170, 5000] {
            for raw in [20i64, 60, 200] {
                assert_eq!(q.evalue_for(raw, n), q.evalue(raw));
                assert_eq!(q.space_lower_bound(n), q.search_space);
            }
        }
    }

    #[test]
    fn per_subject_lower_bound_never_exceeds_the_exact_space() {
        let st = SearchStatistics::with_model(b62_11_1(), 38_000_000, 91_000, EvalueModel::PerSubject);
        for m in [1usize, 5, 24, 25, 60, 170, 1000, 35_000] {
            let q = st.for_query(m);
            for n in [1usize, 3, 24, 25, 26, 40, 100, 420, 3000, 35_000] {
                let (lb, ex) = (q.space_lower_bound(n), q.space_for(n));
                assert!(lb <= ex * (1.0 + 1e-12), "m {m} n {n}: bound {lb} > exact {ex}");
                assert!(ex >= st.db_num_seqs as f64, "m {m} n {n}: at least one residue pair per sequence");
            }
        }
    }

    #[test]
    fn per_subject_evalues_grow_with_subject_length_and_match_the_database_model_on_average() {
        // a database of equal-length subjects: both models agree up to the
        // length adjustment (pairwise vs whole database)
        let (z, n) = (10_000u64, 300usize);
        let db = SearchStatistics::with_model(b62_11_1(), z * n as u64, z, EvalueModel::Database).for_query(250);
        let ps = SearchStatistics::with_model(b62_11_1(), z * n as u64, z, EvalueModel::PerSubject).for_query(250);
        let r = ps.evalue_for(80, n) / db.evalue(80);
        assert!((0.5..2.0).contains(&r), "ratio {r}");
        // a score is less significant in a longer subject
        let mut last = 0.0;
        for len in [50usize, 100, 300, 1000, 10_000] {
            let e = ps.evalue_for(80, len);
            assert!(e > last, "length {len}: {e} <= {last}");
            last = e;
        }
        // bit scores do not depend on the model
        assert_eq!(ps.bitscore(80), db.bitscore(80));
    }

    #[test]
    fn saber_default_maps_to_blast_10_1() {
        let g = GapCosts::from_user(11, 1, GapConvention::Saber).unwrap();
        assert_eq!(g.ncbi(), (10, 1));
        let p = gapped_protein(ScoringSystem::Blosum62, g, false).unwrap();
        assert!((p.lambda - 0.243).abs() < 1e-9);
        assert!((p.k - 0.024).abs() < 1e-9);
        assert!(p.calibrated);
    }

    #[test]
    fn blast_convention_11_1_maps_to_ncbi_11_1() {
        let g = GapCosts::from_user(11, 1, GapConvention::Blast).unwrap();
        assert_eq!(g.open, 12);
        assert_eq!(g.ncbi(), (11, 1));
        let p = gapped_protein(ScoringSystem::Blosum62, g, false).unwrap();
        assert!((p.lambda - 0.267).abs() < 1e-9);
        assert!((p.k - 0.041).abs() < 1e-9);
    }

    #[test]
    fn unsupported_gaps_are_rejected_not_substituted() {
        let g = GapCosts::from_user(5, 1, GapConvention::Saber).unwrap();
        let e = gapped_protein(ScoringSystem::Blosum62, g, false).unwrap_err();
        assert!(e.message.contains("Supported"));
        let p = gapped_protein(ScoringSystem::Blosum62, g, true).unwrap();
        assert!(!p.calibrated);
        assert!(p.source.starts_with("UNCALIBRATED"));
    }

    #[test]
    fn gap_cost_conventions() {
        let s = GapCosts::from_user(11, 1, GapConvention::Saber).unwrap();
        assert_eq!(s.cost(1), 11);
        assert_eq!(s.cost(3), 13);
        let b = GapCosts::from_user(11, 1, GapConvention::Blast).unwrap();
        assert_eq!(b.cost(1), 12);
        assert!(GapCosts::from_user(1, 2, GapConvention::Saber).is_err());
        assert!(GapCosts::from_user(5, 0, GapConvention::Saber).is_err());
    }

    #[test]
    fn length_adjustment_matches_blast_scale() {
        // BLOSUM62 11/1 (alpha 1.9, beta -30), query 300 aa vs a 100M-residue
        // database of 300k sequences: BLAST reports a length adjustment of
        // roughly 100-130 for such searches.
        let p = gapped_protein(ScoringSystem::Blosum62, GapCosts { open: 12, extend: 1 }, false).unwrap();
        let (adj, conv) = compute_length_adjustment(p.k, p.log_k, p.alpha_d_lambda, p.beta, 300, 100_000_000, 300_000);
        assert!(conv);
        assert!(adj > 80 && adj < 160, "adj = {adj}");
        // Larger database -> larger adjustment.
        let (adj2, _) = compute_length_adjustment(p.k, p.log_k, p.alpha_d_lambda, p.beta, 300, 1_000_000_000, 3_000_000);
        assert!(adj2 > adj);
    }

    #[test]
    fn min_score_is_exact_inverse() {
        let p = gapped_protein(ScoringSystem::Blosum62, GapCosts { open: 11, extend: 1 }, false).unwrap();
        let st = SearchStatistics::new(p, 2_000_000, 11_000);
        let q = st.for_query(170);
        for cut in [10.0, 1.0, 1e-3, 1e-10] {
            let s = q.min_score(cut, None);
            assert!(q.evalue(s) <= cut);
            assert!(q.evalue(s - 1) > cut);
        }
        let s = q.min_score(10.0, Some(50.0));
        assert!(q.bitscore(s) >= 50.0 && q.bitscore(s - 1) < 50.0 || q.evalue(s - 1) > 10.0);
    }

    #[test]
    fn nucleotide_even_score_rule() {
        let g = GapCosts::from_user(7, 2, GapConvention::Saber).unwrap(); // BLAST 5/2
        let p = gapped_nucleotide(2, -3, g, false).unwrap();
        assert!(p.even_scores_only);
        let st = SearchStatistics::new(p, 1_000_000, 100);
        let q = st.for_query(500);
        assert_eq!(q.evalue(41), q.evalue(40));
    }
}
