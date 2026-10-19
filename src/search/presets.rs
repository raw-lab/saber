//! Named sensitivity presets (`--preset`). A preset fixes the seeding
//! configuration (stored in the index) and the hit-processing thresholds;
//! every value can still be overridden by its individual flag.
//!
//! The presets were chosen on the *tune* splits (superfamilies with even MD5
//! parity; see bench/README.md) of two benchmarks — SCOP40 all-against-all,
//! and SCOP40 against SCOP40 plus 80,000 reversed protein decoys (the
//! at-scale benchmark) — and evaluated on the held-out *test* splits.
//! Reduced-alphabet spaced seeds with a window identity filter follow the
//! approach DIAMOND uses for its sensitive modes; `default` uses DIAMOND's
//! very-sensitive seed set.

/// One preset.
#[derive(Debug, Clone, Copy)]
pub struct Preset {
    pub name: &'static str,
    pub description: &'static str,
    /// Seed alphabet (`full`, `murphy10`, `custom:...`).
    pub alphabet: &'static str,
    /// Seed shapes (`1` = care position, `0` = wildcard).
    pub shapes: &'static [&'static str],
    /// Window identity filter threshold (0 = off), in units of
    /// seed-alphabet identities + `exact_weight` × exact identities.
    pub hit_filter_ids: u32,
    /// Weight of exact (20-letter) identities in the window filter.
    pub exact_weight: u32,
    /// Minimum ungapped HSP raw score for a candidate region.
    pub hsp_threshold: i32,
    /// Candidate regions kept per query.
    pub max_candidates: usize,
    /// Region threshold as a fraction of the reporting score floor (0 = off);
    /// see `SearchParams::hsp_floor_fraction`.
    pub hsp_floor_fraction: f64,
    /// Score candidate regions inside their HSP diagonals ± this many
    /// diagonals (0 = the whole subject window).
    pub score_band: i64,
    /// Two-stage banded scoring: first pass inside ± this many diagonals
    /// (0 = one pass); see `SearchParams::score_prefilter_band`.
    pub prefilter_band: i64,
    /// Full Smith-Waterman against every subject (no seeding).
    pub exhaustive: bool,
}

/// The 16 weight-7 spaced seeds of DIAMOND's `--very-sensitive` mode (over
/// the Murphy-10 alphabet here), in DIAMOND's order; the faster presets use
/// the first 4 / 8 of them.
pub const DS16: &[&str] = &[
    "11110111",
    "11100100111",
    "110010101011",
    "11010001001011",
    "10101100001101",
    "110100100100011",
    "1010010100010011",
    "1100101000001011",
    "11100000100010101",
    "11000100010010011",
    "11010000001000111",
    "110001001000010011",
    "1010001000100001011",
    "1100010100000010011",
    "1100100000101000011",
    "1101000010000001011",
];

/// Eight weight-6 spaced seeds over Murphy-10, chosen greedily for coverage of
/// the true homologous pairs of the SCOP40 tune split (seed simulator).
pub const W6X8: &[&str] = &["100111011", "1011100011", "1010011011", "11110101", "1101000111", "1111100001", "10101111", "1010010111"];

/// `DS16` followed by `W6X8` (the `more-sensitive` preset).
pub const DS16_W6X8: &[&str] = &[
    "11110111",
    "11100100111",
    "110010101011",
    "11010001001011",
    "10101100001101",
    "110100100100011",
    "1010010100010011",
    "1100101000001011",
    "11100000100010101",
    "11000100010010011",
    "11010000001000111",
    "110001001000010011",
    "1010001000100001011",
    "1100010100000010011",
    "1100100000101000011",
    "1101000010000001011",
    "100111011",
    "1011100011",
    "1010011011",
    "11110101",
    "1101000111",
    "1111100001",
    "10101111",
    "1010010111",
];

pub const PRESETS: &[Preset] = &[
    Preset {
        name: "very-fast",
        description: "Murphy-10, 4 weight-7 spaced seeds, window filter 28 (identities + 2 x exact identities)",
        alphabet: "murphy10",
        shapes: &["11110111", "11100100111", "110010101011", "11010001001011"],
        hit_filter_ids: 28,
        exact_weight: 2,
        hsp_threshold: 40,
        max_candidates: 5000,
        hsp_floor_fraction: 0.7,
        score_band: 48,
        prefilter_band: 0,
        exhaustive: false,
    },
    Preset {
        name: "fast",
        description: "Murphy-10, 8 weight-7 spaced seeds, window filter 27 (identities + 2 x exact identities)",
        alphabet: "murphy10",
        shapes: &["11110111", "11100100111", "110010101011", "11010001001011", "10101100001101", "110100100100011", "1010010100010011", "1100101000001011"],
        hit_filter_ids: 27,
        exact_weight: 2,
        hsp_threshold: 40,
        max_candidates: 5000,
        hsp_floor_fraction: 0.7,
        score_band: 48,
        prefilter_band: 0,
        exhaustive: false,
    },
    Preset {
        name: "default",
        description: "Murphy-10, 16 weight-7 spaced seeds (DIAMOND's very-sensitive set), window filter 27 (identities + 2 x exact identities)",
        alphabet: "murphy10",
        shapes: DS16,
        hit_filter_ids: 27,
        exact_weight: 2,
        hsp_threshold: 40,
        max_candidates: 5000,
        hsp_floor_fraction: 0.65,
        score_band: 48,
        prefilter_band: 16,
        exhaustive: false,
    },
    Preset {
        name: "sensitive",
        description: "Murphy-10, 8 weight-6 spaced seeds, window filter 27 (identities + 2 x exact identities)",
        alphabet: "murphy10",
        shapes: W6X8,
        hit_filter_ids: 27,
        exact_weight: 2,
        hsp_threshold: 40,
        max_candidates: 5000,
        hsp_floor_fraction: 0.65,
        score_band: 48,
        prefilter_band: 16,
        exhaustive: false,
    },
    Preset {
        name: "more-sensitive",
        description: "Murphy-10, the 16 weight-7 and 8 weight-6 spaced seeds, window filter 27 (identities + 2 x exact identities)",
        alphabet: "murphy10",
        shapes: DS16_W6X8,
        hit_filter_ids: 27,
        exact_weight: 2,
        hsp_threshold: 40,
        max_candidates: 8000,
        hsp_floor_fraction: 0.65,
        score_band: 48,
        prefilter_band: 16,
        exhaustive: false,
    },
    Preset {
        name: "legacy",
        description: "SABER 1.x seeding: contiguous full-alphabet 4-mers, no window filter (now bounded, masked and composition-corrected)",
        alphabet: "full",
        shapes: &["1111"],
        hit_filter_ids: 0,
        exact_weight: 0,
        hsp_threshold: 40,
        max_candidates: 5000,
        hsp_floor_fraction: 0.0,
        score_band: 0,
        prefilter_band: 0,
        exhaustive: false,
    },
    Preset {
        name: "exhaustive",
        description: "full Smith-Waterman against every subject (reference mode; slow)",
        alphabet: "full",
        shapes: &["1111"],
        hit_filter_ids: 0,
        exact_weight: 0,
        hsp_threshold: 0,
        max_candidates: usize::MAX,
        hsp_floor_fraction: 0.0,
        score_band: 0,
        prefilter_band: 0,
        exhaustive: true,
    },
];

/// Look up a preset by name (case-insensitive; `very-sensitive` and
/// `ultra-sensitive` are accepted as aliases of `more-sensitive`, `legacy` is
/// SABER 1.x's contiguous full-alphabet 4-mer seeding).
pub fn find(name: &str) -> Result<&'static Preset, String> {
    let n = name.to_ascii_lowercase();
    let n = match n.as_str() {
        "very-sensitive" | "ultra-sensitive" => "more-sensitive",
        "full" | "no-index" => "exhaustive",
        other => other,
    };
    PRESETS.iter().find(|p| p.name == n).ok_or_else(|| {
        format!("unknown preset '{name}' (choose from: {})", PRESETS.iter().map(|p| p.name).collect::<Vec<_>>().join(", "))
    })
}
