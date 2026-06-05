//! SABER: Smith-Waterman homology search with BLAST-like output.
//!
//! See `bin/main.rs` for the CLI. Programmatic users typically want
//! [`algorithm::SmithWaterman`] together with [`scoring::ScoreMatrix`] and
//! [`scoring::EValue`].

pub mod algorithm;
pub mod db;
pub mod doctor;
pub mod io;
pub mod scoring;
pub mod utils;

// GPU module is always present; the CUDA backend inside it is feature-gated.
// The CPU fallback inside this module is what runs when --gpu is requested
// without `--features gpu-cuda` compiled in.
pub mod gpu;

pub use algorithm::{
    AlignmentAlgorithm, AlignmentResult, AlignmentSet, KmerSize, SequenceMode, SmithWaterman,
    TranslationTable,
};
pub use db::{Database, DatabaseIndex};
pub use io::{
    BAMWriter, CSVWriter, Fasta, FastaRecord, JSONWriter, M0Writer, M8Writer, M9Writer,
    OutputFormat, SAMWriter, Writer,
};
pub use scoring::{BitScore, EValue, ScoreMatrix, ScoringSystem};
pub use utils::{ProgressBar, Statistics};

pub const VERSION: &str = "1.0.0";

/// The SABER wordmark: "SABER" in figlet line-art with a sword blade threaded
/// through the centre line (pommel ◉ and crossguard ╪ at the hilt, point ▶ at
/// the tip). Rendered at the top of `--help` and above `--doctor` output.
/// Version-free so it never goes stale; [`banner`] appends the version tagline.
pub const SABER_WORDMARK: &str = r#"    ____    _    ____  _____ ____  
  ╿/ ___|  / \  | __ )| ____|  _ \ 
◉━╪━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━▶
  ╽ ___) / ___ \| |_) | |___|  _ < 
   |____/_/   \_\____/|_____|_| \_\"#;

/// Full banner = wordmark + a centred tagline carrying the live version.
/// Returned as an owned String so the version is always current.
pub fn banner() -> String {
    let width = SABER_WORDMARK.lines().map(|l| l.chars().count()).max().unwrap_or(0);
    let t1 = format!("v{}  \u{b7}  Smith\u{2013}Waterman protein homology", VERSION);
    let t2 = String::from("fast  \u{b7}  sensitive  \u{b7}  low-memory");
    let centre = |s: &str| -> String {
        let n = s.chars().count();
        if n >= width { s.to_string() } else { format!("{}{}", " ".repeat((width - n) / 2), s) }
    };
    format!("\n{}\n\n{}\n{}\n", SABER_WORDMARK, centre(&t1), centre(&t2))
}

pub const AUTHOR: &str = "Richard Allen White III";
