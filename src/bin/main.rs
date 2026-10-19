//! SABER command-line entry point (v1.3).
//!
//! Every search goes through one engine ([`saber_rs::search::engine`]):
//! partitioned seed index → bounded seed processing → ungapped extension →
//! exact SIMD scoring of candidate windows → progressive banded traceback,
//! with gapped Karlin-Altschul statistics. The CLI only translates flags into
//! a [`DriverConfig`]; all v1.0/v1.1 flags are still accepted (see
//! docs/MIGRATION.md for the few whose meaning changed).

use clap::parser::ValueSource;
use clap::{ArgMatches, CommandFactory, FromArgMatches, Parser};
use std::path::PathBuf;

use saber_rs::algorithm::{SequenceMode, TranslationTable};
use saber_rs::alphabet::{SeedAlphabet, SeedShape};
use saber_rs::masking::MaskingMode;
use saber_rs::scoring::stats::{GapConvention, GapCosts};
use saber_rs::scoring::ScoringSystem;
use saber_rs::search::driver::{self, DriverConfig};
use saber_rs::search::engine::SearchParams;
use saber_rs::search::index::SeedConfig;
use saber_rs::search::presets;
use saber_rs::OutputFormat;

#[derive(Parser, Debug)]
#[command(name = "saber", version, author, about = "Smith-Waterman homology search with BLAST-like output")]
struct Args {
    /// Query FASTA. Required for a search; omit only with --makedb.
    #[arg(short, long, value_name = "FILE")]
    query: Option<PathBuf>,

    /// Database FASTA. Required to build an index (a plain search, or
    /// --makedb). Omit when searching against a prebuilt --db.
    #[arg(short = 'd', long, value_name = "FILE")]
    database: Option<PathBuf>,

    /// Search mode: pp (blastp) | nn (blastn) | nx (blastx) | pn (tblastn)
    #[arg(long, default_value = "pp")]
    mode: String,

    /// Sensitivity preset: very-fast | fast | default | sensitive |
    /// more-sensitive | legacy | exhaustive. Individual flags override the
    /// preset's values.
    #[arg(long, value_name = "NAME")]
    preset: Option<String>,

    /// Protein scoring matrix (BLOSUM45/50/62/80/90, PAM30/70/250)
    #[arg(short, long, default_value = "BLOSUM62")]
    matrix: String,

    /// Nucleotide match score (nn only)
    #[arg(long, default_value_t = 2)]
    r#match: i32,

    /// Nucleotide mismatch penalty (nn only)
    #[arg(long, default_value_t = -3, allow_negative_numbers = true)]
    mismatch: i32,

    /// Gap opening penalty (positive). Meaning depends on --gap-convention.
    /// Default 11 (protein) or 7 (nn).
    #[arg(short, long, default_value_t = 11)]
    gap_open: i32,

    /// Gap extension penalty (positive). Default 1 (protein) or 2 (nn).
    #[arg(short = 'e', long, default_value_t = 1)]
    gap_extend: i32,

    /// How -g/-e are interpreted. saber (default, as in v1.0/v1.1): a gap of
    /// length k costs open + (k-1)*extend. blast: open + k*extend (NCBI
    /// BLAST/DIAMOND convention; BLAST 11/1 = SABER 12/1).
    #[arg(long, default_value = "saber", value_name = "saber|blast")]
    gap_convention: String,

    /// Accept scoring parameters without published gapped Karlin-Altschul
    /// constants (E-values are then approximate and marked UNCALIBRATED).
    #[arg(long)]
    allow_uncalibrated_evalues: bool,

    /// Output file (stdout if omitted)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output format: m8|m9|m0|sam|json|csv
    #[arg(long, default_value = "m8")]
    outfmt: String,

    /// E-value threshold (alignments with E > this are dropped)
    #[arg(short = 'E', long, default_value_t = 10.0)]
    evalue: f64,

    /// E-value model. per-subject: the search space of the query and the hit's
    /// subject (BLAST length adjustment for the pair) times the number of
    /// database sequences, so a score is less significant in a long subject
    /// than in a short one (FASTA/SSEARCH, HMMER per-sequence E-values).
    /// database: one search space for the whole database (BLAST; the 1.2
    /// behaviour). Bit scores are the same under both.
    #[arg(long, default_value = "per-subject", value_name = "per-subject|database")]
    evalue_model: String,

    /// Minimum bit-score (applied together with the E-value cutoff)
    #[arg(long)]
    min_bitscore: Option<f64>,

    /// Maximum subjects reported per query
    #[arg(short = 'k', long, default_value_t = 10)]
    max_aligns: usize,

    /// Maximum alignments (HSPs) reported per subject
    #[arg(long, default_value_t = 1)]
    max_hsps: usize,

    /// Threads ("auto" or a number)
    #[arg(short, long, default_value = "auto")]
    threads: String,

    /// Genetic code, by name or NCBI number: standard (1), vertebrate (2), yeast (3),
    /// mold (4), invertebrate (5), ciliate (6), echinoderm (9), bacterial (11),
    /// archaeal (= 11), alt-yeast (12), blepharisma (15)
    #[arg(long, default_value = "standard")]
    genetic_code: String,

    /// Skip queries shorter than this
    #[arg(long)]
    min_length: Option<usize>,

    /// Skip queries longer than this
    #[arg(long)]
    max_length: Option<usize>,

    /// Minimum percent identity (0-100)
    #[arg(long, default_value_t = 0.0)]
    percent_identity: f64,

    /// Minimum query coverage of an alignment (percent of query length)
    #[arg(long, default_value_t = 0.0)]
    query_cover: f64,

    /// Minimum subject coverage of an alignment (percent of subject length)
    #[arg(long, default_value_t = 0.0)]
    subject_cover: f64,

    /// Score candidate windows on a GPU (requires --features gpu-cuda and a
    /// CUDA device; experimental). Falls back to the CPU otherwise.
    #[arg(long)]
    gpu: bool,

    /// [deprecated, ignored] v1.0 k-mer prefilter; the seed index replaces it.
    #[arg(long, hide = true)]
    prefilter: bool,
    /// [deprecated, ignored]
    #[arg(long, default_value_t = 0, hide = true)]
    prefilter_k: usize,
    /// [deprecated, ignored]
    #[arg(long, default_value_t = 0, hide = true)]
    prefilter_min_hits: u32,
    /// [deprecated, ignored]
    #[arg(long, default_value_t = 10000, hide = true)]
    prefilter_top_n: usize,

    /// Full Smith-Waterman against every subject (no seed index); same as
    /// --exhaustive / --preset exhaustive.
    #[arg(long)]
    no_index: bool,

    /// Full Smith-Waterman against every subject (reference mode).
    #[arg(long)]
    exhaustive: bool,

    /// Contiguous seed length (v1.0 --index-k). Sets --seed-shape to k ones
    /// and, unless --seed-alphabet is given, the full 20-letter alphabet.
    /// Nucleotide default 11.
    #[arg(long, default_value_t = 0)]
    index_k: usize,

    /// Seed alphabet: full | murphy10 | custom:<groups> (e.g.
    /// "custom:A KR EDNQ C G H ILVM FYW P ST"). Stored in the index.
    #[arg(long, value_name = "ALPHABET")]
    seed_alphabet: Option<String>,

    /// Comma-separated seed shapes, 1 = care position, 0 = wildcard
    /// (e.g. 110111 or 1101011,111011). Stored in the index.
    #[arg(long, value_name = "SHAPES")]
    seed_shape: Option<String>,

    /// Skip query seeds occurring more than N times in the database
    /// (0 = automatic: 64x the mean bucket size, at least 256).
    #[arg(long, default_value_t = 0)]
    seed_max_occ: usize,

    /// Window identity filter: a seed hit is extended only if a 32-residue
    /// window around it reaches N = seed-alphabet identities + the exact
    /// identity weight x exact identities (0 = off). Preset-dependent
    /// default (27 or 28 with weight 2; 1.2's presets: 11-13 with weight 0).
    #[arg(long)]
    hit_filter_ids: Option<u32>,

    /// Weight of exact (20-letter) identities in the window filter value.
    /// Preset-dependent (2; 0 for legacy).
    #[arg(long)]
    hit_filter_exact_weight: Option<u32>,

    /// Minimum seed hits on one diagonal (v1.0). 1 = single-hit
    /// triggering; >= 2 enables two-hit triggering within --two-hit-window.
    #[arg(long, default_value_t = 1)]
    index_min_hits: u32,

    /// [accepted, no effect] v1.0 diagonal bucket size; v1.2 tracks exact
    /// diagonals.
    #[arg(long, default_value_t = 16, hide = true)]
    index_diag_bucket: i32,

    /// [accepted, no effect]
    #[arg(long, hide = true)]
    keep_ascii: bool,

    /// Database minimizer window (1 = index every seed position; larger =
    /// smaller index, lower sensitivity). 0 = automatic (1).
    #[arg(long, default_value_t = 0)]
    minimizer_window: usize,

    /// Build a persistent on-disk index from --database, write it to this
    /// path and exit. The index is written partition by partition (bounded
    /// memory) and atomically renamed into place when complete.
    #[arg(long, value_name = "PATH")]
    makedb: Option<String>,

    /// Search a persistent index built with --makedb (memory-mapped; the
    /// --database FASTA is not needed).
    #[arg(long, value_name = "PATH")]
    db: Option<String>,

    /// Verify index checksums and residue codes when opening --db (reads
    /// the whole file once).
    #[arg(long)]
    verify_index: bool,

    /// Run environment and self-integrity diagnostics, print a report, and
    /// exit (non-zero if any check fails).
    #[arg(long, default_value_t = false)]
    doctor: bool,

    /// Candidate regions kept per query (heuristic limit; drops are counted
    /// and reported). Preset-dependent default (5000).
    #[arg(long)]
    index_max_candidates: Option<usize>,

    /// Neighborhood-word threshold T for full-alphabet seeds (BLAST-style:
    /// each query word is expanded to all words scoring >= T). 0 = off.
    /// Implies the full alphabet, contiguous --index-k seeds (default 4) and
    /// two-hit triggering within --two-hit-window, unless set explicitly.
    #[arg(long, default_value_t = 0)]
    neighborhood_t: i32,

    /// Alias for --preset sensitive. (v1.0's neighborhood seeding is
    /// available as --neighborhood-t 13.)
    #[arg(long, default_value_t = false)]
    sensitive: bool,

    /// Two-hit window (used with --index-min-hits >= 2 or --neighborhood-t).
    #[arg(long, default_value_t = 40)]
    two_hit_window: u32,

    /// Subject window slack when --score-band is 0: each candidate is
    /// scored against the subject segment spanned by its seed diagonals plus
    /// this many residues on each side. 0 = whole subject.
    #[arg(long, default_value_t = 64)]
    band: u32,

    /// Minimum (composition-corrected) ungapped HSP score for a candidate.
    #[arg(long)]
    hsp_threshold: Option<i32>,

    /// Region threshold: a candidate region is kept only if one of its
    /// ungapped segments, or the sum of two or more, reaches this fraction of
    /// the reporting raw-score floor at -E or at E = 10, whichever is lower
    /// (0 = off; preset-dependent, advanced).
    #[arg(long, hide = true)]
    hsp_floor_fraction: Option<f64>,

    /// Score each candidate region inside its seeded diagonals ± N (the
    /// diagonals of the subject's strong regions within 1,024 of each other
    /// are joined, so a long indel between two seeded parts stays inside).
    /// 0 = the subject window of --band. Preset-dependent (48; 0 for legacy).
    #[arg(long)]
    score_band: Option<u32>,

    /// Two-stage banded scoring: narrow band half-width of the first pass
    /// (0 = off; preset-dependent) and the fraction of the reporting floor a
    /// region must reach in it to be scored in the full band.
    #[arg(long, hide = true)]
    score_prefilter_band: Option<u32>,
    #[arg(long, default_value_t = 0.75, hide = true)]
    score_prefilter_frac: f64,

    /// X-drop for ungapped extension.
    #[arg(long, default_value_t = 20)]
    x_drop: i32,

    /// Seed masking (database at index build, query at search):
    /// tantan | none. Masking affects seeding only, never alignment scores.
    #[arg(long, default_value = "tantan", value_name = "tantan|none")]
    masking: String,

    /// Composition-based score correction (protein): on | off.
    #[arg(long, default_value = "on", value_name = "on|off")]
    comp_bias: String,

    /// Approximate memory budget, e.g. 4G or 512M (default: 3/4 of RAM).
    /// Bounds query batches and makes an in-memory index build spill to a
    /// temporary on-disk index (--tmpdir) when it would not fit.
    #[arg(long, value_name = "SIZE")]
    memory_budget: Option<String>,

    /// Directory for temporary on-disk indexes.
    #[arg(long, value_name = "DIR")]
    tmpdir: Option<PathBuf>,

    /// Queries per batch (0 = automatic). Results are written in input order.
    #[arg(long, default_value_t = 0)]
    query_batch: usize,

    /// Batch scheduling: auto | scan (seed tables of a whole query batch,
    /// the database streamed once per batch; protein searches) | query (each
    /// query looks its seeds up in the index) | partition (each index
    /// partition is scanned once per query batch). Results are identical.
    #[arg(long, default_value = "auto", value_name = "auto|scan|query|partition")]
    batch_order: String,

    /// Also store seed postings in the index (--makedb, or an in-memory
    /// build), so that it can be searched with --batch-order query|partition.
    /// The default scan schedule needs none.
    #[arg(long)]
    index_postings: bool,

    /// Query residues per scan batch (0 = automatic; advanced).
    #[arg(long, default_value_t = 0, hide = true)]
    scan_batch_residues: usize,

    /// Survivors above which one (query, subject) pair is scanned and visited
    /// on its own, in query-position slices (0 = automatic; advanced).
    #[arg(long, default_value_t = 0, hide = true)]
    scan_pair_survivors: usize,

    /// Write a JSON run summary (timings, peak RSS, work counters,
    /// heuristic limits that were hit) to this file.
    #[arg(long, value_name = "FILE")]
    stats: Option<PathBuf>,

    /// Residues per index partition (advanced).
    #[arg(long, default_value_t = 64 << 20, hide = true)]
    partition_residues: usize,

    /// Seed-hit buffer per worker thread, in hits (8 bytes each; advanced).
    /// Default: derived from --memory-budget (1/16 of it across threads),
    /// clamped to [65536, 262144].
    #[arg(long, hide = true)]
    hit_buffer: Option<usize>,

    /// Memory cap per traceback direction matrix, in MB (advanced). Beyond it
    /// the best banded alignment found so far is reported with exact=false
    /// and counted as band_capped.
    #[arg(long, default_value_t = 32, hide = true)]
    traceback_max_mb: usize,

    /// Initial traceback band pad (advanced; the band widens automatically).
    #[arg(long, default_value_t = 16, hide = true)]
    traceback_band: i64,

    #[arg(short, long)]
    verbose: bool,

    /// Quiet mode (errors only)
    #[arg(long)]
    quiet: bool,
}

fn explicit(m: &ArgMatches, id: &str) -> bool {
    m.value_source(id) == Some(ValueSource::CommandLine)
}

fn main() {
    let cmd = Args::command().before_help(saber_rs::banner());
    let matches = cmd.get_matches();
    let args = match Args::from_arg_matches(&matches) {
        Ok(a) => a,
        Err(e) => e.exit(),
    };
    if let Err(e) = real_main(args, &matches) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn real_main(args: Args, m: &ArgMatches) -> Result<(), String> {
    let log_level = if args.verbose {
        log::LevelFilter::Debug
    } else if args.quiet {
        log::LevelFilter::Error
    } else {
        log::LevelFilter::Warn
    };
    env_logger::Builder::from_default_env().filter_level(log_level).init();

    let threads = if args.threads.eq_ignore_ascii_case("auto") {
        num_cpus::get()
    } else {
        args.threads.parse::<usize>().map_err(|_| "invalid --threads")?.max(1)
    };
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().ok();

    if args.doctor {
        use std::path::Path;
        let opts = saber_rs::doctor::DoctorOptions {
            query: args.query.as_deref(),
            database: args.database.as_deref(),
            db_index: args.db.as_deref().map(Path::new),
            gpu_requested: args.gpu,
            cuda_compiled: cfg!(feature = "gpu-cuda"),
            threads,
        };
        let report = saber_rs::doctor::run(&opts);
        println!("{}", saber_rs::banner());
        print!("{}", report.render());
        if report.healthy() {
            return Ok(());
        }
        std::process::exit(1);
    }

    let mode: SequenceMode = args.mode.parse().map_err(|e: String| format!("--mode: {e}"))?;
    let outfmt: OutputFormat = args.outfmt.parse().map_err(|e: String| format!("--outfmt: {e}"))?;
    let genetic_code: TranslationTable = args.genetic_code.parse().map_err(|e: String| format!("--genetic-code: {e}"))?;
    let matrix: ScoringSystem = args.matrix.parse().map_err(|e: String| format!("--matrix: {e}"))?;
    let quiet = args.quiet;
    let note = |s: String| {
        if !quiet {
            eprintln!("[saber] {s}");
        }
    };

    // ---- deprecated / no-op flags ----
    if outfmt == OutputFormat::BAM {
        note("--outfmt bam writes SAM text (BAM encoding is not implemented); convert with `samtools view -b`".to_string());
    }
    for (id, name) in [("prefilter", "--prefilter"), ("prefilter_k", "--prefilter-k"), ("prefilter_min_hits", "--prefilter-min-hits"), ("prefilter_top_n", "--prefilter-top-n")] {
        if explicit(m, id) {
            note(format!("{name} is deprecated and ignored: the seed index replaces the v1.0 prefilter"));
        }
    }
    if explicit(m, "index_diag_bucket") {
        note("--index-diag-bucket is accepted for compatibility and has no effect (exact diagonals are tracked)".into());
    }

    // ---- preset ----
    let preset_name = match (&args.preset, args.sensitive, args.no_index || args.exhaustive) {
        (_, _, true) => "exhaustive".to_string(),
        (Some(p), _, _) => p.clone(),
        (None, true, _) => "sensitive".to_string(),
        (None, false, _) => "default".to_string(),
    };
    let preset = presets::find(&preset_name)?;

    // ---- seeding configuration ----
    let neighborhood = args.neighborhood_t > 0;
    let masking: MaskingMode = args.masking.parse().map_err(|e: String| format!("--masking: {e}"))?;
    let mut seed = if mode == SequenceMode::NucleotideToNucleotide {
        SeedConfig::nucleotide(if args.index_k > 0 { args.index_k } else { 11 })
    } else {
        let mut s = SeedConfig::protein_default();
        let alpha = match &args.seed_alphabet {
            Some(a) => a.clone(),
            None if args.index_k > 0 || neighborhood => "full".into(),
            None => preset.alphabet.into(),
        };
        s.alphabet = SeedAlphabet::parse(&alpha)?;
        s.shapes = match &args.seed_shape {
            Some(sh) => sh.split(',').map(|x| SeedShape::parse(x.trim())).collect::<Result<Vec<_>, _>>()?,
            None if args.index_k > 0 || neighborhood => vec![SeedShape::contiguous(if args.index_k > 0 { args.index_k } else { 4 })],
            None => preset.shapes.iter().map(|x| SeedShape::parse(x)).collect::<Result<Vec<_>, _>>()?,
        };
        s.masking = masking;
        s
    };
    seed.minimizer_window = args.minimizer_window.max(1);
    seed.validate()?;

    // ---- search parameters ----
    let mut p = SearchParams::default();
    p.evalue = args.evalue;
    p.min_bitscore = args.min_bitscore;
    p.max_target_seqs = args.max_aligns.max(1);
    p.max_hsps = args.max_hsps.max(1);
    p.min_identity = args.percent_identity;
    p.min_query_cover = args.query_cover;
    p.min_subject_cover = args.subject_cover;
    p.hsp_threshold = args.hsp_threshold.unwrap_or(preset.hsp_threshold);
    p.x_drop_ungapped = args.x_drop;
    p.hsp_floor_fraction = args.hsp_floor_fraction.unwrap_or(preset.hsp_floor_fraction).clamp(0.0, 1.0);
    p.hit_filter_min_ids = match args.hit_filter_ids {
        Some(v) => v,
        None if args.seed_alphabet.is_some() || args.seed_shape.is_some() || args.index_k > 0 || neighborhood => 0,
        None => preset.hit_filter_ids,
    };
    if mode == SequenceMode::NucleotideToNucleotide && args.hit_filter_ids.is_none() {
        p.hit_filter_min_ids = 0;
    }
    p.hit_filter_exact_weight = args.hit_filter_exact_weight.unwrap_or(preset.exact_weight);
    p.max_candidates = args.index_max_candidates.unwrap_or(preset.max_candidates).max(1);
    p.seed_max_occ = args.seed_max_occ;
    p.query_masking = masking;
    p.comp_bias = match args.comp_bias.to_ascii_lowercase().as_str() {
        "on" | "1" | "true" | "yes" => true,
        "off" | "0" | "false" | "no" => false,
        o => return Err(format!("--comp-bias: expected on|off, got '{o}'")),
    };
    p.exhaustive = preset.exhaustive;
    p.window_slack = if args.band == 0 { i64::MAX / 4 } else { args.band as i64 };
    p.score_band = args.score_band.map_or(preset.score_band, |v| v as i64);
    p.score_prefilter_band = args.score_prefilter_band.map_or(preset.prefilter_band, |v| v as i64);
    p.score_prefilter_frac = args.score_prefilter_frac.clamp(0.0, 1.0);
    p.scan_pair_survivors = args.scan_pair_survivors;
    p.band_pad = args.traceback_band.max(1);
    p.traceback_max_bytes = args.traceback_max_mb.max(1) << 20;
    p.neighborhood_t = args.neighborhood_t.max(0);
    p.alignment_strings = matches!(outfmt, OutputFormat::BlastM0 | OutputFormat::SAM | OutputFormat::BAM);
    if args.index_min_hits >= 2 || (neighborhood && !explicit(m, "index_min_hits")) {
        p.two_hit_window = args.two_hit_window.max(1);
    }

    // ---- gap costs (nn defaults to 7/2 = BLAST 5/2 for 2/-3 scoring) ----
    let conv: GapConvention = args.gap_convention.parse().map_err(|e: String| format!("--gap-convention: {e}"))?;
    let (go, ge) = if mode == SequenceMode::NucleotideToNucleotide && !explicit(m, "gap_open") && !explicit(m, "gap_extend") {
        match conv {
            GapConvention::Saber => (7, 2),
            GapConvention::Blast => (5, 2),
        }
    } else {
        (args.gap_open, args.gap_extend)
    };
    let gaps = GapCosts::from_user(go, ge, conv)?;

    let memory_budget = match &args.memory_budget {
        Some(s) => Some(driver::parse_size(s)?),
        None => driver::default_memory_budget(),
    };
    p.hit_buffer = match args.hit_buffer {
        Some(h) => h.max(1024),
        // 256 K hits (2 MB) per thread: enough that a typical query needs one
        // or two subject-range passes; larger buffers did not run faster and
        // cost resident memory.
        None => memory_budget.map_or(1 << 18, |b| ((b / 16 / 8 / threads as u64) as usize).clamp(1 << 16, 1 << 18)),
    };
    let cfg = DriverConfig {
        mode,
        query: args.query.clone(),
        database: args.database.clone(),
        db_index: args.db.as_ref().map(PathBuf::from),
        output: args.output.clone(),
        outfmt,
        matrix,
        gaps,
        nt_match: args.r#match,
        nt_mismatch: args.mismatch,
        allow_uncalibrated: args.allow_uncalibrated_evalues,
        evalue_model: saber_rs::scoring::EvalueModel::parse(&args.evalue_model)
            .ok_or_else(|| format!("unknown --evalue-model '{}' (per-subject | database)", args.evalue_model))?,
        genetic_code,
        threads,
        memory_budget,
        params: p,
        seed,
        partition_residues: args.partition_residues.max(1 << 16),
        min_length: args.min_length,
        max_length: args.max_length,
        stats_path: args.stats.clone(),
        verify_index: args.verify_index,
        tmpdir: args.tmpdir.clone().unwrap_or_else(std::env::temp_dir),
        quiet,
        query_batch: args.query_batch,
        preset: preset.name.to_string(),
        gpu: args.gpu,
        batch_order: args.batch_order.parse()?,
        index_postings: args.index_postings,
        scan_batch_residues: args.scan_batch_residues,
    };

    // ---- --makedb ----
    if let Some(out) = &args.makedb {
        if mode == SequenceMode::NucleotideToProtein {
            note("nx searches use a protein index; building a protein (pp) index".into());
        }
        let (bytes, secs) = driver::makedb(&cfg, std::path::Path::new(out))?;
        note(format!("wrote index '{out}' ({:.1} MB) in {:.2}s", bytes as f64 / 1e6, secs));
        return Ok(());
    }

    // ---- search ----
    if cfg.db_index.is_some() && !driver::scan_eligible(&cfg) && (explicit(m, "seed_alphabet") || explicit(m, "seed_shape") || explicit(m, "index_k") || explicit(m, "minimizer_window")) {
        note("seeding flags are fixed when the index is built; the --db index's own seed configuration is used".into());
    }
    let s = driver::run(&cfg)?;
    if !quiet {
        eprintln!(
            "[saber] {} queries, {} alignments | preset {} | index {} ({:.2}s) | search {:.2}s | total {:.2}s | peak RSS {:.1} MB",
            s.queries,
            s.alignments,
            s.preset,
            s.index_source,
            s.index_build_s,
            s.search_s,
            s.total_s,
            s.peak_rss_kb as f64 / 1024.0
        );
        for l in &s.heuristic_limits {
            if args.verbose || !l.starts_with("seed-and-extend heuristic") {
                eprintln!("[saber] limit: {l}");
            }
        }
    }
    Ok(())
}
