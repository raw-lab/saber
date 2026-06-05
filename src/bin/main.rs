//! SABER command-line entry point.
//!
//! Compared to v2.0 this version actually plumbs the user's `--mode` selection
//! through to the aligner, loads the requested scoring matrix, computes real
//! E-values/bit-scores, and runs the query loop in parallel via rayon.

use clap::{CommandFactory, FromArgMatches, Parser};
use log::info;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use saber_rs::{
    algorithm::{
        index::{KmerIndex, DEFAULT_K_NUCLEOTIDE, DEFAULT_K_PROTEIN, DEFAULT_MINIMIZER_WINDOW_PROTEIN},
        neighborhood::NeighborhoodGenerator,
        prefilter::{prefilter as run_prefilter, PrefilterConfig},
        translate_six_frames, AlignmentResult, SequenceMode, SmithWaterman, TranslationTable,
    },
    gpu::{cpu_fallback::CpuExecutor, try_init_gpu, GpuExecutor},
    io::{
        CSVWriter, Fasta, FastaRecord, JSONWriter, M0Writer, M8Writer, M9Writer, SAMWriter, Writer,
    },
    scoring::{EValue, ScoreMatrix, ScoringSystem},
    utils::Timer,
    OutputFormat,
};

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

    /// Protein scoring matrix (BLOSUM45/50/62/80/90, PAM30/70/250)
    #[arg(short, long, default_value = "BLOSUM62")]
    matrix: String,

    /// Nucleotide match score (nn only)
    #[arg(long, default_value_t = 2)]
    r#match: i32,

    /// Nucleotide mismatch penalty (nn only)
    #[arg(long, default_value_t = -3)]
    mismatch: i32,

    /// Gap opening penalty (positive)
    #[arg(short, long, default_value_t = 11)]
    gap_open: i32,

    /// Gap extension penalty (positive)
    #[arg(short = 'e', long, default_value_t = 1)]
    gap_extend: i32,

    /// Output file (stdout if omitted)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output format: m8|m9|m0|sam|json|csv
    #[arg(long, default_value = "m8")]
    outfmt: String,

    /// E-value threshold (alignments with E > this are dropped)
    #[arg(short = 'E', long, default_value_t = 10.0)]
    evalue: f64,

    /// Minimum bit-score; overrides E-value cutoff when set
    #[arg(long)]
    min_bitscore: Option<f64>,

    /// Maximum alignments reported per query
    #[arg(short = 'k', long, default_value_t = 10)]
    max_aligns: usize,

    /// Threads ("auto" or a number)
    #[arg(short, long, default_value = "auto")]
    threads: String,

    /// Genetic code: standard|bacterial|vertebrate|...|11
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

    /// Try GPU if compiled with --features gpu-cuda
    #[arg(long)]
    gpu: bool,

    /// Enable k-mer prefilter to skip obvious non-hits (recommended for large DBs).
    /// Default: off (every subject gets full SW). Turn on with `--prefilter`.
    #[arg(long)]
    prefilter: bool,

    /// Prefilter k-mer length (3-5 for protein, 11-15 for nucleotide).
    #[arg(long, default_value_t = 0)]
    prefilter_k: usize,

    /// Minimum query k-mer hits a subject must have to survive prefiltering.
    #[arg(long, default_value_t = 0)]
    prefilter_min_hits: u32,

    /// Top-N subjects to keep after prefiltering (most-hits-first).
    #[arg(long, default_value_t = 10000)]
    prefilter_top_n: usize,

    /// Disable the inverted k-mer index (default on for pp/nn modes).
    /// With the index off, every subject runs through full SW — much slower
    /// for large DBs but useful for benchmarking the index's contribution.
    #[arg(long)]
    no_index: bool,

    /// Index k-mer length. Defaults: 4 for protein, 11 for nucleotide.
    /// Smaller k → more sensitive at low identity but bigger index and more
    /// candidates to filter. Larger k → faster but less sensitive at low ID.
    #[arg(long, default_value_t = 0)]
    index_k: usize,

    /// Minimum clustered seed hits on a single diagonal for a subject to
    /// become a candidate. Default 1 (preserves full sensitivity — matches
    /// `--no-index` recall). Higher = faster but may miss remote homologs.
    /// Used when neighborhood expansion is OFF; for neighborhood mode, the
    /// `--two-hit-window` parameter controls selectivity instead.
    #[arg(long, default_value_t = 1)]
    index_min_hits: u32,

    /// Diagonal bucket size for seed clustering. Hits with diagonals within
    /// this many cells are considered on the same alignment. Default 16.
    #[arg(long, default_value_t = 16)]
    index_diag_bucket: i32,

    /// (Debug) Skip the post-index ASCII-drop step. Default off; the drop
    /// saves ~9 MB on bundled bench but causes the traceback to reconstruct
    /// ASCII from the index. Useful for isolating the regression cost.
    #[arg(long, hide = true)]
    keep_ascii: bool,

    /// Minimizer window for sparse subject indexing. The subject side indexes
    /// only the minimizer k-mer of each window-length window, cutting index
    /// memory by ~2/(window+1). 1 = dense (index every k-mer, lowest memory
    /// savings, highest sensitivity). Larger windows = less memory, slightly
    /// lower sensitivity. The query side always looks up every k-mer, so the
    /// sensitivity cost is small. Default 1 for nucleotide; protein default
    /// is set in code after CLI parse.
    #[arg(long, default_value_t = 0)]
    minimizer_window: usize,

    /// Build a persistent on-disk index from the --database FASTA, write it
    /// to this path, and exit (no search). Reuse with --db to skip the
    /// in-memory index build on every run and to memory-map the index --
    /// essential for UniRef90-scale databases that do not fit in RAM.
    /// Protein mode; honors --index-k and --minimizer-window.
    #[arg(long, value_name = "PATH")]
    makedb: Option<String>,

    /// Use a persistent index previously built with --makedb. The index is
    /// memory-mapped (not copied into RAM), so peak memory includes only the
    /// pages actually touched. The --database FASTA is not needed with --db
    /// (subject IDs and lengths are stored in the index).
    #[arg(long, value_name = "PATH")]
    db: Option<String>,

    /// Run environment and self-integrity diagnostics, print a report, and
    /// exit. Verifies CPU SIMD features, the Smith-Waterman kernel, the
    /// indexed search path, persistent-index save/open, scoring matrices,
    /// GPU availability, and temp-dir writability. If --query / --database /
    /// --db are also given, those files are validated too. Exits non-zero if
    /// any check fails.
    #[arg(long, default_value_t = false)]
    doctor: bool,


    /// Maximum candidates to forward from index to SIMD SW. Default 5000.
    /// On a very large DB you may want to raise this; on a small DB lowering
    /// it makes no difference.
    #[arg(long, default_value_t = 5000)]
    index_max_candidates: usize,

    /// Neighborhood-word score threshold T. 0 = exact-match-only k-mer
    /// lookup (default — on most databases the neighborhood approach helps
    /// less than the index lookup overhead costs). BLAST blastp default
    /// is 11 for BLOSUM62, SWORD's is 13. Use --neighborhood-t 11 with
    /// `--index-min-hits 2` and `--two-hit-window 40` for BLAST-style
    /// seeding on extremely large, biologically-diverse databases where
    /// remote-homology recall is critical.
    #[arg(long, default_value_t = 0)]
    neighborhood_t: i32,

    /// Convenience preset for maximum sensitivity: enables BLAST-style
    /// neighborhood-word seeding (equivalent to --neighborhood-t 13
    /// --two-hit-window 40). Recovers remote homologs (matches/exceeds SWORD
    /// recall on the bio benchmark) at a substantial memory and time cost —
    /// the expanded seed list is large. Prefer the default seeding unless you
    /// specifically need deep remote-homology recall.
    #[arg(long, default_value_t = false)]
    sensitive: bool,

    /// Window size for BLAST's two-hit heuristic (max distance between two
    /// hits on the same diagonal to count as a "double hit"). Default 40,
    /// matching BLAST's blastp `-window_size`. Only used when neighborhood
    /// expansion is enabled.
    #[arg(long, default_value_t = 40)]
    two_hit_window: u32,

    /// Banded SW slack on each side of the best seed diagonal, in residues.
    /// Subjects in the indexed path are sliced to the diagonal band before
    /// SW, shrinking the DP matrix from `q × subj_len` to `q × 2*band`.
    /// Default 64 — comfortably wider than typical biological indels
    /// (~20 aa) while still cutting work substantially on long subjects.
    /// Set to 0 to disable banding (run full SW on every candidate).
    #[arg(long, default_value_t = 64)]
    band: u32,

    /// Minimum ungapped HSP score (X-drop extension along seed diagonals)
    /// for a candidate to survive to gapped SW. This is BLAST's real filter
    /// between seeding and gapped DP. Default 40 — approximately matches
    /// BLAST's S_min for BLOSUM62 with typical 100-300aa queries, preserves
    /// every significant hit (E ≤ 1) on bundled benchmarks vs `--no-index`.
    /// Tighten (--hsp-threshold 45-50) for ~2× speedup on workloads where
    /// you don't care about very weak distant homology. Set to 0 to disable.
    #[arg(long, default_value_t = 40)]
    hsp_threshold: i32,

    /// X-drop value for ungapped extension (BLOSUM62-tuned default 20,
    /// matching BLAST's `-xdrop_ungap`). Lower = ungapped stops sooner.
    #[arg(long, default_value_t = 20)]
    x_drop: i32,

    #[arg(short, long)]
    verbose: bool,

    /// Quiet mode (errors only)
    #[arg(long)]
    quiet: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Build the clap command and attach the SABER sword banner so it appears
    // at the top of both `-h` and `--help`. We go through CommandFactory rather
    // than Args::parse() purely so the multi-line Unicode banner can live in a
    // normal const. get_matches() still handles --help/--version (printing the
    // banner) and exits on its own.
    let cmd = Args::command().before_help(saber_rs::banner());
    let matches = cmd.get_matches();
    let mut args = Args::from_arg_matches(&matches).map_err(|e| e.to_string())?;

    let log_level = if args.verbose {
        log::LevelFilter::Debug
    } else if args.quiet {
        log::LevelFilter::Error
    } else {
        log::LevelFilter::Info
    };
    env_logger::Builder::from_default_env().filter_level(log_level).init();

    let threads = if args.threads.eq_ignore_ascii_case("auto") {
        num_cpus::get()
    } else {
        args.threads.parse::<usize>().map_err(|_| "invalid --threads")?
    };
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().ok();

    // ---- saber doctor: diagnostics + self-test, then exit ----
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
        } else {
            std::process::exit(1);
        }
    }

    let mode: SequenceMode = args.mode.parse().map_err(|e: String| format!("mode: {}", e))?;

    // --sensitive is a convenience preset: enable neighborhood-word seeding
    // (T=13) unless the user already set an explicit --neighborhood-t.
    if args.sensitive && args.neighborhood_t == 0 {
        args.neighborhood_t = 13;
        info!("--sensitive: enabling neighborhood-word seeding (T=13)");
    }

    let outfmt: OutputFormat = args.outfmt.parse().map_err(|e: String| format!("outfmt: {}", e))?;
    let genetic_code: TranslationTable =
        args.genetic_code.parse().map_err(|e: String| format!("genetic_code: {}", e))?;
    let matrix_system: ScoringSystem = args.matrix.parse().map_err(|e: String| format!("matrix: {}", e))?;

    info!("SABER v{}  mode={}  outfmt={}  threads={}", saber_rs::VERSION, mode, outfmt, threads);

    // ---- Choose execution backend ----
    // - If --gpu is set and we have a CUDA executor compiled in, try to init it.
    // - Otherwise (default), use the CPU fallback which runs SW across all cores via rayon.
    let executor: Arc<dyn GpuExecutor> = if args.gpu {
        match try_init_gpu() {
            Some(gpu) => {
                info!("Using GPU backend: {}", gpu.backend_name());
                gpu
            }
            None => {
                if cfg!(feature = "gpu-cuda") {
                    eprintln!(
                        "[saber] --gpu requested but no CUDA device available; falling back to CPU."
                    );
                } else {
                    eprintln!(
                        "[saber] --gpu requested but binary was not built with --features gpu-cuda; \
                         falling back to multithreaded CPU. Rebuild with `cargo build --release \
                         --features gpu-cuda` on a host with the CUDA toolkit to enable GPU."
                    );
                }
                Arc::new(CpuExecutor::new())
            }
        }
    } else {
        Arc::new(CpuExecutor::new())
    };
    info!("Alignment backend: {}", executor.backend_name());

    // ---- Load sequences ----
    // ============================================================
    //  Index acquisition: three routes.
    //   (A) --makedb PATH : build from --database, save, exit.
    //   (B) --db PATH      : memory-map a prebuilt index; no FASTA DB needed.
    //   (C) plain search   : build the index in memory from --database.
    // ============================================================

    // ---- Route A: --makedb (build + save + exit) ----
    if let Some(makedb_path) = args.makedb.clone() {
        if !matches!(mode, SequenceMode::ProteinToProtein | SequenceMode::NucleotideToNucleotide) {
            return Err("--makedb supports only pp (protein) and nn (nucleotide) modes".into());
        }
        let db_path = args.database.clone()
            .ok_or("--makedb requires --database <FILE> to build the index from")?;
        let mut t = Timer::new("makedb");
        t.start();
        let database = Fasta::from_file(&db_path).map_err(|e| format!("db load: {}", e))?;

        // Resolve the minimizer window the same way the search path does, so a
        // prebuilt index matches what a plain run would have produced.
        if args.minimizer_window == 0 {
            args.minimizer_window = match mode {
                SequenceMode::ProteinToProtein => {
                    let r: usize = database.records().iter().map(|x| x.length).sum();
                    if r < 2_000_000 { 1 } else { DEFAULT_MINIMIZER_WINDOW_PROTEIN }
                }
                _ => 1,
            };
        }

        let subjects: Vec<&[u8]> = database.records().iter().map(|r| r.sequence.as_bytes()).collect();
        let t_idx = Instant::now();
        let mut idx = match mode {
            SequenceMode::NucleotideToNucleotide => {
                let k = if args.index_k > 0 { args.index_k } else { DEFAULT_K_NUCLEOTIDE };
                KmerIndex::build_nucleotide_windowed(&subjects, k, args.minimizer_window)
            }
            _ => {
                let k = if args.index_k > 0 { args.index_k } else { DEFAULT_K_PROTEIN };
                KmerIndex::build_protein_windowed(&subjects, k, args.minimizer_window)
            }
        };
        // Attach subject metadata so the saved index is self-contained
        // (search via --db needs no FASTA).
        let ids: Vec<String> = database.records().iter().map(|r| r.id.clone()).collect();
        let lens: Vec<u32> = database.records().iter().map(|r| r.length as u32).collect();
        idx.set_subject_meta(&ids, &lens);
        let build_dt = t_idx.elapsed();

        idx.save(std::path::Path::new(&makedb_path)).map_err(|e| format!("index save: {}", e))?;
        let sz = std::fs::metadata(&makedb_path).map(|m| m.len()).unwrap_or(0);
        info!(
            "Wrote index '{}': k={}, window={}, {} subjects, {} postings, {:.1} MB on disk, built in {:.3}s",
            makedb_path, idx.k(), args.minimizer_window, idx.num_subjects(),
            idx.num_postings(), sz as f64 / 1e6, build_dt.as_secs_f64(),
        );
        t.stop_and_print();
        return Ok(());
    }

    // ---- Load queries (needed for both search routes) ----
    let mut t = Timer::new("load");
    t.start();
    let query_path = args.query.clone().ok_or("--query <FILE> is required for a search")?;
    let queries = Fasta::from_file(&query_path).map_err(|e| format!("query load: {}", e))?;

    // ---- Route B vs C: prebuilt mmap index, or build in memory ----
    // `database` is Some only on route C (plain search); on route B the FASTA
    // DB is not loaded at all.
    let mut database: Option<Fasta> = None;

    let (kmer_index, db_residues): (Option<Arc<KmerIndex>>, usize) = if let Some(db_path) = args.db.clone() {
        // Route B: memory-map the prebuilt index.
        if !matches!(mode, SequenceMode::ProteinToProtein | SequenceMode::NucleotideToNucleotide) {
            return Err("--db supports only pp (protein) and nn (nucleotide) modes".into());
        }
        let idx = KmerIndex::open(std::path::Path::new(&db_path))
            .map_err(|e| format!("opening --db '{}': {}. (Rebuild with --makedb if the format/endianness differs.)", db_path, e))?;
        if !idx.has_subject_meta() {
            return Err(format!("index '{}' lacks subject metadata; rebuild it with this SABER version", db_path).into());
        }
        let db_residues = idx.total_residues().max(1);
        info!(
            "Mapped index '{}': k={}, {} subjects, {} postings, {:.1} MB virtual (paged on demand)",
            db_path, idx.k(), idx.num_subjects(), idx.num_postings(),
            idx.memory_bytes() as f64 / 1e6,
        );
        (Some(Arc::new(idx)), db_residues)
    } else {
        // Route C: build the index in memory from the FASTA DB.
        let db_path = args.database.clone()
            .ok_or("a search needs either --db <INDEX> or --database <FASTA>")?;
        let db = Fasta::from_file(&db_path).map_err(|e| format!("db load: {}", e))?;

        // Adaptive minimizer-window resolution (see route A for rationale).
        if args.minimizer_window == 0 {
            args.minimizer_window = match mode {
                SequenceMode::ProteinToProtein => {
                    let r: usize = db.records().iter().map(|x| x.length).sum();
                    if r < 2_000_000 { 1 } else { DEFAULT_MINIMIZER_WINDOW_PROTEIN }
                }
                _ => 1,
            };
            info!("Minimizer window (auto): {} ({})", args.minimizer_window,
                  if args.minimizer_window == 1 { "dense" } else { "sparse" });
        }

        let dstats = db.stats();
        let db_residues = dstats.total_residues.max(1);

        let index = if !args.no_index
            && matches!(mode, SequenceMode::ProteinToProtein | SequenceMode::NucleotideToNucleotide)
        {
            let t_idx = Instant::now();
            let subjects: Vec<&[u8]> = db.records().iter().map(|r| r.sequence.as_bytes()).collect();
            let mut idx = match mode {
                SequenceMode::NucleotideToNucleotide => {
                    let k = if args.index_k > 0 { args.index_k } else { DEFAULT_K_NUCLEOTIDE };
                    KmerIndex::build_nucleotide_windowed(&subjects, k, args.minimizer_window)
                }
                _ => {
                    let k = if args.index_k > 0 { args.index_k } else { DEFAULT_K_PROTEIN };
                    KmerIndex::build_protein_windowed(&subjects, k, args.minimizer_window)
                }
            };
            // Attach subject metadata so the search reads id/length from the
            // index and we can drop the FASTA ASCII below.
            let ids: Vec<String> = db.records().iter().map(|r| r.id.clone()).collect();
            let lens: Vec<u32> = db.records().iter().map(|r| r.length as u32).collect();
            idx.set_subject_meta(&ids, &lens);
            let dt = t_idx.elapsed();
            info!("Index built: k={}, {} postings, {:.2} MB, {:.3}s",
                  idx.k(), idx.num_postings(), idx.memory_bytes() as f64 / 1e6, dt.as_secs_f64());
            Some(Arc::new(idx))
        } else {
            None
        };

        // Keep the FASTA only when we have no index (the --no-index / linear
        // path needs the ASCII sequences). With an index, metadata lives in
        // the index, so drop the FASTA entirely to save memory.
        if index.is_none() {
            database = Some(db);
        }
        (index, db_residues)
    };

    let qstats = queries.stats();
    info!("loaded {} queries ({} residues)", qstats.total_sequences, qstats.total_residues);
    t.stop_and_print();

    // ---- Build the aligner and an E-value calculator for this matrix ----
    let aligner = Arc::new(build_aligner(mode, matrix_system, args.gap_open, args.gap_extend,
                                          args.r#match, args.mismatch));
    let ev_params = evalue_params_for(mode, matrix_system);


    // (Memory note: when an index is built, the FASTA DB is dropped at the end
    // of route C above — subject IDs/lengths live in the index, and traceback
    // reconstructs ASCII from the index's encoded buffer. With --db, the FASTA
    // is never loaded at all. Only the --no-index path keeps `database`.)


    // ---- Build neighborhood generator (protein modes only) ----
    // T > 0 enables BLAST-style neighborhood word expansion. Each query
    // k-mer expands to all k-mers scoring ≥ T against it; this lets
    // remote homologs produce dense seed clusters on real diagonals
    // and is what makes index-based search both fast and sensitive.
    let neighborhood: Option<Arc<NeighborhoodGenerator>> = if args.neighborhood_t > 0
        && kmer_index.is_some()
        && matches!(mode, SequenceMode::ProteinToProtein)
    {
        let m = ScoreMatrix::new(matrix_system);
        let k = kmer_index.as_ref().unwrap().k();
        let gen = NeighborhoodGenerator::new(&m, k, args.neighborhood_t);
        info!("Neighborhood expansion: T={}, k={}", args.neighborhood_t, k);
        Some(Arc::new(gen))
    } else {
        None
    };

    // ---- Open output ----
    let out_file: Box<dyn Write + Send> = if let Some(p) = &args.output {
        Box::new(BufWriter::new(File::create(p)?))
    } else {
        Box::new(BufWriter::new(std::io::stdout()))
    };
    let mut writer: Box<dyn Writer> = match outfmt {
        OutputFormat::BlastM8 => Box::new(M8Writer::new(out_file)),
        OutputFormat::BlastM9 => Box::new(M9Writer::new(out_file)),
        OutputFormat::BlastM0 => Box::new(M0Writer::new(out_file)),
        OutputFormat::SAM | OutputFormat::BAM => Box::new(SAMWriter::new(out_file)),
        OutputFormat::JSON => Box::new(JSONWriter::new(out_file)),
        OutputFormat::CSV => Box::new(CSVWriter::new(out_file)),
    };
    writer.write_header()?;

    // ---- Run search ----
    let mut search = Timer::new("search");
    search.start();

    // For pp / nn we just iterate. For nx we 6-frame the query. For pn we 6-frame the database
    // up-front (once) and then treat the search as protein-protein per frame.
    // Subject records for paths that still need FASTA ASCII (translated modes,
    // --no-index). When an index is in use, `database` is None and the indexed
    // search reads subject id/length from the index itself.
    let empty_records: Vec<FastaRecord> = Vec::new();
    let db_records: &[FastaRecord] = database.as_ref().map(|d| d.records()).unwrap_or(&empty_records);

    let translated_db: Option<Vec<TranslatedSubject>> = match mode {
        SequenceMode::ProteinToNucleotide => {
            Some(db_records.iter().flat_map(|rec| {
                translate_six_frames(rec.sequence.as_bytes(), genetic_code)
                    .into_iter()
                    .map(move |tf| TranslatedSubject {
                        id: format!("{}|frame{:+}", rec.id, signed_frame(tf.frame)),
                        seq: tf.protein,
                    })
            }).collect())
        }
        _ => None,
    };

    let mut total_alignments_written = 0usize;
    // For protein-on-protein and nucleotide-on-nucleotide we always route
    // through the executor pipeline (SIMD score-only + CPU traceback on top
    // hits). This is identical in shape to the GPU path; switching executor
    // CUDA-vs-CPU is the only difference. Wall-time speedup vs the v2.1
    // scalar pipeline is ~50-65× on protein data with AVX2.
    let use_simd_pipeline = matches!(
        mode,
        SequenceMode::ProteinToProtein | SequenceMode::NucleotideToNucleotide
    );
    let scoring = if matches!(mode, SequenceMode::NucleotideToNucleotide) {
        saber_rs::algorithm::Scoring::NucleotideSimple {
            match_score: args.r#match,
            mismatch: args.mismatch,
        }
    } else {
        let m = Arc::new(ScoreMatrix::new(matrix_system));
        saber_rs::algorithm::Scoring::Matrix(m)
    };

    for query in queries.records() {
        if let Some(min) = args.min_length { if query.length < min { continue; } }
        if let Some(max) = args.max_length { if query.length > max { continue; } }

        // Compute and emit alignments for this query. Each branch builds a Vec<AlignmentResult>
        // sorted by score descending, capped at max_aligns.
        let aligns: Vec<AlignmentResult> = match mode {
            SequenceMode::ProteinToProtein | SequenceMode::NucleotideToNucleotide => {
                let q_bytes = query.sequence.as_bytes();
                if let Some(idx) = kmer_index.as_ref() {
                    // Index-based fast path: seed lookup (with optional neighborhood
                    // expansion) → SIMD SW on candidates → CPU traceback for top hits.
                    run_one_to_many_indexed(
                        idx,
                        neighborhood.as_deref(),
                        &executor,
                        &aligner,
                        query,
                        q_bytes,
                        db_records,
                        db_residues,
                        ev_params,
                        &scoring,
                        &args,
                        &mode,
                    )
                } else if use_simd_pipeline {
                    // SIMD without index — falls back to score-batch over every subject.
                    run_one_to_many_gpu(&executor, &aligner, query, q_bytes,
                                         db_records, db_residues, ev_params,
                                         &scoring, &args, &mode)
                } else {
                    run_one_to_many(&aligner, query, q_bytes, db_records, db_residues,
                                    ev_params, &args, &mode)
                }
            }
            SequenceMode::NucleotideToProtein => {
                // 6-frame query against protein DB
                let frames = translate_six_frames(query.sequence.as_bytes(), genetic_code);
                let mut all: Vec<AlignmentResult> = frames.par_iter().flat_map(|tf| {
                    let qid = format!("{}|frame{:+}", query.id, signed_frame(tf.frame));
                    let qbytes = tf.protein.as_bytes();
                    let pseudo = FastaRecord::new(qid.clone(), qid, tf.protein.clone());
                    run_one_to_many(&aligner, &pseudo, qbytes, db_records, db_residues,
                                    ev_params, &args, &mode)
                }).collect();
                all.sort_by(|a, b| b.score.cmp(&a.score));
                all.truncate(args.max_aligns);
                all
            }
            SequenceMode::ProteinToNucleotide => {
                // Protein query against translated DB (built up-front above).
                let dbtr = translated_db.as_ref().unwrap();
                let q_bytes = query.sequence.as_bytes();
                let mut all: Vec<AlignmentResult> = dbtr.par_iter().filter_map(|ts| {
                    let mut a = aligner.align(q_bytes, ts.seq.as_bytes());
                    if a.score <= 0 { return None; }
                    a.query_id = query.id.clone();
                    a.subject_id = ts.id.clone();
                    finalize(&mut a, query.length, ts.seq.len(), db_residues, ev_params, &args);
                    if a.evalue > args.evalue { return None; }
                    if let Some(b) = args.min_bitscore { if a.bitscore < b { return None; } }
                    if a.identity < args.percent_identity { return None; }
                    Some(a)
                }).collect();
                all.sort_by(|a, b| b.score.cmp(&a.score));
                all.truncate(args.max_aligns);
                all
            }
        };

        for a in &aligns {
            writer.write_alignment(a)?;
        }
        total_alignments_written += aligns.len();
    }

    writer.write_footer()?;
    writer.flush()?;
    search.stop_and_print();

    if !args.quiet {
        eprintln!("[saber] wrote {} alignments", total_alignments_written);
        if let Some(p) = &args.output {
            eprintln!("[saber] results: {}", p.display());
        }
    }
    Ok(())
}

/// Build the aligner appropriate for the mode and parameters.
fn build_aligner(
    mode: SequenceMode,
    matrix_system: ScoringSystem,
    gap_open: i32,
    gap_extend: i32,
    nt_match: i32,
    nt_mismatch: i32,
) -> SmithWaterman {
    match mode {
        SequenceMode::NucleotideToNucleotide => {
            SmithWaterman::nucleotide(gap_open, gap_extend, nt_match, nt_mismatch)
        }
        // pp, nx, pn all align in protein space → use a substitution matrix.
        _ => {
            let m = Arc::new(ScoreMatrix::new(matrix_system));
            SmithWaterman::protein(gap_open, gap_extend, m)
        }
    }
}

/// Pick lambda/K for the mode. For nucleotide alignment with +2/-3 the classic
/// blastn defaults are λ≈0.625, K≈0.41; for protein we ask the matrix.
fn evalue_params_for(mode: SequenceMode, matrix: ScoringSystem) -> (f64, f64) {
    match mode {
        SequenceMode::NucleotideToNucleotide => (0.625, 0.41),
        _ => {
            let (l, k, _h) = ScoreMatrix::new(matrix).karlin_altschul_params();
            (l, k)
        }
    }
}

/// Set evalue / bitscore / coverages on an `AlignmentResult` once it has been computed.
#[allow(clippy::too_many_arguments)]
fn finalize(
    a: &mut AlignmentResult,
    query_len: usize,
    subject_len: usize,
    db_residues: usize,
    (lambda, kappa): (f64, f64),
    _args: &Args,
) {
    let ev = EValue::new(lambda, kappa, db_residues, query_len.max(1));
    a.bitscore = ev.score_to_bitscore(a.score);
    a.evalue = ev.score_to_evalue(a.score);
    a.query_coverage = if query_len > 0 { 100.0 * a.align_length as f64 / query_len as f64 } else { 0.0 };
    a.subject_coverage = if subject_len > 0 { 100.0 * a.align_length as f64 / subject_len as f64 } else { 0.0 };
}

/// Align one query against every record in `database`, filter, sort, truncate.
///
/// When `args.prefilter` is set, first applies a k-mer prefilter to the subject
/// list and only runs full SW on the survivors. This is the path that closes
/// the throughput gap to DIAMOND / MMseqs2 — we skip the obvious non-hits and
/// only pay for full alignment on plausible candidates.
fn run_one_to_many(
    aligner: &Arc<SmithWaterman>,
    query: &FastaRecord,
    query_bytes: &[u8],
    subjects: &[FastaRecord],
    db_residues: usize,
    ev_params: (f64, f64),
    args: &Args,
    mode: &SequenceMode,
) -> Vec<AlignmentResult> {
    // Apply prefilter if requested. Otherwise, every subject is a candidate.
    let candidates: Vec<usize> = if args.prefilter {
        let cfg = prefilter_config_for(mode, args);
        let subj_bytes: Vec<Vec<u8>> = subjects
            .iter()
            .map(|s| s.sequence.as_bytes().to_vec())
            .collect();
        let result = run_prefilter(query_bytes, &subj_bytes, &cfg);
        if !args.quiet {
            log::debug!(
                "prefilter: {} → {} candidates ({}% retained)",
                subjects.len(),
                result.kept.len(),
                if subjects.is_empty() {
                    0.0
                } else {
                    100.0 * result.kept.len() as f64 / subjects.len() as f64
                }
            );
        }
        result.kept
    } else {
        (0..subjects.len()).collect()
    };

    let mut out: Vec<AlignmentResult> = candidates
        .par_iter()
        .filter_map(|&idx| {
            let s = &subjects[idx];
            let mut a = aligner.align(query_bytes, s.sequence.as_bytes());
            if a.score <= 0 {
                return None;
            }
            a.query_id = query.id.clone();
            a.subject_id = s.id.clone();
            finalize(&mut a, query.length, s.length, db_residues, ev_params, args);
            if a.evalue > args.evalue {
                return None;
            }
            if let Some(b) = args.min_bitscore {
                if a.bitscore < b {
                    return None;
                }
            }
            if a.identity < args.percent_identity {
                return None;
            }
            Some(a)
        })
        .collect();
    out.sort_by(|a, b| b.score.cmp(&a.score));
    out.truncate(args.max_aligns);
    out
}

/// GPU-accelerated alignment pipeline (CUDASW++-style).
///
/// 1. Optional k-mer prefilter on host (rayon).
/// 2. Score-only batch through the GPU executor — one CUDA thread per
///    (query, subject) pair.
/// 3. Take the top-N candidates by raw score (overshoot N=4×max_aligns to
///    give E-value filtering room).
/// 4. Run CPU traceback only on those candidates to get full
///    [`AlignmentResult`] objects with CIGAR strings, %-identity, etc.
///
/// On a 100k-subject DB, this typically converts a multi-second CPU search
/// into a sub-second one: only 40-100 CPU tracebacks instead of 100k.
#[allow(clippy::too_many_arguments)]
fn run_one_to_many_gpu(
    executor: &Arc<dyn GpuExecutor>,
    aligner: &Arc<SmithWaterman>,
    query: &FastaRecord,
    query_bytes: &[u8],
    subjects: &[FastaRecord],
    db_residues: usize,
    ev_params: (f64, f64),
    scoring: &saber_rs::algorithm::Scoring,
    args: &Args,
    mode: &SequenceMode,
) -> Vec<AlignmentResult> {
    // Step 1: prefilter (optional).
    let candidates: Vec<usize> = if args.prefilter {
        let cfg = prefilter_config_for(mode, args);
        let subj_bytes: Vec<Vec<u8>> = subjects
            .iter()
            .map(|s| s.sequence.as_bytes().to_vec())
            .collect();
        run_prefilter(query_bytes, &subj_bytes, &cfg).kept
    } else {
        (0..subjects.len()).collect()
    };

    if candidates.is_empty() {
        return Vec::new();
    }

    // Step 2: GPU score-only.
    let candidate_seqs: Vec<Vec<u8>> = candidates
        .iter()
        .map(|&i| subjects[i].sequence.as_bytes().to_vec())
        .collect();
    let scores = match executor.score_batch(
        query_bytes,
        &candidate_seqs,
        scoring,
        args.gap_open,
        args.gap_extend,
    ) {
        Ok(s) => s,
        Err(e) => {
            log::error!("GPU score_batch failed: {e}. Falling back to CPU for this query.");
            return run_one_to_many(
                aligner, query, query_bytes, subjects, db_residues, ev_params, args, mode,
            );
        }
    };

    // Step 3: rank by score, take overshoot of max_aligns to allow E-value drops.
    let overshoot = (args.max_aligns + args.max_aligns / 4 + 50).max(50);
    let mut ranked: Vec<(usize, i32)> = scores
        .iter()
        .map(|gs| (candidates[gs.subject_idx], gs.score))
        .filter(|(_, sc)| *sc > 0)
        .collect();
    ranked.sort_unstable_by(|a, b| b.1.cmp(&a.1));
    ranked.truncate(overshoot);

    // Step 4: CPU traceback on the survivors → full AlignmentResult.
    let mut out: Vec<AlignmentResult> = ranked
        .par_iter()
        .filter_map(|&(idx, _gpu_score)| {
            let s = &subjects[idx];
            let mut a = aligner.align(query_bytes, s.sequence.as_bytes());
            if a.score <= 0 {
                return None;
            }
            a.query_id = query.id.clone();
            a.subject_id = s.id.clone();
            finalize(&mut a, query.length, s.length, db_residues, ev_params, args);
            if a.evalue > args.evalue {
                return None;
            }
            if let Some(b) = args.min_bitscore {
                if a.bitscore < b {
                    return None;
                }
            }
            if a.identity < args.percent_identity {
                return None;
            }
            Some(a)
        })
        .collect();
    out.sort_by(|a, b| b.score.cmp(&a.score));
    out.truncate(args.max_aligns);
    out
}

/// Index-accelerated alignment pipeline.
///
/// 1. Look up query k-mers in the inverted index → seed hits.
///    If `neighborhood` is supplied, each query k-mer expands to all k-mers
///    scoring ≥ T against it (BLAST T-score expansion) — dramatically more
///    seed hits per real homolog, while random subjects don't benefit.
/// 2. Cluster hits by (subject, diagonal-bucket); rank subjects by best
///    bucket's hit count. Top-N become candidates.
/// 3. Run SIMD SW on candidate subjects via the executor.
/// 4. Rank by SW score, take overshoot of `max_aligns × 4`.
/// 5. CPU traceback only for those → full `AlignmentResult`.
/// 6. Apply E-value / identity / bitscore filters, sort, truncate.
#[allow(clippy::too_many_arguments)]
fn run_one_to_many_indexed(
    index: &KmerIndex,
    neighborhood: Option<&NeighborhoodGenerator>,
    executor: &Arc<dyn GpuExecutor>,
    aligner: &Arc<SmithWaterman>,
    query: &FastaRecord,
    query_bytes: &[u8],
    subjects: &[FastaRecord],
    db_residues: usize,
    ev_params: (f64, f64),
    scoring: &saber_rs::algorithm::Scoring,
    args: &Args,
    mode: &SequenceMode,
) -> Vec<AlignmentResult> {
    // ---- Seed lookup ----
    let hits = match mode {
        SequenceMode::NucleotideToNucleotide => index.query_nucleotide(query_bytes),
        _ => {
            if let Some(nb) = neighborhood {
                index.query_protein_with_neighborhood(query_bytes, nb)
            } else {
                index.query_protein(query_bytes)
            }
        }
    };
    if hits.is_empty() {
        return Vec::new();
    }

    // ---- Candidate clustering ----
    // With neighborhood expansion (T>0), use BLAST's two-hit heuristic:
    // pairs of close hits on the same exact diagonal. Without neighborhood,
    // use the looser bucket-based clustering.
    let candidates = if neighborhood.is_some() {
        index.candidates_two_hit(&hits, args.two_hit_window, args.index_max_candidates)
    } else {
        index.candidates(
            &hits,
            args.index_min_hits,
            args.index_diag_bucket,
            args.index_max_candidates,
        )
    };
    log::debug!(
        "query={}: {} seed hits → {} candidates (of {} subjects)",
        query.id,
        hits.len(),
        candidates.len(),
        subjects.len()
    );
    if candidates.is_empty() {
        return Vec::new();
    }

    // ---- Ungapped X-drop filter (BLAST's HSP stage) ----
    // For each candidate, find the maximum ungapped HSP score across
    // *all* its seed hits — not just one anchor. The issue with single-
    // anchor extension is that when a seed lands outside the alignment's
    // high-scoring core, X-drop terminates before reaching it, and the
    // candidate gets rejected even though it's a real homolog.
    //
    // We also track the HSP's query range for each candidate. That range
    // tells the traceback step which slice of the query+subject is worth
    // running scalar SW on — typically <100 residues even when the query
    // is 1500 residues long.
    let mut hsp_range: std::collections::HashMap<u32, (i32, usize, usize, i32)> =
        std::collections::HashMap::new();
    // Maps subject_idx → (best_hsp_score, best_q_start, best_q_end, diagonal)

    let candidates = if matches!(mode, SequenceMode::ProteinToProtein) && args.hsp_threshold > 0 {
        let q_enc = saber_rs::gpu::encoding::encode_protein(query_bytes);
        let matrix = match scoring {
            saber_rs::algorithm::Scoring::Matrix(m) => ScoreMatrix::new(m.system),
            _ => {
                log::error!("indexed protein path got non-matrix scoring; this is a bug");
                return Vec::new();
            }
        };

        // Build the candidate-subject set for fast lookup during hit iteration.
        let candidate_subjects: std::collections::HashSet<u32> =
            candidates.iter().map(|c| c.subject_idx).collect();

        for hit in &hits {
            if !candidate_subjects.contains(&hit.subject_idx) {
                continue;
            }
            let s_enc = index.encoded_subject(hit.subject_idx as usize);
            let q_anchor = hit.q_pos as usize;
            let s_anchor = hit.s_pos as usize;
            let (hsp, q_start, q_end) = saber_rs::algorithm::ungapped::extend_ungapped_with_range(
                &q_enc, s_enc, q_anchor, s_anchor, &matrix, args.x_drop,
            );
            let e = hsp_range.entry(hit.subject_idx).or_insert((0, 0, 0, 0));
            if hsp > e.0 {
                *e = (hsp, q_start, q_end, hit.diagonal);
            }
        }

        let n_before = candidates.len();
        let filtered: Vec<_> = candidates
            .into_iter()
            .filter(|c| hsp_range.get(&c.subject_idx).map(|r| r.0).unwrap_or(0) >= args.hsp_threshold)
            .collect();
        log::debug!(
            "query={}: ungapped filter {} → {} candidates (threshold S={})",
            query.id, n_before, filtered.len(), args.hsp_threshold,
        );
        filtered
    } else {
        candidates
    };

    if candidates.is_empty() {
        return Vec::new();
    }

    // ---- SIMD SW on candidates ----
    // Two important optimizations vs the executor.score_batch() path:
    //
    // 1. **Skip re-encoding**: the KmerIndex already encoded every subject
    //    when it was built. We pass those encoded slices directly into the
    //    SIMD kernel, avoiding the Vec<u8>→Vec<u8>→encoded triple-allocation
    //    per candidate. On the 100-query × 5003-subject bench this is ~25%
    //    of total wall time.
    //
    // 2. **Banded slicing (protein only)**: each candidate has a
    //    `best_diagonal` from the seed clustering — the diagonal where its
    //    seeds concentrated. The real alignment is almost certainly within
    //    ±band residues of that diagonal, so we restrict the subject we
    //    feed into SW to the slice `[s_start, s_end]` corresponding to
    //    `query × subject_band`. This shrinks the DP matrix from
    //    `q × subj_len` to `q × 2band`. Pure win when `subj_len > 2band`.
    //    Band default 64 chosen to be comfortably wider than typical
    //    biological indels (~20 aa) while still cutting work substantially
    //    on long subjects.
    //
    // Falls back to the executor path on nucleotide mode (which doesn't
    // use the protein matrix-encoded path) or if anything goes wrong.

    let use_encoded_path = matches!(mode, SequenceMode::ProteinToProtein);

    let scores: Vec<i32> = if use_encoded_path {
        // Encode the query once, reuse it across all candidates.
        let q_enc = saber_rs::gpu::encoding::encode_protein(query_bytes);

        // Build sliced subject view per candidate. Owned Vec only when the
        // slice differs from the full subject — for short subjects we just
        // borrow the index's encoded buffer.
        let band: i32 = args.band as i32;
        let q_len = query_bytes.len() as i32;
        let mut sliced_owned: Vec<Vec<u8>> = Vec::with_capacity(candidates.len());
        let mut sliced_refs: Vec<&[u8]> = Vec::with_capacity(candidates.len());

        for c in &candidates {
            let full = index.encoded_subject(c.subject_idx as usize);
            let slen = full.len() as i32;
            let d = c.best_diagonal;
            let s_lo = (d - band).max(0).min(slen) as usize;
            let s_hi = (d + q_len + band).max(0).min(slen) as usize;
            // Worth slicing only when the slice is meaningfully smaller than
            // the full subject. SIMD SW cost scales with the longest sequence
            // in the lane batch, so trimming 10% of a 500-aa subject doesn't
            // help if another lane has a 1500-aa subject. The owned Vec
            // allocation, on the other hand, has a fixed cost we eat
            // regardless. Threshold: slice must be ≤ 75% of full length
            // AND save at least 64 residues. Otherwise borrow the full
            // encoded subject — zero allocation.
            let want_slice = band > 0
                && s_hi > s_lo
                && (s_hi - s_lo) * 4 < full.len() * 3      // ≥ 25% reduction
                && full.len() - (s_hi - s_lo) >= 64;       // and ≥ 64 saved
            if !want_slice {
                sliced_refs.push(full);
            } else {
                sliced_owned.push(full[s_lo..s_hi].to_vec());
                sliced_refs.push(&[]);
            }
        }
        // Second pass: fix up the borrows now that sliced_owned is stable.
        let mut owned_iter = sliced_owned.iter();
        for r in sliced_refs.iter_mut() {
            if r.is_empty() {
                if let Some(o) = owned_iter.next() {
                    *r = o.as_slice();
                }
            }
        }

        // Matrix needed for the encoded SW path. In ProteinToProtein mode
        // the scoring is always Scoring::Matrix, so this match is total.
        // We rebuild via ScoreMatrix::new() since ScoreMatrix isn't Clone.
        let matrix: ScoreMatrix = match scoring {
            saber_rs::algorithm::Scoring::Matrix(m) => ScoreMatrix::new(m.system),
            _ => {
                log::error!("indexed protein path got non-matrix scoring; this is a bug");
                return Vec::new();
            }
        };

        // GPU FAST PATH:
        // When a CUDA executor is loaded (--features gpu-cuda, GPU present),
        // we'd lose more than we'd gain doing the slicing dance here — GPU
        // wants large contiguous batches and pays a fixed launch overhead.
        // So for GPU we feed it the full (banded but un-sliced) subjects via
        // score_batch, accepting the small extra DP work in exchange for
        // amortized launch cost. The encoded-path CPU SIMD is the right
        // choice for CPU; GPU executor is the right choice when present.
        if executor.backend_name() == "cuda" {
            let candidate_seqs: Vec<Vec<u8>> = candidates
                .iter()
                .map(|c| {
                    let full = index.encoded_subject(c.subject_idx as usize);
                    let slen = full.len() as i32;
                    let band: i32 = args.band as i32;
                    let q_len = query_bytes.len() as i32;
                    if band <= 0 {
                        return full.to_vec();
                    }
                    let s_lo = (c.best_diagonal - band).max(0).min(slen) as usize;
                    let s_hi = (c.best_diagonal + q_len + band).max(0).min(slen) as usize;
                    if s_hi > s_lo {
                        full[s_lo..s_hi].to_vec()
                    } else {
                        full.to_vec()
                    }
                })
                .collect();
            match executor.score_batch(
                &q_enc,
                &candidate_seqs,
                scoring,
                args.gap_open,
                args.gap_extend,
            ) {
                Ok(gs) => return_gpu_scores(gs),
                Err(e) => {
                    log::warn!(
                        "GPU score_batch failed: {e}. Falling back to CPU SIMD for this query."
                    );
                    saber_rs::algorithm::sw_simd::sw_simd_protein_batch_encoded(
                        &q_enc,
                        &sliced_refs,
                        &matrix,
                        args.gap_open,
                        args.gap_extend,
                    )
                }
            }
        } else {
            saber_rs::algorithm::sw_simd::sw_simd_protein_batch_encoded(
                &q_enc,
                &sliced_refs,
                &matrix,
                args.gap_open,
                args.gap_extend,
            )
        }
    } else {
        // Nucleotide path — go through the executor as before. The database's
        // ASCII strings are gone post-drop, so reconstruct from the index's
        // encoded buffer (which is what we want anyway — same source of truth
        // the protein path uses).
        let candidate_seqs: Vec<Vec<u8>> = candidates
            .iter()
            .map(|c| index.ascii_subject(c.subject_idx as usize))
            .collect();
        match executor.score_batch(
            query_bytes,
            &candidate_seqs,
            scoring,
            args.gap_open,
            args.gap_extend,
        ) {
            Ok(gs) => gs.into_iter().map(|x| x.score).collect(),
            Err(e) => {
                // Database ASCII has been dropped by main() before search begins
                // (it's redundant with the index's encoded buffer). Calling the
                // legacy run_one_to_many here would dereference empty strings,
                // so we degrade gracefully: log the failure, return no hits for
                // this query, and let the rest of the run proceed.
                log::error!(
                    "score_batch failed on indexed nucleotide path: {e}. \
                     Skipping this query (re-run with --no-index to get scalar fallback)."
                );
                return Vec::new();
            }
        }
    };

    // ---- Rank, traceback top-K, filter ----
    // Overshoot accounts for hits that pass SIMD scoring but fail the
    // E-value/identity/bitscore filters. SIMD scores are exact SW scores
    // so the ordering is essentially correct; we just need a buffer for
    // post-filter dropouts. Empirically 25% extra is plenty; capped to
    // ≥ 50 to handle very small max_aligns (e.g., user asks for top 1).
    let overshoot = (args.max_aligns + args.max_aligns / 4 + 50).max(50);
    let mut ranked: Vec<(u32, i32)> = scores
        .iter()
        .enumerate()
        .map(|(i, &sc)| (candidates[i].subject_idx, sc))
        .filter(|(_, sc)| *sc > 0)
        .collect();
    ranked.sort_unstable_by(|a, b| b.1.cmp(&a.1));
    ranked.truncate(overshoot);

    // Build a candidate-idx → diagonal lookup for banded traceback below.
    // We use the same band slack as the SIMD scoring step.
    let cand_diag: std::collections::HashMap<u32, i32> = candidates
        .iter()
        .map(|c| (c.subject_idx, c.best_diagonal))
        .collect();

    let mut out: Vec<AlignmentResult> = ranked
        .par_iter()
        .filter_map(|&(sidx, _)| {
            // Subject ASCII bytes for traceback. The DB ASCII has been dropped
            // by main() (memory diet); reconstruct from the index's encoded
            // buffer. Under --keep-ascii (and an in-memory DB) we can borrow
            // the original sequence directly.
            let has_meta = index.has_subject_meta();
            let s_opt = if has_meta { None } else { subjects.get(sidx as usize) };
            let s_ascii: Vec<u8> = match s_opt {
                Some(s) if args.keep_ascii && !s.sequence.is_empty() => s.sequence.as_bytes().to_vec(),
                _ => index.ascii_subject(sidx as usize),
            };
            let s_bytes: &[u8] = &s_ascii;

            // WINDOWED + BANDED TRACEBACK:
            //
            // 1. If we have an HSP range from the ungapped step, use it to
            //    slice BOTH the query and the subject to a window around
            //    the alignment. The full DP matrix shrinks from
            //    (query_len × subject_len) to (window × window).
            //
            // 2. If we don't have an HSP range (e.g. ungapped filter was
            //    disabled), fall back to slicing only the subject, using
            //    the candidate's best_diagonal as the band center.
            //
            // 3. The slack values are deliberately generous — scalar SW
            //    needs room to find the optimal endpoints, and a 100-cell
            //    matrix is still ~1000× cheaper than a 1500×1500 matrix.
            //    A real alignment of length L has at most L gaps total,
            //    so window_slack = 2×L is the theoretical max needed.
            //    We use a fixed slack to keep things bounded.
            let (q_lo, q_hi, s_lo, s_hi) = if let Some(&(_, q_s, q_e, diag)) = hsp_range.get(&sidx) {
                // Window-based slicing.
                //
                // Slack must be generous enough that gapped SW can extend
                // beyond the X-drop ungapped HSP boundary, but tight enough
                // that the DP matrix doesn't blow up memory on huge queries.
                //
                // We use a fixed slack of max(band*4, 256) and HARD-CAP the
                // window size so the DP matrix never exceeds ~50 MB. For
                // very long HSPs (> 1500 aa) the cap will trim the window
                // around the HSP center, which is a sensible heuristic:
                // most of the alignment information is at the high-scoring
                // central residues, and the cap matches BLAST's behavior
                // of breaking very long HSPs into pieces.
                let q_len = query_bytes.len() as i32;
                let s_len = s_bytes.len() as i32;
                let slack = (args.band as i32 * 4).max(256);
                // Hard cap on window dimension. 2048 × 2048 × 12 bytes ≈ 50 MB.
                const MAX_WINDOW: i32 = 2048;
                let q_center = (q_s as i32 + q_e as i32) / 2;
                let s_center = q_center + diag;
                // Half-width around the HSP center.
                let q_lo_want = (q_s as i32) - slack;
                let q_hi_want = (q_e as i32) + slack;
                let q_lo = q_lo_want.max(q_center - MAX_WINDOW / 2).max(0).min(q_len) as usize;
                let q_hi = q_hi_want.min(q_center + MAX_WINDOW / 2).min(q_len).max(0) as usize;
                let s_lo_want = s_center - (q_e as i32 - q_s as i32) / 2 - slack;
                let s_hi_want = s_center + (q_e as i32 - q_s as i32) / 2 + slack;
                let s_lo = s_lo_want.max(s_center - MAX_WINDOW / 2).max(0).min(s_len) as usize;
                let s_hi = s_hi_want.min(s_center + MAX_WINDOW / 2).min(s_len).max(0) as usize;
                (q_lo, q_hi, s_lo, s_hi)
            } else if let Some(&d) = cand_diag.get(&sidx) {
                let q_len = query_bytes.len() as i32;
                let s_len = s_bytes.len() as i32;
                let tb_slack = (args.band as i32 * 2).max(96);
                let lo = (d - tb_slack).max(0).min(s_len) as usize;
                let hi = (d + q_len + tb_slack).max(0).min(s_len) as usize;
                (0, query_bytes.len(), lo, hi)
            } else {
                (0, query_bytes.len(), 0, s_bytes.len())
            };

            if q_hi <= q_lo || s_hi <= s_lo {
                return None;
            }

            let q_slice = &query_bytes[q_lo..q_hi];
            let s_slice = &s_bytes[s_lo..s_hi];
            let mut a = aligner.align(q_slice, s_slice);
            if a.score <= 0 {
                return None;
            }
            // Fix up coords to refer to full sequences, not the slices.
            a.query_start += q_lo;
            a.query_end += q_lo;
            a.subject_start += s_lo;
            a.subject_end += s_lo;

            a.query_id = query.id.clone();
            let (subj_id, subj_len) = if has_meta {
                (index.subject_id(sidx as usize).to_string(), index.subject_len(sidx as usize))
            } else {
                let s = &subjects[sidx as usize];
                (s.id.clone(), s.length)
            };
            a.subject_id = subj_id;
            finalize(&mut a, query.length, subj_len, db_residues, ev_params, args);
            if a.evalue > args.evalue {
                return None;
            }
            if let Some(b) = args.min_bitscore {
                if a.bitscore < b {
                    return None;
                }
            }
            if a.identity < args.percent_identity {
                return None;
            }
            Some(a)
        })
        .collect();
    out.sort_by(|a, b| b.score.cmp(&a.score));
    out.truncate(args.max_aligns);
    out
}

struct TranslatedSubject {
    id: String,
    seq: String,
}

fn signed_frame(frame: u8) -> i8 {
    if frame < 3 { (frame as i8) + 1 } else { -((frame as i8) - 2) }
}

/// Translate a Vec<GpuScore> from the executor into the parallel Vec<i32>
/// that the rest of the indexed pipeline expects. `GpuScore` carries the
/// score plus the subject index it refers to — but since the executor
/// preserves caller-input order (per its contract), we can just take
/// scores in order.
fn return_gpu_scores(gs: Vec<saber_rs::gpu::GpuScore>) -> Vec<i32> {
    gs.into_iter().map(|s| s.score).collect()
}

/// Build a [`PrefilterConfig`] from CLI args, with sensible defaults per mode.
fn prefilter_config_for(mode: &SequenceMode, args: &Args) -> PrefilterConfig {
    let base = match mode {
        SequenceMode::NucleotideToNucleotide => PrefilterConfig::nucleotide_default(),
        _ => PrefilterConfig::protein_default(),
    };
    PrefilterConfig {
        k: if args.prefilter_k > 0 { args.prefilter_k } else { base.k },
        min_hits: if args.prefilter_min_hits > 0 { args.prefilter_min_hits } else { base.min_hits },
        top_n: args.prefilter_top_n,
        alphabet_size: base.alphabet_size,
    }
}
