//! End-to-end search driver shared by the CLI and library users: obtains an
//! index (opened, built in memory, or spilled to a temporary file under a
//! memory budget), streams queries in bounded batches, runs them through the
//! [`Engine`] on a fixed-size thread pool and writes results in input order.
//!
//! Scheduling policy (no nested parallelism):
//! * a batch with at least `2 × threads` queries is processed query-parallel
//!   (one engine call per query, each with a thread-local [`Worker`]);
//! * smaller batches (few, typically expensive queries) are processed one
//!   query at a time with the *database* split into subject ranges searched
//!   in parallel, so a single query can use every core.
//! Output order is the input order in both cases; at most one batch of
//! results is buffered (bounded memory, natural back-pressure).

use std::cell::RefCell;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use rayon::prelude::*;

use crate::algorithm::translation::{frame_label, protein_to_nt, CodonTable, TranslationTable};
use crate::algorithm::{AlignmentResult, SequenceMode};
use crate::alphabet::{encode_nucleotide, EncodedProtein};
use crate::io::stream::{FastaStream, OwnedRecord};
use crate::io::{CSVWriter, JSONWriter, M0Writer, M8Writer, M9Writer, OutputFormat, SAMWriter, Writer};
use crate::scoring::matrices::{ScoreMatrix, ScoringSystem};
use crate::scoring::stats::{gapped_nucleotide, gapped_protein, GapCosts, KarlinParams, SearchStatistics};

use super::engine::{Counters, Engine, Hit, SearchParams, Worker};
use super::index::{BuildOptions, IndexBuilder, SeedConfig, SeedIndex, SeqKind};

/// Everything the driver needs (filled from CLI flags or by library users).
#[derive(Debug, Clone)]
pub struct DriverConfig {
    pub mode: SequenceMode,
    pub query: Option<PathBuf>,
    pub database: Option<PathBuf>,
    pub db_index: Option<PathBuf>,
    pub output: Option<PathBuf>,
    pub outfmt: OutputFormat,
    pub matrix: ScoringSystem,
    pub gaps: GapCosts,
    pub nt_match: i32,
    pub nt_mismatch: i32,
    pub allow_uncalibrated: bool,
    /// E-value model (per-subject search spaces by default).
    pub evalue_model: crate::scoring::EvalueModel,
    pub genetic_code: TranslationTable,
    pub threads: usize,
    /// Approximate memory budget in bytes (None = unlimited).
    pub memory_budget: Option<u64>,
    pub params: SearchParams,
    pub seed: SeedConfig,
    pub partition_residues: usize,
    pub min_length: Option<usize>,
    pub max_length: Option<usize>,
    pub stats_path: Option<PathBuf>,
    pub verify_index: bool,
    pub tmpdir: PathBuf,
    pub quiet: bool,
    /// Queries per batch (0 = automatic).
    pub query_batch: usize,
    pub preset: String,
    /// Use a GPU score backend when one is available (`--gpu`).
    pub gpu: bool,
    /// Scheduling of query batches against the index.
    pub batch_order: BatchOrder,
    /// Build seed postings even when the scan schedule will be used
    /// (`--index-postings`; needed to search the index with the query or
    /// partition schedules later).
    pub index_postings: bool,
    /// Query residues per scan batch (0 = automatic).
    pub scan_batch_residues: usize,
}

/// How a batch of queries is scheduled against the index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BatchOrder {
    /// The scan schedule for protein searches it supports; otherwise
    /// query-major for in-memory indexes and partition-major for
    /// multi-partition memory-mapped indexes larger than half the memory budget.
    Auto,
    /// Seed tables of a whole query batch, database streamed once per batch
    /// ([`Engine::search_scan`]); needs no seed postings.
    Scan,
    /// Each query scans every partition (one engine call per query).
    Query,
    /// Each partition is scanned once per batch ([`Engine::search_batch`]).
    Partition,
}

impl std::str::FromStr for BatchOrder {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s.to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "scan" => Ok(Self::Scan),
            "query" | "query-major" => Ok(Self::Query),
            "partition" | "partition-major" => Ok(Self::Partition),
            o => Err(format!("unknown batch order '{o}' (expected auto|scan|query|partition)")),
        }
    }
}

impl DriverConfig {
    /// CLI-equivalent defaults for `mode` with the `default` preset:
    /// BLOSUM62, gaps 11/1 (protein) or 7/2 with +2/-3 (nn), E <= 10,
    /// 10 subjects per query, tantan seed masking, composition correction.
    pub fn new(mode: SequenceMode) -> Self {
        let nn = mode == SequenceMode::NucleotideToNucleotide;
        let mut c = Self {
            mode,
            query: None,
            database: None,
            db_index: None,
            output: None,
            outfmt: OutputFormat::BlastM8,
            matrix: ScoringSystem::Blosum62,
            gaps: if nn { GapCosts { open: 7, extend: 2 } } else { GapCosts { open: 11, extend: 1 } },
            nt_match: 2,
            nt_mismatch: -3,
            allow_uncalibrated: false,
            evalue_model: crate::scoring::EvalueModel::PerSubject,
            genetic_code: TranslationTable::Standard,
            threads: rayon::current_num_threads(),
            memory_budget: default_memory_budget(),
            params: SearchParams::default(),
            seed: if nn { SeedConfig::nucleotide(11) } else { SeedConfig::protein_default() },
            partition_residues: 64 << 20,
            min_length: None,
            max_length: None,
            stats_path: None,
            verify_index: false,
            tmpdir: std::env::temp_dir(),
            quiet: true,
            query_batch: 0,
            preset: "default".into(),
            gpu: false,
            batch_order: BatchOrder::Auto,
            index_postings: false,
            scan_batch_residues: 0,
        };
        if !nn {
            c.apply_preset("default").expect("built-in preset");
        }
        c
    }

    /// Apply a named preset (seed alphabet/shapes, hit filter, HSP threshold,
    /// candidate cap, exhaustive flag) — see [`super::presets`].
    pub fn apply_preset(&mut self, name: &str) -> Result<(), String> {
        let p = super::presets::find(name)?;
        self.preset = p.name.to_string();
        self.params.exhaustive = p.exhaustive;
        self.params.max_candidates = p.max_candidates;
        self.params.hsp_threshold = p.hsp_threshold;
        self.params.hsp_floor_fraction = p.hsp_floor_fraction;
        self.params.score_band = p.score_band;
        self.params.score_prefilter_band = p.prefilter_band;
        if self.mode != SequenceMode::NucleotideToNucleotide {
            self.params.hit_filter_min_ids = p.hit_filter_ids;
            self.params.hit_filter_exact_weight = p.exact_weight;
            self.seed.alphabet = crate::alphabet::SeedAlphabet::parse(p.alphabet)?;
            self.seed.shapes = p.shapes.iter().map(|s| crate::alphabet::SeedShape::parse(s)).collect::<Result<Vec<_>, _>>()?;
        }
        Ok(())
    }
}

/// Summary of a run (also written as JSON with `--stats`).
#[derive(Debug, Clone, serde::Serialize)]
pub struct RunSummary {
    pub saber_version: String,
    pub mode: String,
    pub preset: String,
    pub threads: usize,
    pub queries: u64,
    pub queries_skipped_length: u64,
    pub alignments: u64,
    pub db_subjects: u64,
    pub db_residues: u64,
    pub index_bytes: u64,
    pub index_source: String,
    pub index_build_s: f64,
    pub search_s: f64,
    pub total_s: f64,
    pub peak_rss_kb: u64,
    pub evalue_params: String,
    pub evalue_calibrated: bool,
    pub seed_alphabet: String,
    pub seed_shapes: Vec<String>,
    pub masking: String,
    pub comp_bias: bool,
    pub heuristic_limits: Vec<String>,
    pub score_backend: String,
    pub memory_budget_bytes: Option<u64>,
    /// Batches run query-parallel / partition-major / intra-query parallel.
    pub batches_query_major: u64,
    pub batches_partition_major: u64,
    pub batches_intra_query: u64,
    /// Batches run with the scan schedule.
    pub batches_scan: u64,
    /// Database seed mask used by the scan schedule (index | recomputed | none).
    pub scan_mask: String,
    /// Seconds spent preparing the scan schedule's database data (once per run).
    pub scan_db_s: f64,
    pub counters: Counters,
}

/// Peak resident set size of this process (kB), from /proc (Linux only).
pub fn peak_rss_kb() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| s.lines().find(|l| l.starts_with("VmHWM:")).map(|l| l.to_string()))
        .and_then(|l| l.split_whitespace().nth(1).and_then(|v| v.parse().ok()))
        .unwrap_or(0)
}

/// Default memory budget: 3/4 of physical RAM (Linux), else unlimited.
pub fn default_memory_budget() -> Option<u64> {
    std::fs::read_to_string("/proc/meminfo")
        .ok()
        .and_then(|s| s.lines().find(|l| l.starts_with("MemTotal:")).map(|l| l.to_string()))
        .and_then(|l| l.split_whitespace().nth(1).and_then(|v| v.parse::<u64>().ok()))
        .map(|kb| kb * 1024 / 4 * 3)
}

/// Parse sizes like `512M`, `2G`, `1.5g`, `4096` (bytes).
pub fn parse_size(s: &str) -> Result<u64, String> {
    let t = s.trim().to_ascii_uppercase();
    let (num, mul) = match t.chars().last() {
        Some('K') => (&t[..t.len() - 1], 1u64 << 10),
        Some('M') => (&t[..t.len() - 1], 1u64 << 20),
        Some('G') => (&t[..t.len() - 1], 1u64 << 30),
        Some('T') => (&t[..t.len() - 1], 1u64 << 40),
        _ => (&t[..], 1u64),
    };
    let v: f64 = num.trim().parse().map_err(|_| format!("invalid size '{s}'"))?;
    if v <= 0.0 {
        return Err(format!("size must be positive: '{s}'"));
    }
    Ok((v * mul as f64) as u64)
}

/// True when the configuration is searched with the scan schedule (protein
/// databases; single-hit seeding without neighborhood words or database
/// minimizers; not the exhaustive reference mode; no GPU backend).
pub fn scan_eligible(cfg: &DriverConfig) -> bool {
    cfg.mode != SequenceMode::NucleotideToNucleotide
        && matches!(cfg.batch_order, BatchOrder::Auto | BatchOrder::Scan)
        && !cfg.params.exhaustive
        && cfg.params.two_hit_window == 0
        && cfg.params.seed_min_score == 0
        && cfg.params.neighborhood_t == 0
        && cfg.seed.minimizer_window <= 1
        && !cfg.seed.shapes.is_empty()
        && !cfg.gpu
}

fn kind_for(mode: SequenceMode) -> SeqKind {
    match mode {
        SequenceMode::NucleotideToNucleotide => SeqKind::Nucleotide,
        _ => SeqKind::Protein,
    }
}

/// Build an index from a FASTA file (streaming). `to_file` = external-memory
/// build into that path; otherwise in memory.
pub fn build_index(cfg: &DriverConfig, fasta: &Path, to_file: Option<&Path>) -> Result<(SeedIndex, f64), String> {
    let t0 = Instant::now();
    let mut seed = cfg.seed.clone();
    seed.kind = kind_for(cfg.mode);
    if seed.kind == SeqKind::Nucleotide && seed.alphabet.size != 4 {
        let k = seed.shapes.first().map(|s| s.span).unwrap_or(11).max(8);
        seed = SeedConfig::nucleotide(k);
    }
    let stats_bytes: u64 = (0..seed.shapes.len()).map(|i| seed.num_codes(i)).sum();
    let mut opt = BuildOptions::new(seed);
    opt.postings = cfg.index_postings || !scan_eligible(cfg);
    opt.scan_stats = scan_eligible(cfg) && kind_for(cfg.mode) == SeqKind::Protein && stats_bytes <= 256 << 20;
    opt.max_partition_residues = cfg.partition_residues;
    opt.masking_matrix = ScoringSystem::Blosum62;
    opt.source = fasta.display().to_string();
    let mut b = match to_file {
        Some(p) => IndexBuilder::to_file(opt, p).map_err(|e| e.to_string())?,
        None => IndexBuilder::in_memory(opt).map_err(|e| e.to_string())?,
    };
    if cfg.mode == SequenceMode::ProteinToNucleotide {
        b.set_translated(cfg.genetic_code, to_file).map_err(|e| e.to_string())?;
    }
    let mut fs = FastaStream::open(fasta).map_err(|e| format!("database: {e}"))?;
    let mut buf = Vec::new();
    let mut noncanon = 0usize;
    while let Some(r) = fs.next_record().map_err(|e| format!("database {}: {e}", fasta.display()))? {
        let id = r.id.to_string();
        match cfg.mode {
            SequenceMode::ProteinToNucleotide => b.add_nucleotide_translated(&id, r.seq).map_err(|e| e.to_string())?,
            SequenceMode::NucleotideToNucleotide => b.add_ascii(&id, r.seq, &mut buf).map_err(|e| e.to_string())?,
            _ => {
                noncanon += crate::alphabet::count_noncanonical_protein(r.seq);
                b.add_ascii(&id, r.seq, &mut buf).map_err(|e| e.to_string())?
            }
        }
    }
    if noncanon > 0 && !cfg.quiet {
        eprintln!("[saber] note: {noncanon} database residues were not in ARNDCQEGHILKMFPSTWYVBZX* and are scored as X");
    }
    let idx = match to_file {
        Some(p) => {
            b.finish_to_file().map_err(|e| e.to_string())?;
            SeedIndex::open(p).map_err(|e| e.to_string())?
        }
        None => b.finish_in_memory().map_err(|e| e.to_string())?,
    };
    Ok((idx, t0.elapsed().as_secs_f64()))
}

/// Karlin-Altschul parameters for the configured scoring system.
pub fn karlin_for(cfg: &DriverConfig) -> Result<KarlinParams, String> {
    match cfg.mode {
        SequenceMode::NucleotideToNucleotide => {
            gapped_nucleotide(cfg.nt_match, cfg.nt_mismatch, cfg.gaps, cfg.allow_uncalibrated).map_err(|e| e.to_string())
        }
        _ => gapped_protein(cfg.matrix, cfg.gaps, cfg.allow_uncalibrated).map_err(|e| e.to_string()),
    }
}

/// Database length/sequence count used for statistics.
fn db_stats_size(idx: &SeedIndex) -> (u64, u64) {
    if idx.is_translated() {
        let m = idx.meta();
        ((m.nt_residues / 3).max(1), m.nt_sequences.max(1))
    } else {
        (idx.total_residues(), idx.num_subjects() as u64)
    }
}

thread_local! {
    static WORKER: RefCell<Worker> = RefCell::new(Worker::new());
}

/// Make an index-compatible engine.
pub fn make_engine<'a>(cfg: &DriverConfig, idx: &'a SeedIndex) -> Result<Engine<'a>, String> {
    let kp = karlin_for(cfg)?;
    let (dl, dn) = db_stats_size(idx);
    let st = SearchStatistics::with_model(kp, dl, dn, cfg.evalue_model);
    let matrix = match idx.kind() {
        SeqKind::Nucleotide => ScoreMatrix::nucleotide(cfg.nt_match, cfg.nt_mismatch),
        SeqKind::Protein => ScoreMatrix::new(cfg.matrix),
    };
    let mut p = cfg.params.clone();
    p.subject_group = if idx.is_translated() { 6 } else { 1 };
    if idx.kind() == SeqKind::Nucleotide {
        p.comp_bias = false;
        p.query_masking = crate::masking::MaskingMode::None;
    }
    let mut eng = Engine::new(idx, matrix, cfg.gaps, st, p);
    if idx.kind() == SeqKind::Protein && cfg.mode != SequenceMode::NucleotideToNucleotide && cfg.seed.minimizer_window <= 1 {
        eng = eng.with_scan_seeds(super::engine::ScanSeeds { alphabet: cfg.seed.alphabet.clone(), shapes: cfg.seed.shapes.clone() });
    }
    if cfg.gpu {
        match crate::gpu::try_init_backend() {
            Some(b) => {
                if !cfg.quiet {
                    eprintln!("[saber] region scoring on GPU backend '{}' (experimental: not validated on hardware in this release)", b.name());
                }
                eng = eng.with_backend(b);
            }
            None => {
                if !cfg.quiet {
                    eprintln!(
                        "[saber] --gpu requested but no GPU backend is available ({}); using the CPU kernels",
                        if cfg!(feature = "gpu-cuda") { "no CUDA device initialised" } else { "binary built without --features gpu-cuda" }
                    );
                }
            }
        }
    }
    if cfg.params.neighborhood_t > 0 && !eng.neighborhood_active() && !cfg.quiet {
        eprintln!("[saber] note: --neighborhood-t needs the full 20-letter seed alphabet with seed weight 2..8; neighborhood expansion is off");
    }
    Ok(eng)
}

/// Encoded query frames for one record: (frame index or None, codes).
pub(crate) fn query_frames(mode: SequenceMode, rec: &OwnedRecord, ct: &CodonTable) -> Vec<(Option<usize>, Vec<u8>)> {
    match mode {
        SequenceMode::NucleotideToProtein => {
            let nt = encode_nucleotide(&rec.seq);
            (0..6)
                .map(|f| {
                    let mut out = Vec::new();
                    ct.translate_codes(&nt, f, &mut out);
                    (Some(f), out)
                })
                .collect()
        }
        SequenceMode::NucleotideToNucleotide => vec![(None, encode_nucleotide(&rec.seq))],
        _ => vec![(None, EncodedProtein::from_ascii(&rec.seq).as_slice().to_vec())],
    }
}

/// Convert an engine hit to the public result type (coordinates mapped back
/// to nucleotides for translated searches).
fn to_result(eng: &Engine<'_>, idx: &SeedIndex, rec: &OwnedRecord, frame: Option<usize>, h: &Hit, calibrated: bool) -> AlignmentResult {
    let (mut qs, mut qe) = (h.q_start, h.q_end);
    let mut qframe = 0i8;
    if let Some(f) = frame {
        let (a, b) = protein_to_nt(f, h.q_start, h.q_end, rec.seq.len());
        qs = a;
        qe = b;
        qframe = frame_label(f);
    }
    let (mut ss, mut se) = (h.s_start, h.s_end);
    let mut sframe = 0i8;
    let mut slen = h.s_len;
    if idx.is_translated() {
        let seq = (h.subject / 6) as usize;
        let f = (h.subject % 6) as usize;
        let nl = idx.nt_len(seq);
        let (a, b) = protein_to_nt(f, h.s_start, h.s_end, nl);
        ss = a;
        se = b;
        sframe = frame_label(f);
        slen = nl;
    }
    let qlen_report = rec.seq.len();
    let _ = eng;
    // Translated sequences keep SABER 1.x's "|frame±N" ID suffix so frame
    // and strand stay visible in every output format (coordinates are now
    // nucleotide positions, BLAST-style: start > end on the minus strand).
    let query_id = if frame.is_some() { format!("{}|frame{:+}", rec.id, qframe) } else { rec.id.clone() };
    let subject_id = if idx.is_translated() {
        format!("{}|frame{:+}", idx.subject_id(h.subject as usize), sframe)
    } else {
        idx.subject_id(h.subject as usize).to_string()
    };
    AlignmentResult {
        query_id,
        subject_id,
        score: h.score,
        evalue: h.evalue,
        bitscore: h.bitscore,
        identity: h.identity(),
        query_coverage: h.query_cover(),
        subject_coverage: h.subject_cover(),
        query_start: qs,
        query_end: qe,
        subject_start: ss,
        subject_end: se,
        align_length: h.columns,
        matches: h.matches,
        mismatches: h.mismatches,
        gaps: h.gap_columns,
        cigar: h.cigar(),
        query_seq: h.q_aln.clone(),
        subject_seq: h.s_aln.clone(),
        match_string: h.midline.clone(),
        gap_opens: h.gap_opens,
        positives: h.positives,
        query_length: qlen_report,
        subject_length: slen,
        query_frame: qframe,
        subject_frame: sframe,
        evalue_calibrated: calibrated,
        exact: h.exact_in_window,
    }
}

/// Merge per-frame hits of one record (translated queries: best alignments
/// per subject across frames, then -k / --max-hsps) and convert them.
pub(crate) fn finish_record(
    eng: &Engine<'_>,
    idx: &SeedIndex,
    cfg: &DriverConfig,
    rec: &OwnedRecord,
    mut all: Vec<(Option<usize>, Hit)>,
    calibrated: bool,
) -> Vec<AlignmentResult> {
    if cfg.mode == SequenceMode::NucleotideToProtein {
        all.sort_by(|a, b| b.1.score.cmp(&a.1.score).then(a.1.subject.cmp(&b.1.subject)).then(a.0.cmp(&b.0)));
        let mut per: std::collections::HashMap<u64, usize> = std::collections::HashMap::new();
        let mut kept = Vec::new();
        let mut subjects = 0usize;
        for (f, h) in all {
            let c = per.entry(h.subject).or_insert(0);
            if *c >= eng.params.max_hsps {
                continue;
            }
            if *c == 0 {
                if subjects >= eng.params.max_target_seqs {
                    continue;
                }
                subjects += 1;
            }
            *c += 1;
            kept.push((f, h));
        }
        all = kept;
    }
    all.iter().map(|(f, h)| to_result(eng, idx, rec, *f, h, calibrated)).collect()
}

/// Run one query (frame) on the calling thread with its thread-local worker.
pub(crate) fn search_local(eng: &Engine<'_>, codes: &[u8]) -> (Vec<Hit>, Counters) {
    WORKER.with(|w| {
        let mut w = w.borrow_mut();
        let before = w.counters;
        let h = eng.search(codes, &mut w);
        let mut d = w.counters;
        d.sub(&before);
        (h, d)
    })
}

fn open_writer(cfg: &DriverConfig) -> Result<Box<dyn Writer>, String> {
    let sink: Box<dyn Write + Send> = match &cfg.output {
        Some(p) => Box::new(std::io::BufWriter::new(std::fs::File::create(p).map_err(|e| format!("{}: {e}", p.display()))?)),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };
    Ok(match cfg.outfmt {
        OutputFormat::BlastM8 => Box::new(M8Writer::new(sink)),
        OutputFormat::BlastM9 => Box::new(M9Writer::new(sink)),
        OutputFormat::BlastM0 => Box::new(M0Writer::new(sink)),
        OutputFormat::SAM | OutputFormat::BAM => Box::new(SAMWriter::new(sink)),
        OutputFormat::JSON => Box::new(JSONWriter::new(sink)),
        OutputFormat::CSV => Box::new(CSVWriter::new(sink)),
    })
}

/// Obtain the index for a search: open `--db`, or build from `--database`
/// (in memory, or in `tmpdir` when the estimated in-memory size exceeds the
/// memory budget). Returns (index, build seconds, source description, temp path).
pub fn obtain_index(cfg: &DriverConfig) -> Result<(SeedIndex, f64, String, Option<PathBuf>), String> {
    if let Some(p) = &cfg.db_index {
        let t0 = Instant::now();
        let idx = SeedIndex::open(p).map_err(|e| e.to_string())?;
        if cfg.verify_index {
            idx.verify_checksums().map_err(|e| e.to_string())?;
            idx.verify_alphabet().map_err(|e| e.to_string())?;
        }
        let want_kind = kind_for(cfg.mode);
        if idx.kind() != want_kind {
            return Err(format!("index '{}' holds {:?} sequences but mode {} needs {:?}", p.display(), idx.kind(), cfg.mode, want_kind));
        }
        let want_tr = cfg.mode == SequenceMode::ProteinToNucleotide;
        if idx.is_translated() != want_tr {
            return Err(if want_tr {
                format!("mode pn needs a translated index (build it with --mode pn --makedb); '{}' is not translated", p.display())
            } else {
                format!("'{}' is a translated (pn) index; use --mode pn with it", p.display())
            });
        }
        return Ok((idx, t0.elapsed().as_secs_f64(), format!("mapped {}", p.display()), None));
    }
    let fasta = cfg.database.as_ref().ok_or("a search needs --db <INDEX> or --database <FASTA>")?;
    let fbytes = std::fs::metadata(fasta).map(|m| m.len()).unwrap_or(0);
    // In-memory index ≈ residues × (1 + 4·shapes·(6 for translated)) bytes
    // (residues only for the scan schedule).
    let shapes = if scan_eligible(cfg) && !cfg.index_postings { 0 } else { cfg.seed.shapes.len().max(1) as u64 };
    let factor = if cfg.mode == SequenceMode::ProteinToNucleotide { 2 } else { 1 };
    let est = fbytes * (1 + 4 * shapes) * factor;
    let spill = cfg.memory_budget.map_or(false, |b| est > b / 2);
    if spill {
        let tmp = cfg.tmpdir.join(format!("saber-{}-{}.sdx", std::process::id(), Instant::now().elapsed().as_nanos()));
        let (idx, t) = build_index(cfg, fasta, Some(&tmp))?;
        release_free_memory();
        if !cfg.quiet {
            eprintln!(
                "[saber] estimated in-memory index {:.1} MB exceeds half the memory budget; built a temporary on-disk index at {} (memory-mapped)",
                est as f64 / 1e6,
                tmp.display()
            );
        }
        return Ok((idx, t, format!("temporary on-disk index {}", tmp.display()), Some(tmp)));
    }
    let (idx, t) = build_index(cfg, fasta, None)?;
    release_free_memory();
    Ok((idx, t, format!("in-memory from {}", fasta.display()), None))
}

/// Return freed heap pages to the OS. glibc keeps memory freed by the index
/// build (per-partition masks, counting tables) in its arenas, where the
/// search threads (which allocate from their own arenas) cannot reuse it; it
/// would otherwise stay resident for the whole search.
fn release_free_memory() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    // SAFETY: malloc_trim only releases free memory; no pointers are invalidated.
    unsafe {
        libc::malloc_trim(0);
    }
}

/// Run a complete search.
pub fn run(cfg: &DriverConfig) -> Result<RunSummary, String> {
    let t_all = Instant::now();
    let (idx, build_s, source, tmp) = obtain_index(cfg)?;
    let res = run_with_index(cfg, &idx, build_s, &source);
    drop(idx);
    if let Some(t) = tmp {
        let _ = std::fs::remove_file(t);
    }
    let mut s = res?;
    s.total_s = t_all.elapsed().as_secs_f64();
    s.peak_rss_kb = peak_rss_kb();
    if let Some(p) = &cfg.stats_path {
        let j = serde_json::to_string_pretty(&s).map_err(|e| e.to_string())?;
        std::fs::write(p, j).map_err(|e| format!("{}: {e}", p.display()))?;
    }
    Ok(s)
}

/// Run a search against an already obtained index.
pub fn run_with_index(cfg: &DriverConfig, idx: &SeedIndex, build_s: f64, source: &str) -> Result<RunSummary, String> {
    let eng = make_engine(cfg, idx)?;
    let calibrated = eng.stats.params.calibrated;
    if !calibrated && !cfg.quiet {
        eprintln!("[saber] WARNING: {} — E-values and bit scores are NOT calibrated for these settings.", eng.stats.params.source);
    }
    let qpath = cfg.query.as_ref().ok_or("--query <FILE> is required for a search")?;
    let mut qs = FastaStream::open(qpath).map_err(|e| format!("query: {e}"))?;
    let mut writer = open_writer(cfg)?;
    writer.write_header().map_err(|e| e.to_string())?;
    let ct = CodonTable::new(cfg.genetic_code);
    let threads = cfg.threads.max(1);
    let use_scan = match cfg.batch_order {
        BatchOrder::Scan => {
            if !eng.scan_available() {
                return Err("--batch-order scan is not available for this mode/configuration (protein databases with single-hit seeding only)".into());
            }
            true
        }
        BatchOrder::Auto => scan_eligible(cfg) && eng.scan_available(),
        _ => false,
    };
    if !use_scan && !idx.has_postings() {
        return Err(format!(
            "the index '{}' has no seed postings (it was built for the scan schedule); search it with the default schedule, or rebuild it with --index-postings",
            source
        ));
    }
    let batch = if cfg.query_batch > 0 {
        cfg.query_batch
    } else if use_scan {
        1 << 22
    } else {
        (threads * 64).max(256)
    };
    // Bound the residues per batch by the memory budget (queries are small;
    // results dominate): ~1/16 of the budget, at least 1 M residues. The scan
    // schedule streams the database once per batch, so its batches are as
    // large as its seed tables allow (≈ 4 bytes per residue and shape) —
    // except against a small database, where a pass over it costs little and
    // smaller batches keep the seed tables (and the memory) small.
    let max_res = if use_scan {
        if cfg.scan_batch_residues > 0 {
            cfg.scan_batch_residues
        } else {
            let per_res = 16 + 4 * cfg.seed.shapes.len() as u64;
            let by_budget = cfg.memory_budget.map_or(4usize << 20, |b| ((b / 8 / per_res) as usize).clamp(1 << 18, 8 << 20));
            let by_db = (idx.total_residues() as usize / 4).clamp(1 << 19, 8 << 20);
            by_budget.min(by_db)
        }
    } else {
        cfg.memory_budget.map_or(64usize << 20, |b| ((b / 16) as usize).max(1 << 20))
    };
    let mut batches_scan = 0u64;
    let mut total = Counters::default();
    let (mut batches_query, mut batches_part, mut batches_intra) = (0u64, 0u64, 0u64);
    let partition_major_auto = idx.num_partitions() > 1
        && idx.is_mapped()
        && cfg.memory_budget.map_or(false, |b| idx.data_bytes() > b / 2);
    let mut nq = 0u64;
    let mut nskip = 0u64;
    let mut nal = 0u64;
    let t0 = Instant::now();
    loop {
        let recs = qs.next_batch(batch, max_res).map_err(|e| format!("query {}: {e}", qpath.display()))?;
        if recs.is_empty() {
            break;
        }
        let keep: Vec<&OwnedRecord> = recs
            .iter()
            .filter(|r| {
                let l = r.seq.len();
                let ok = cfg.min_length.map_or(true, |m| l >= m) && cfg.max_length.map_or(true, |m| l <= m);
                if !ok {
                    nskip += 1;
                }
                ok
            })
            .collect();
        if keep.is_empty() {
            continue;
        }
        // One job per query frame.
        let mut jobs: Vec<(usize, Option<usize>, Vec<u8>)> = Vec::new();
        for (ri, r) in keep.iter().enumerate() {
            for (f, codes) in query_frames(cfg.mode, r, &ct) {
                jobs.push((ri, f, codes));
            }
        }
        if use_scan {
            batches_scan += 1;
            let views: Vec<&[u8]> = jobs.iter().map(|(_, _, c)| c.as_slice()).collect();
            let per_job = eng.search_scan(&views);
            let mut grouped: Vec<Vec<(Option<usize>, Hit)>> = (0..keep.len()).map(|_| Vec::new()).collect();
            for ((ri, f, _), (hits, c)) in jobs.iter().zip(per_job) {
                total.add(&c);
                grouped[*ri].extend(hits.into_iter().map(|h| (*f, h)));
            }
            let results: Vec<Vec<AlignmentResult>> = grouped
                .into_par_iter()
                .zip(keep.par_iter())
                .map(|(all, r)| finish_record(&eng, idx, cfg, r, all, calibrated))
                .collect();
            for res in &results {
                for a in res {
                    writer.write_alignment(a).map_err(|e| e.to_string())?;
                }
                nal += res.len() as u64;
            }
            nq += keep.len() as u64;
            continue;
        }
        let order = if jobs.len() < 2 * threads && threads > 1 {
            batches_intra += 1;
            None
        } else if cfg.batch_order == BatchOrder::Partition || (cfg.batch_order == BatchOrder::Auto && partition_major_auto) {
            batches_part += 1;
            Some(BatchOrder::Partition)
        } else {
            batches_query += 1;
            Some(BatchOrder::Query)
        };
        let per_job: Vec<(Vec<Hit>, Counters)> = match order {
            None => jobs.iter().map(|(_, _, c)| eng.search_parallel(c)).collect(),
            Some(BatchOrder::Partition) => {
                let views: Vec<&[u8]> = jobs.iter().map(|(_, _, c)| c.as_slice()).collect();
                eng.search_batch(&views)
            }
            _ => {
                // Longest jobs first (LPT): work stealing then balances the
                // short ones around them instead of a long query finishing last.
                let mut order: Vec<usize> = (0..jobs.len()).collect();
                order.sort_by_key(|&i| std::cmp::Reverse(jobs[i].2.len()));
                let done: Vec<(usize, (Vec<Hit>, Counters))> =
                    order.par_iter().map(|&i| (i, search_local(&eng, &jobs[i].2))).collect();
                let mut slots: Vec<Option<(Vec<Hit>, Counters)>> = (0..jobs.len()).map(|_| None).collect();
                for (i, r) in done {
                    slots[i] = Some(r);
                }
                slots.into_iter().map(|r| r.expect("every job ran")).collect()
            }
        };
        let mut grouped: Vec<Vec<(Option<usize>, Hit)>> = (0..keep.len()).map(|_| Vec::new()).collect();
        for ((ri, f, _), (hits, c)) in jobs.iter().zip(per_job) {
            total.add(&c);
            grouped[*ri].extend(hits.into_iter().map(|h| (*f, h)));
        }
        let results: Vec<Vec<AlignmentResult>> = grouped
            .into_par_iter()
            .zip(keep.par_iter())
            .map(|(all, r)| finish_record(&eng, idx, cfg, r, all, calibrated))
            .collect();
        for res in &results {
            for a in res {
                writer.write_alignment(a).map_err(|e| e.to_string())?;
            }
            nal += res.len() as u64;
        }
        nq += keep.len() as u64;
    }
    writer.write_footer().map_err(|e| e.to_string())?;
    writer.flush().map_err(|e| e.to_string())?;
    let search_s = t0.elapsed().as_secs_f64();
    let mut limits = Vec::new();
    if total.regions_dropped_by_cap > 0 {
        limits.push(format!(
            "{} candidate regions were dropped by --index-max-candidates {} (heuristic limit; the search is not exhaustive)",
            total.regions_dropped_by_cap, eng.params.max_candidates
        ));
    }
    if total.seeds_frequent_skipped > 0 {
        limits.push(format!(
            "{} query seeds were skipped as over-frequent in the database (--seed-max-occ caps {:?})",
            total.seeds_frequent_skipped,
            if use_scan { eng.scan_caps() } else { eng.occ_caps().to_vec() }
        ));
    }
    if total.band_capped > 0 {
        limits.push(format!(
            "{} alignments hit the traceback memory cap before the band reproduced the window score (reported with exact=false)",
            total.band_capped
        ));
    }
    if !cfg.params.exhaustive {
        limits.push("seed-and-extend heuristic: alignments without a qualifying seed/HSP are not examined (use --exhaustive for full Smith-Waterman)".into());
    }
    let sc = if use_scan { cfg.seed.clone() } else { idx.seed_config().clone() };
    let (scan_mask, scan_db_s) = if use_scan { eng.scan_info() } else { (String::new(), 0.0) };
    Ok(RunSummary {
        saber_version: crate::VERSION.to_string(),
        mode: format!("{}", cfg.mode),
        preset: cfg.preset.clone(),
        threads,
        queries: nq,
        queries_skipped_length: nskip,
        alignments: nal,
        db_subjects: idx.num_subjects() as u64,
        db_residues: idx.total_residues(),
        index_bytes: idx.data_bytes(),
        index_source: source.to_string(),
        index_build_s: build_s,
        search_s,
        total_s: 0.0,
        peak_rss_kb: 0,
        evalue_params: eng.stats.params.source.clone(),
        evalue_calibrated: calibrated,
        seed_alphabet: sc.alphabet.name.clone(),
        seed_shapes: sc.shapes.iter().map(|s| s.pattern.clone()).collect(),
        masking: format!("db={} query={}", idx.seed_config().masking, if eng.params.query_masking == crate::masking::MaskingMode::Tantan { "tantan" } else { "none" }),
        comp_bias: eng.params.comp_bias,
        heuristic_limits: limits,
        score_backend: if total.backend_windows > 0 { "gpu".into() } else { "cpu".into() },
        memory_budget_bytes: cfg.memory_budget,
        batches_query_major: batches_query,
        batches_partition_major: batches_part,
        batches_intra_query: batches_intra,
        batches_scan,
        scan_mask,
        scan_db_s,
        counters: total,
    })
}

/// Build and save a persistent index (`--makedb`).
pub fn makedb(cfg: &DriverConfig, out: &Path) -> Result<(u64, f64), String> {
    let fasta = cfg.database.as_ref().ok_or("--makedb requires --database <FASTA>")?;
    let (idx, t) = build_index(cfg, fasta, Some(out))?;
    if cfg.verify_index {
        idx.verify_checksums().map_err(|e| e.to_string())?;
    }
    let sz = std::fs::metadata(out).map(|m| m.len()).unwrap_or(0);
    Ok((sz, t))
}
