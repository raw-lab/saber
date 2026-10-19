//! `saber doctor` — environment and self-integrity diagnostics.
//!
//! Run via `saber --doctor`. The doctor performs a battery of independent
//! checks and prints a readable report, then exits non-zero if any check
//! FAILED (warnings alone do not fail the run). It is meant to answer, in one
//! command, "is this SABER build healthy on this machine, and will a search
//! actually work?" — covering the things that bite people in practice:
//!
//!   * Is this the binary/version I think it is?
//!   * Does the CPU have the SIMD features the fast path needs (AVX2)?
//!   * Does the Smith-Waterman kernel produce the known-correct score on a
//!     hand-verified pair? (catches a miscompiled or mis-tuned build)
//!   * Does the indexed search agree with the brute-force path? (catches
//!     seeding/index regressions)
//!   * Does the persistent index save→open→query round-trip faithfully?
//!     (catches serialization / mmap / endianness bugs)
//!   * Are the scoring matrices loadable and symmetric?
//!   * Is the GPU backend available, and if requested, usable?
//!   * Can we actually write a temp file (needed for --makedb)?
//!   * Optionally: are a user-supplied --query / --database / --db readable
//!     and well-formed?
//!
//! Each check is deliberately small and self-contained so a failure points at
//! exactly one subsystem.

use std::fmt;
use std::path::Path;

use crate::algorithm::reference::sw_reference;
use crate::algorithm::sw_simd::{score_many, KernelCounters, SimdScratch};
use crate::algorithm::SmithWaterman;
use crate::alphabet::{EncodedProtein, ProteinCodes, SeedAlphabet, SeedShape};
use crate::masking::MaskingMode;
use crate::scoring::stats::{gapped_protein, GapConvention, GapCosts, SearchStatistics};
use crate::scoring::{ScoreMatrix, ScoringSystem};
use crate::search::engine::{Engine, SearchParams, Worker};
use crate::search::index::{BuildOptions, IndexBuilder, SeedConfig, SeedIndex};

/// Outcome of a single diagnostic check.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Pass,
    Warn,
    Fail,
    /// Skipped — not applicable on this build/host (e.g. GPU check without a
    /// device). Never counts against the exit code.
    Skip,
}

impl Status {
    fn glyph(self) -> &'static str {
        match self {
            Status::Pass => "[ OK ]",
            Status::Warn => "[WARN]",
            Status::Fail => "[FAIL]",
            Status::Skip => "[SKIP]",
        }
    }
}

impl fmt::Display for Status {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.glyph())
    }
}

/// One line of the report.
pub struct Check {
    pub name: String,
    pub status: Status,
    pub detail: String,
}

impl Check {
    fn new(name: impl Into<String>, status: Status, detail: impl Into<String>) -> Self {
        Check { name: name.into(), status, detail: detail.into() }
    }
}

/// Accumulates checks and tracks whether any hard-failed.
#[derive(Default)]
pub struct Report {
    pub checks: Vec<Check>,
}

impl Report {
    fn push(&mut self, name: impl Into<String>, status: Status, detail: impl Into<String>) {
        self.checks.push(Check::new(name, status, detail));
    }

    pub fn n_pass(&self) -> usize { self.checks.iter().filter(|c| c.status == Status::Pass).count() }
    pub fn n_warn(&self) -> usize { self.checks.iter().filter(|c| c.status == Status::Warn).count() }
    pub fn n_fail(&self) -> usize { self.checks.iter().filter(|c| c.status == Status::Fail).count() }
    pub fn n_skip(&self) -> usize { self.checks.iter().filter(|c| c.status == Status::Skip).count() }

    /// True if every check passed (warnings allowed, no failures).
    pub fn healthy(&self) -> bool {
        self.n_fail() == 0
    }

    /// Render the full report to a String (so the caller controls stdout vs
    /// a buffer; tests render to a buffer).
    pub fn render(&self) -> String {
        let mut out = String::new();
        out.push_str("SABER doctor — self-diagnostics\n");
        out.push_str("================================\n\n");
        // Align the check names for readability.
        let width = self.checks.iter().map(|c| c.name.len()).max().unwrap_or(0);
        for c in &self.checks {
            out.push_str(&format!("{}  {:<width$}  {}\n", c.status, c.name, c.detail, width = width));
        }
        out.push('\n');
        out.push_str(&format!(
            "Summary: {} passed, {} warning(s), {} failed, {} skipped.\n",
            self.n_pass(), self.n_warn(), self.n_fail(), self.n_skip(),
        ));
        if self.healthy() {
            if self.n_warn() == 0 {
                out.push_str("Result: HEALTHY — everything checks out.\n");
            } else {
                out.push_str("Result: HEALTHY (with warnings) — searches will work; see warnings above.\n");
            }
        } else {
            out.push_str("Result: PROBLEMS DETECTED — see [FAIL] lines above.\n");
        }
        out
    }
}

/// Options controlling which optional checks run.
#[derive(Default)]
pub struct DoctorOptions<'a> {
    /// If set, validate this query FASTA is readable and parses.
    pub query: Option<&'a Path>,
    /// If set, validate this database FASTA is readable and parses.
    pub database: Option<&'a Path>,
    /// If set, validate this persistent index opens and is self-consistent.
    pub db_index: Option<&'a Path>,
    /// Whether the user asked for the GPU backend (`--gpu`).
    pub gpu_requested: bool,
    /// Whether the binary was compiled with a CUDA backend.
    pub cuda_compiled: bool,
    /// Number of worker threads the run would use.
    pub threads: usize,
}

/// Run all diagnostics and return the assembled report.
pub fn run(opts: &DoctorOptions) -> Report {
    let mut r = Report::default();

    check_build(&mut r, opts);
    check_cpu_features(&mut r);
    check_scoring_matrices(&mut r);
    check_statistics(&mut r);
    check_sw_kernel(&mut r);
    check_simd_vs_reference(&mut r);
    check_indexed_vs_exhaustive(&mut r);
    check_persistent_index(&mut r);
    check_temp_writable(&mut r);
    check_gpu(&mut r, opts);
    check_user_files(&mut r, opts);

    r
}

// ---------------------------------------------------------------------------
// Individual checks
// ---------------------------------------------------------------------------

fn check_build(r: &mut Report, opts: &DoctorOptions) {
    r.push(
        "build / version",
        Status::Pass,
        format!(
            "SABER v{} ({}-bit, {} target, {} thread(s) configured, index format v{})",
            crate::VERSION,
            (std::mem::size_of::<usize>() * 8),
            std::env::consts::ARCH,
            opts.threads.max(1),
            crate::search::index::FORMAT_VERSION,
        ),
    );
}

fn check_cpu_features(r: &mut Report) {
    #[cfg(target_arch = "x86_64")]
    {
        let avx2 = is_x86_feature_detected!("avx2");
        let popcnt = is_x86_feature_detected!("popcnt");
        let bmi1 = is_x86_feature_detected!("bmi1");
        if avx2 && popcnt && bmi1 {
            r.push("CPU SIMD (AVX2)", Status::Pass, "AVX2/POPCNT/BMI1 present — SIMD seed filter, scoring and banded traceback active");
        } else {
            r.push(
                "CPU SIMD (AVX2)",
                Status::Warn,
                format!("AVX2={avx2} POPCNT={popcnt} BMI1={bmi1}: scalar fallbacks are used (identical results, much slower)"),
            );
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        r.push("CPU SIMD (AVX2)", Status::Skip, format!("non-x86_64 architecture ({}); scalar kernels are used", std::env::consts::ARCH));
    }
}

fn check_scoring_matrices(r: &mut Report) {
    let systems = [
        ("BLOSUM62", ScoringSystem::Blosum62),
        ("BLOSUM45", ScoringSystem::Blosum45),
        ("BLOSUM50", ScoringSystem::Blosum50),
        ("BLOSUM80", ScoringSystem::Blosum80),
        ("BLOSUM90", ScoringSystem::Blosum90),
        ("PAM30", ScoringSystem::Pam30),
        ("PAM70", ScoringSystem::Pam70),
        ("PAM250", ScoringSystem::Pam250),
    ];
    let mut bad: Vec<String> = Vec::new();
    for (name, sys) in systems {
        let m = ScoreMatrix::new(sys);
        let symmetric = (0..24).all(|i| (0..24).all(|j| m.matrix[i][j] == m.matrix[j][i]));
        if !symmetric {
            bad.push(name.to_string());
        }
    }
    if bad.is_empty() {
        r.push("scoring matrices", Status::Pass, format!("{} built-in matrices load and are symmetric", systems.len()));
    } else {
        r.push("scoring matrices", Status::Fail, format!("non-symmetric matrix/matrices: {}", bad.join(", ")));
    }
}

fn check_statistics(r: &mut Report) {
    // BLOSUM62 with BLAST gap costs 11/1 (= SABER 12/1) has published gapped
    // parameters lambda=0.267, K=0.041 (NCBI blast_stat.c).
    let g = GapCosts::from_user(11, 1, GapConvention::Blast).expect("valid");
    match gapped_protein(ScoringSystem::Blosum62, g, false) {
        Ok(k) if (k.lambda - 0.267).abs() < 1e-9 && (k.k - 0.041).abs() < 1e-9 && k.calibrated => {
            r.push("gapped statistics", Status::Pass, "BLOSUM62 11/1 (BLAST convention): lambda=0.267 K=0.041 from NCBI tables");
        }
        Ok(k) => r.push("gapped statistics", Status::Fail, format!("unexpected parameters lambda={} K={}", k.lambda, k.k)),
        Err(e) => r.push("gapped statistics", Status::Fail, e.to_string()),
    }
}

fn check_sw_kernel(r: &mut Report) {
    // Known-answer test on the scalar aligner.
    let sw = SmithWaterman::nucleotide(10, 1, 2, -3);
    let res = sw.align(b"ACGTACGT", b"ACGTACGT");
    let unrel = sw.align(b"ACGTACGT", b"TTTTTTTT");
    if res.score == 16 && res.identity > 99.9 && res.cigar == "8M" && res.score > unrel.score {
        r.push("Smith-Waterman kernel", Status::Pass, "known-answer test passes (score=16, id=100%, CIGAR=8M)");
    } else {
        r.push(
            "Smith-Waterman kernel",
            Status::Fail,
            format!("known-answer MISMATCH: score={} id={:.1}% cigar={:?}", res.score, res.identity, res.cigar),
        );
    }
}

fn rng_seq(state: &mut u64, len: usize, alpha: u64) -> Vec<u8> {
    (0..len)
        .map(|_| {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            (*state % alpha) as u8
        })
        .collect()
}

fn check_simd_vs_reference(r: &mut Report) {
    // The production kernel (16-bit SIMD with exact 32-bit repair) against
    // the independent i64 full-matrix reference, including pairs whose score
    // exceeds the 16-bit range.
    let m = ScoreMatrix::new(ScoringSystem::Blosum62);
    let gaps = GapCosts { open: 11, extend: 1 };
    let mut st = 0x5eed_u64;
    let q = rng_seq(&mut st, 180, 20);
    let mut subs: Vec<Vec<u8>> = (0..40).map(|i| rng_seq(&mut st, 20 + i * 7, 20)).collect();
    let long = vec![17u8; 3500]; // W×3500: score 38500 > i16::MAX
    subs.push(long.clone());
    let views: Vec<ProteinCodes<'_>> = subs.iter().map(|s| ProteinCodes::new_unchecked(s)).collect();
    let mut sc = SimdScratch::new();
    let mut kc = KernelCounters::default();
    let fast = score_many(ProteinCodes::new_unchecked(&q), None, &views[..40], &m.matrix, gaps, &mut sc, &mut kc);
    let lq = vec![17u8; 3500];
    let big = score_many(ProteinCodes::new_unchecked(&lq), None, &views[40..], &m.matrix, gaps, &mut sc, &mut kc);
    let mut bad = 0;
    for (i, s) in subs[..40].iter().enumerate() {
        if sw_reference(&q, s, &m.matrix, None, 11, 1).0 != fast[i] as i64 {
            bad += 1;
        }
    }
    let big_ok = big[0] == 38500;
    if bad == 0 && big_ok {
        r.push(
            "SIMD kernel vs reference",
            Status::Pass,
            format!("41 pairs identical to the i64 reference ({} saturated lane(s) repaired exactly)", kc.saturation_retries),
        );
    } else {
        r.push("SIMD kernel vs reference", Status::Fail, format!("{bad} mismatching pair(s); long pair score {} (want 38500)", big[0]));
    }
}

fn tiny_db() -> Vec<(String, Vec<u8>)> {
    let mut st = 0xd0c7_u64;
    let aa = b"ARNDCQEGHILKMFPSTWYV";
    let rand = |st: &mut u64, n: usize| -> Vec<u8> { rng_seq(st, n, 20).into_iter().map(|c| aa[c as usize]).collect() };
    let q = b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR".to_vec();
    let mut db = Vec::new();
    db.push(("homolog_identical".to_string(), q.clone()));
    let mut mutated = q.clone();
    for (i, c) in mutated.iter_mut().enumerate() {
        if i % 4 == 1 {
            *c = aa[(i * 7) % 20];
        }
    }
    db.push(("homolog_75pct".to_string(), mutated));
    for i in 0..30 {
        db.push((format!("decoy_{i}"), rand(&mut st, 60 + i * 9)));
    }
    db.push(("query".to_string(), q));
    db
}

fn build_tiny(db: &[(String, Vec<u8>)], to: Option<&Path>) -> Result<SeedIndex, String> {
    let mut cfg = SeedConfig::protein_default();
    cfg.alphabet = SeedAlphabet::murphy10();
    cfg.shapes = vec![SeedShape::parse("110111")?];
    cfg.masking = MaskingMode::Tantan;
    let opt = BuildOptions::new(cfg);
    let mut b = match to {
        Some(p) => IndexBuilder::to_file(opt, p).map_err(|e| e.to_string())?,
        None => IndexBuilder::in_memory(opt).map_err(|e| e.to_string())?,
    };
    let mut scratch = Vec::new();
    for (id, s) in db {
        b.add_ascii(id, s, &mut scratch).map_err(|e| e.to_string())?;
    }
    match to {
        Some(p) => {
            b.finish_to_file().map_err(|e| e.to_string())?;
            SeedIndex::open(p).map_err(|e| e.to_string())
        }
        None => b.finish_in_memory().map_err(|e| e.to_string()),
    }
}

fn tiny_engine(idx: &SeedIndex, exhaustive: bool) -> Engine<'_> {
    let gaps = GapCosts::from_user(11, 1, GapConvention::Saber).expect("valid");
    let kp = gapped_protein(ScoringSystem::Blosum62, gaps, false).expect("calibrated");
    let st = SearchStatistics::new(kp, idx.total_residues(), idx.num_subjects() as u64);
    let mut p = SearchParams::default();
    p.max_target_seqs = 100;
    p.evalue = 1e-3;
    p.exhaustive = exhaustive;
    p.hit_filter_min_ids = 12;
    Engine::new(idx, ScoreMatrix::new(ScoringSystem::Blosum62), gaps, st, p)
}

fn check_indexed_vs_exhaustive(r: &mut Report) {
    let db = tiny_db();
    let idx = match build_tiny(&db, None) {
        Ok(i) => i,
        Err(e) => {
            r.push("indexed search vs exhaustive", Status::Fail, format!("index build failed: {e}"));
            return;
        }
    };
    let q = EncodedProtein::from_ascii(&db[db.len() - 1].1);
    let mut w = Worker::new();
    let fast = tiny_engine(&idx, false).search(q.as_slice(), &mut w);
    let full = tiny_engine(&idx, true).search(q.as_slice(), &mut w);
    let key = |h: &crate::search::engine::Hit| (h.subject, h.score, h.q_start, h.s_start);
    let same = fast.iter().map(key).collect::<Vec<_>>() == full.iter().map(key).collect::<Vec<_>>();
    let names: Vec<&str> = fast.iter().map(|h| idx.subject_id(h.subject as usize)).collect();
    let found = names.contains(&"homolog_identical") && names.contains(&"homolog_75pct") && names.contains(&"query");
    let decoy = names.iter().any(|n| n.starts_with("decoy"));
    if same && found && !decoy {
        r.push("indexed search vs exhaustive", Status::Pass, "seeded search = exhaustive SW on a planted-homolog test (3 hits, no decoys at E<=1e-3)");
    } else {
        r.push(
            "indexed search vs exhaustive",
            Status::Fail,
            format!("identical={same} homologs_found={found} decoy_reported={decoy} (hits: {names:?})"),
        );
    }
}

fn check_persistent_index(r: &mut Report) {
    let db = tiny_db();
    let mut path = std::env::temp_dir();
    path.push(format!("saber_doctor_{}_{}.sdx", std::process::id(), now_nanos()));
    let res = (|| -> Result<(bool, u64), String> {
        let mem = build_tiny(&db, None)?;
        let mapped = build_tiny(&db, Some(&path))?;
        mapped.verify_checksums().map_err(|e| e.to_string())?;
        mapped.verify_alphabet().map_err(|e| e.to_string())?;
        let q = EncodedProtein::from_ascii(&db[db.len() - 1].1);
        let mut w = Worker::new();
        let a = tiny_engine(&mem, false).search(q.as_slice(), &mut w);
        let b = tiny_engine(&mapped, false).search(q.as_slice(), &mut w);
        let same = a.len() == b.len()
            && a.iter().zip(&b).all(|(x, y)| x.subject == y.subject && x.score == y.score && x.cigar() == y.cigar())
            && (0..mem.num_subjects()).all(|i| mem.subject_id(i) == mapped.subject_id(i) && mem.subject_codes(i) == mapped.subject_codes(i));
        Ok((same, std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0)))
    })();
    let _ = std::fs::remove_file(&path);
    match res {
        Ok((true, sz)) => r.push("persistent index round-trip", Status::Pass, format!("build→write→mmap→verify→search is faithful ({sz} bytes)")),
        Ok((false, _)) => r.push("persistent index round-trip", Status::Fail, "mapped index search differs from the in-memory index"),
        Err(e) => r.push("persistent index round-trip", Status::Fail, e),
    }
}

fn check_temp_writable(r: &mut Report) {
    // --makedb needs to write an index file. Confirm the temp dir is writable
    // so that failure mode is reported here rather than mid-build.
    let dir = std::env::temp_dir();
    let mut path = dir.clone();
    path.push(format!("saber_doctor_wtest_{}_{}", std::process::id(), now_nanos()));
    match std::fs::write(&path, b"saber-doctor-write-test") {
        Ok(()) => {
            let _ = std::fs::remove_file(&path);
            r.push("temp dir writable", Status::Pass, format!("{} is writable", dir.display()));
        }
        Err(e) => {
            r.push(
                "temp dir writable",
                Status::Warn,
                format!("cannot write to {}: {} (--makedb may fail; set TMPDIR)", dir.display(), e),
            );
        }
    }
}

fn check_gpu(r: &mut Report, opts: &DoctorOptions) {
    if !opts.cuda_compiled {
        let status = if opts.gpu_requested { Status::Warn } else { Status::Skip };
        let detail = if opts.gpu_requested {
            "--gpu requested but this binary was built WITHOUT --features gpu-cuda; searches run on the CPU"
        } else {
            "binary built without --features gpu-cuda (CPU-only build)"
        };
        r.push("GPU backend", status, detail.to_string());
        return;
    }
    let Some(be) = crate::gpu::try_init_backend() else {
        let status = if opts.gpu_requested { Status::Warn } else { Status::Skip };
        r.push("GPU backend", status, "CUDA compiled in but no device initialised; searches run on the CPU".to_string());
        return;
    };
    // Score parity: the GPU must reproduce the CPU kernel exactly.
    let m = ScoreMatrix::new(ScoringSystem::Blosum62);
    let gaps = GapCosts { open: 11, extend: 1 };
    let mut st = 0x6e0_u64;
    let q = rng_seq(&mut st, 250, 24);
    let bias: Vec<i8> = (0..q.len()).map(|i| (i % 5) as i8 - 2).collect();
    let wins: Vec<Vec<u8>> = (0..300).map(|i| rng_seq(&mut st, 1 + (i * 13) % 700, 24)).collect();
    let views: Vec<&[u8]> = wins.iter().map(|w| w.as_slice()).collect();
    let pv: Vec<ProteinCodes<'_>> = wins.iter().map(|w| ProteinCodes::new_unchecked(w)).collect();
    let mut sc = SimdScratch::new();
    let mut kc = KernelCounters::default();
    let cpu = score_many(ProteinCodes::new_unchecked(&q), Some(&bias), &pv, &m.matrix, gaps, &mut sc, &mut kc);
    match be.score_windows(&q, Some(&bias), &views, &m.matrix, gaps) {
        Ok(g) if g == cpu => r.push("GPU backend", Status::Pass, format!("{}: 300 windows score-identical to the CPU kernel", be.name())),
        Ok(g) => {
            let n = g.iter().zip(&cpu).filter(|(a, b)| a != b).count();
            r.push("GPU backend", Status::Fail, format!("{}: {n} of 300 scores differ from the CPU kernel — do not use --gpu", be.name()))
        }
        Err(e) => r.push("GPU backend", Status::Fail, format!("{}: scoring failed: {e}", be.name())),
    }
}

fn check_user_files(r: &mut Report, opts: &DoctorOptions) {
    use crate::io::Fasta;

    if let Some(q) = opts.query {
        match Fasta::from_file(q) {
            Ok(f) => {
                let n = f.records().len();
                if n == 0 {
                    r.push("query FASTA", Status::Fail, format!("{}: parsed but contains 0 sequences", q.display()));
                } else {
                    let res: usize = f.records().iter().map(|x| x.length).sum();
                    r.push("query FASTA", Status::Pass, format!("{}: {} sequence(s), {} residues", q.display(), n, res));
                }
            }
            Err(e) => r.push("query FASTA", Status::Fail, format!("{}: {}", q.display(), e)),
        }
    }

    if let Some(d) = opts.database {
        match Fasta::from_file(d) {
            Ok(f) => {
                let n = f.records().len();
                if n == 0 {
                    r.push("database FASTA", Status::Fail, format!("{}: parsed but contains 0 sequences", d.display()));
                } else {
                    let res: usize = f.records().iter().map(|x| x.length).sum();
                    r.push("database FASTA", Status::Pass, format!("{}: {} sequence(s), {} residues", d.display(), n, res));
                }
            }
            Err(e) => r.push("database FASTA", Status::Fail, format!("{}: {}", d.display(), e)),
        }
    }

    if let Some(idx_path) = opts.db_index {
        match SeedIndex::open(idx_path).and_then(|idx| {
            idx.verify_checksums()?;
            idx.verify_alphabet()?;
            Ok(idx)
        }) {
            Ok(idx) => r.push(
                "persistent --db index",
                Status::Pass,
                format!(
                    "{}: valid v{} index ({:?}, {} subjects, {} residues, {} partition(s), seeds {} {:?}, checksums OK)",
                    idx_path.display(),
                    idx.meta().format_version,
                    idx.kind(),
                    idx.num_subjects(),
                    idx.total_residues(),
                    idx.num_partitions(),
                    idx.seed_config().alphabet.name,
                    idx.seed_config().shapes.iter().map(|s| s.pattern.clone()).collect::<Vec<_>>()
                ),
            ),
            Err(e) => r.push("persistent --db index", Status::Fail, format!("{}: {e}", idx_path.display())),
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Monotonic-ish nanosecond tag for unique temp filenames. Falls back to 0 if
/// the clock is unavailable (never happens in practice).
fn now_nanos() -> u128 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_nanos()).unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn doctor_core_checks_pass_in_this_build() {
        // The build-independent core checks (kernel, matrices, index,
        // persistence) must all pass in a healthy build. We exclude the
        // host-dependent ones (CPU features, GPU, user files) from the
        // assertion since those legitimately vary.
        let mut r = Report::default();
        check_scoring_matrices(&mut r);
        check_statistics(&mut r);
        check_sw_kernel(&mut r);
        check_simd_vs_reference(&mut r);
        check_indexed_vs_exhaustive(&mut r);
        check_persistent_index(&mut r);
        check_temp_writable(&mut r);
        assert_eq!(r.n_fail(), 0, "core doctor checks must not fail:\n{}", r.render());
    }

    #[test]
    fn report_renders_and_summarizes() {
        let mut r = Report::default();
        r.push("a", Status::Pass, "ok");
        r.push("b", Status::Warn, "careful");
        let s = r.render();
        assert!(s.contains("[ OK ]"));
        assert!(s.contains("[WARN]"));
        assert!(s.contains("1 passed"));
        assert!(r.healthy()); // warnings don't fail
        r.push("c", Status::Fail, "broken");
        assert!(!r.healthy());
    }
}
