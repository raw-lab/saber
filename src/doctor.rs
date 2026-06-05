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

use crate::algorithm::index::KmerIndex;
use crate::algorithm::SmithWaterman;
use crate::scoring::{ScoreMatrix, ScoringSystem};

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
    check_sw_kernel(&mut r);
    check_indexed_vs_bruteforce(&mut r);
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
            "SABER v{} ({}-bit, {} target, {} thread(s) configured)",
            crate::VERSION,
            (std::mem::size_of::<usize>() * 8),
            std::env::consts::ARCH,
            opts.threads.max(1),
        ),
    );
}

fn check_cpu_features(r: &mut Report) {
    // AVX2 drives the inter-sequence SIMD Smith-Waterman; without it the code
    // still runs (scalar fallback) but ~10-30× slower. This is the single
    // most common "why is SABER slow on this box" cause, so we surface it
    // explicitly. Non-x86 architectures skip (the scalar path is correct).
    #[cfg(target_arch = "x86_64")]
    {
        let avx2 = is_x86_feature_detected!("avx2");
        let avx512 = is_x86_feature_detected!("avx512f");
        let sse41 = is_x86_feature_detected!("sse4.1");
        if avx2 {
            r.push(
                "CPU SIMD (AVX2)",
                Status::Pass,
                format!(
                    "AVX2 present (sse4.1={}, avx512f={}) — fast SIMD SW path active",
                    sse41, avx512
                ),
            );
        } else {
            r.push(
                "CPU SIMD (AVX2)",
                Status::Warn,
                "AVX2 NOT detected — SABER will use the scalar SW fallback (much slower). \
                 Check that the binary was not built with a restrictive target-cpu, and that \
                 the host actually supports AVX2."
                    .to_string(),
            );
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        r.push(
            "CPU SIMD (AVX2)",
            Status::Skip,
            format!("non-x86_64 architecture ({}); scalar SW path is used", std::env::consts::ARCH),
        );
    }
}

fn check_scoring_matrices(r: &mut Report) {
    // Load each built-in matrix and verify it is square-symmetric over the
    // 20 standard residues (a corrupted/edited matrix would break scoring).
    let systems = [
        ("BLOSUM62", ScoringSystem::Blosum62),
        ("BLOSUM45", ScoringSystem::Blosum45),
        ("BLOSUM80", ScoringSystem::Blosum80),
        ("PAM30", ScoringSystem::Pam30),
        ("PAM70", ScoringSystem::Pam70),
        ("PAM250", ScoringSystem::Pam250),
    ];
    let mut bad: Vec<String> = Vec::new();
    let mut loaded = 0usize;
    for (name, sys) in systems {
        let m = ScoreMatrix::new(sys);
        loaded += 1;
        // Symmetry over the 20 standard amino acids (indices 0..20).
        let mut symmetric = true;
        'outer: for i in 0..20 {
            for j in 0..20 {
                if m.matrix[i][j] != m.matrix[j][i] {
                    symmetric = false;
                    break 'outer;
                }
            }
        }
        if !symmetric {
            bad.push(name.to_string());
        }
    }
    if bad.is_empty() {
        r.push(
            "scoring matrices",
            Status::Pass,
            format!("{} built-in matrices load and are symmetric", loaded),
        );
    } else {
        r.push(
            "scoring matrices",
            Status::Fail,
            format!("non-symmetric matrix/matrices: {}", bad.join(", ")),
        );
    }
}

fn check_sw_kernel(r: &mut Report) {
    // Known-answer test. Two identical 8-residue nucleotide sequences with
    // match=+2 must score exactly 8*2 = 16 with 100% identity and an "8M"
    // CIGAR. A miscompiled DP, a broken traceback, or a corrupted match score
    // all show up here.
    let sw = SmithWaterman::nucleotide(10, 1, 2, -3);
    let res = sw.align(b"ACGTACGT", b"ACGTACGT");
    let score_ok = res.score == 16;
    let id_ok = res.identity > 99.9;
    let cigar_ok = res.cigar == "8M";

    // Self-vs-shifted: an identical pair must out-score an unrelated pair.
    let unrel = sw.align(b"ACGTACGT", b"TTTTTTTT");
    let order_ok = res.score > unrel.score;

    if score_ok && id_ok && cigar_ok && order_ok {
        r.push(
            "Smith-Waterman kernel",
            Status::Pass,
            "known-answer test passes (score=16, id=100%, CIGAR=8M, ordering correct)",
        );
    } else {
        r.push(
            "Smith-Waterman kernel",
            Status::Fail,
            format!(
                "known-answer MISMATCH: score={} (want 16), id={:.1}% (want 100), cigar={:?} (want 8M), ordering_ok={}",
                res.score, res.identity, res.cigar, order_ok
            ),
        );
    }
}

fn check_indexed_vs_bruteforce(r: &mut Report) {
    // The indexed seed→extend path must find the obvious homolog that a
    // brute-force all-vs-all SW would. We build a tiny protein DB where two
    // subjects are clear matches to the query and one is junk, then confirm
    // the index produces candidates for the matches and not the junk.
    let subjects: Vec<&[u8]> = vec![
        b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF",   // 0: matches query
        b"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW",   // 1: junk
        b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF",   // 2: dup of 0
    ];
    let idx = KmerIndex::build_protein(&subjects, 4);
    let query = b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF";
    let hits = idx.query_protein(query);
    let cands = idx.candidates(&hits, 2, 16, 100);
    let has0 = cands.iter().any(|c| c.subject_idx == 0);
    let has2 = cands.iter().any(|c| c.subject_idx == 2);
    let has_junk = cands.iter().any(|c| c.subject_idx == 1);

    if has0 && has2 && !has_junk {
        r.push(
            "indexed search vs brute force",
            Status::Pass,
            "seed/candidate stage finds both true homologs and rejects the decoy",
        );
    } else {
        r.push(
            "indexed search vs brute force",
            Status::Fail,
            format!(
                "candidate mismatch: found subj0={}, subj2={}, junk_wrongly_included={}",
                has0, has2, has_junk
            ),
        );
    }
}

fn check_persistent_index(r: &mut Report) {
    // Build → save → mmap-open → query, and confirm the mapped index returns
    // exactly the same seed hits and subject metadata as the in-memory one.
    // This exercises the serialization layout, the 8-byte alignment, and the
    // mmap reinterpretation end-to-end.
    let subjects: Vec<&[u8]> = vec![
        b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF",
        b"GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHA",
        b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
    ];
    let ids = vec!["sp|A".to_string(), "sp|B".to_string(), "tr|C".to_string()];
    let lens: Vec<u32> = subjects.iter().map(|s| s.len() as u32).collect();

    let mut mem = KmerIndex::build_protein_windowed(&subjects, 4, 2);
    mem.set_subject_meta(&ids, &lens);

    let mut path = std::env::temp_dir();
    path.push(format!("saber_doctor_{}_{}.sdx", std::process::id(), now_nanos()));

    let save_res = mem.save(&path);
    if let Err(e) = save_res {
        r.push("persistent index round-trip", Status::Fail, format!("save failed: {}", e));
        return;
    }

    let mapped = match KmerIndex::open(&path) {
        Ok(m) => m,
        Err(e) => {
            let _ = std::fs::remove_file(&path);
            r.push("persistent index round-trip", Status::Fail, format!("open failed: {}", e));
            return;
        }
    };

    let q = b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF";
    let h_mem = mem.query_protein(q);
    let h_map = mapped.query_protein(q);

    let mut ok = h_mem.len() == h_map.len()
        && mem.num_subjects() == mapped.num_subjects()
        && mem.num_postings() == mapped.num_postings()
        && mapped.has_subject_meta();
    if ok {
        for i in 0..subjects.len() {
            if mem.subject_id(i) != mapped.subject_id(i)
                || mem.subject_len(i) != mapped.subject_len(i)
                || mem.encoded_subject(i) != mapped.encoded_subject(i)
            {
                ok = false;
                break;
            }
        }
    }
    if ok {
        for (a, b) in h_mem.iter().zip(h_map.iter()) {
            if a.subject_idx != b.subject_idx || a.q_pos != b.q_pos || a.s_pos != b.s_pos {
                ok = false;
                break;
            }
        }
    }

    let sz = std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0);
    let _ = std::fs::remove_file(&path);

    if ok {
        r.push(
            "persistent index round-trip",
            Status::Pass,
            format!("save→mmap-open→query is faithful ({} bytes, {} subjects)", sz, mapped.num_subjects()),
        );
    } else {
        r.push(
            "persistent index round-trip",
            Status::Fail,
            "mmap'd index disagrees with in-memory index (serialization/mmap bug)",
        );
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
    // Report the GPU situation honestly: whether CUDA was compiled in, and
    // (only if requested) whether a device actually initializes.
    if !opts.cuda_compiled {
        let status = if opts.gpu_requested { Status::Warn } else { Status::Skip };
        let detail = if opts.gpu_requested {
            "--gpu requested but this binary was built WITHOUT --features gpu-cuda; \
             rebuild with the CUDA toolkit present to enable GPU. Running on CPU."
                .to_string()
        } else {
            "binary built without --features gpu-cuda (CPU-only build)".to_string()
        };
        r.push("GPU backend", status, detail);
        return;
    }

    // CUDA was compiled in — try to bring up a device.
    match crate::gpu::try_init_gpu() {
        Some(exec) => r.push(
            "GPU backend",
            Status::Pass,
            format!("CUDA device available (backend: {})", exec.backend_name()),
        ),
        None => {
            let status = if opts.gpu_requested { Status::Warn } else { Status::Skip };
            r.push(
                "GPU backend",
                status,
                "CUDA compiled in but no usable device initialized; SABER will run on CPU".to_string(),
            );
        }
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
        match KmerIndex::open(idx_path) {
            Ok(idx) => {
                if !idx.has_subject_meta() {
                    r.push(
                        "persistent --db index",
                        Status::Warn,
                        format!("{}: opens but lacks subject metadata; rebuild with current SABER", idx_path.display()),
                    );
                } else {
                    r.push(
                        "persistent --db index",
                        Status::Pass,
                        format!(
                            "{}: valid (k={}, {} subjects, {} postings)",
                            idx_path.display(), idx.k(), idx.num_subjects(), idx.num_postings()
                        ),
                    );
                }
            }
            Err(e) => r.push("persistent --db index", Status::Fail, format!("{}: {}", idx_path.display(), e)),
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
        check_sw_kernel(&mut r);
        check_indexed_vs_bruteforce(&mut r);
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
