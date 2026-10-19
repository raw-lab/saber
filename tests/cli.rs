//! End-to-end CLI regression tests: pathological inputs, determinism across
//! thread counts and schedules, long alignments, translated coordinates and
//! index validation.

use std::path::{Path, PathBuf};
use std::process::Command;

fn saber() -> Command {
    Command::new(env!("CARGO_BIN_EXE_saber"))
}

struct Tmp(PathBuf);
impl Tmp {
    fn new(tag: &str) -> Self {
        let d = std::env::temp_dir().join(format!("saber_cli_{tag}_{}_{}", std::process::id(), nanos()));
        std::fs::create_dir_all(&d).unwrap();
        Tmp(d)
    }
    fn file(&self, name: &str, content: &str) -> PathBuf {
        let p = self.0.join(name);
        std::fs::write(&p, content).unwrap();
        p
    }
    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}
impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn nanos() -> u128 {
    std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn protein(&mut self, n: usize) -> String {
        let aa = b"ARNDCQEGHILKMFPSTWYV";
        (0..n).map(|_| aa[(self.next() % 20) as usize] as char).collect()
    }
    fn dna(&mut self, n: usize) -> String {
        (0..n).map(|_| b"ACGT"[(self.next() % 4) as usize] as char).collect()
    }
}

fn run_ok(args: &[&str]) -> String {
    let out = saber().args(args).output().expect("run saber");
    assert!(out.status.success(), "saber {:?} failed:\n{}", args, String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

fn run_err(args: &[&str]) -> String {
    let out = saber().args(args).output().expect("run saber");
    assert!(!out.status.success(), "saber {:?} unexpectedly succeeded", args);
    String::from_utf8_lossy(&out.stderr).to_string()
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

fn json_stats(p: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

/// A small protein DB with planted homologs of the queries.
fn protein_case(t: &Tmp, nq: usize) -> (PathBuf, PathBuf) {
    let mut r = Rng(0xfeed);
    let mut q = String::new();
    let mut db = String::new();
    for i in 0..nq {
        let n = 80 + (r.next() % 300) as usize;
        let p = r.protein(n);
        q += &format!(">q{i}\n{p}\n");
        // homolog: every 3rd residue substituted
        let h: String = p.chars().enumerate().map(|(k, c)| if k % 3 == 1 { 'A' } else { c }).collect();
        db += &format!(">h{i}\n{}{}{}\n", r.protein(40), h, r.protein(25));
    }
    for i in 0..200 {
        let n = 50 + (r.next() % 500) as usize;
        db += &format!(">r{i}\n{}\n", r.protein(n));
    }
    (t.file("q.fa", &q), t.file("db.fa", &db))
}

#[test]
fn empty_inputs_are_handled() {
    let t = Tmp::new("empty");
    let (q, db) = protein_case(&t, 3);
    let empty = t.file("empty.fa", "");
    let out = run_ok(&["-q", s(&empty), "-d", s(&db), "--quiet"]);
    assert!(out.is_empty());
    let out = run_ok(&["-q", s(&q), "-d", s(&empty), "--quiet"]);
    assert!(out.is_empty());
}

#[test]
fn unusual_residues_do_not_crash_or_misencode() {
    let t = Tmp::new("weird");
    let q = t.file("q.fa", ">q lower case, ambiguity codes, stop, digits\nmktayiakqrqisfvkshfsrqXXXXBZJUO*12 3mktayiakqrq\n>x only X\nXXXXXXXXXXXXXXXXXXXX\n>n1\n\n");
    let db = t.file("db.fa", ">s1\nMKTAYIAKQRQISFVKSHFSRQBZJUOXMKTAYIAKQRQ\n>s2\nXXXXXXXXXXXXXXXXXXXXXXXX\n");
    let out = run_ok(&["-q", s(&q), "-d", s(&db), "--quiet", "--outfmt", "json"]);
    assert!(out.contains("\"subject_id\": \"s1\"") || out.contains("\"subject_id\":\"s1\""), "{out}");
    assert!(!out.contains("s2"), "an all-X subject must not be reported:\n{out}");
}

#[test]
fn low_complexity_and_repeats_stay_bounded() {
    // Adversarial input: long homopolymers and short tandem repeats. With
    // masking on, seeds are suppressed; with masking off, a tiny hit buffer
    // forces subject-range passes and streaming instead of materialising all
    // hits at once. Both must finish and report their limits.
    let t = Tmp::new("lowc");
    let mut db = String::new();
    for i in 0..40 {
        db += &format!(">polyq{i}\n{}\n", "Q".repeat(400 + i * 10));
        db += &format!(">rep{i}\n{}\n", "PEST".repeat(150));
    }
    let q = t.file("q.fa", &format!(">polyq\n{}\n>rep\n{}\n", "Q".repeat(500), "PESTPEST".repeat(60)));
    let dbf = t.file("db.fa", &db);
    let st = t.path("masked.json");
    run_ok(&["-q", s(&q), "-d", s(&dbf), "--quiet", "--stats", s(&st), "-o", s(&t.path("m.m8"))]);
    let m = json_stats(&st);
    assert!(m["counters"]["query_residues_masked"].as_u64().unwrap() > 900, "tantan must mask the repeats: {m}");
    // Without masking, the default frequency cap skips the over-represented
    // words and says so.
    let st1 = t.path("cap.json");
    run_ok(&["-q", s(&q), "-d", s(&dbf), "--quiet", "--masking", "none", "--stats", s(&st1), "-o", s(&t.path("c.m8"))]);
    let c1 = json_stats(&st1);
    assert!(c1["counters"]["seeds_frequent_skipped"].as_u64().unwrap() > 0, "{c1}");
    assert!(c1["heuristic_limits"].as_array().unwrap().iter().any(|l| l.as_str().unwrap().contains("over-frequent")), "{c1}");
    // Index schedule: a tiny hit buffer forces bounded subject-range passes
    // and streaming of oversized subjects.
    let st2 = t.path("unmasked.json");
    run_ok(&[
        "-q", s(&q), "-d", s(&dbf), "--quiet", "--masking", "none", "--preset", "fast", "--batch-order", "query", "--hit-buffer", "1024", "--seed-max-occ", "100000000", "--stats",
        s(&st2), "-o", s(&t.path("u.m8")),
    ]);
    let u = json_stats(&st2);
    let c = &u["counters"];
    assert!(c["hits"].as_u64().unwrap() > 100_000, "unmasked repeats generate many hits: {c}");
    assert!(c["subject_range_passes"].as_u64().unwrap() > 2, "hits are processed in bounded passes: {c}");
    assert!(c["streamed_subjects"].as_u64().unwrap() > 0, "oversized subjects are streamed: {c}");
    // Scan schedule (the default): work is cut by chunks, subjects and
    // queries so that the window-filter survivors held stay bounded; the
    // output is the index schedule's.
    let st3 = t.path("scan.json");
    run_ok(&["-q", s(&q), "-d", s(&dbf), "--quiet", "--masking", "none", "--preset", "fast", "--seed-max-occ", "100000000", "--stats", s(&st3), "-o", s(&t.path("s.m8"))]);
    let c3 = json_stats(&st3)["counters"].clone();
    assert!(c3["scan_splits"].as_u64().unwrap() > 0, "survivor-bounded cuts: {c3}");
    assert_eq!(c3["hits"], c["hits"], "both schedules see every seed hit");
    assert_eq!(std::fs::read(t.path("s.m8")).unwrap(), std::fs::read(t.path("u.m8")).unwrap(), "scan and index schedules must agree");
    // A tiny pair limit: shape scans stop at subject boundaries and every
    // (query, subject) pair above the limit is scanned and visited on its
    // own in query-position slices; the output does not change.
    let st4 = t.path("pairs.json");
    run_ok(&[
        "-q", s(&q), "-d", s(&dbf), "--quiet", "--masking", "none", "--preset", "fast", "--seed-max-occ", "100000000", "--scan-pair-survivors", "500", "--stats", s(&st4), "-o",
        s(&t.path("p.m8")),
    ]);
    let c4 = json_stats(&st4)["counters"].clone();
    assert!(c4["scan_pairs_streamed"].as_u64().unwrap() > 0, "oversized pairs streamed: {c4}");
    assert!(c4["scan_splits"].as_u64().unwrap() > c3["scan_splits"].as_u64().unwrap(), "more cuts with the small limit: {c4}");
    assert_eq!(std::fs::read(t.path("p.m8")).unwrap(), std::fs::read(t.path("u.m8")).unwrap(), "streamed pairs must give the same output");
}

#[test]
fn results_identical_across_threads_batches_and_schedules() {
    let t = Tmp::new("det");
    let (q, db) = protein_case(&t, 12);
    let idx = t.path("db.sdx");
    // postings stored, so that the index schedules can search it too
    run_ok(&["--makedb", s(&idx), "-d", s(&db), "--partition-residues", "20000", "--index-postings", "--quiet"]);
    let base = run_ok(&["-q", s(&q), "-d", s(&db), "-t", "1", "--quiet", "-k", "50"]);
    assert!(base.lines().count() >= 12);
    for extra in [
        vec!["-t", "2"],
        vec!["-t", "2", "--query-batch", "1"],
        vec!["-t", "2", "--scan-batch-residues", "700"],
        vec!["-t", "2", "--db", s(&idx)],
        vec!["-t", "2", "--db", s(&idx), "--batch-order", "partition"],
        vec!["-t", "2", "--db", s(&idx), "--batch-order", "query"],
    ] {
        let mut a = vec!["-q", s(&q), "--quiet", "-k", "50"];
        if !extra.contains(&"--db") {
            a.extend(["-d", s(&db)]);
        }
        a.extend(extra.iter().copied());
        assert_eq!(run_ok(&a), base, "output differs for {extra:?}");
    }
}

#[test]
fn scan_schedule_matches_index_schedules_with_multi_shape_seeds() {
    // Several spaced seeds, the exact-identity-weighted window filter, the
    // floor-relative region filter and banded scoring: the scan schedule
    // must report exactly what the index schedules report, for both E-value
    // models and thread counts.
    let t = Tmp::new("scaneq");
    let (q, db) = protein_case(&t, 16);
    let idx = t.path("db.sdx");
    let seeds = ["--seed-alphabet", "murphy10", "--seed-shape", "1101011,111011,1100111,1011011"];
    let mut mk = vec!["--makedb", s(&idx), "-d", s(&db), "--partition-residues", "30000", "--index-postings", "--quiet"];
    mk.extend(seeds);
    run_ok(&mk);
    for model in ["per-subject", "database"] {
        let mut base_args = vec!["-q", s(&q), "--db", s(&idx), "--quiet", "-k", "50", "-E", "10"];
        base_args.extend(seeds);
        base_args.extend(["--hit-filter-ids", "14", "--hit-filter-exact-weight", "2", "--hsp-floor-fraction", "0.5", "--score-band", "24", "--evalue-model", model]);
        let mut a = base_args.clone();
        a.extend(["-t", "1", "--batch-order", "query"]);
        let base = run_ok(&a);
        assert!(base.lines().count() >= 16, "planted homologs must be found ({model})");
        for extra in [vec!["-t", "1", "--batch-order", "scan"], vec!["-t", "2", "--batch-order", "scan"], vec!["-t", "2", "--batch-order", "partition"]] {
            let mut a = base_args.clone();
            a.extend(extra.iter().copied());
            assert_eq!(run_ok(&a), base, "output differs for {extra:?} ({model})");
        }
    }
}

#[test]
fn stringent_evalue_only_filters_an_evalue_10_search() {
    // The region threshold and the two-stage scoring threshold are relative
    // to the floor at E = 10 (or a more permissive -E), so a stringent -E or
    // --min-bitscore reports exactly the -E 10 alignments that pass it.
    let t = Tmp::new("ecut");
    // divergent homologs (about 50 % substitutions, short indels): strong
    // gapped alignments made of weak ungapped segments
    let mut r = Rng(0xec07);
    let (mut qs, mut dbs) = (String::new(), String::new());
    for i in 0..24 {
        let n = 250 + (r.next() % 150) as usize;
        let p = r.protein(n);
        qs += &format!(">q{i}\n{p}\n");
        for k in 0..3 {
            let mut h = String::new();
            for c in p.chars() {
                match r.next() % 100 {
                    0..=2 => {}
                    3..=5 => {
                        h.push(c);
                        h.push_str(&r.protein(1));
                    }
                    x if x < 50 + 4 * k as u64 => h.push_str(&r.protein(1)),
                    _ => h.push(c),
                }
            }
            let (pre, suf) = (r.protein(30), r.protein(30));
            dbs += &format!(">h{i}_{k}\n{pre}{h}{suf}\n");
        }
    }
    for i in 0..300 {
        let n = 80 + (r.next() % 400) as usize;
        dbs += &format!(">r{i}\n{}\n", r.protein(n));
    }
    let (q, db) = (t.file("q.fa", &qs), t.file("db.fa", &dbs));
    let run = |extra: &[&str]| -> Vec<String> {
        let mut a = vec!["-q", s(&q), "-d", s(&db), "--quiet", "-k", "1000", "--outfmt", "json"];
        a.extend_from_slice(extra);
        let v: serde_json::Value = serde_json::from_str(&run_ok(&a)).unwrap();
        v.as_array().unwrap().iter().map(|h| h.to_string()).collect()
    };
    let parse = |h: &String| -> serde_json::Value { serde_json::from_str(h).unwrap() };
    let all = run(&["-E", "10"]);
    assert!(all.len() >= 24);
    for (extra, keep) in [
        (vec!["-E", "1e-10"], Box::new(|h: &serde_json::Value| h["evalue"].as_f64().unwrap() <= 1e-10) as Box<dyn Fn(&serde_json::Value) -> bool>),
        (vec!["-E", "10", "--min-bitscore", "60"], Box::new(|h: &serde_json::Value| h["bitscore"].as_f64().unwrap() >= 60.0)),
    ] {
        let want: Vec<String> = all.iter().filter(|h| keep(&parse(h))).cloned().collect();
        assert_eq!(run(&extra), want, "{extra:?}");
    }
}

#[test]
fn evalue_models_share_alignments_and_bit_scores() {
    let t = Tmp::new("evm");
    let (q, db) = protein_case(&t, 10);
    let run = |model: &str| -> Vec<Vec<String>> {
        let out = run_ok(&["-q", s(&q), "-d", s(&db), "--quiet", "-k", "50", "-E", "1e-3", "--evalue-model", model]);
        let mut v: Vec<Vec<String>> = out.lines().map(|l| l.split('\t').map(String::from).collect()).collect();
        v.sort();
        v
    };
    let (ps, db_) = (run("per-subject"), run("database"));
    assert!(!ps.is_empty() && !db_.is_empty());
    // the planted homologs are far below the cut-off under both models:
    // same pairs, coordinates and bit scores; only the E-values differ
    let key = |r: &Vec<String>| (r[0].clone(), r[1].clone(), r[6].clone(), r[7].clone(), r[8].clone(), r[9].clone(), r[11].clone());
    let a: Vec<_> = ps.iter().map(key).collect();
    let b: Vec<_> = db_.iter().map(key).collect();
    assert_eq!(a, b);
    let err = run_err(&["-q", s(&q), "-d", s(&db), "--evalue-model", "bogus"]);
    assert!(err.contains("evalue-model"), "{err}");
}

#[test]
fn long_alignment_is_not_truncated() {
    // 6000-residue protein against a copy with a 300-residue deletion and a
    // 40-residue insertion: v1.1 capped traceback windows at 2048 residues.
    let t = Tmp::new("long");
    let mut r = Rng(42);
    let a = r.protein(6000);
    let b = format!("{}{}{}{}", &a[..2500], r.protein(40), &a[2800..5000], &a[5000..]);
    let q = t.file("q.fa", &format!(">a\n{a}\n"));
    let db = t.file("db.fa", &format!(">b\n{b}\n>c\n{}\n", r.protein(300)));
    let st = t.path("long.json");
    let out = run_ok(&["-q", s(&q), "-d", s(&db), "--quiet", "--outfmt", "json", "--stats", s(&st)]);
    let c = json_stats(&st)["counters"].clone();
    assert!(c["band_scored_regions"].as_u64().unwrap() > 0, "long pair must take the banded scoring path: {c}");
    let v: serde_json::Value = serde_json::from_str(&out).unwrap();
    let hit = v.as_array().unwrap().iter().find(|h| h["subject_id"] == "b").expect("hit to b").clone();
    let (qs, qe) = (hit["query_range"][0].as_u64().unwrap(), hit["query_range"][1].as_u64().unwrap());
    let (ss, se) = (hit["subject_range"][0].as_u64().unwrap(), hit["subject_range"][1].as_u64().unwrap());
    assert!(qe - qs + 1 > 5900 && se - ss + 1 > 5600, "alignment must span the whole protein: q {qs}-{qe} s {ss}-{se}");
    assert!(hit["score"].as_i64().unwrap() > 25_000, "score above the 16-bit range");
    // CIGAR consumes exactly the reported ranges
    let cigar = hit["cigar"].as_str().unwrap();
    let (mut ql, mut sl, mut n) = (0u64, 0u64, 0u64);
    for ch in cigar.chars() {
        if let Some(d) = ch.to_digit(10) {
            n = n * 10 + d as u64;
        } else {
            match ch {
                'M' | '=' | 'X' => {
                    ql += n;
                    sl += n
                }
                'I' => ql += n,
                'D' => sl += n,
                _ => panic!("unexpected op {ch}"),
            }
            n = 0;
        }
    }
    assert_eq!((ql, sl), (qe - qs + 1, se - ss + 1), "CIGAR {cigar}");
    assert_eq!(hit["exact"].as_bool(), Some(true));
}

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            o => o,
        })
        .collect()
}

#[test]
fn translated_searches_report_frames_and_nucleotide_coordinates() {
    let t = Tmp::new("tr");
    let mut r = Rng(7);
    let prot = r.protein(120);
    // back-translate with one fixed codon per residue
    let codon = |c: char| match c {
        'A' => "GCT", 'R' => "CGT", 'N' => "AAC", 'D' => "GAT", 'C' => "TGC", 'Q' => "CAG", 'E' => "GAA", 'G' => "GGT", 'H' => "CAC", 'I' => "ATC",
        'L' => "CTG", 'K' => "AAA", 'M' => "ATG", 'F' => "TTC", 'P' => "CCG", 'S' => "TCT", 'T' => "ACC", 'W' => "TGG", 'Y' => "TAC", _ => "GTT",
    };
    let cds: String = prot.chars().map(codon).collect();
    let left = r.dna(100);
    let right = r.dna(77);
    let plus = format!("{left}{cds}{right}");
    let minus = revcomp(&plus);
    let nq = t.file("nq.fa", &format!(">plus\n{plus}\n>minus\n{minus}\n"));
    let pdb = t.file("pdb.fa", &format!(">P\n{prot}\n>decoy\n{}\n", r.protein(200)));
    let out = run_ok(&["--mode", "nx", "-q", s(&nq), "-d", s(&pdb), "--quiet", "-E", "1e-5"]);
    let rows: Vec<Vec<&str>> = out.lines().map(|l| l.split('\t').collect()).collect();
    let p = rows.iter().find(|r| r[0].starts_with("plus|frame+")).expect("plus-strand hit with |frame+N suffix");
    assert_eq!((p[6], p[7]), ("101", &*(100 + cds.len()).to_string()), "plus-strand nt coordinates");
    let m = rows.iter().find(|r| r[0].starts_with("minus|frame-")).expect("minus-strand hit with |frame-N suffix");
    let total = plus.len();
    assert_eq!(m[6], (total - 100).to_string(), "minus-strand start (BLAST style: start > end)");
    assert_eq!(m[7], (total - 100 - cds.len() + 1).to_string());
    // pn: protein query against the translated nucleotide database
    let pq = t.file("pq.fa", &format!(">P\n{prot}\n"));
    let ndb = t.file("ndb.fa", &format!(">minus\n{minus}\n>noise\n{}\n", r.dna(600)));
    let out = run_ok(&["--mode", "pn", "-q", s(&pq), "-d", s(&ndb), "--quiet", "-E", "1e-5"]);
    let row: Vec<&str> = out.lines().next().expect("pn hit").split('\t').collect();
    assert!(row[1].starts_with("minus|frame-"), "{row:?}");
    assert_eq!((row[8], row[9]), (&*(total - 100).to_string(), &*(total - 100 - cds.len() + 1).to_string()));
    // persistent indexes for both translated modes
    let pidx = t.path("p.sdx");
    run_ok(&["--mode", "nx", "--makedb", s(&pidx), "-d", s(&pdb), "--quiet"]);
    let out2 = run_ok(&["--mode", "nx", "-q", s(&nq), "--db", s(&pidx), "--quiet", "-E", "1e-5"]);
    assert_eq!(out2, run_ok(&["--mode", "nx", "-q", s(&nq), "-d", s(&pdb), "--quiet", "-E", "1e-5"]));
    let tidx = t.path("t.sdx");
    run_ok(&["--mode", "pn", "--makedb", s(&tidx), "-d", s(&ndb), "--quiet"]);
    let out3 = run_ok(&["--mode", "pn", "-q", s(&pq), "--db", s(&tidx), "--quiet", "-E", "1e-5"]);
    assert_eq!(out3, out, "pn search of the translated index = in-memory pn search");
    let err = run_err(&["--mode", "pp", "-q", s(&pq), "--db", s(&tidx), "--quiet"]);
    assert!(err.contains("translated"), "{err}");
}

#[test]
fn corrupt_truncated_and_foreign_indexes_are_rejected() {
    let t = Tmp::new("idx");
    let (q, db) = protein_case(&t, 2);
    let idx = t.path("db.sdx");
    run_ok(&["--makedb", s(&idx), "-d", s(&db), "--quiet"]);
    let bytes = std::fs::read(&idx).unwrap();
    // truncated
    let tr = t.path("trunc.sdx");
    std::fs::write(&tr, &bytes[..bytes.len() / 2]).unwrap();
    let e = run_err(&["-q", s(&q), "--db", s(&tr), "--quiet"]);
    assert!(e.to_lowercase().contains("trunc") || e.contains("size"), "{e}");
    // flipped residue byte: detected by --verify-index checksums
    let mut fl = bytes.clone();
    let mid = fl.len() / 2;
    fl[mid] ^= 0x5a;
    let fp = t.path("flip.sdx");
    std::fs::write(&fp, &fl).unwrap();
    let e = run_err(&["-q", s(&q), "--db", s(&fp), "--verify-index", "--quiet"]);
    assert!(e.to_lowercase().contains("checksum") || e.contains("corrupt"), "{e}");
    // not an index at all
    let junk = t.file("junk.sdx", &"x".repeat(4096));
    let e = run_err(&["-q", s(&q), "--db", s(&junk), "--quiet"]);
    assert!(e.contains("not a SABER") || e.to_lowercase().contains("magic") || e.contains("index"), "{e}");
    // no temporary or partial file left behind by --makedb
    let leftovers: Vec<_> = std::fs::read_dir(&t.0).unwrap().filter_map(|e| e.ok()).map(|e| e.file_name().to_string_lossy().to_string()).filter(|n| n.contains("tmp")).collect();
    assert!(leftovers.is_empty(), "{leftovers:?}");
}

#[test]
fn uncalibrated_scoring_is_rejected_unless_allowed() {
    let t = Tmp::new("cal");
    let (q, db) = protein_case(&t, 2);
    let e = run_err(&["-q", s(&q), "-d", s(&db), "-g", "3", "-e", "3", "--quiet"]);
    assert!(e.contains("allow-uncalibrated"), "{e}");
    let out = run_ok(&["-q", s(&q), "-d", s(&db), "-g", "3", "-e", "3", "--allow-uncalibrated-evalues", "--quiet", "--outfmt", "json"]);
    assert!(out.contains("\"evalue_calibrated\": false") || out.contains("\"evalue_calibrated\":false"), "{out}");
}

#[test]
fn negative_mismatch_and_bam_alias_are_handled() {
    // 1.0/1.1 rejected `--mismatch -3` (clap read -3 as a flag) although the
    // README documented that form; both spellings must now give the same run.
    let t = Tmp::new("nnargs");
    let mut r = Rng(0x5eed_1234);
    let subj = r.dna(600);
    let q = t.file("q.fna", &format!(">q\n{}\n", &subj[100..400]));
    let db = t.file("db.fna", &format!(">s\n{subj}\n>other\n{}\n", r.dna(600)));
    let a = run_ok(&["--mode", "nn", "-q", s(&q), "-d", s(&db), "--match", "1", "--mismatch", "-2", "-g", "2", "-e", "2", "--gap-convention", "blast", "--quiet"]);
    let b = run_ok(&["--mode", "nn", "-q", s(&q), "-d", s(&db), "--match", "1", "--mismatch=-2", "-g", "2", "-e", "2", "--gap-convention", "blast", "--quiet"]);
    assert_eq!(a, b);
    assert!(a.lines().next().unwrap_or("").starts_with("q\ts\t100.00\t300\t"), "{a}");
    // `--outfmt bam` has always written SAM text; it now says so on stderr.
    let out = saber().args(["--mode", "nn", "-q", s(&q), "-d", s(&db), "--outfmt", "bam"]).output().unwrap();
    assert!(out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("writes SAM text"));
    assert!(String::from_utf8_lossy(&out.stdout).starts_with("@HD"));
}

#[test]
fn every_output_format_is_well_formed() {
    let t = Tmp::new("fmt");
    let (q, db) = protein_case(&t, 3);
    // protein_case's homologs carry an A at every third position, a period-3
    // pattern that tantan masks for seeding; search them unmasked here.
    let run = |fmt: &str| run_ok(&["-q", s(&q), "-d", s(&db), "--quiet", "-E", "1e-5", "--masking", "none", "--outfmt", fmt]);
    let m8 = run("m8");
    let rows: Vec<&str> = m8.lines().collect();
    assert!(rows.len() >= 3, "{m8}");
    assert!(rows.iter().all(|r| r.split('\t').count() == 12), "{m8}");
    // m9 = comment lines + the same rows
    let m9 = run("m9");
    assert!(m9.lines().next().unwrap().starts_with("# SABER"), "{m9}");
    assert!(m9.contains("# Fields: query id, subject id"), "{m9}");
    assert_eq!(m9.lines().filter(|l| !l.starts_with('#')).collect::<Vec<_>>(), rows);
    // csv = the same 12 fields, comma-separated
    let csv = run("csv");
    assert!(csv.starts_with("query_id,subject_id,"), "{csv}");
    let crow: Vec<&str> = csv.lines().skip(1).collect();
    assert_eq!(crow.len(), rows.len(), "{csv}");
    assert!(crow.iter().all(|r| r.split(',').count() == 12), "{csv}");
    assert_eq!(crow[0].split(',').take(2).collect::<Vec<_>>(), rows[0].split('\t').take(2).collect::<Vec<_>>());
    // m0 = one pairwise block per alignment
    let m0 = run("m0");
    assert_eq!(m0.matches(" Score = ").count(), rows.len(), "{m0}");
    assert!(m0.contains("Query ") && m0.contains("Sbjct "), "{m0}");
    // SAM = header + one record per alignment, CIGAR consistent with the query range
    let sam = run("sam");
    assert!(sam.starts_with("@HD\t"), "{sam}");
    let recs: Vec<Vec<&str>> = sam.lines().filter(|l| !l.starts_with('@')).map(|l| l.split('\t').collect()).collect();
    assert_eq!(recs.len(), rows.len());
    for (r, m) in recs.iter().zip(&rows) {
        let m: Vec<&str> = m.split('\t').collect();
        assert!(r.len() >= 11 && r[11].starts_with("AS:i:"), "{r:?}");
        let (mut qlen, mut n) = (0usize, 0usize);
        for c in r[5].chars() {
            if let Some(d) = c.to_digit(10) {
                n = n * 10 + d as usize;
            } else {
                if c == 'M' || c == 'I' {
                    qlen += n;
                }
                n = 0;
            }
        }
        let (qs, qe): (usize, usize) = (m[6].parse().unwrap(), m[7].parse().unwrap());
        assert_eq!(qlen, qe - qs + 1, "CIGAR {} vs query range {qs}-{qe}", r[5]);
        assert_eq!(r[9].len(), qlen, "SEQ holds the aligned query residues");
    }
    // JSON = an array of objects with the same pairs
    let js: serde_json::Value = serde_json::from_str(&run("json")).unwrap();
    let arr = js.as_array().unwrap();
    assert_eq!(arr.len(), rows.len());
    assert_eq!(arr[0]["subject_id"].as_str().unwrap(), rows[0].split('\t').nth(1).unwrap());
}

#[test]
fn capped_traceback_of_a_repetitive_long_pair_stays_exact() {
    // A 3,000-residue protein made of diverged copies of one domain, searched
    // against itself: its candidate region spans thousands of diagonals (the
    // repeats), so a band over the whole region exceeds a small traceback
    // memory cap. The self-alignment must still be found exactly.
    let t = Tmp::new("cap");
    let mut r = Rng(0xabc);
    let unit = r.protein(120);
    let aa = b"ARNDCQEGHILKMFPSTWYV";
    let mut long = String::new();
    for _ in 0..25 {
        for c in unit.chars() {
            long.push(if r.next() % 4 == 0 { aa[(r.next() % 20) as usize] as char } else { c });
        }
    }
    let q = t.file("q.fa", &format!(">rep\n{long}\n"));
    let db = t.file("db.fa", &format!(">rep\n{long}\n>other\n{}\n", r.protein(500)));
    let out = run_ok(&["-q", s(&q), "-d", s(&db), "--quiet", "--traceback-max-mb", "1", "-E", "1e-3", "--outfmt", "json"]);
    let js: serde_json::Value = serde_json::from_str(&out).unwrap();
    let hit = js.as_array().unwrap().iter().find(|h| h["subject_id"] == "rep").expect("self hit");
    assert_eq!(hit["query_range"], serde_json::json!([1, 3000]), "{hit}");
    assert_eq!(hit["subject_range"], serde_json::json!([1, 3000]), "{hit}");
    assert_eq!(hit["exact"], serde_json::json!(true), "{hit}");
}
