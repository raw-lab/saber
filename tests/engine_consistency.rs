//! Engine-level regression tests:
//! * intra-query parallel search returns exactly the sequential results
//!   (several partitions, sub-partition ranges, both seed alphabets);
//! * a tiny hit buffer (forcing subject-range passes and streamed subjects)
//!   returns exactly the unbounded results;
//! * exhaustive mode reproduces the naive full-matrix Smith-Waterman score
//!   for every subject, and every reported alignment re-scores to its score.

use saber_rs::algorithm::reference::{score_alignment, sw_reference};
use saber_rs::alphabet::{EncodedProtein, SeedAlphabet, SeedShape};
use saber_rs::masking::MaskingMode;
use saber_rs::scoring::stats::{gapped_protein, GapConvention, GapCosts, SearchStatistics};
use saber_rs::scoring::{ScoreMatrix, ScoringSystem};
use saber_rs::search::engine::{Engine, Hit, SearchParams, Worker};
use saber_rs::search::index::{BuildOptions, IndexBuilder, SeedConfig, SeedIndex};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

const AA: &[u8] = b"ARNDCQEGHILKMFPSTWYV";

fn random_protein(r: &mut Rng, len: usize) -> Vec<u8> {
    (0..len).map(|_| AA[r.below(20) as usize]).collect()
}

fn mutate(r: &mut Rng, seq: &[u8], sub_pct: u64, indel_pct: u64) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len() + 8);
    for &c in seq {
        let x = r.below(100);
        if x < indel_pct / 2 {
            continue; // deletion
        }
        if x < indel_pct {
            out.push(AA[r.below(20) as usize]); // insertion
        }
        if r.below(100) < sub_pct {
            out.push(AA[r.below(20) as usize]);
        } else {
            out.push(c);
        }
    }
    out
}

/// Database with planted homologs of each query (various divergences) among
/// random sequences.
fn dataset(seed: u64, nq: usize, nrand: usize) -> (Vec<Vec<u8>>, Vec<(String, Vec<u8>)>) {
    let mut r = Rng(seed);
    let queries: Vec<Vec<u8>> = (0..nq).map(|_| {
        let l = 60 + r.below(400) as usize;
        random_protein(&mut r, l)
    }).collect();
    let mut db = Vec::new();
    for (qi, q) in queries.iter().enumerate() {
        for (k, &(sub, indel)) in [(10u64, 2u64), (35, 6), (55, 8)].iter().enumerate() {
            let lf = r.below(150) as usize;
            let rf = r.below(150) as usize;
            let mut s = random_protein(&mut r, lf);
            let lo = r.below(q.len() as u64 / 3) as usize;
            let hi = q.len() - r.below(q.len() as u64 / 3) as usize;
            s.extend(mutate(&mut r, &q[lo..hi], sub, indel));
            s.extend(random_protein(&mut r, rf));
            db.push((format!("hom_{qi}_{k}"), s));
        }
    }
    for i in 0..nrand {
        let l = 30 + r.below(700) as usize;
        db.push((format!("rand_{i}"), random_protein(&mut r, l)));
    }
    // interleave deterministically so homologs land in different partitions
    let n = db.len();
    let mut order: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        let j = r.below(i as u64 + 1) as usize;
        order.swap(i, j);
    }
    let db = order.into_iter().map(|i| db[i].clone()).collect();
    (queries, db)
}

fn build(db: &[(String, Vec<u8>)], alphabet: &str, shapes: &[&str], masking: MaskingMode, part_res: usize) -> SeedIndex {
    let mut cfg = SeedConfig::protein_default();
    cfg.alphabet = SeedAlphabet::parse(alphabet).unwrap();
    cfg.shapes = shapes.iter().map(|s| SeedShape::parse(s).unwrap()).collect();
    cfg.masking = masking;
    let mut opt = BuildOptions::new(cfg);
    opt.max_partition_residues = part_res;
    let mut b = IndexBuilder::in_memory(opt).unwrap();
    let mut scratch = Vec::new();
    for (id, s) in db {
        b.add_ascii(id, s, &mut scratch).unwrap();
    }
    b.finish_in_memory().unwrap()
}

fn engine<'a>(idx: &'a SeedIndex, p: SearchParams) -> Engine<'a> {
    let gaps = GapCosts::from_user(11, 1, GapConvention::Saber).unwrap();
    let kp = gapped_protein(ScoringSystem::Blosum62, gaps, false).unwrap();
    let st = SearchStatistics::new(kp, idx.total_residues(), idx.num_subjects() as u64);
    Engine::new(idx, ScoreMatrix::new(ScoringSystem::Blosum62), gaps, st, p)
}

type Key = (u64, i32, usize, usize, usize, usize, String, bool);

fn keys(h: &[Hit]) -> Vec<Key> {
    h.iter().map(|h| (h.subject, h.score, h.q_start, h.q_end, h.s_start, h.s_end, h.cigar(), h.exact_in_window)).collect()
}

fn check_parallel(alphabet: &str, shapes: &[&str], hf: u32, part_res: usize) {
    check_parallel_with(alphabet, shapes, hf, part_res, 0.0, 10.0)
}

fn check_parallel_with(alphabet: &str, shapes: &[&str], hf: u32, part_res: usize, floor_frac: f64, evalue: f64) {
    let (queries, db) = dataset(0x5eed_1234 ^ part_res as u64, 12, 500);
    let idx = build(&db, alphabet, shapes, MaskingMode::Tantan, part_res);
    let mut p = SearchParams::default();
    p.max_target_seqs = 50;
    p.max_hsps = 2;
    p.hit_filter_min_ids = hf;
    p.max_candidates = 40; // small cap: exercises the deterministic merge
    p.hsp_floor_fraction = floor_frac;
    p.evalue = evalue;
    let eng = engine(&idx, p);
    let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
    let mut w = Worker::new();
    let mut total = 0usize;
    let encoded: Vec<EncodedProtein> = queries.iter().map(|q| EncodedProtein::from_ascii(q)).collect();
    let views: Vec<&[u8]> = encoded.iter().map(|e| e.as_slice()).collect();
    let batch = pool.install(|| eng.search_batch(&views));
    let mut seq_regions = 0u64;
    let mut batch_regions = 0u64;
    for (qi, codes) in encoded.iter().enumerate() {
        let before = w.counters;
        let seq = eng.search(codes.as_slice(), &mut w);
        seq_regions += w.counters.regions - before.regions;
        let (par, _) = pool.install(|| eng.search_parallel(codes.as_slice()));
        assert_eq!(keys(&seq), keys(&par), "parallel result differs ({alphabet} {shapes:?}, partitions {})", idx.num_partitions());
        assert_eq!(keys(&seq), keys(&batch[qi].0), "partition-major batch result differs ({alphabet} {shapes:?})");
        batch_regions += batch[qi].1.regions;
        total += seq.len();
    }
    assert_eq!(seq_regions, batch_regions, "same candidate regions in both schedules");
    assert!(total >= queries.len(), "planted homologs should be found (got {total})");
}

#[test]
fn parallel_equals_sequential_many_partitions_full_alphabet() {
    check_parallel("full", &["1111"], 0, 8_000);
}

#[test]
fn parallel_equals_sequential_one_partition_murphy10() {
    check_parallel("murphy10", &["110111"], 12, 64 << 20);
}

#[test]
fn parallel_equals_sequential_with_floor_relative_region_filter() {
    check_parallel_with("murphy10", &["110111"], 12, 30_000, 0.6, 1e-5);
}

#[test]
fn parallel_equals_sequential_few_partitions_two_shapes() {
    check_parallel("murphy10", &["1101011", "111011"], 11, 40_000);
}

#[test]
fn tiny_hit_buffer_and_streaming_equal_unbounded() {
    // Repetitive subjects produce far more hits than a 1024-hit buffer.
    let mut r = Rng(77);
    let motif = random_protein(&mut r, 23);
    let mut db: Vec<(String, Vec<u8>)> = Vec::new();
    for i in 0..30 {
        let mut s = Vec::new();
        let copies = 5 + r.below(120) as usize;
        for _ in 0..copies {
            s.extend(mutate(&mut r, &motif, 5, 0));
        }
        db.push((format!("rep_{i}"), s));
    }
    let (queries, rnd) = dataset(91, 4, 80);
    db.extend(rnd);
    let mut queries = queries;
    let mut q = Vec::new();
    for _ in 0..6 {
        q.extend_from_slice(&motif);
    }
    queries.push(q);
    let idx = build(&db, "full", &["1111"], MaskingMode::None, 64 << 20);
    let mk = |buf: usize| {
        let mut p = SearchParams::default();
        p.max_target_seqs = 100;
        p.max_hsps = 3;
        p.query_masking = MaskingMode::None;
        p.comp_bias = false;
        p.hit_buffer = buf;
        // the repeats must produce hits (the automatic frequency cap would
        // skip their seeds)
        p.seed_max_occ = 1 << 30;
        p
    };
    let big = engine(&idx, mk(1 << 22));
    let small = engine(&idx, mk(1024));
    let mut wb = Worker::new();
    let mut ws = Worker::new();
    for q in &queries {
        let c = EncodedProtein::from_ascii(q);
        let a = big.search(c.as_slice(), &mut wb);
        let b = small.search(c.as_slice(), &mut ws);
        assert_eq!(keys(&a), keys(&b));
    }
    assert!(ws.counters.streamed_subjects > 0, "test must exercise the streaming path");
    assert!(ws.counters.subject_range_passes > wb.counters.subject_range_passes);
}

#[test]
fn exhaustive_mode_matches_reference_smith_waterman() {
    let (queries, db) = dataset(4242, 5, 40);
    let idx = build(&db, "full", &["1111"], MaskingMode::None, 64 << 20);
    let mut p = SearchParams::default();
    p.exhaustive = true;
    p.evalue = 1e12;
    p.max_target_seqs = db.len();
    p.comp_bias = false;
    p.query_masking = MaskingMode::None;
    let eng = engine(&idx, p);
    let m = ScoreMatrix::new(ScoringSystem::Blosum62);
    let mut w = Worker::new();
    for q in &queries {
        let c = EncodedProtein::from_ascii(q);
        let hits = eng.search(c.as_slice(), &mut w);
        for h in &hits {
            let s = idx.subject_codes(h.subject as usize);
            let (best, _, _) = sw_reference(c.as_slice(), s, &m.matrix, None, 11, 1);
            assert_eq!(h.score as i64, best, "exhaustive score must equal full Smith-Waterman");
            // the reported alignment re-scores to the reported score
            let ops = parse_cigar(&h.cigar());
            let (rs, qe, se) = score_alignment(c.as_slice(), s, h.q_start - 1, h.s_start - 1, &ops, &m.matrix, None, 11, 1)
                .expect("CIGAR must fit the sequences");
            assert_eq!(rs, h.score as i64, "alignment re-scores to its score");
            assert_eq!((qe, se), (h.q_end, h.s_end), "CIGAR consumes exactly the reported ranges");
        }
        // every subject with a positive reference score is reported
        let positive = (0..idx.num_subjects())
            .filter(|&i| sw_reference(c.as_slice(), idx.subject_codes(i), &m.matrix, None, 11, 1).0 > 0)
            .count();
        assert_eq!(hits.len(), positive);
    }
}

fn parse_cigar(c: &str) -> Vec<(u8, u32)> {
    let mut out = Vec::new();
    let mut n = 0u32;
    for ch in c.bytes() {
        if ch.is_ascii_digit() {
            n = n * 10 + (ch - b'0') as u32;
        } else {
            out.push((ch, n));
            n = 0;
        }
    }
    out
}

#[test]
fn scalar_and_simd_hit_paths_agree() {
    let (queries, db) = dataset(0xabcdef, 10, 300);
    for (alpha, shape, hf) in [("murphy10", "110111", 12u32), ("full", "1111", 0u32)] {
        let idx = build(&db, alpha, &[shape], MaskingMode::Tantan, 64 << 20);
        let mut p = SearchParams::default();
        p.max_target_seqs = 50;
        p.hit_filter_min_ids = hf;
        let fast = engine(&idx, p.clone());
        let mut slow = engine(&idx, p);
        slow.force_scalar_hit_path();
        let (mut w1, mut w2) = (Worker::new(), Worker::new());
        for q in &queries {
            let c = EncodedProtein::from_ascii(q);
            assert_eq!(keys(&fast.search(c.as_slice(), &mut w1)), keys(&slow.search(c.as_slice(), &mut w2)), "{alpha} {shape}");
        }
        assert_eq!(w1.counters.ungapped_extensions, w2.counters.ungapped_extensions);
        assert_eq!(w1.counters.hits_filtered, w2.counters.hits_filtered);
    }
}

#[test]
fn scan_schedule_scalar_and_avx2_paths_match_the_index_schedule() {
    // Dense (weight 3) and blocked (weight 7) seed tables, the exact-identity
    // weighted window filter and banded scoring: the scan schedule (AVX2 and
    // forced-scalar paths) must report exactly the query-major results.
    use saber_rs::search::engine::ScanSeeds;
    let (queries, db) = dataset(0x5ca9, 10, 300);
    let shapes = ["11110111", "110010101011", "111"];
    let idx = build(&db, "murphy10", &shapes, MaskingMode::Tantan, 64 << 20);
    let mut p = SearchParams::default();
    p.max_target_seqs = 50;
    p.hit_filter_min_ids = 22;
    p.hit_filter_exact_weight = 2;
    p.score_band = 32;
    p.score_prefilter_band = 8;
    let seeds = || ScanSeeds { alphabet: SeedAlphabet::parse("murphy10").unwrap(), shapes: shapes.iter().map(|s| SeedShape::parse(s).unwrap()).collect() };
    let index_sched = engine(&idx, p.clone());
    let fast = engine(&idx, p.clone()).with_scan_seeds(seeds());
    let mut slow = engine(&idx, p).with_scan_seeds(seeds());
    slow.force_scalar_hit_path();
    assert!(fast.scan_available() && slow.scan_available());
    let encoded: Vec<EncodedProtein> = queries.iter().map(|q| EncodedProtein::from_ascii(q)).collect();
    let views: Vec<&[u8]> = encoded.iter().map(|e| e.as_slice()).collect();
    let mut w = Worker::new();
    let want: Vec<Vec<Key>> = views.iter().map(|q| keys(&index_sched.search(q, &mut w))).collect();
    assert!(want.iter().map(|k| k.len()).sum::<usize>() >= 20, "planted homologs must be found");
    for (name, eng) in [("avx2", &fast), ("scalar", &slow)] {
        let got: Vec<Vec<Key>> = eng.search_scan(&views).iter().map(|(h, _)| keys(h)).collect();
        assert_eq!(got, want, "scan schedule ({name}) differs from the query-major schedule");
    }
}

#[test]
fn neighborhood_seeding_finds_a_superset_of_exact_seed_hits() {
    let (queries, db) = dataset(0x1234, 6, 150);
    let idx = build(&db, "full", &["1111"], MaskingMode::Tantan, 64 << 20);
    let mut p = SearchParams::default();
    p.max_target_seqs = 100;
    let exact = engine(&idx, p.clone());
    p.neighborhood_t = 13;
    let nb = engine(&idx, p);
    assert!(nb.neighborhood_active());
    let (mut w1, mut w2) = (Worker::new(), Worker::new());
    for q in &queries {
        let c = EncodedProtein::from_ascii(q);
        let a: std::collections::HashSet<u64> = exact.search(c.as_slice(), &mut w1).iter().map(|h| h.subject).collect();
        let b: std::collections::HashSet<u64> = nb.search(c.as_slice(), &mut w2).iter().map(|h| h.subject).collect();
        assert!(a.is_subset(&b), "neighborhood seeding lost subjects");
    }
    assert!(w2.counters.neighbor_seeds > 0 && w2.counters.hits > w1.counters.hits);
}
