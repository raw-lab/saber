//! Development probe for the v1.2 engine: builds an in-memory index and runs
//! a search, printing m8 and a counters JSON on stderr. Parameters via env:
//!   ALPHA=full|murphy10  SHAPES=1111,11011  MASK=tantan|none  BIAS=1|0
//!   HSP=40  XDROP=20  TWOHIT=0  K=1000  E=10  EXH=0  MINWIN=1  OCC=0
//!   GO=11 GE=1 SLACK=64 BANDPAD=16 MAXCAND=5000
use rayon::prelude::*;
use saber_rs::alphabet::{EncodedProtein, SeedAlphabet, SeedShape};
use saber_rs::io::stream::FastaStream;
use saber_rs::masking::MaskingMode;
use saber_rs::scoring::stats::{gapped_protein, GapConvention, GapCosts, SearchStatistics};
use saber_rs::scoring::{ScoreMatrix, ScoringSystem};
use saber_rs::search::engine::{Counters, Engine, SearchParams, Worker};
use saber_rs::search::index::{BuildOptions, IndexBuilder, SeedConfig};
use std::io::Write;
use std::time::Instant;

fn env<T: std::str::FromStr>(k: &str, d: T) -> T {
    std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (dbp, qp, outp) = (&a[1], &a[2], &a[3]);
    let t0 = Instant::now();
    let mut cfg = SeedConfig::protein_default();
    cfg.alphabet = SeedAlphabet::parse(&std::env::var("ALPHA").unwrap_or("full".into())).unwrap();
    cfg.shapes = std::env::var("SHAPES").unwrap_or("1111".into()).split(',').map(|s| SeedShape::parse(s).unwrap()).collect();
    cfg.masking = std::env::var("MASK").unwrap_or("tantan".into()).parse().unwrap();
    cfg.minimizer_window = env("MINWIN", 1usize);
    let mut b = IndexBuilder::in_memory(BuildOptions::new(cfg)).unwrap();
    let mut fs = FastaStream::open(dbp).unwrap();
    let mut buf = Vec::new();
    while let Some(r) = fs.next_record().unwrap() {
        let id = r.id.to_string();
        b.add_ascii(&id, r.seq, &mut buf).unwrap();
    }
    let idx = b.finish_in_memory().unwrap();
    let t_build = t0.elapsed().as_secs_f64();
    let m = ScoreMatrix::new(ScoringSystem::Blosum62);
    let gaps = GapCosts::from_user(env("GO", 11), env("GE", 1), GapConvention::Saber).unwrap();
    let kp = gapped_protein(ScoringSystem::Blosum62, gaps, false).unwrap();
    let st = SearchStatistics::new(kp, idx.total_residues(), idx.num_subjects() as u64);
    let mut p = SearchParams::default();
    p.max_target_seqs = env("K", 1000);
    p.evalue = env("E", 10.0);
    p.hsp_threshold = env("HSP", 40);
    p.x_drop_ungapped = env("XDROP", 20);
    p.two_hit_window = env("TWOHIT", 0);
    p.seed_min_score = env("SEEDMIN", 0);
    p.single_hit_score = env("SINGLE", 0);
    p.hit_filter_min_ids = env("HFIDS", 0);
    p.exhaustive = env("EXH", 0) == 1;
    p.comp_bias = env("BIAS", 1) == 1;
    p.query_masking = if std::env::var("QMASK").unwrap_or("tantan".into()) == "none" { MaskingMode::None } else { MaskingMode::Tantan };
    p.seed_max_occ = env("OCC", 0);
    p.window_slack = env("SLACK", 64);
    p.band_pad = env("BANDPAD", 16);
    p.max_candidates = env("MAXCAND", 5000);
    let eng = Engine::new(&idx, ScoreMatrix::new(ScoringSystem::Blosum62), gaps, st, p);
    let _ = &m;
    // queries
    let mut qs = Vec::new();
    let mut fq = FastaStream::open(qp).unwrap();
    while let Some(r) = fq.next_record().unwrap() {
        qs.push((r.id.to_string(), EncodedProtein::from_ascii(r.seq)));
    }
    let t1 = Instant::now();
    let results: Vec<(Vec<String>, Counters)> = qs
        .par_iter()
        .map_init(Worker::new, |w, (id, q)| {
            let before = w.counters;
            let hits = eng.search(q.as_slice(), w);
            let lines = hits
                .iter()
                .map(|h| {
                    format!(
                        "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{:.1}",
                        id, idx.subject_id(h.subject as usize), h.identity(), h.columns, h.mismatches, h.gap_opens,
                        h.q_start, h.q_end, h.s_start, h.s_end, h.evalue, h.bitscore
                    )
                })
                .collect();
            let mut d = w.counters;
            // per-query delta
            macro_rules! sub { ($($f:ident),*) => { $( d.$f -= before.$f; )* } }
            sub!(queries, query_residues, query_residues_masked, seed_positions, seeds_frequent_skipped, postings_visited, hits,
                 subject_range_passes, streamed_subjects, ungapped_extensions, extensions_avoided, two_hit_waits, seeds_low_score, hits_filtered, hsps, regions,
                 regions_dropped_by_cap, regions_below_evalue, simd_cells, scalar_cells, saturation_retries, banded_cells,
                 band_widenings, band_capped, tracebacks, filtered_identity_coverage, alignments, t_prep_ns, t_seed_ns,
                 t_extend_ns, t_score_ns, t_traceback_ns);
            (lines, d)
        })
        .collect();
    let t_search = t1.elapsed().as_secs_f64();
    let mut tot = Counters::default();
    let mut f = std::io::BufWriter::new(std::fs::File::create(outp).unwrap());
    for (lines, c) in &results {
        tot.add(c);
        for l in lines {
            writeln!(f, "{l}").unwrap();
        }
    }
    eprintln!("build {:.2}s search {:.2}s", t_build, t_search);
    eprintln!("{}", serde_json::to_string_pretty(&tot).unwrap());
}
