//! Micro-benchmarks of the alignment kernels (cells/second).
use saber_rs::algorithm::banded::{banded_align, banded_score, Band, BandScratch, Scheme};
use saber_rs::algorithm::sw_simd::{score_many, KernelCounters, SimdScratch};
use saber_rs::alphabet::ProteinCodes;
use saber_rs::scoring::stats::GapCosts;
use saber_rs::scoring::{ScoreMatrix, ScoringSystem};
use std::time::Instant;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self, m: u64) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 33) % m
    }
}

fn main() {
    let m = ScoreMatrix::new(ScoringSystem::Blosum62);
    let scheme = Scheme { matrix: &m.matrix, gap_open: 11, gap_extend: 1 };
    let mut r = Lcg(1);
    let n_pairs = 2000;
    let len = 300;
    let q: Vec<u8> = (0..len).map(|_| r.next(20) as u8).collect();
    let subs: Vec<Vec<u8>> = (0..n_pairs).map(|_| (0..len).map(|_| r.next(20) as u8).collect()).collect();
    let mut bs = BandScratch::new();
    for half in [16i64, 64] {
        let t = Instant::now();
        let mut cells = 0u64;
        let mut acc = 0i64;
        for s in &subs {
            let x = banded_score(&q, s, None, scheme, Band { lo: -half, hi: half }, &mut bs);
            cells += x.cells;
            acc += x.score as i64;
        }
        let dt = t.elapsed().as_secs_f64();
        println!("banded_score  half={half:3}  {:.2} ns/cell  ({:.3} GCUPS) [{acc}]", dt * 1e9 / cells as f64, cells as f64 / dt / 1e9);
        let t = Instant::now();
        let mut cells = 0u64;
        for s in &subs {
            if let Ok(Some(a)) = banded_align(&q, s, None, scheme, Band { lo: -half, hi: half }, &mut bs, usize::MAX) {
                cells += a.cells;
            }
        }
        let dt = t.elapsed().as_secs_f64();
        println!("banded_align  half={half:3}  {:.2} ns/cell  ({:.3} GCUPS)", dt * 1e9 / cells as f64, cells as f64 / dt / 1e9);
    }
    let views: Vec<ProteinCodes> = subs.iter().map(|s| ProteinCodes::new(s).unwrap()).collect();
    let mut sc = SimdScratch::new();
    let mut kc = KernelCounters::default();
    let t = Instant::now();
    let sc_out = score_many(ProteinCodes::new(&q).unwrap(), None, &views, &m.matrix, GapCosts { open: 11, extend: 1 }, &mut sc, &mut kc);
    let dt = t.elapsed().as_secs_f64();
    println!("score_many (AVX2 16 lanes) {:.3} GCUPS  [{}]", kc.simd_cells as f64 / dt / 1e9, sc_out.iter().map(|&x| x as i64).sum::<i64>());
}
