//! Micro-benchmark of the banded affine-gap kernels (traceback and
//! score-only) on related protein pairs, reported as ns per band cell.
//! `cargo run --release --example bench_banded [len] [half_band] [pairs]`
use saber_rs::algorithm::banded::{banded_align, banded_score, Band, BandScratch, Scheme};
use saber_rs::scoring::{ScoreMatrix, ScoringSystem};
use std::time::Instant;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self, m: u64) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 33) % m
    }
}

fn related(r: &mut Lcg, a: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < a.len() {
        let x = r.next(100);
        if x < 3 {
            for _ in 0..1 + r.next(6) {
                out.push(r.next(20) as u8);
            }
        } else if x < 6 {
            i += 1 + r.next(6) as usize;
            continue;
        } else if x < 50 {
            out.push(r.next(20) as u8);
        } else {
            out.push(a[i]);
        }
        i += 1;
    }
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let len: usize = args.get(1).map_or(300, |s| s.parse().unwrap());
    let half: i64 = args.get(2).map_or(64, |s| s.parse().unwrap());
    let pairs: usize = args.get(3).map_or(3000, |s| s.parse().unwrap());
    let m = ScoreMatrix::new(ScoringSystem::Blosum62);
    let scheme = Scheme { matrix: &m.matrix, gap_open: 11, gap_extend: 1 };
    let mut r = Lcg(7);
    let data: Vec<(Vec<u8>, Vec<u8>)> = (0..pairs)
        .map(|_| {
            let a: Vec<u8> = (0..len).map(|_| r.next(20) as u8).collect();
            let b = related(&mut r, &a);
            (a, b)
        })
        .collect();
    let band = Band { lo: -half, hi: half };
    let mut sc = BandScratch::new();
    for round in 0..3 {
        let t0 = Instant::now();
        let mut cells = 0u64;
        let mut chk = 0i64;
        for (a, b) in &data {
            if let Ok(Some(al)) = banded_align(a, b, None, scheme, band, &mut sc, usize::MAX) {
                cells += al.cells;
                chk += al.score as i64;
            }
        }
        let tb = t0.elapsed().as_secs_f64();
        let t1 = Instant::now();
        let mut cells2 = 0u64;
        let mut chk2 = 0i64;
        for (a, b) in &data {
            let s = banded_score(a, b, None, scheme, band, &mut sc);
            cells2 += s.cells;
            chk2 += s.score as i64;
        }
        let so = t1.elapsed().as_secs_f64();
        println!(
            "round {round}: traceback {:.3} ns/cell ({:.1} M cells, checksum {chk}); score-only {:.3} ns/cell (checksum {chk2})",
            tb * 1e9 / cells as f64,
            cells as f64 / 1e6,
            so * 1e9 / cells2 as f64
        );
    }
}
