//! Time the tantan masker alone, single-threaded, and print a checksum of the
//! masks so two builds can be compared. Usage: mask_bench <in.fa> [max_residues]
use saber_rs::alphabet::EncodedProtein;
use saber_rs::io::Fasta;
use saber_rs::masking::{Tantan, TantanScratch};
use saber_rs::scoring::{ScoreMatrix, ScoringSystem};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let fa = Fasta::from_file(&args[1]).expect("fasta");
    let cap: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);
    let m = ScoreMatrix::new(ScoringSystem::Blosum62);
    let t = Tantan::protein(&m);
    let mut seqs = Vec::new();
    let mut tot = 0usize;
    for r in fa.records() {
        let e = EncodedProtein::from_ascii(r.sequence.as_bytes());
        tot += e.as_slice().len();
        seqs.push(e.as_slice().to_vec());
        if tot >= cap {
            break;
        }
    }
    let mut sc = TantanScratch::default();
    let mut mask = Vec::new();
    let mut h: u64 = 0xcbf29ce484222325;
    let mut masked = 0usize;
    let t0 = Instant::now();
    for s in &seqs {
        masked += t.mask(s, &mut sc, &mut mask);
        for (i, &b) in mask.iter().enumerate() {
            if b {
                h = (h ^ i as u64).wrapping_mul(0x100000001b3);
            }
        }
        h = (h ^ 0xff).wrapping_mul(0x100000001b3);
    }
    let dt = t0.elapsed().as_secs_f64();
    println!(
        "sequences {} residues {} masked {} time {:.3} s ({:.1} ns/residue) checksum {:016x}",
        seqs.len(),
        tot,
        masked,
        dt,
        dt * 1e9 / tot as f64,
        h
    );
}
