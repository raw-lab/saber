//! Print tantan-style per-residue repeat probabilities (same layout as
//! `tantan -p -f 1`) so SABER's masker can be validated against the reference
//! tantan binary. Usage: mask_probs <in.fa> [lambda]
use saber_rs::alphabet::EncodedProtein;
use saber_rs::io::Fasta;
use saber_rs::masking::{Tantan, TantanParams, TantanScratch};
use saber_rs::scoring::{ScoreMatrix, ScoringSystem};
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let fa = Fasta::from_file(&args[1]).expect("fasta");
    let m = ScoreMatrix::new(ScoringSystem::Blosum62);
    let lambda = args.get(2).map(|s| s.parse::<f64>().unwrap());
    let t = Tantan::with_params(&m, TantanParams::default(), lambda);
    eprintln!("lambda = {}", t.lambda());
    let mut sc = TantanScratch::default();
    let out = std::io::stdout();
    let mut w = std::io::BufWriter::new(out.lock());
    for r in fa.records() {
        let e = EncodedProtein::from_ascii(r.sequence.as_bytes());
        let p = t.repeat_probs(e.as_slice(), &mut sc);
        writeln!(w, ">{}", r.id).unwrap();
        for x in p {
            writeln!(w, "{:.3}", x).unwrap();
        }
    }
}
