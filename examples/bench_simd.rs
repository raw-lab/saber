//! Quick benchmark: SIMD vs scalar SW throughput on protein data.
use std::time::Instant;
use saber_rs::algorithm::{sw_simd_protein_batch, SmithWaterman};
use saber_rs::gpu::encoding::encode_protein;
use saber_rs::scoring::{ScoreMatrix, ScoringSystem};
use std::sync::Arc;

fn random_seq(len: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    let alpha = b"ARNDCQEGHILKMFPSTWYV";
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        out.push(alpha[(state >> 33) as usize % alpha.len()]);
    }
    out
}

fn main() {
    let matrix = Arc::new(ScoreMatrix::new(ScoringSystem::Blosum62));
    let query_ascii = random_seq(300, 42);          // 300-aa query
    let query_enc = encode_protein(&query_ascii);

    // 2000 subjects, lengths 200-400, similar distribution to a small DB shard
    let n_subjects = 2000;
    let subjects_ascii: Vec<Vec<u8>> = (0..n_subjects)
        .map(|i| random_seq(200 + (i % 200), 1000 + i as u64))
        .collect();
    let subjects_enc: Vec<Vec<u8>> = subjects_ascii.iter().map(|s| encode_protein(s)).collect();

    let total_cells: u64 = subjects_ascii.iter()
        .map(|s| (query_ascii.len() as u64) * (s.len() as u64))
        .sum();
    println!("Benchmark: 1 query × {} subjects, query_len={}, total DP cells={}M",
        n_subjects, query_ascii.len(), total_cells / 1_000_000);

    // --- Scalar (current SABER SmithWaterman) ---
    let aligner = SmithWaterman::protein(11, 1, Arc::clone(&matrix));
    let t0 = Instant::now();
    let mut scalar_total = 0i64;
    for (q_ascii, subj_ascii) in std::iter::repeat(&query_ascii).zip(subjects_ascii.iter()) {
        let res = aligner.align(q_ascii, subj_ascii);
        scalar_total += res.score as i64;
    }
    let scalar_dt = t0.elapsed();
    let scalar_gcups = total_cells as f64 / scalar_dt.as_secs_f64() / 1e9;
    println!("  Scalar (full traceback):  {:.3}s  ({:.2} GCUPS)  checksum={}",
        scalar_dt.as_secs_f64(), scalar_gcups, scalar_total);

    // --- SIMD score-only ---
    let t0 = Instant::now();
    let simd_scores = sw_simd_protein_batch(&query_enc, &subjects_enc, &matrix, 11, 1);
    let simd_dt = t0.elapsed();
    let simd_gcups = total_cells as f64 / simd_dt.as_secs_f64() / 1e9;
    let simd_total: i64 = simd_scores.iter().map(|&s| s as i64).sum();
    println!("  SIMD AVX2 (score-only):   {:.3}s  ({:.2} GCUPS)  checksum={}",
        simd_dt.as_secs_f64(), simd_gcups, simd_total);

    println!("  Speedup: {:.1}x", scalar_dt.as_secs_f64() / simd_dt.as_secs_f64());
}
