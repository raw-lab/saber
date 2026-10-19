//! Head-to-head comparison of SABER's Smith-Waterman kernel against two
//! established Rust alignment libraries, on identical (query, subject) pairs.
//!
//! * **rust-bio** (`bio::alignment::pairwise`) — the reference scalar affine-gap
//!   Smith-Waterman implementation most Rust bioinformatics code reaches for.
//! * **block-aligner** — a SIMD adaptive-banded aligner, i.e. the closest
//!   published Rust analogue to what SABER's kernel is doing.
//!
//! Gap model note: rust-bio and block-aligner both charge `open + k*extend`
//! for a gap of length k, while SABER (and BLAST) charge
//! `open + (k-1)*extend`. Passing them `open = -(11-1) = -10, extend = -1`
//! makes all three price a length-k gap identically at `-(10+k)`, so the
//! scores are directly comparable.

use bio::alignment::pairwise::Aligner as BioAligner;
use bio::scores::blosum62;
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use saber_rs::algorithm::index::KmerIndex;
use saber_rs::io::fasta::Fasta;
use saber_rs::scoring::matrices::ScoreMatrix;
use saber_rs::scoring::ScoringSystem;
use std::time::Instant;

/// Total DP cells for a set of pairs — the denominator for GCUPS.
fn cells(pairs: &[(Vec<u8>, Vec<u8>)]) -> f64 {
    pairs.iter().map(|(q, s)| (q.len() * s.len()) as f64).sum()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let db_path = &args[1];
    let q_path = &args[2];
    let max_pairs: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(2000);

    let db = Fasta::from_file(db_path).unwrap();
    let queries = Fasta::from_file(q_path).unwrap();
    let subjects: Vec<&[u8]> = db.records().iter().map(|r| r.sequence.as_bytes()).collect();

    // Build the pair set the way the real pipeline does: index → seed → cluster.
    // Benchmarking on the pairs SABER actually aligns is the honest comparison;
    // random pairs would understate everyone's cache behaviour equally but
    // wouldn't reflect the real length distribution.
    let idx = KmerIndex::build_protein(&subjects, 4);
    let mut pairs: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    'outer: for q in queries.records() {
        let hits = idx.query_protein(q.sequence.as_bytes());
        let cands = idx.candidates(&hits, 2, 16, 40);
        for c in &cands {
            pairs.push((
                q.sequence.as_bytes().to_vec(),
                subjects[c.subject_idx as usize].to_vec(),
            ));
            if pairs.len() >= max_pairs {
                break 'outer;
            }
        }
    }

    let total_cells = cells(&pairs);
    println!(
        "pairs={}  total_cells={:.2e}  mean_q={:.0}aa  mean_s={:.0}aa",
        pairs.len(),
        total_cells,
        pairs.iter().map(|p| p.0.len()).sum::<usize>() as f64 / pairs.len() as f64,
        pairs.iter().map(|p| p.1.len()).sum::<usize>() as f64 / pairs.len() as f64,
    );
    println!();

    // ---------------- SABER: SIMD inter-sequence kernel ---------------------
    // Batched one query against many subjects, which is the shape it is built
    // for and the shape the pipeline uses.
    let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
    let t = Instant::now();
    let mut saber_scores: Vec<i32> = Vec::with_capacity(pairs.len());
    {
        let mut i = 0;
        while i < pairs.len() {
            // Group consecutive pairs that share a query.
            let mut j = i;
            while j < pairs.len() && pairs[j].0 == pairs[i].0 {
                j += 1;
            }
            let q_enc = saber_rs::gpu::encoding::encode_protein(&pairs[i].0);
            let subs: Vec<Vec<u8>> = pairs[i..j]
                .iter()
                .map(|p| saber_rs::gpu::encoding::encode_protein(&p.1))
                .collect();
            let refs: Vec<&[u8]> = subs.iter().map(|s| s.as_slice()).collect();
            let sc = saber_rs::algorithm::sw_simd::sw_simd_protein_batch_encoded(
                &q_enc, &refs, &matrix, 11, 1,
            );
            saber_scores.extend(sc);
            i = j;
        }
    }
    let t_saber = t.elapsed().as_secs_f64();

    // ---------------- rust-bio: scalar affine SW ----------------------------
    let t = Instant::now();
    let mut bio_scores: Vec<i32> = Vec::with_capacity(pairs.len());
    for (q, s) in &pairs {
        let mut a = BioAligner::with_capacity(q.len(), s.len(), -10, -1, &blosum62);
        bio_scores.push(a.local(q, s).score);
    }
    let t_bio = t.elapsed().as_secs_f64();

    // ---------------- block-aligner: SIMD adaptive banded -------------------
    // X-drop + local-start is its closest analogue to Smith-Waterman. It is an
    // adaptive-band heuristic, so it is allowed to under-report on alignments
    // that wander outside the band; we report how often that happens rather
    // than hiding it.
    let block_size = 32usize;
    let max_block = 256usize;
    let gaps = Gaps { open: -11, extend: -1 };
    let t = Instant::now();
    let mut ba_scores: Vec<i32> = Vec::with_capacity(pairs.len());
    for (q, s) in &pairs {
        let pq = PaddedBytes::from_bytes::<AAMatrix>(q, max_block);
        let ps = PaddedBytes::from_bytes::<AAMatrix>(s, max_block);
        let mut blk = Block::<false, true, true>::new(q.len(), s.len(), max_block);
        blk.align(&pq, &ps, &BLOSUM62, gaps, block_size..=max_block, 50);
        ba_scores.push(blk.res().score);
    }
    let t_ba = t.elapsed().as_secs_f64();

    // ---------------- Report -------------------------------------------------
    let gcups = |t: f64| total_cells / t / 1e9;
    println!("{:<24} {:>9} {:>9} {:>12}", "ENGINE", "sec", "GCUPS", "vs SABER");
    println!(
        "{:<24} {:>9.3} {:>9.2} {:>12}",
        "SABER (SIMD AVX2)", t_saber, gcups(t_saber), "1.00x"
    );
    println!(
        "{:<24} {:>9.3} {:>9.2} {:>11.2}x",
        "rust-bio (scalar SW)", t_bio, gcups(t_bio), t_bio / t_saber
    );
    println!(
        "{:<24} {:>9.3} {:>9.2} {:>11.2}x",
        "block-aligner (SIMD)", t_ba, gcups(t_ba), t_ba / t_saber
    );
    println!();

    // Correctness: rust-bio's full DP is the ground truth for exact SW.
    let mut exact = 0usize;
    let mut saber_lower = 0usize;
    let mut saber_higher = 0usize;
    let mut ba_exact = 0usize;
    let mut ba_lower = 0usize;
    for i in 0..pairs.len() {
        if saber_scores[i] == bio_scores[i] {
            exact += 1;
        } else if saber_scores[i] < bio_scores[i] {
            saber_lower += 1;
        } else {
            saber_higher += 1;
        }
        if ba_scores[i] == bio_scores[i] {
            ba_exact += 1;
        } else if ba_scores[i] < bio_scores[i] {
            ba_lower += 1;
        }
    }
    let n = pairs.len() as f64;
    println!("SCORE AGREEMENT vs rust-bio full DP (exact Smith-Waterman):");
    println!(
        "  SABER          exact {:>6} ({:5.1}%)   lower {:>5}   higher {:>5}",
        exact,
        100.0 * exact as f64 / n,
        saber_lower,
        saber_higher
    );
    println!(
        "  block-aligner  exact {:>6} ({:5.1}%)   lower {:>5}   higher {:>5}",
        ba_exact,
        100.0 * ba_exact as f64 / n,
        ba_lower,
        pairs.len() - ba_exact - ba_lower
    );
}
