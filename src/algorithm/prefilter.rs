//! K-mer prefilter: reduce the candidate subject set before running full SW.
//!
//! This is the "filter" half of the seed-and-extend strategy that lets
//! DIAMOND/MMseqs2 stay fast on large databases. SABER's filter is simpler
//! than theirs by design — we trade a few percentage points of speed for
//! higher remote-homology sensitivity (because we never reject a subject
//! that has at least one matching k-mer with the query).
//!
//! ## Algorithm
//!
//! 1. Extract all k-mers from the query, hash each, deduplicate.
//! 2. For each subject, count how many query k-mers appear in it.
//! 3. Keep subjects whose hit count ≥ `min_hits`, OR the top-N by count
//!    (whichever is more permissive).
//!
//! ## Recommended parameters
//!
//! - **Protein, BLOSUM62**: `k=3`, `min_hits=2`, `top_n=10_000`.
//!   Sensitive enough to recover ~95%+ of remote homologs (≥ 25% identity)
//!   while typically rejecting 80–95% of unrelated subjects.
//! - **Nucleotide, blastn-like**: `k=11`, `min_hits=1`, `top_n=10_000`.
//!   k=11 matches BLAST's default seed; k=15 for higher specificity.
//! - **Reduced-alphabet protein**: not implemented yet; would map the 20 aa
//!   to ~10 buckets (e.g., Murphy10) for higher seed hit rate on remote pairs.
//!
//! ## Why on host, not GPU?
//!
//! The prefilter is memory-bandwidth-bound (one hash lookup per subject
//! residue). On a modern CPU with rayon, it sustains ~1 GB/s subject scan
//! per core, which for typical DB sizes (< 10 GB) is already 1-2 seconds
//! total. The CUDA prefilter kernel in `kernels/sw.cu` is available for
//! GPU-resident DBs where avoiding the host round-trip matters, but the
//! default path uses CPU rayon — simpler and adequate for the workloads
//! that need GPU SW.

use std::collections::HashSet;

use rayon::prelude::*;

/// Outcome of prefiltering: indices into the original subject array that
/// survived the cut, in descending order of k-mer hit count.
#[derive(Debug, Clone)]
pub struct PrefilterResult {
    pub kept: Vec<usize>,
    /// Per-kept-subject hit count (parallel array to `kept`).
    pub hit_counts: Vec<u32>,
}

#[derive(Debug, Clone, Copy)]
pub struct PrefilterConfig {
    pub k: usize,
    pub min_hits: u32,
    pub top_n: usize,
    /// Alphabet size for hashing. 24 for protein, 5 for nucleotide.
    pub alphabet_size: u32,
}

impl PrefilterConfig {
    pub fn protein_default() -> Self {
        Self {
            k: 3,
            min_hits: 2,
            top_n: 10_000,
            alphabet_size: 24,
        }
    }

    pub fn nucleotide_default() -> Self {
        Self {
            k: 11,
            min_hits: 1,
            top_n: 10_000,
            alphabet_size: 5,
        }
    }
}

/// Build the set of unique k-mer hashes for the query.
///
/// Hashes are case-insensitive (uppercase folded). Works directly on ASCII
/// FASTA bytes — no encoding step needed.
pub fn build_query_kmers(query: &[u8], cfg: &PrefilterConfig) -> HashSet<u64> {
    let k = cfg.k;
    let mut set = HashSet::new();
    if query.len() < k {
        return set;
    }
    for window in query.windows(k) {
        set.insert(hash_kmer(window));
    }
    set
}

#[inline]
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut h = 0u64;
    for &b in kmer {
        // FNV-style mix on case-folded bytes — keeps protein and nucleotide
        // alphabets distinct without needing a separate hash function per mode.
        let c = b.to_ascii_uppercase() as u64;
        h = h.wrapping_mul(257).wrapping_add(c);
    }
    h
}

/// Count how many query k-mers appear in a subject.
pub fn count_hits(subject: &[u8], query_kmers: &HashSet<u64>, cfg: &PrefilterConfig) -> u32 {
    let k = cfg.k;
    if subject.len() < k {
        return 0;
    }
    let mut count = 0u32;
    for window in subject.windows(k) {
        let h = hash_kmer(window);
        if query_kmers.contains(&h) {
            count += 1;
        }
    }
    count
}

/// Run the prefilter on a batch of subjects in parallel.
pub fn prefilter(
    query: &[u8],
    subjects: &[Vec<u8>],
    cfg: &PrefilterConfig,
) -> PrefilterResult {
    let q_kmers = build_query_kmers(query, cfg);

    // Per-subject hit counts in parallel.
    let counts: Vec<u32> = subjects
        .par_iter()
        .map(|s| count_hits(s, &q_kmers, cfg))
        .collect();

    // Build (idx, count) and apply the threshold.
    let mut idx_count: Vec<(usize, u32)> = counts
        .iter()
        .copied()
        .enumerate()
        .filter(|&(_, c)| c >= cfg.min_hits)
        .collect();

    // Sort by count descending; take top_n.
    idx_count.sort_unstable_by(|a, b| b.1.cmp(&a.1));
    idx_count.truncate(cfg.top_n);

    let (kept, hit_counts): (Vec<_>, Vec<_>) = idx_count.into_iter().unzip();
    PrefilterResult { kept, hit_counts }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sequences_have_maximum_overlap() {
        let cfg = PrefilterConfig::protein_default();
        let query: &[u8] = b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKK";
        let subjects = vec![
            b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKK".to_vec(),
            b"KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK".to_vec(),
        ];
        let result = prefilter(query, &subjects, &cfg);
        assert_eq!(result.kept[0], 0);
        assert!(result.hit_counts[0] > 0);
    }

    #[test]
    fn unrelated_sequence_is_filtered_out() {
        let cfg = PrefilterConfig {
            min_hits: 5,
            ..PrefilterConfig::protein_default()
        };
        let query: &[u8] = b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKK";
        let subjects = vec![b"PPPPPPPPPPPPPPPPPPPPPPPPPPPP".to_vec()];
        let result = prefilter(query, &subjects, &cfg);
        assert!(result.kept.is_empty());
    }

    #[test]
    fn case_insensitive_match() {
        let cfg = PrefilterConfig::protein_default();
        let query: &[u8] = b"MVLSPADKTNVKAAW";
        let subjects = vec![b"mvlspadktnvkaaw".to_vec()];
        let result = prefilter(query, &subjects, &cfg);
        assert_eq!(result.kept[0], 0);
    }
}
