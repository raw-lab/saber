//! K-mer index for candidate filtering. Used by `db::Database` callers that
//! want to pre-filter the subject set before running full Smith-Waterman.
//!
//! The default CLI runs full S-W against every subject (which is what BLAST's
//! `-task blastp-short` does too) and skips this index. Kept here so external
//! consumers of the library API can opt into k-mer prefiltering.

use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum KmerSize {
    Three,
    Four,
    Five,
}

impl KmerSize {
    pub fn as_usize(&self) -> usize {
        match self {
            KmerSize::Three => 3,
            KmerSize::Four => 4,
            KmerSize::Five => 5,
        }
    }
}

pub struct KmerIndex {
    kmer_positions: HashMap<u64, Vec<(usize, usize)>>, // (seq_id, position)
    kmer_size: usize,
    alphabet_size: usize,
}

impl KmerIndex {
    pub fn new(kmer_size: KmerSize, alphabet_size: usize) -> Self {
        Self {
            kmer_positions: HashMap::new(),
            kmer_size: kmer_size.as_usize(),
            alphabet_size,
        }
    }

    pub fn index_sequence(&mut self, sequence: &[u8], seq_id: usize) {
        let k = self.kmer_size;
        if sequence.len() < k {
            return;
        }
        for pos in 0..=(sequence.len() - k) {
            let kmer = &sequence[pos..pos + k];
            let hash = self.hash_kmer(kmer);
            self.kmer_positions
                .entry(hash)
                .or_insert_with(Vec::new)
                .push((seq_id, pos));
        }
    }

    /// Score candidate database sequences by how many query k-mers hit them.
    /// Returned in descending order of hit count.
    pub fn find_candidates(&self, query: &[u8], max_candidates: usize) -> Vec<(usize, i32)> {
        let k = self.kmer_size;
        if query.len() < k {
            return Vec::new();
        }
        let mut scores: HashMap<usize, i32> = HashMap::new();
        for pos in 0..=(query.len() - k) {
            let hash = self.hash_kmer(&query[pos..pos + k]);
            if let Some(positions) = self.kmer_positions.get(&hash) {
                for (seq_id, _) in positions.iter() {
                    *scores.entry(*seq_id).or_insert(0) += 1;
                }
            }
        }
        let mut v: Vec<(usize, i32)> = scores.into_iter().collect();
        v.sort_by(|a, b| b.1.cmp(&a.1));
        v.truncate(max_candidates);
        v
    }

    fn hash_kmer(&self, kmer: &[u8]) -> u64 {
        let mut h = 0u64;
        for &b in kmer {
            h = h
                .wrapping_mul(self.alphabet_size as u64)
                .wrapping_add(b.to_ascii_uppercase() as u64);
        }
        h
    }

    pub fn stats(&self) -> KmerIndexStats {
        let mut total = 0usize;
        let mut max_p = 0usize;
        for (_, v) in self.kmer_positions.iter() {
            let n = v.len();
            total += n;
            if n > max_p {
                max_p = n;
            }
        }
        KmerIndexStats {
            unique_kmers: self.kmer_positions.len(),
            total_kmer_positions: total,
            max_positions_per_kmer: max_p,
            kmer_size: self.kmer_size,
        }
    }
}

#[derive(Debug)]
pub struct KmerIndexStats {
    pub unique_kmers: usize,
    pub total_kmer_positions: usize,
    pub max_positions_per_kmer: usize,
    pub kmer_size: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_then_find_ranks_by_hits() {
        let mut idx = KmerIndex::new(KmerSize::Three, 256);
        idx.index_sequence(b"ACGTACGTACGT", 0);
        idx.index_sequence(b"TTTTTTTTTTTT", 1);
        let cand = idx.find_candidates(b"ACGTACGT", 10);
        assert_eq!(cand[0].0, 0); // seq 0 should rank highest
    }

    #[test]
    fn empty_query_returns_empty() {
        let mut idx = KmerIndex::new(KmerSize::Three, 256);
        idx.index_sequence(b"ACGT", 0);
        assert!(idx.find_candidates(b"", 10).is_empty());
        assert!(idx.find_candidates(b"A", 10).is_empty());
    }

    #[test]
    fn case_insensitive_hash() {
        let mut idx = KmerIndex::new(KmerSize::Three, 256);
        idx.index_sequence(b"acgtacgt", 0);
        let cand = idx.find_candidates(b"ACGTACGT", 10);
        assert!(!cand.is_empty());
        assert_eq!(cand[0].0, 0);
    }
}
