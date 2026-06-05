//! Neighborhood word generation — BLAST's T-score k-mer expansion.
//!
//! For each query k-mer, generate all k-mers that score ≥ T against it
//! under the chosen substitution matrix. These "neighbor" k-mers become
//! additional seeds in the index lookup. With T=11 for BLOSUM62 (BLAST's
//! blastp default), a typical 4-mer has ~10-50 neighbors.
//!
//! ## Why this matters
//!
//! Exact-match-only k-mer indexing requires `min_hits=1` to find remote
//! homologs (because 25%-identity proteins share very few exact 4-mers).
//! That makes ~95% of subjects into candidates and negates the index's
//! filtering value.
//!
//! With neighborhood expansion, real homologs have multiple "near match"
//! seeds clustering on the same diagonal, while random subjects don't.
//! `min_hits=2` (BLAST's two-hit heuristic) then becomes selective enough
//! to filter ≥99% of subjects without losing sensitivity.
//!
//! ## Algorithm
//!
//! Depth-first enumeration over the k positions. At each depth `d`, try
//! each candidate residue in score-descending order. Prune the branch if
//! `score_so_far + best_possible_remaining < threshold`. Since each row
//! is pre-sorted by score, once a branch is pruned all subsequent residues
//! at that depth are also pruneable → `break` instead of `continue`.

use crate::scoring::matrices::ScoreMatrix;

/// Number of standard amino acids in the protein neighborhood (we exclude
/// B/Z/X/J/*). Same as `PROTEIN_ALPHABET_FOR_INDEX` in `index.rs`.
const PROTEIN_AA_COUNT: usize = 20;

pub struct NeighborhoodGenerator {
    k: usize,
    threshold: i32,
    /// For each query residue q ∈ 0..20, the substitution row sorted by
    /// score descending: `(target_residue, score)`. Pre-sorted to enable
    /// early-termination pruning.
    sorted_by_score: Vec<Vec<(u8, i16)>>,
    /// Best (highest) score per query residue — used for the
    /// best-possible-remaining bound during pruning.
    best_score: Vec<i16>,
}

impl NeighborhoodGenerator {
    /// Construct a generator for k-mer length `k` and score threshold T.
    /// T=11 is the BLAST blastp default for BLOSUM62; T=13 matches SWORD's
    /// default; lower T = more neighbors (more sensitive, slower).
    pub fn new(matrix: &ScoreMatrix, k: usize, threshold: i32) -> Self {
        assert!(k >= 2 && k <= 8, "k must be in [2, 8]");
        let mut sorted = Vec::with_capacity(PROTEIN_AA_COUNT);
        let mut best = Vec::with_capacity(PROTEIN_AA_COUNT);
        for q in 0..PROTEIN_AA_COUNT {
            let mut row: Vec<(u8, i16)> = (0..PROTEIN_AA_COUNT)
                .map(|d| (d as u8, matrix.score_by_index(q, d) as i16))
                .collect();
            row.sort_by(|a, b| b.1.cmp(&a.1));
            best.push(row[0].1);
            sorted.push(row);
        }
        Self {
            k,
            threshold,
            sorted_by_score: sorted,
            best_score: best,
        }
    }

    pub fn k(&self) -> usize { self.k }
    pub fn threshold(&self) -> i32 { self.threshold }

    /// Generate all k-mer hashes whose alignment score against `query_kmer`
    /// is ≥ threshold. Hashes use `alphabet_size`-radix encoding (same as
    /// [`crate::algorithm::index::KmerIndex`]). Output is appended to `out`.
    ///
    /// `query_kmer` must contain encoded residues 0..20 (ambiguous residues
    /// like X cause the k-mer to be skipped entirely — returns without
    /// pushing anything).
    pub fn expand(&self, query_kmer: &[u8], alphabet_size: u32, out: &mut Vec<u32>) {
        debug_assert_eq!(query_kmer.len(), self.k);
        // Skip ambiguous k-mers entirely (caller may want to do this before
        // calling, but we double-check for safety).
        if query_kmer.iter().any(|&b| (b as usize) >= PROTEIN_AA_COUNT) {
            return;
        }

        // best_remaining[i] = sum of best_score[query_kmer[i..]]; used as the
        // optimistic upper bound on the score that can still be added.
        let mut best_remaining = vec![0i32; self.k + 1];
        for i in (0..self.k).rev() {
            let q = query_kmer[i] as usize;
            best_remaining[i] = best_remaining[i + 1] + self.best_score[q] as i32;
        }

        self.expand_rec(query_kmer, 0, 0, 0, &best_remaining, alphabet_size, out);
    }

    fn expand_rec(
        &self,
        query_kmer: &[u8],
        depth: usize,
        score_so_far: i32,
        hash_so_far: u32,
        best_remaining: &[i32],
        alphabet_size: u32,
        out: &mut Vec<u32>,
    ) {
        if depth == self.k {
            if score_so_far >= self.threshold {
                out.push(hash_so_far);
            }
            return;
        }
        let q = query_kmer[depth] as usize;
        let row = &self.sorted_by_score[q];
        let needed_remaining = self.threshold - score_so_far;

        for &(b, s) in row.iter() {
            // Best case from here: this residue's score + best of all later positions.
            let best_total_from_here = s as i32 + best_remaining[depth + 1];
            if best_total_from_here < needed_remaining {
                // Row is sorted descending → no later residue can reach threshold either.
                break;
            }
            let new_score = score_so_far + s as i32;
            let new_hash = hash_so_far * alphabet_size + b as u32;
            self.expand_rec(
                query_kmer,
                depth + 1,
                new_score,
                new_hash,
                best_remaining,
                alphabet_size,
                out,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::ScoringSystem;

    #[test]
    fn exact_kmer_is_always_in_its_own_neighborhood() {
        let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
        let gen = NeighborhoodGenerator::new(&matrix, 4, 11);
        // Encoded "MVLS" → indices 12, 19, 10, 15 (M, V, L, S in 0..19 order)
        let kmer = [12u8, 19, 10, 15];
        let mut nbrs = Vec::new();
        gen.expand(&kmer, 20, &mut nbrs);

        // The exact k-mer's own hash must be present.
        let self_hash = kmer.iter().fold(0u32, |h, &b| h * 20 + b as u32);
        assert!(nbrs.contains(&self_hash), "self hash missing from neighborhood");
    }

    #[test]
    fn higher_threshold_yields_fewer_neighbors() {
        let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
        let kmer = [12u8, 19, 10, 15]; // MVLS

        let mut low = Vec::new();
        NeighborhoodGenerator::new(&matrix, 4, 5).expand(&kmer, 20, &mut low);

        let mut mid = Vec::new();
        NeighborhoodGenerator::new(&matrix, 4, 11).expand(&kmer, 20, &mut mid);

        let mut high = Vec::new();
        NeighborhoodGenerator::new(&matrix, 4, 20).expand(&kmer, 20, &mut high);

        assert!(low.len() >= mid.len(), "T=5 should yield at least as many as T=11");
        assert!(mid.len() >= high.len(), "T=11 should yield at least as many as T=20");

        // T=5 should be substantially bigger than T=20 for any reasonable k-mer.
        assert!(low.len() > high.len(), "got T=5 → {}, T=20 → {}", low.len(), high.len());
    }

    #[test]
    fn unique_neighbors() {
        let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
        let gen = NeighborhoodGenerator::new(&matrix, 4, 11);
        let kmer = [12u8, 19, 10, 15];
        let mut nbrs = Vec::new();
        gen.expand(&kmer, 20, &mut nbrs);
        let n_before = nbrs.len();
        nbrs.sort_unstable();
        nbrs.dedup();
        assert_eq!(nbrs.len(), n_before, "neighbors should already be unique");
    }

    #[test]
    fn neighborhood_size_is_reasonable_for_t11() {
        // BLAST blastp default T=11. Typical k=3-mers have 50-200 neighbors,
        // k=4-mers have 10-100 depending on residue rarity.
        let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
        let gen = NeighborhoodGenerator::new(&matrix, 4, 11);
        let kmer = [12u8, 19, 10, 15]; // MVLS — common-ish residues
        let mut nbrs = Vec::new();
        gen.expand(&kmer, 20, &mut nbrs);
        assert!(nbrs.len() >= 1, "should have at least the self-hash");
        assert!(nbrs.len() < 5000, "way too many neighbors: {}", nbrs.len());
    }

    #[test]
    fn empty_for_ambiguous_kmer() {
        let matrix = ScoreMatrix::new(ScoringSystem::Blosum62);
        let gen = NeighborhoodGenerator::new(&matrix, 4, 11);
        // 255 means "ambiguous residue from encoding"
        let kmer = [12u8, 19, 255, 15];
        let mut nbrs = Vec::new();
        gen.expand(&kmer, 20, &mut nbrs);
        assert!(nbrs.is_empty());
    }
}
