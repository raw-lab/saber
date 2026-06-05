//! Database loader + k-mer index wrapper.

use crate::algorithm::{LegacyKmerIndex as KmerIndex, KmerSize};
use crate::io::Fasta;
use std::path::Path;

pub struct Database {
    pub sequences: Fasta,
    pub index: DatabaseIndex,
}

impl Database {
    pub fn from_file<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let sequences = Fasta::from_file(path)?;
        let index = DatabaseIndex::new(&sequences);
        Ok(Self { sequences, index })
    }

    pub fn stats(&self) -> DatabaseStats {
        let fs = self.sequences.stats();
        let is = self.index.stats();
        DatabaseStats {
            num_sequences: fs.total_sequences,
            num_residues: fs.total_residues,
            min_seq_length: fs.min_length,
            max_seq_length: fs.max_length,
            mean_seq_length: fs.mean_length,
            total_kmer_positions: is.total_kmer_positions,
            unique_kmers: is.unique_kmers,
            kmer_size: is.kmer_size,
        }
    }
}

pub struct DatabaseIndex {
    pub kmer_index: KmerIndex,
    pub num_sequences: usize,
    pub total_residues: usize,
}

impl DatabaseIndex {
    pub fn new(sequences: &Fasta) -> Self {
        let mut kmer_index = KmerIndex::new(KmerSize::Three, 256);
        for (i, record) in sequences.records().iter().enumerate() {
            kmer_index.index_sequence(record.sequence.as_bytes(), i);
        }
        Self { kmer_index, num_sequences: sequences.len(), total_residues: sequences.total_length() }
    }

    pub fn stats(&self) -> IndexStats {
        let ks = self.kmer_index.stats();
        IndexStats {
            num_sequences: self.num_sequences,
            total_residues: self.total_residues,
            kmer_size: ks.kmer_size,
            unique_kmers: ks.unique_kmers,
            total_kmer_positions: ks.total_kmer_positions,
            max_positions_per_kmer: ks.max_positions_per_kmer,
        }
    }
}

#[derive(Debug)]
pub struct DatabaseStats {
    pub num_sequences: usize,
    pub num_residues: usize,
    pub min_seq_length: usize,
    pub max_seq_length: usize,
    pub mean_seq_length: f64,
    pub total_kmer_positions: usize,
    pub unique_kmers: usize,
    pub kmer_size: usize,
}

#[derive(Debug)]
pub struct IndexStats {
    pub num_sequences: usize,
    pub total_residues: usize,
    pub kmer_size: usize,
    pub unique_kmers: usize,
    pub total_kmer_positions: usize,
    pub max_positions_per_kmer: usize,
}
