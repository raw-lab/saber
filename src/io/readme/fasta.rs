/// FASTA file reader and utilities
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use serde::{Deserialize, Serialize};

/// FASTA sequence record
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastaRecord {
    /// Sequence ID (header without '>')
    pub id: String,
    /// Full header line
    pub description: String,
    /// Sequence data
    pub sequence: String,
    /// Sequence length
    pub length: usize,
}

impl FastaRecord {
    /// Create a new FASTA record
    pub fn new(id: String, description: String, sequence: String) -> Self {
        let length = sequence.len();
        Self {
            id,
            description,
            sequence,
            length,
        }
    }

    /// Get the sequence as bytes
    pub fn as_bytes(&self) -> Vec<u8> {
        self.sequence.as_bytes().to_vec()
    }

    /// Get length of sequence
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if sequence is empty
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// FASTA file reader
pub struct FastaReader<R: Read> {
    reader: BufReader<R>,
}

impl<R: Read> FastaReader<R> {
    /// Create a new FASTA reader from any Read source
    pub fn new(reader: R) -> Self {
        Self {
            reader: BufReader::new(reader),
        }
    }

    /// Read all records from the file
    pub fn read_all(&mut self) -> std::io::Result<Vec<FastaRecord>> {
        let mut records = Vec::new();
        let mut current_id = String::new();
        let mut current_desc = String::new();
        let mut current_seq = String::new();

        let mut line = String::new();
        loop {
            line.clear();
            let bytes_read = self.reader.read_line(&mut line)?;
            if bytes_read == 0 {
                break;
            }
            let line = line.trim_end();

            if line.starts_with('>') {
                // Save previous record if exists
                if !current_id.is_empty() {
                    records.push(FastaRecord::new(
                        current_id.clone(),
                        current_desc.clone(),
                        current_seq.clone(),
                    ));
                    current_seq.clear();
                }

                // Parse new header
                let header = &line[1..];
                let parts: Vec<&str> = header.splitn(2, ' ').collect();
                current_id = parts[0].to_string();
                current_desc = header.to_string();
            } else if !line.is_empty() {
                // Append to sequence
                current_seq.push_str(&line.to_uppercase());
            }
        }

        // Save last record
        if !current_id.is_empty() {
            records.push(FastaRecord::new(current_id, current_desc, current_seq));
        }

        Ok(records)
    }
}

/// Simple FASTA reader wrapper
pub struct Fasta {
    records: Vec<FastaRecord>,
    _index: usize,
}

impl Fasta {
    /// Load FASTA file from path
    pub fn from_file<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let mut reader = FastaReader::new(file);
        let records = reader.read_all()?;
        Ok(Self {
            records,
            _index: 0,
        })
    }

    /// Load FASTA from string
    pub fn from_string(content: &str) -> std::io::Result<Self> {
        let mut reader = FastaReader::new(content.as_bytes());
        let records = reader.read_all()?;
        Ok(Self {
            records,
            _index: 0,
        })
    }

    /// Get number of records
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Get record by index
    pub fn get(&self, index: usize) -> Option<&FastaRecord> {
        self.records.get(index)
    }

    /// Get all records
    pub fn records(&self) -> &[FastaRecord] {
        &self.records
    }

    /// Free the ASCII sequence strings of every record while keeping their
    /// IDs, descriptions, and lengths intact. Used after an inverted k-mer
    /// index has been built — the index already holds a compact encoded
    /// copy of every subject, so the original ASCII strings can be dropped
    /// to free memory. Saves ~1 byte per residue (typically several MB on
    /// any non-trivial DB).
    ///
    /// After this call, `records()[i].sequence` will be an empty string and
    /// callers must use the index's encoded subjects instead. `id`, `desc`,
    /// and `length` remain valid (the length pre-clear is preserved).
    pub fn drop_sequence_strings(&mut self) {
        for r in &mut self.records {
            // r.length was set at parse time; keep it. r.sequence: free its
            // backing buffer by replacing with an empty String.
            r.sequence = String::new();
        }
    }

    /// Total sequence length
    pub fn total_length(&self) -> usize {
        self.records.iter().map(|r| r.length).sum()
    }

    /// Get statistics
    pub fn stats(&self) -> FastaStats {
        let total_seqs = self.records.len();
        let total_residues: usize = self.records.iter().map(|r| r.length).sum();
        let min_length = self.records.iter().map(|r| r.length).min().unwrap_or(0);
        let max_length = self.records.iter().map(|r| r.length).max().unwrap_or(0);
        let mean_length = if total_seqs > 0 {
            total_residues as f64 / total_seqs as f64
        } else {
            0.0
        };

        FastaStats {
            total_sequences: total_seqs,
            total_residues,
            min_length,
            max_length,
            mean_length,
        }
    }
}

/// Statistics about FASTA file
#[derive(Debug)]
pub struct FastaStats {
    pub total_sequences: usize,
    pub total_residues: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub mean_length: f64,
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_parsing() {
        let content = ">seq1\nACGT\n>seq2\nTGCA\n";
        let fasta = Fasta::from_string(content).unwrap();
        assert_eq!(fasta.len(), 2);
        assert_eq!(fasta.get(0).unwrap().id, "seq1");
    }

    #[test]
    fn test_fasta_stats() {
        let content = ">seq1\nACGTACGT\n>seq2\nTGCATGCA\n";
        let fasta = Fasta::from_string(content).unwrap();
        let stats = fasta.stats();
        assert_eq!(stats.total_sequences, 2);
        assert_eq!(stats.total_residues, 16);
    }
}
