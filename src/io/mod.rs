/// I/O utilities for reading and writing sequence alignments
pub mod fasta;
pub mod output;

pub use fasta::{Fasta, FastaRecord, FastaReader};
pub use output::{Writer, M8Writer, M9Writer, M0Writer, SAMWriter, BAMWriter, JSONWriter, CSVWriter};

use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum OutputFormat {
    BlastM0,  // Formatted alignment blocks
    BlastM8,  // Tab-delimited
    BlastM9,  // Tab-delimited with comments
    SAM,      // Sequence Alignment Map
    BAM,      // Binary SAM
    JSON,     // JSON format
    CSV,      // Comma-separated
}

impl fmt::Display for OutputFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::BlastM0 => write!(f, "m0"),
            Self::BlastM8 => write!(f, "m8"),
            Self::BlastM9 => write!(f, "m9"),
            Self::SAM => write!(f, "sam"),
            Self::BAM => write!(f, "bam"),
            Self::JSON => write!(f, "json"),
            Self::CSV => write!(f, "csv"),
        }
    }
}

impl std::str::FromStr for OutputFormat {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "m0" | "bm0" => Ok(Self::BlastM0),
            "m8" | "bm8" => Ok(Self::BlastM8),
            "m9" | "bm9" => Ok(Self::BlastM9),
            "sam" => Ok(Self::SAM),
            "bam" => Ok(Self::BAM),
            "json" => Ok(Self::JSON),
            "csv" => Ok(Self::CSV),
            _ => Err(format!("Unknown output format: {}", s)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn test_output_format_parsing() {
        assert_eq!(OutputFormat::from_str("m8").unwrap(), OutputFormat::BlastM8);
        assert_eq!(OutputFormat::from_str("sam").unwrap(), OutputFormat::SAM);
    }
}
