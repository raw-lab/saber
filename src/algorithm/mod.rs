//! Sequence alignment algorithms.

pub mod smith_waterman;
pub mod kmer;
pub mod translation;
pub mod prefilter;
pub mod sw_simd;
pub mod index;
pub mod neighborhood;
pub mod ungapped;

pub use smith_waterman::{SmithWaterman, AlignmentResult, AlignmentSet, Scoring};
pub use kmer::{KmerIndex as LegacyKmerIndex, KmerSize};
pub use translation::{
    TranslationTable, TranslatedFrame,
    translate_frame, translate_six_frames, translate_dna_to_protein, reverse_complement,
};
pub use prefilter::{PrefilterConfig, PrefilterResult, prefilter as run_prefilter};
pub use sw_simd::sw_simd_protein_batch;
pub use index::{KmerIndex, SeedHit, Candidate, DEFAULT_K_PROTEIN, DEFAULT_K_NUCLEOTIDE};
pub use neighborhood::NeighborhoodGenerator;
pub use ungapped::{extend_ungapped, extend_ungapped_with_range, DEFAULT_X_DROP, DEFAULT_HSP_THRESHOLD};

use std::fmt;

/// What we're aligning, and which scoring backend to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SequenceMode {
    /// Protein query against protein database (blastp).
    ProteinToProtein,
    /// Nucleotide query against nucleotide database (blastn).
    NucleotideToNucleotide,
    /// 6-frame translated nucleotide query against protein database (blastx).
    NucleotideToProtein,
    /// Protein query against 6-frame translated nucleotide database (tblastn).
    ProteinToNucleotide,
}

impl fmt::Display for SequenceMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ProteinToProtein => write!(f, "blastp (P->P)"),
            Self::NucleotideToNucleotide => write!(f, "blastn (N->N)"),
            Self::NucleotideToProtein => write!(f, "blastx (N->P)"),
            Self::ProteinToNucleotide => write!(f, "tblastn (P->N)"),
        }
    }
}

impl std::str::FromStr for SequenceMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "pp" | "blastp" | "protein-protein" => Ok(Self::ProteinToProtein),
            "nn" | "blastn" | "nucleotide-nucleotide" => Ok(Self::NucleotideToNucleotide),
            "nx" | "blastx" | "nucleotide-protein" => Ok(Self::NucleotideToProtein),
            "pn" | "tblastn" | "protein-nucleotide" => Ok(Self::ProteinToNucleotide),
            other => Err(format!("Unknown sequence mode: {} (expected pp/nn/nx/pn)", other)),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AlignmentAlgorithm {
    SmithWaterman,
    NeedlemanWunsch,
    SemiGlobal,
    Overlap,
}

impl fmt::Display for AlignmentAlgorithm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SmithWaterman => write!(f, "SW"),
            Self::NeedlemanWunsch => write!(f, "NW"),
            Self::SemiGlobal => write!(f, "HW"),
            Self::Overlap => write!(f, "OV"),
        }
    }
}

impl std::str::FromStr for AlignmentAlgorithm {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "SW" => Ok(Self::SmithWaterman),
            "NW" => Ok(Self::NeedlemanWunsch),
            "HW" => Ok(Self::SemiGlobal),
            "OV" => Ok(Self::Overlap),
            other => Err(format!("Unknown alignment algorithm: {}", other)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn sequence_mode_aliases_work() {
        assert_eq!("pp".parse::<SequenceMode>().unwrap(), SequenceMode::ProteinToProtein);
        assert_eq!("blastp".parse::<SequenceMode>().unwrap(), SequenceMode::ProteinToProtein);
        assert_eq!("blastn".parse::<SequenceMode>().unwrap(), SequenceMode::NucleotideToNucleotide);
        assert_eq!("nx".parse::<SequenceMode>().unwrap(), SequenceMode::NucleotideToProtein);
    }

    #[test]
    fn unknown_mode_errors() {
        assert!("xx".parse::<SequenceMode>().is_err());
    }

    #[test]
    fn algorithm_parsing() {
        assert_eq!(AlignmentAlgorithm::from_str("SW").unwrap(), AlignmentAlgorithm::SmithWaterman);
    }
}
