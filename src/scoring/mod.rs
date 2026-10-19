/// Scoring matrices and statistical calculations for alignments
pub mod matrices;
pub mod evalue;
pub mod stats;
pub mod bias;

pub use matrices::{ScoreMatrix, ScoringSystem};
pub use evalue::{EValue, BitScore};
pub use stats::{EvalueModel, GapConvention, GapCosts, KarlinParams, SearchStatistics, QueryStatistics};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_module_structure() {
        // Ensure all public items are accessible
        let _matrix_type = ScoringSystem::Blosum62;
        let _evalue = EValue::default();
    }
}
