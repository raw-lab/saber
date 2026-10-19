//! Legacy library helper kept for API compatibility. Since v1.2 it runs the
//! same engine as the CLI ([`crate::search::api::Searcher`]) in exhaustive
//! mode (full Smith-Waterman against every subject, as before), with SIMD
//! scoring, banded traceback, gapped Karlin-Altschul statistics and
//! consumed-residue coverage.
//!
//! Behaviour changes versus v1.1 (see docs/MIGRATION.md):
//! * E-values/bit scores come from gapped parameters for the aligner's
//!   matrix and gap costs. The `ev` argument is kept for source
//!   compatibility but no longer used; scoring systems without published
//!   gapped parameters produce results flagged `evalue_calibrated = false`.
//! * Coverage counts residues consumed on each sequence (v1.1 used the
//!   alignment column count, which includes gaps).

use crate::algorithm::{AlignmentSet, Scoring, SequenceMode, SmithWaterman};
use crate::db::Database;
use crate::io::FastaRecord;
use crate::scoring::stats::GapCosts;
use crate::scoring::EValue;
use crate::search::api::Searcher;
use crate::search::driver::DriverConfig;

pub struct DatabaseSearch {
    pub database: Database,
    pub aligner: SmithWaterman,
    pub ev: EValue,
    searcher: Result<Searcher, String>,
}

impl DatabaseSearch {
    pub fn new(database: Database, aligner: SmithWaterman, ev: EValue) -> Self {
        let mode = match aligner.scoring {
            Scoring::Matrix(_) => SequenceMode::ProteinToProtein,
            Scoring::NucleotideSimple { .. } => SequenceMode::NucleotideToNucleotide,
        };
        let mut cfg = DriverConfig::new(mode);
        cfg.apply_preset("exhaustive").expect("built-in preset");
        cfg.allow_uncalibrated = true;
        cfg.gaps = GapCosts { open: aligner.gap_open.max(aligner.gap_extend), extend: aligner.gap_extend.max(1) };
        match &aligner.scoring {
            Scoring::Matrix(m) => cfg.matrix = m.system,
            Scoring::NucleotideSimple { match_score, mismatch } => {
                cfg.nt_match = *match_score;
                cfg.nt_mismatch = *mismatch;
            }
        }
        cfg.params.comp_bias = false; // v1.1 semantics: plain substitution scores
        let recs: Vec<(&str, &[u8])> =
            database.sequences.records().iter().map(|r| (r.id.as_str(), r.sequence.as_bytes())).collect();
        let searcher = Searcher::from_records(recs, cfg);
        Self { database, aligner, ev, searcher }
    }

    /// Search one query against every subject. Returns up to `max_alignments`
    /// hits sorted by score descending, filtered by `evalue_threshold`.
    pub fn search(&self, query: &FastaRecord, evalue_threshold: f64, max_alignments: usize) -> AlignmentSet {
        let mut alignments = Vec::new();
        if let Ok(s) = &self.searcher {
            let mut cfg = s.config().clone();
            cfg.params.evalue = evalue_threshold;
            cfg.params.max_target_seqs = max_alignments.max(1);
            // Re-run with per-call thresholds on the shared index.
            if let Ok(eng) = crate::search::driver::make_engine(&cfg, s.index()) {
                let calibrated = eng.stats.params.calibrated;
                let ct = crate::algorithm::translation::CodonTable::new(cfg.genetic_code);
                let rec = crate::io::stream::OwnedRecord {
                    id: query.id.clone(),
                    header: query.description.clone(),
                    seq: query.sequence.as_bytes().to_vec(),
                    ordinal: 0,
                };
                let mut all = Vec::new();
                for (f, codes) in crate::search::driver::query_frames(cfg.mode, &rec, &ct) {
                    let (h, _) = crate::search::driver::search_local(&eng, &codes);
                    all.extend(h.into_iter().map(|h| (f, h)));
                }
                alignments = crate::search::driver::finish_record(&eng, s.index(), &cfg, &rec, all, calibrated);
            }
        }
        alignments.truncate(max_alignments);
        AlignmentSet { query_id: query.id.clone(), query_length: query.length, alignments }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{ScoreMatrix, ScoringSystem};

    #[test]
    fn legacy_database_search_runs_the_engine() {
        let dir = std::env::temp_dir().join(format!("saber_dbsearch_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("db.fa");
        std::fs::write(&fa, ">a\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF\n>b\nWWWWPPPPGGGG\n>c\nMVLSPADKTNVKAAWGKVGAHAGE\n").unwrap();
        let db = Database::from_file(&fa).unwrap();
        let m = std::sync::Arc::new(ScoreMatrix::new(ScoringSystem::Blosum62));
        let ds = DatabaseSearch::new(db, SmithWaterman::protein(11, 1, m), EValue::new(0.3176, 0.134, 1000, 1000));
        let q = FastaRecord {
            id: "q".into(),
            description: "q".into(),
            sequence: "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF".into(),
            length: 47,
        };
        let set = ds.search(&q, 1e-3, 10);
        let ids: Vec<&str> = set.alignments.iter().map(|a| a.subject_id.as_str()).collect();
        assert_eq!(ids, vec!["a", "c"]);
        assert!(set.alignments[0].evalue_calibrated);
        // coverage counts consumed residues
        assert!((set.alignments[1].subject_coverage - 100.0).abs() < 1e-9);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
