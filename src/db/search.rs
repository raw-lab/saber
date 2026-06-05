//! Database search wrapper. Composes `SmithWaterman` and `Database` for callers
//! who want a library-level helper rather than the CLI.

use rayon::prelude::*;

use crate::algorithm::{AlignmentResult, AlignmentSet, SmithWaterman};
use crate::db::Database;
use crate::io::FastaRecord;
use crate::scoring::EValue;

pub struct DatabaseSearch {
    pub database: Database,
    pub aligner: SmithWaterman,
    pub ev: EValue,
}

impl DatabaseSearch {
    pub fn new(database: Database, aligner: SmithWaterman, ev: EValue) -> Self {
        Self { database, aligner, ev }
    }

    /// Search one query against every subject. Returns up to `max_alignments`
    /// hits sorted by score descending, filtered by `evalue_threshold`.
    pub fn search(
        &self,
        query: &FastaRecord,
        evalue_threshold: f64,
        max_alignments: usize,
    ) -> AlignmentSet {
        let qbytes = query.sequence.as_bytes();
        let mut hits: Vec<AlignmentResult> = self.database.sequences.records().par_iter().filter_map(|s| {
            let mut a = self.aligner.align(qbytes, s.sequence.as_bytes());
            if a.score <= 0 { return None; }
            a.query_id = query.id.clone();
            a.subject_id = s.id.clone();
            a.bitscore = self.ev.score_to_bitscore(a.score);
            a.evalue = self.ev.score_to_evalue(a.score);
            a.query_coverage = if query.length > 0 { 100.0 * a.align_length as f64 / query.length as f64 } else { 0.0 };
            a.subject_coverage = if s.length > 0 { 100.0 * a.align_length as f64 / s.length as f64 } else { 0.0 };
            if a.evalue > evalue_threshold { return None; }
            Some(a)
        }).collect();
        hits.sort_by(|a, b| b.score.cmp(&a.score));
        hits.truncate(max_alignments);

        AlignmentSet { query_id: query.id.clone(), query_length: query.length, alignments: hits }
    }
}
