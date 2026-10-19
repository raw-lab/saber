//! Library API. [`Searcher`] runs exactly the same engine, statistics and
//! result conversion as the `saber` CLI (it is built on the same driver
//! functions), so library and command-line results are identical for the
//! same [`DriverConfig`].
//!
//! ```no_run
//! use saber_rs::algorithm::SequenceMode;
//! use saber_rs::search::api::Searcher;
//! use saber_rs::search::driver::DriverConfig;
//!
//! let mut cfg = DriverConfig::new(SequenceMode::ProteinToProtein);
//! cfg.params.evalue = 1e-3;
//! let s = Searcher::from_records([("subj1", b"MKTAYIAKQRQISFVKSHFSRQ".as_slice())], cfg).unwrap();
//! let hits = s.search("q1", b"MKTAYIAKQRQISFVKSHF").unwrap();
//! # let _ = hits;
//! ```

use std::path::Path;

use rayon::prelude::*;

use crate::algorithm::translation::CodonTable;
use crate::algorithm::{AlignmentResult, SequenceMode};
use crate::io::stream::OwnedRecord;

use super::driver::{self, DriverConfig};
use super::engine::{Counters, Hit};
use super::index::{BuildOptions, IndexBuilder, SeedIndex, SeqKind};

/// A reusable search handle: one index plus one configuration.
pub struct Searcher {
    index: SeedIndex,
    cfg: DriverConfig,
}

impl Searcher {
    /// Build an in-memory index from `(id, sequence)` records (ASCII
    /// residues; nucleotide for `nn`/`pn`, protein otherwise).
    pub fn from_records<'a, I>(records: I, cfg: DriverConfig) -> Result<Self, String>
    where
        I: IntoIterator<Item = (&'a str, &'a [u8])>,
    {
        let mut seed = cfg.seed.clone();
        seed.kind = if cfg.mode == SequenceMode::NucleotideToNucleotide { SeqKind::Nucleotide } else { SeqKind::Protein };
        let mut opt = BuildOptions::new(seed);
        opt.max_partition_residues = cfg.partition_residues;
        let mut b = IndexBuilder::in_memory(opt).map_err(|e| e.to_string())?;
        if cfg.mode == SequenceMode::ProteinToNucleotide {
            b.set_translated(cfg.genetic_code, None).map_err(|e| e.to_string())?;
        }
        let mut scratch = Vec::new();
        for (id, seq) in records {
            if cfg.mode == SequenceMode::ProteinToNucleotide {
                b.add_nucleotide_translated(id, seq).map_err(|e| e.to_string())?;
            } else {
                b.add_ascii(id, seq, &mut scratch).map_err(|e| e.to_string())?;
            }
        }
        let index = b.finish_in_memory().map_err(|e| e.to_string())?;
        Ok(Self { index, cfg })
    }

    /// Build an in-memory index from a FASTA file (streamed).
    pub fn from_fasta(path: &Path, cfg: DriverConfig) -> Result<Self, String> {
        let (index, _) = driver::build_index(&cfg, path, None)?;
        Ok(Self { index, cfg })
    }

    /// Open a persistent index written by `saber --makedb` (memory-mapped).
    pub fn open(path: &Path, cfg: DriverConfig) -> Result<Self, String> {
        let index = SeedIndex::open(path).map_err(|e| e.to_string())?;
        if cfg.verify_index {
            index.verify_checksums().map_err(|e| e.to_string())?;
            index.verify_alphabet().map_err(|e| e.to_string())?;
        }
        Ok(Self { index, cfg })
    }

    pub fn index(&self) -> &SeedIndex {
        &self.index
    }

    pub fn config(&self) -> &DriverConfig {
        &self.cfg
    }

    /// Search one query (ASCII residues).
    pub fn search(&self, id: &str, seq: &[u8]) -> Result<Vec<AlignmentResult>, String> {
        Ok(self.search_many(&[(id.to_string(), seq.to_vec())])?.0.pop().unwrap_or_default())
    }

    /// Search many queries (parallel over queries on the current rayon pool).
    /// Returns per-query results in input order plus the summed counters.
    pub fn search_many(&self, queries: &[(String, Vec<u8>)]) -> Result<(Vec<Vec<AlignmentResult>>, Counters), String> {
        let eng = driver::make_engine(&self.cfg, &self.index)?;
        let calibrated = eng.stats.params.calibrated;
        let ct = CodonTable::new(self.cfg.genetic_code);
        let recs: Vec<OwnedRecord> = queries
            .iter()
            .enumerate()
            .map(|(i, (id, s))| OwnedRecord { id: id.clone(), header: id.clone(), seq: s.clone(), ordinal: i as u64 })
            .collect();
        let out: Vec<(Vec<AlignmentResult>, Counters)> = recs
            .par_iter()
            .map(|r| {
                let mut all: Vec<(Option<usize>, Hit)> = Vec::new();
                let mut c = Counters::default();
                for (f, codes) in driver::query_frames(self.cfg.mode, r, &ct) {
                    let (h, d) = driver::search_local(&eng, &codes);
                    c.add(&d);
                    all.extend(h.into_iter().map(|h| (f, h)));
                }
                (driver::finish_record(&eng, &self.index, &self.cfg, r, all, calibrated), c)
            })
            .collect();
        let mut total = Counters::default();
        let mut res = Vec::with_capacity(out.len());
        for (r, c) in out {
            total.add(&c);
            res.push(r);
        }
        Ok((res, total))
    }
}
