//! Partitioned seed index, on-disk format v4.
//!
//! ## What is stored
//!
//! * every subject in the canonical *alignment* encoding (protein codes
//!   0..=23 or nucleotide codes 0..=4) — lossless for the 24-letter protein
//!   alphabet, unlike format v3 which collapsed `B/Z/X/*` to one sentinel;
//! * subject identifiers and 64-bit global offsets;
//! * for every **partition** (≤ 65 535 subjects and a bounded number of
//!   residues) and every **seed shape**: bucket offsets and the postings of
//!   that partition. Offsets are a dense `u32` table (one entry per seed code)
//!   or, when the code space is much larger than the number of distinct seeds
//!   in the partition (e.g. weight-7 Murphy-10 seeds: 10⁷ codes), a sparse
//!   table: an occupancy bitmap with per-64-code ranks plus one `u32` start
//!   per non-empty bucket. Postings are 4 bytes (`local subject | pos<<16`)
//!   when every subject of *that partition* is shorter than 65 536 residues,
//!   otherwise 8 bytes — one long sequence no longer doubles the storage of
//!   the whole database;
//! * the seeding configuration (seed alphabet, shapes, minimizer window,
//!   database-side masking) so a search can never silently use a different
//!   seeding scheme than the one the index was built with.
//!
//! Within a bucket, postings are ordered by (local subject, position), which
//! lets the search restrict a lookup to a subject range by binary search —
//! the basis of the bounded, partitioned seed processing in
//! [`super::engine`].
//!
//! ## File layout
//!
//! ```text
//! [header: 256 bytes][section 0][section 1]...[metadata JSON]
//! ```
//! The header records the file size, the metadata offset/length and CRC-32s
//! of itself and of the metadata. Every section is 8-byte aligned and has a
//! CRC-32 in the metadata. Files are written to a temporary name and renamed
//! into place only after a successful `fsync`, so an interrupted build never
//! leaves a truncated index under the final name; truncated or foreign files
//! are rejected on open with an explanatory error. Format v3 files (SABER ≤
//! 1.1.x) are rejected with migration instructions (rebuild with `--makedb`).
//!
//! ## Bounded-memory construction
//!
//! [`IndexBuilder`] consumes subjects one at a time. Only the subjects of the
//! partition being built are held in memory together with that partition's
//! offset table and postings; finished partitions are either kept in memory
//! (plain in-memory search) or written straight to the output file
//! (`--makedb`, external-memory build). Buckets are counted in place in the
//! final offset table and filled back to front, so no separate count or
//! cursor array is allocated. A partition is limited to fewer than 2³² seed
//! positions (checked), so counts are provably bounded — no saturating
//! arithmetic is involved.

use std::fs::File;
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::crc32::{crc32, Crc32};
use crate::alphabet::{SeedAlphabet, SeedShape, AA_X, NT_N};
use crate::masking::{MaskingMode, Tantan, TantanScratch};
use crate::scoring::matrices::ScoreMatrix;
use crate::scoring::ScoringSystem;

/// Format version written by this SABER.
pub const FORMAT_VERSION: u32 = 4;
const MAGIC: &[u8; 8] = b"SABERIDX";
const HEADER_BYTES: usize = 256;
const ENDIAN_MARK: u32 = 0x0102_0304;
/// Largest number of seed codes (alphabet^weight) per shape.
pub const MAX_CODES: u64 = 1 << 26;
/// Subjects per partition (local ids are 16-bit).
pub const MAX_PARTITION_SUBJECTS: usize = 65_535;

/// Kind of residues stored in the index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SeqKind {
    Protein,
    Nucleotide,
}

/// Seeding configuration (stored in the index).
#[derive(Debug, Clone)]
pub struct SeedConfig {
    pub kind: SeqKind,
    pub alphabet: SeedAlphabet,
    pub shapes: Vec<SeedShape>,
    /// Database-side minimizer window (1 = index every seed position).
    pub minimizer_window: usize,
    /// Database-side masking applied when choosing seed positions.
    pub masking: MaskingMode,
}

impl SeedConfig {
    /// SABER's historical protein default: contiguous 4-mers, full alphabet.
    pub fn protein_default() -> Self {
        Self {
            kind: SeqKind::Protein,
            alphabet: SeedAlphabet::full(),
            shapes: vec![SeedShape::contiguous(4)],
            minimizer_window: 1,
            masking: MaskingMode::Tantan,
        }
    }

    pub fn nucleotide(k: usize) -> Self {
        let mut map = [255u8; 24];
        for (i, m) in map.iter_mut().enumerate().take(4) {
            *m = i as u8;
        }
        Self {
            kind: SeqKind::Nucleotide,
            alphabet: SeedAlphabet { name: "acgt".into(), size: 4, map },
            shapes: vec![SeedShape::contiguous(k)],
            minimizer_window: 1,
            masking: MaskingMode::None,
        }
    }

    /// Number of codes for shape `s`.
    pub fn num_codes(&self, s: usize) -> u64 {
        (self.alphabet.size as u64).pow(self.shapes[s].weight() as u32)
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.shapes.is_empty() {
            return Err("at least one seed shape is required".into());
        }
        for (i, _) in self.shapes.iter().enumerate() {
            let n = self.num_codes(i);
            if n > MAX_CODES {
                return Err(format!(
                    "seed shape '{}' with a {}-letter alphabet has {} codes (> {}); use a smaller weight or a reduced alphabet",
                    self.shapes[i].pattern, self.alphabet.size, n, MAX_CODES
                ));
            }
        }
        if self.minimizer_window == 0 || self.minimizer_window > 64 {
            return Err("minimizer window must be in [1, 64]".into());
        }
        Ok(())
    }

    /// Seed code at window start `p`, or `None` if any care position is not
    /// seedable (ambiguous residue) or masked.
    #[inline]
    pub fn code_at(&self, shape: &SeedShape, seq: &[u8], mask: Option<&[bool]>, p: usize) -> Option<u32> {
        if p + shape.span > seq.len() {
            return None;
        }
        let a = self.alphabet.size as u32;
        let mut c: u32 = 0;
        for &o in &shape.offsets {
            let i = p + o as usize;
            if let Some(m) = mask {
                if m[i] {
                    return None;
                }
            }
            let r = self.alphabet.map[(seq[i] as usize).min(23)];
            if r == 255 {
                return None;
            }
            c = c * a + r as u32;
        }
        Some(c)
    }
}

#[inline]
fn mix32(mut x: u32) -> u32 {
    x ^= x >> 16;
    x = x.wrapping_mul(0x7feb_352d);
    x ^= x >> 15;
    x = x.wrapping_mul(0x846c_a68b);
    x ^= x >> 16;
    x
}

/// Call `f(pos, code)` for every *selected* seed position of `seq`:
/// all valid positions when `window <= 1`, otherwise (w, k)-minimizers over
/// the valid positions (smallest mixed hash, earliest on ties; each selected
/// position reported once).
pub fn for_each_seed<F: FnMut(usize, u32)>(cfg: &SeedConfig, shape: &SeedShape, seq: &[u8], mask: Option<&[bool]>, mut f: F) {
    if seq.len() < shape.span {
        return;
    }
    let last = seq.len() - shape.span;
    let w = cfg.minimizer_window;
    if w <= 1 {
        for p in 0..=last {
            if let Some(c) = cfg.code_at(shape, seq, mask, p) {
                f(p, c);
            }
        }
        return;
    }
    // Monotone deque over window [p-w+1, p] of valid positions.
    let mut dq: std::collections::VecDeque<(usize, u32, u32)> = std::collections::VecDeque::with_capacity(w + 1);
    let mut last_emit: Option<usize> = None;
    for p in 0..=last {
        if let Some(c) = cfg.code_at(shape, seq, mask, p) {
            let r = mix32(c);
            while let Some(&(_, br, _)) = dq.back() {
                if br > r {
                    dq.pop_back();
                } else {
                    break;
                }
            }
            dq.push_back((p, r, c));
        }
        if p + 1 >= w {
            let lo = p + 1 - w;
            while let Some(&(fp, _, _)) = dq.front() {
                if fp < lo {
                    dq.pop_front();
                } else {
                    break;
                }
            }
            if let Some(&(fp, _, fc)) = dq.front() {
                if last_emit != Some(fp) {
                    f(fp, fc);
                    last_emit = Some(fp);
                }
            }
        }
    }
    // Sequences shorter than one window: emit the overall minimizer.
    if last + 1 < w {
        if let Some(&(fp, _, fc)) = dq.front() {
            if last_emit != Some(fp) {
                f(fp, fc);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Metadata
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SectionMeta {
    pub name: String,
    pub offset: u64,
    pub len: u64,
    pub crc32: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PartitionMeta {
    pub first_subject: u64,
    pub num_subjects: u32,
    pub residues: u64,
    pub max_len: u64,
    pub compact: bool,
    /// Postings count per shape.
    pub postings: Vec<u64>,
    /// Per shape: offsets stored sparsely (section `obits.P.S` holds the
    /// occupancy bitmap and ranks, `off.P.S` one start per non-empty bucket).
    /// Absent/false = dense offset table.
    #[serde(default)]
    pub sparse: Vec<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndexMeta {
    pub format_version: u32,
    pub saber_version: String,
    pub kind: SeqKind,
    pub alphabet_name: String,
    pub alphabet_size: u8,
    pub alphabet_map: Vec<u8>,
    pub shapes: Vec<String>,
    pub minimizer_window: u32,
    pub masking: String,
    pub masking_matrix: String,
    pub num_subjects: u64,
    pub total_residues: u64,
    pub max_subject_len: u64,
    pub partitions: Vec<PartitionMeta>,
    pub sections: Vec<SectionMeta>,
    pub source: String,
    /// Per-shape total postings and a high-frequency threshold suggestion.
    pub total_postings: Vec<u64>,
    /// True for a translated nucleotide database (tblastn / `pn` mode):
    /// subject `i` is frame `i % 6` of nucleotide sequence `i / 6`.
    #[serde(default)]
    pub translated: bool,
    #[serde(default)]
    pub genetic_code: String,
    #[serde(default)]
    pub nt_sequences: u64,
    #[serde(default)]
    pub nt_residues: u64,
    /// False for an index without seed postings (residues, identifiers and
    /// the database seed mask only; searched with the scan schedule).
    #[serde(default = "default_true")]
    pub postings: bool,
    /// Seed statistics for the scan schedule (built with the index's seed
    /// alphabet and shapes): per shape, the number of seedable database
    /// positions; the codes occurring more than `scan_freq_min` times are in
    /// sections `scanfreq.<shape>` as `code << 32 | min(count, u32::MAX)`.
    #[serde(default)]
    pub scan_totals: Vec<u64>,
    /// Listing threshold of the `scanfreq` sections (4096 in 1.2 indexes).
    #[serde(default = "default_scan_freq_min")]
    pub scan_freq_min: u64,
}

/// Codes occurring at most this often are not listed in `scanfreq` sections
/// (the automatic over-frequency cap is never below it).
pub const SCAN_FREQ_MIN: u64 = 128;

/// Automatic over-frequency cap of a seed shape: 64x the mean number of
/// database positions per code, at least `OCC_CAP_FLOOR`.
pub const OCC_CAP_FLOOR: u64 = 256;

/// The automatic cap for `total` seed positions over `ncodes` codes.
pub fn auto_occ_cap(total: u64, ncodes: u64) -> u64 {
    ((total as f64 / ncodes.max(1) as f64 * 64.0).ceil() as u64).max(OCC_CAP_FLOOR)
}

fn default_scan_freq_min() -> u64 {
    4096
}

fn default_true() -> bool {
    true
}

impl IndexMeta {
    pub fn seed_config(&self) -> Result<SeedConfig, String> {
        let mut map = [255u8; 24];
        if self.alphabet_map.len() != 24 {
            return Err("index metadata: alphabet map must have 24 entries".into());
        }
        map.copy_from_slice(&self.alphabet_map);
        let shapes = self.shapes.iter().map(|s| SeedShape::parse(s)).collect::<Result<Vec<_>, _>>()?;
        let cfg = SeedConfig {
            kind: self.kind,
            alphabet: SeedAlphabet { name: self.alphabet_name.clone(), size: self.alphabet_size, map },
            shapes,
            minimizer_window: self.minimizer_window as usize,
            masking: self.masking.parse()?,
        };
        cfg.validate()?;
        Ok(cfg)
    }
}

// ---------------------------------------------------------------------------
// Storage
// ---------------------------------------------------------------------------

/// Readable zero bytes kept before and after the residues of an in-memory
/// index (the scan schedule's vector loads may run past a subject's ends).
pub(crate) const RESIDUE_SLACK: usize = 64;

/// 8-byte aligned owned buffer.
#[derive(Default, Clone)]
pub struct AlignedBuf {
    words: Vec<u64>,
    /// Byte offset of the data in `words` (a multiple of 8; zero bytes of
    /// slack before it — see [`AlignedBuf::with_head`]).
    off: usize,
    len: usize,
}

impl AlignedBuf {
    pub fn from_bytes(b: &[u8]) -> Self {
        let mut a = Self::with_len(b.len());
        a.bytes_mut().copy_from_slice(b);
        a
    }
    pub fn with_len(len: usize) -> Self {
        Self { words: vec![0u64; (len + 7) / 8], off: 0, len }
    }
    /// Empty buffer with at least `head` readable zero bytes before its data.
    pub fn with_head(head: usize) -> Self {
        let head = (head + 7) & !7;
        Self { words: vec![0u64; head / 8], off: head, len: 0 }
    }
    /// Make at least `n` zero bytes after the data readable.
    pub fn ensure_tail(&mut self, n: usize) {
        let need = (self.off + self.len + n + 7) / 8;
        if need > self.words.len() {
            self.words.resize(need, 0);
        }
    }
    /// The data with `slack` bytes before and after it, if the buffer has them.
    pub fn padded(&self, slack: usize) -> Option<&[u8]> {
        if self.off >= slack && self.words.len() * 8 >= self.off + self.len + slack {
            Some(&bytemuck::cast_slice::<u64, u8>(&self.words)[self.off - slack..self.off + self.len + slack])
        } else {
            None
        }
    }
    pub fn bytes(&self) -> &[u8] {
        &bytemuck::cast_slice::<u64, u8>(&self.words)[self.off..self.off + self.len]
    }
    pub fn bytes_mut(&mut self) -> &mut [u8] {
        let (off, len) = (self.off, self.len);
        &mut bytemuck::cast_slice_mut::<u64, u8>(&mut self.words)[off..off + len]
    }
    pub fn from_u32(v: Vec<u32>) -> Self {
        Self::from_bytes(bytemuck::cast_slice(&v))
    }
    /// Zeroed buffer holding `n` u32 values (8-byte aligned, no copy later).
    pub fn zeroed_u32(n: usize) -> Self {
        Self::with_len(n * 4)
    }
    /// The first `len / 4` u32 values.
    pub fn as_u32_mut(&mut self) -> &mut [u32] {
        let (o, n) = (self.off / 4, self.len / 4);
        &mut bytemuck::cast_slice_mut::<u64, u32>(&mut self.words)[o..o + n]
    }
    pub fn from_u64(v: Vec<u64>) -> Self {
        let len = v.len() * 8;
        Self { words: v, off: 0, len }
    }
    pub fn capacity_bytes(&self) -> usize {
        self.words.capacity() * 8
    }
    pub fn len(&self) -> usize {
        self.len
    }
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
    /// Append bytes (amortised growth; capacity beyond the written words is
    /// never touched, so it does not become resident).
    pub fn extend_from_slice(&mut self, b: &[u8]) {
        let new_len = self.len + b.len();
        let need = (self.off + new_len + 7) / 8;
        if need > self.words.len() {
            self.words.resize(need, 0);
        }
        let (a, e) = (self.off + self.len, self.off + new_len);
        bytemuck::cast_slice_mut::<u64, u8>(&mut self.words)[a..e].copy_from_slice(b);
        self.len = new_len;
    }
}

enum Storage {
    Owned(Vec<AlignedBuf>),
    Mapped(memmap2::Mmap),
}

/// Postings of one bucket.
#[derive(Clone, Copy)]
pub enum Postings<'a> {
    Compact(&'a [u32]),
    Wide(&'a [u64]),
}

impl<'a> Postings<'a> {
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            Postings::Compact(v) => v.len(),
            Postings::Wide(v) => v.len(),
        }
    }
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// (local subject, position) of element `i`.
    #[inline]
    pub fn get(&self, i: usize) -> (u32, u32) {
        match self {
            Postings::Compact(v) => {
                let p = v[i];
                (p & 0xFFFF, p >> 16)
            }
            Postings::Wide(v) => {
                let p = v[i];
                (p as u32, (p >> 32) as u32)
            }
        }
    }
    /// First element whose local subject is >= `s` (postings are sorted by
    /// subject within a bucket).
    #[inline]
    pub fn lower_bound(&self, s: u32) -> usize {
        match self {
            Postings::Compact(v) => v.partition_point(|&p| (p & 0xFFFF) < s),
            Postings::Wide(v) => v.partition_point(|&p| (p as u32) < s),
        }
    }
}

/// Section handles of one partition.
#[derive(Debug, Clone)]
struct PartSections {
    offsets: Vec<usize>,  // per shape: section index (dense table or sparse starts)
    postings: Vec<usize>, // per shape: section index
    obits: Vec<Option<usize>>, // per shape: sparse bitmap/rank section
}

/// An opened (mapped) or freshly built (owned) v4 index.
pub struct SeedIndex {
    meta: IndexMeta,
    cfg: SeedConfig,
    storage: Storage,
    /// (byte offset, len) of every section within the storage.
    ranges: Vec<(usize, usize)>,
    sec_enc: usize,
    sec_enc_off: usize,
    sec_ids: usize,
    sec_id_off: usize,
    sec_nt_len: Option<usize>,
    /// Database seed-mask runs (global [start, end) pairs), when stored.
    sec_mask: Option<usize>,
    /// Per index shape: section of its frequent-code list (scan statistics).
    sec_scanfreq: Vec<usize>,
    parts: Vec<PartSections>,
    path: Option<PathBuf>,
}

/// Error type for index operations.
#[derive(Debug)]
pub struct IndexError(pub String);
impl std::fmt::Display for IndexError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for IndexError {}
impl From<std::io::Error> for IndexError {
    fn from(e: std::io::Error) -> Self {
        IndexError(e.to_string())
    }
}
fn ierr<T>(m: impl Into<String>) -> Result<T, IndexError> {
    Err(IndexError(m.into()))
}

impl SeedIndex {
    fn section_bytes(&self, s: usize) -> &[u8] {
        let (o, l) = self.ranges[s];
        match &self.storage {
            Storage::Owned(v) => v[s].bytes(),
            Storage::Mapped(m) => &m[o..o + l],
        }
    }

    fn section_u32(&self, s: usize) -> &[u32] {
        bytemuck::cast_slice(self.section_bytes(s))
    }

    fn section_u64(&self, s: usize) -> &[u64] {
        bytemuck::cast_slice(self.section_bytes(s))
    }

    pub fn meta(&self) -> &IndexMeta {
        &self.meta
    }
    pub fn seed_config(&self) -> &SeedConfig {
        &self.cfg
    }
    pub fn kind(&self) -> SeqKind {
        self.cfg.kind
    }
    pub fn num_subjects(&self) -> usize {
        self.meta.num_subjects as usize
    }
    pub fn total_residues(&self) -> u64 {
        self.meta.total_residues
    }
    pub fn num_partitions(&self) -> usize {
        self.meta.partitions.len()
    }
    pub fn partition(&self, p: usize) -> &PartitionMeta {
        &self.meta.partitions[p]
    }
    pub fn path(&self) -> Option<&Path> {
        self.path.as_deref()
    }
    pub fn is_mapped(&self) -> bool {
        matches!(self.storage, Storage::Mapped(_))
    }

    /// True for a translated (six-frame) nucleotide database.
    pub fn is_translated(&self) -> bool {
        self.meta.translated
    }

    /// True when the index stores seed postings (the per-query index
    /// schedules need them; the scan schedule does not).
    pub fn has_postings(&self) -> bool {
        self.meta.postings
    }

    /// Database-side seed mask as global half-open residue intervals
    /// (`[s0, e0, s1, e1, ...]`), when the index stores it (built by SABER
    /// >= 1.3 with masking on). `None` = not stored.
    pub fn mask_runs(&self) -> Option<&[u64]> {
        self.sec_mask.map(|s| self.section_u64(s))
    }

    /// Scan-schedule seed statistics of index shape `si`: (seedable database
    /// positions, codes occurring more than [`SCAN_FREQ_MIN`] times as
    /// `code << 32 | count`), when the index stores them.
    pub fn scan_freq(&self, si: usize) -> Option<(u64, &[u64])> {
        let tot = *self.meta.scan_totals.get(si)?;
        let sec = *self.sec_scanfreq.get(si)?;
        Some((tot, self.section_u64(sec)))
    }

    /// Listing threshold of [`Self::scan_freq`]: caps below it cannot use
    /// the stored statistics.
    pub fn scan_freq_min(&self) -> u64 {
        self.meta.scan_freq_min
    }

    /// Nucleotide length of original sequence `seq` (translated DBs only).
    pub fn nt_len(&self, seq: usize) -> usize {
        match self.sec_nt_len {
            Some(sidx) => self.section_u64(sidx)[seq] as usize,
            None => 0,
        }
    }

    /// All encoded residues (every subject, concatenated). Subject slices
    /// returned by [`Self::subject_codes`] point into this buffer, so reads
    /// slightly before/after a subject stay in bounds when they stay inside it.
    #[inline]
    pub fn residues(&self) -> &[u8] {
        self.section_bytes(self.sec_enc)
    }

    /// [`Self::residues`] with `slack` readable bytes before and after (their
    /// values are unspecified), when the storage has them: a mapped index
    /// (the residues lie between the header and the metadata) or one built
    /// in memory (see [`RESIDUE_SLACK`]).
    pub(crate) fn residues_padded(&self, slack: usize) -> Option<&[u8]> {
        let (o, l) = self.ranges[self.sec_enc];
        match &self.storage {
            Storage::Owned(v) => v[self.sec_enc].padded(slack),
            Storage::Mapped(m) => {
                if o >= slack && o + l + slack <= m.len() {
                    Some(&m[o - slack..o + l + slack])
                } else {
                    None
                }
            }
        }
    }

    /// Start offset of every subject in [`Self::residues`] (`num_subjects + 1`
    /// entries; subject `i` is `residues()[off[i]..off[i + 1]]`).
    #[inline]
    pub fn subject_offsets(&self) -> &[u64] {
        self.section_u64(self.sec_enc_off)
    }

    /// Encoded residues of subject `i` (alignment encoding).
    #[inline]
    pub fn subject_codes(&self, i: usize) -> &[u8] {
        let off = self.section_u64(self.sec_enc_off);
        let (lo, hi) = (off[i] as usize, off[i + 1] as usize);
        &self.section_bytes(self.sec_enc)[lo..hi]
    }

    #[inline]
    pub fn subject_len(&self, i: usize) -> usize {
        let off = self.section_u64(self.sec_enc_off);
        (off[i + 1] - off[i]) as usize
    }

    pub fn subject_id(&self, i: usize) -> &str {
        let off = self.section_u64(self.sec_id_off);
        let b = self.section_bytes(self.sec_ids);
        std::str::from_utf8(&b[off[i] as usize..off[i + 1] as usize]).unwrap_or("?")
    }

    /// Postings range `[lo, hi)` of `code` for shape `shape` in partition `part`.
    #[inline]
    fn bucket_range(&self, part: usize, shape: usize, code: u32) -> (usize, usize) {
        let ps = &self.parts[part];
        let off = self.section_u32(ps.offsets[shape]);
        match ps.obits[shape] {
            None => (off[code as usize] as usize, off[code as usize + 1] as usize),
            Some(sb) => match sparse_rank(self.section_u32(sb), code) {
                Some(r) => (off[r] as usize, off[r + 1] as usize),
                None => (0, 0),
            },
        }
    }

    /// Bucket of `code` for shape `shape` in partition `part`.
    #[inline]
    pub fn bucket(&self, part: usize, shape: usize, code: u32) -> Postings<'_> {
        let (lo, hi) = self.bucket_range(part, shape, code);
        let ps = &self.parts[part];
        if self.meta.partitions[part].compact {
            Postings::Compact(&self.section_u32(ps.postings[shape])[lo..hi])
        } else {
            Postings::Wide(&self.section_u64(ps.postings[shape])[lo..hi])
        }
    }

    /// Size of a bucket without materialising it.
    #[inline]
    pub fn bucket_len(&self, part: usize, shape: usize, code: u32) -> usize {
        let (lo, hi) = self.bucket_range(part, shape, code);
        hi - lo
    }

    /// True if any partition stores shape offsets sparsely.
    pub fn has_sparse_offsets(&self) -> bool {
        self.parts.iter().any(|p| p.obits.iter().any(|o| o.is_some()))
    }

    /// Total bytes of all sections (index size).
    pub fn data_bytes(&self) -> u64 {
        self.ranges.iter().map(|r| r.1 as u64).sum()
    }

    /// Bytes held in RAM by an owned index (0 for a mapped one: mapped pages
    /// are managed by the OS page cache and are not a strict RSS bound).
    pub fn owned_bytes(&self) -> u64 {
        match &self.storage {
            Storage::Owned(v) => v.iter().map(|b| b.capacity_bytes() as u64).sum(),
            Storage::Mapped(_) => 0,
        }
    }

    fn from_parts(meta: IndexMeta, storage: Storage, ranges: Vec<(usize, usize)>, path: Option<PathBuf>) -> Result<Self, IndexError> {
        let cfg = meta.seed_config().map_err(IndexError)?;
        let find = |name: &str| -> Result<usize, IndexError> {
            meta.sections
                .iter()
                .position(|s| s.name == name)
                .ok_or_else(|| IndexError(format!("index is missing section '{name}'")))
        };
        let sec_enc = find("enc")?;
        let sec_enc_off = find("enc_off")?;
        let sec_ids = find("ids")?;
        let sec_id_off = find("id_off")?;
        let sec_nt_len = meta.sections.iter().position(|s| s.name == "nt_len");
        let sec_mask = meta.sections.iter().position(|s| s.name == "maskruns");
        let sec_scanfreq: Vec<usize> = if meta.scan_totals.len() == cfg.shapes.len() {
            (0..cfg.shapes.len()).filter_map(|i| meta.sections.iter().position(|s| s.name == format!("scanfreq.{i}"))).collect()
        } else {
            Vec::new()
        };
        if !meta.scan_totals.is_empty() && sec_scanfreq.len() != meta.scan_totals.len() {
            return ierr("index is missing scan statistics sections");
        }
        if meta.translated && sec_nt_len.is_none() {
            return ierr("translated index is missing its nucleotide-length section");
        }
        let mut parts = Vec::with_capacity(meta.partitions.len());
        for p in 0..meta.partitions.len() {
            let mut offsets = Vec::new();
            let mut postings = Vec::new();
            let mut obits = Vec::new();
            let nshapes = if meta.postings { cfg.shapes.len() } else { 0 };
            for s in 0..nshapes {
                offsets.push(find(&format!("off.{p}.{s}"))?);
                postings.push(find(&format!("post.{p}.{s}"))?);
                let sparse = meta.partitions[p].sparse.get(s).copied().unwrap_or(false);
                obits.push(if sparse { Some(find(&format!("obits.{p}.{s}"))?) } else { None });
            }
            parts.push(PartSections { offsets, postings, obits });
        }
        let idx = Self { meta, cfg, storage, ranges, sec_enc, sec_enc_off, sec_ids, sec_id_off, sec_nt_len, sec_mask, sec_scanfreq, parts, path };
        idx.validate_structure()?;
        Ok(idx)
    }

    /// Structural validation (cheap, always run on open/build).
    fn validate_structure(&self) -> Result<(), IndexError> {
        let n = self.num_subjects();
        let enc_off = self.section_bytes(self.sec_enc_off);
        if enc_off.len() != (n + 1) * 8 {
            return ierr("corrupt index: enc_off section has the wrong length");
        }
        let enc_off = self.section_u64(self.sec_enc_off);
        if enc_off[0] != 0 || enc_off[n] as usize != self.section_bytes(self.sec_enc).len() || enc_off[n] != self.meta.total_residues {
            return ierr("corrupt index: subject offsets do not match the residue section (corrupt or truncated index)");
        }
        if enc_off.windows(2).any(|w| w[1] < w[0]) {
            return ierr("corrupt index: subject offsets are not monotonic (corrupt index)");
        }
        let id_off = self.section_u64(self.sec_id_off);
        if id_off.len() != n + 1 || id_off[n] as usize != self.section_bytes(self.sec_ids).len() {
            return ierr("corrupt index: identifier offsets do not match the identifier section");
        }
        if id_off.windows(2).any(|w| w[1] < w[0]) {
            return ierr("corrupt index: identifier offsets are not monotonic");
        }
        let mut next_subject = 0u64;
        for (pi, pm) in self.meta.partitions.iter().enumerate() {
            if pm.first_subject != next_subject || pm.num_subjects as usize > MAX_PARTITION_SUBJECTS {
                return ierr(format!("corrupt index: partition {pi} has inconsistent subject range"));
            }
            next_subject += pm.num_subjects as u64;
            if pm.compact && pm.max_len >= 65_536 {
                return ierr(format!("corrupt index: partition {pi} is compact but holds a subject of length {}", pm.max_len));
            }
            if !self.meta.postings {
                continue;
            }
            if !pm.sparse.is_empty() && pm.sparse.len() != self.cfg.shapes.len() {
                return ierr(format!("corrupt index: partition {pi} has {} offset-layout flags for {} shapes", pm.sparse.len(), self.cfg.shapes.len()));
            }
            for s in 0..self.cfg.shapes.len() {
                let ncodes = self.cfg.num_codes(s) as usize;
                let off = self.section_u32(self.parts[pi].offsets[s]);
                let entries = match self.parts[pi].obits[s] {
                    None => ncodes + 1,
                    Some(sb) => {
                        if !pm.compact {
                            return ierr(format!("corrupt index: partition {pi} shape {s} has sparse offsets but wide postings"));
                        }
                        let ob = self.section_u32(sb);
                        let nblk = (ncodes + 63) / 64;
                        if ob.len() != nblk * 4 {
                            return ierr(format!("corrupt index: partition {pi} shape {s} sparse bitmap has the wrong size"));
                        }
                        // Ranks must be the running popcount; no bit past the last code.
                        let mut r = 0u64;
                        for blk in 0..nblk {
                            if ob[blk * 4 + 2] as u64 != r || ob[blk * 4 + 3] != 0 {
                                return ierr(format!("corrupt index: partition {pi} shape {s} sparse ranks are inconsistent"));
                            }
                            let word = ob[blk * 4] as u64 | (ob[blk * 4 + 1] as u64) << 32;
                            let valid = (ncodes - blk * 64).min(64);
                            if valid < 64 && word >> valid != 0 {
                                return ierr(format!("corrupt index: partition {pi} shape {s} sparse bitmap marks codes out of range"));
                            }
                            r += word.count_ones() as u64;
                        }
                        r as usize + 1
                    }
                };
                if off.len() != entries || off[0] != 0 {
                    return ierr(format!("corrupt index: partition {pi} shape {s} offset table has the wrong size"));
                }
                let width = if pm.compact { 4 } else { 8 };
                let plen = self.section_bytes(self.parts[pi].postings[s]).len() / width;
                if off[entries - 1] as usize != plen || plen as u64 != pm.postings[s] {
                    return ierr(format!("corrupt index: partition {pi} shape {s} postings length mismatch (truncated?)"));
                }
                // Spot-check monotonicity (full check is O(codes) and cheap).
                if off.windows(2).any(|w| w[1] < w[0]) {
                    return ierr(format!("corrupt index: partition {pi} shape {s} offsets not monotonic"));
                }
            }
        }
        if next_subject != n as u64 {
            return ierr("corrupt index: partitions do not cover all subjects");
        }
        if let Some(ms) = self.sec_mask {
            let runs = self.section_u64(ms);
            if runs.len() % 2 != 0 || runs.chunks(2).any(|r| r[0] >= r[1] || r[1] > self.meta.total_residues) || runs.windows(2).any(|w| w[1] < w[0]) {
                return ierr("corrupt index: seed-mask runs are inconsistent");
            }
        }
        if self.meta.translated {
            let ns = self.meta.nt_sequences as usize;
            if n != ns * 6 {
                return ierr("translated index: subject count is not six frames per sequence");
            }
            if let Some(sidx) = self.sec_nt_len {
                if self.section_bytes(sidx).len() != ns * 8 {
                    return ierr("translated index: nucleotide-length section has the wrong size");
                }
            }
        }
        Ok(())
    }

    /// Verify every section CRC-32 (reads the whole file).
    pub fn verify_checksums(&self) -> Result<(), IndexError> {
        for (i, s) in self.meta.sections.iter().enumerate() {
            let got = crc32(self.section_bytes(i));
            if got != s.crc32 {
                return ierr(format!("index section '{}' checksum mismatch (corrupt file)", s.name));
            }
        }
        Ok(())
    }

    /// Verify that every stored residue is a valid code (reads all residues).
    pub fn verify_alphabet(&self) -> Result<(), IndexError> {
        let max = match self.cfg.kind {
            SeqKind::Protein => 23u8,
            SeqKind::Nucleotide => NT_N,
        };
        if let Some(p) = self.section_bytes(self.sec_enc).iter().position(|&c| c > max) {
            return ierr(format!("index residue section holds invalid code at byte {p}"));
        }
        Ok(())
    }

    /// Open (memory-map) a v4 index file.
    pub fn open(path: &Path) -> Result<Self, IndexError> {
        let file = File::open(path).map_err(|e| IndexError(format!("{}: {e}", path.display())))?;
        let flen = file.metadata()?.len();
        // SAFETY: read-only map of a file we never modify while mapped.
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        if mmap.len() < 16 || &mmap[0..8] != MAGIC {
            return ierr(format!("{}: not a SABER index file", path.display()));
        }
        let rd32 = |o: usize| u32::from_le_bytes(mmap[o..o + 4].try_into().unwrap());
        let rd64 = |o: usize| u64::from_le_bytes(mmap[o..o + 8].try_into().unwrap());
        let version = rd32(8);
        if version != FORMAT_VERSION {
            let hint = match version {
                2 | 3 => " — this index was written by SABER <= 1.1.x (formats v2/v3 store ambiguous residues lossily and \
                          predate seed-scheme metadata). Rebuild it with:  saber --makedb NEW.sdx -d DB.fasta",
                _ => " — rebuild with this SABER version (saber --makedb)",
            };
            return ierr(format!("{}: SABER index format v{version}, this build reads v{FORMAT_VERSION}{hint}", path.display()));
        }
        if mmap.len() < HEADER_BYTES {
            return ierr(format!("{}: truncated index (header incomplete)", path.display()));
        }
        if rd32(12) != ENDIAN_MARK {
            return ierr(format!("{}: index byte order differs from this machine; rebuild it here", path.display()));
        }
        let hdr_crc = rd32(252);
        if crc32(&mmap[0..252]) != hdr_crc {
            return ierr(format!("{}: index header checksum mismatch (corrupt file)", path.display()));
        }
        let file_size = rd64(16);
        if file_size != flen {
            return ierr(format!(
                "{}: index is {} bytes but its header says {} — the file is truncated or was modified; rebuild it",
                path.display(), flen, file_size
            ));
        }
        let meta_off = rd64(24) as usize;
        let meta_len = rd64(32) as usize;
        let meta_crc = rd32(40);
        if meta_off.checked_add(meta_len).map_or(true, |e| e > mmap.len()) {
            return ierr(format!("{}: index metadata out of bounds", path.display()));
        }
        let mbytes = &mmap[meta_off..meta_off + meta_len];
        if crc32(mbytes) != meta_crc {
            return ierr(format!("{}: index metadata checksum mismatch", path.display()));
        }
        let meta: IndexMeta = serde_json::from_slice(mbytes)
            .map_err(|e| IndexError(format!("{}: index metadata unreadable: {e}", path.display())))?;
        let mut ranges = Vec::with_capacity(meta.sections.len());
        for s in &meta.sections {
            let o = s.offset as usize;
            let l = s.len as usize;
            if o % 8 != 0 || o < HEADER_BYTES || o.checked_add(l).map_or(true, |e| e > meta_off) {
                return ierr(format!("{}: index section '{}' out of bounds or misaligned", path.display(), s.name));
            }
            ranges.push((o, l));
        }
        Self::from_parts(meta, Storage::Mapped(mmap), ranges, Some(path.to_path_buf()))
    }
}

// ---------------------------------------------------------------------------
// Builder
// ---------------------------------------------------------------------------

/// Where finished sections go.
enum Sink {
    Memory(Vec<(String, AlignedBuf)>),
    File { f: std::io::BufWriter<File>, tmp: PathBuf, final_path: PathBuf, pos: u64, sections: Vec<SectionMeta> },
}

impl Sink {
    /// Copy a side file into the sink as one section (file sinks only).
    fn push_file(&mut self, name: String, path: &Path, len: u64) -> Result<(), IndexError> {
        match self {
            Sink::Memory(_) => {
                let b = std::fs::read(path)?;
                self.push(name, AlignedBuf::from_bytes(&b))
            }
            Sink::File { f, pos, sections, .. } => {
                let pad = ((*pos + 7) & !7) - *pos;
                f.write_all(&[0u8; 8][..pad as usize])?;
                *pos += pad;
                let mut src = std::io::BufReader::with_capacity(1 << 20, File::open(path)?);
                let mut crc = Crc32::new();
                let mut buf = vec![0u8; 1 << 20];
                let mut copied = 0u64;
                loop {
                    let n = std::io::Read::read(&mut src, &mut buf)?;
                    if n == 0 {
                        break;
                    }
                    crc.update(&buf[..n]);
                    f.write_all(&buf[..n])?;
                    copied += n as u64;
                }
                if copied != len {
                    return ierr(format!("side file {} has {copied} bytes, expected {len}", path.display()));
                }
                sections.push(SectionMeta { name, offset: *pos, len, crc32: crc.finish() });
                *pos += len;
                Ok(())
            }
        }
    }

    fn push(&mut self, name: String, data: AlignedBuf) -> Result<(), IndexError> {
        match self {
            Sink::Memory(v) => v.push((name, data)),
            Sink::File { f, pos, sections, .. } => {
                let pad = ((*pos + 7) & !7) - *pos;
                f.write_all(&[0u8; 8][..pad as usize])?;
                *pos += pad;
                let b = data.bytes();
                f.write_all(b)?;
                sections.push(SectionMeta { name, offset: *pos, len: b.len() as u64, crc32: crc32(b) });
                *pos += b.len() as u64;
            }
        }
        Ok(())
    }
}

/// A growing byte stream kept in memory or spilled to a side file
/// (external-memory build), moved into the index as one section at finish.
enum ByteStream {
    Mem(AlignedBuf),
    File { w: std::io::BufWriter<File>, path: PathBuf, len: u64 },
}

impl ByteStream {
    fn mem() -> Self {
        ByteStream::Mem(AlignedBuf::default())
    }
    fn file(path: PathBuf) -> Result<Self, IndexError> {
        let w = std::io::BufWriter::with_capacity(1 << 20, File::create(&path).map_err(|e| IndexError(format!("{}: {e}", path.display())))?);
        Ok(ByteStream::File { w, path, len: 0 })
    }
    fn write(&mut self, b: &[u8]) -> Result<(), IndexError> {
        match self {
            ByteStream::Mem(v) => v.extend_from_slice(b),
            ByteStream::File { w, len, .. } => {
                w.write_all(b)?;
                *len += b.len() as u64;
            }
        }
        Ok(())
    }
    /// Append a whole buffer; an empty in-memory stream takes it without a copy.
    fn append_owned(&mut self, buf: AlignedBuf) -> Result<(), IndexError> {
        match self {
            // (a buffer with head slack keeps it: the bytes are copied in)
            ByteStream::Mem(v) if v.is_empty() && v.off == 0 => {
                *v = buf;
                Ok(())
            }
            _ => self.write(buf.bytes()),
        }
    }
    fn write_u64(&mut self, x: u64) -> Result<(), IndexError> {
        self.write(&x.to_le_bytes())
    }
    /// Move the stream into the sink as section `name`.
    fn into_section(self, name: &str, sink: &mut Sink) -> Result<(), IndexError> {
        match self {
            ByteStream::Mem(v) => sink.push(name.to_string(), v),
            ByteStream::File { w, path, len } => {
                let mut f = w.into_inner().map_err(|e| IndexError(e.to_string()))?;
                f.flush()?;
                drop(f);
                let r = sink.push_file(name.to_string(), &path, len);
                let _ = std::fs::remove_file(&path);
                r
            }
        }
    }
    fn cleanup(&self) {
        if let ByteStream::File { path, .. } = self {
            let _ = std::fs::remove_file(path);
        }
    }
}

/// How bucket offsets are stored (see the module documentation).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OffsetsLayout {
    /// Sparse when the code space is large and the sparse table is at least
    /// 20 % smaller than the dense one (compact partitions only).
    Auto,
    Dense,
    /// Sparse for every compact partition (tests).
    Sparse,
}

/// Seed-code spaces up to this size always use a dense offset table.
const SPARSE_MIN_CODES: usize = 1 << 20;

/// Options controlling index construction.
#[derive(Debug, Clone)]
pub struct BuildOptions {
    pub cfg: SeedConfig,
    /// Maximum residues per partition (bounds build memory).
    pub max_partition_residues: usize,
    /// Matrix used by the masker (protein only).
    pub masking_matrix: ScoringSystem,
    pub source: String,
    pub offsets_layout: OffsetsLayout,
    /// Build seed postings (needed by the per-query index schedules; the
    /// scan schedule only needs residues and the seed mask).
    pub postings: bool,
    /// Store scan-schedule seed statistics (frequent codes per shape).
    pub scan_stats: bool,
}

impl BuildOptions {
    pub fn new(cfg: SeedConfig) -> Self {
        Self {
            cfg,
            max_partition_residues: 64 << 20,
            masking_matrix: ScoringSystem::Blosum62,
            source: String::new(),
            offsets_layout: OffsetsLayout::Auto,
            postings: true,
            scan_stats: false,
        }
    }
}

/// Maximal runs of `true` as half-open intervals.
fn bool_runs(m: &[bool]) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < m.len() {
        if m[i] {
            let s = i;
            while i < m.len() && m[i] {
                i += 1;
            }
            out.push((s as u32, i as u32));
        } else {
            i += 1;
        }
    }
    out
}

/// Expand masked runs of one subject into `buf` (None when nothing is masked).
fn expand_runs<'a>(runs: &[(u32, u32)], len: usize, buf: &'a mut Vec<bool>) -> Option<&'a [bool]> {
    if runs.is_empty() {
        return None;
    }
    buf.clear();
    buf.resize(len, false);
    for &(a, b) in runs {
        buf[a as usize..b as usize].iter_mut().for_each(|x| *x = true);
    }
    Some(&buf[..])
}

/// Rank of a non-empty code in a sparse offset table: `ob` holds, per block
/// of 64 codes, `[bits_lo, bits_hi, rank, 0]`. `None` = empty bucket.
#[inline]
fn sparse_rank(ob: &[u32], code: u32) -> Option<usize> {
    let b = (code >> 6) as usize * 4;
    let word = ob[b] as u64 | (ob[b + 1] as u64) << 32;
    let k = code & 63;
    if (word >> k) & 1 == 0 {
        return None;
    }
    Some(ob[b + 2] as usize + (word & ((1u64 << k) - 1)).count_ones() as usize)
}

/// Streaming, partition-at-a-time index builder.
pub struct IndexBuilder {
    opt: BuildOptions,
    sink: Sink,
    // whole-database streams (memory or side files)
    enc_all: ByteStream,
    enc_off: ByteStream,
    ids: ByteStream,
    id_off: ByteStream,
    ids_len: u64,
    // current partition
    part_first: u64,
    part_codes: AlignedBuf,
    part_off: Vec<usize>,
    part_max_len: usize,
    partitions: Vec<PartitionMeta>,
    n_subjects: u64,
    total_residues: u64,
    max_len: u64,
    total_postings: Vec<u64>,
    masker: Option<Tantan>,
    /// Residues skipped from seeding by masking (statistics).
    pub masked_residues: u64,
    /// Global seed-mask runs `[start, end)` of all finished partitions.
    mask_runs: Vec<u64>,
    /// Scan statistics: per shape, code counts (u8, with an overflow map for
    /// codes seen 255 times or more) and totals.
    scan_counts: Vec<Vec<u8>>,
    /// Statistics lists computed after the last partition (in-memory
    /// residue stream): one shape at a time instead of every shape's
    /// counters for the whole build.
    scan_lists: Option<Vec<Vec<u64>>>,
    scan_over: Vec<std::collections::HashMap<u32, u64>>,
    scan_totals: Vec<u64>,
    translated: Option<(crate::algorithm::translation::CodonTable, ByteStream)>,
    nt_seqs: u64,
    nt_residues: u64,
    frame_buf: Vec<u8>,
}

impl IndexBuilder {
    /// Builder that keeps everything in memory (for an in-memory search).
    pub fn in_memory(opt: BuildOptions) -> Result<Self, IndexError> {
        Self::new(opt, Sink::Memory(Vec::new()), None)
    }

    /// Builder that writes partitions to `path` as they complete
    /// (external-memory build). The file appears under `path` only after a
    /// successful [`IndexBuilder::finish`].
    pub fn to_file(opt: BuildOptions, path: &Path) -> Result<Self, IndexError> {
        let tmp = path.with_extension(format!("tmp-{}", std::process::id()));
        let mut f = std::io::BufWriter::with_capacity(1 << 20, File::create(&tmp).map_err(|e| IndexError(format!("{}: {e}", tmp.display())))?);
        f.write_all(&[0u8; HEADER_BYTES])?;
        let side = |tag: &str| path.with_extension(format!("{tag}-tmp-{}", std::process::id()));
        let streams = (
            ByteStream::file(side("enc"))?,
            ByteStream::file(side("encoff"))?,
            ByteStream::file(side("ids"))?,
            ByteStream::file(side("idoff"))?,
        );
        Self::new(
            opt,
            Sink::File { f, tmp, final_path: path.to_path_buf(), pos: HEADER_BYTES as u64, sections: Vec::new() },
            Some(streams),
        )
    }

    fn new(
        opt: BuildOptions,
        sink: Sink,
        streams: Option<(ByteStream, ByteStream, ByteStream, ByteStream)>,
    ) -> Result<Self, IndexError> {
        opt.cfg.validate().map_err(IndexError)?;
        let masker = match (opt.cfg.kind, opt.cfg.masking) {
            (SeqKind::Protein, MaskingMode::Tantan) => Some(Tantan::protein(&ScoreMatrix::new(opt.masking_matrix))),
            _ => None,
        };
        let nshapes = opt.cfg.shapes.len();
        let (enc_all, mut enc_off, ids, mut id_off) = streams.unwrap_or_else(|| {
            (ByteStream::Mem(AlignedBuf::with_head(RESIDUE_SLACK)), ByteStream::mem(), ByteStream::mem(), ByteStream::mem())
        });
        enc_off.write_u64(0)?;
        id_off.write_u64(0)?;
        Ok(Self {
            opt,
            sink,
            enc_all,
            enc_off,
            ids,
            id_off,
            ids_len: 0,
            part_first: 0,
            part_codes: AlignedBuf::default(),
            part_off: vec![0],
            part_max_len: 0,
            partitions: Vec::new(),
            n_subjects: 0,
            total_residues: 0,
            max_len: 0,
            total_postings: vec![0; nshapes],
            masker,
            masked_residues: 0,
            mask_runs: Vec::new(),
            scan_counts: Vec::new(),
            scan_lists: None,
            scan_over: Vec::new(),
            scan_totals: Vec::new(),
            translated: None,
            nt_seqs: 0,
            nt_residues: 0,
            frame_buf: Vec::new(),
        })
    }

    /// Turn this builder into a translated (six-frame) database builder.
    /// Must be called before any subject is added.
    pub fn set_translated(&mut self, table: crate::algorithm::translation::TranslationTable, spill: Option<&Path>) -> Result<(), IndexError> {
        if self.n_subjects > 0 {
            return ierr("set_translated must be called before adding subjects");
        }
        if self.opt.cfg.kind != SeqKind::Protein {
            return ierr("translated databases use a protein seed configuration");
        }
        let st = match spill {
            Some(p) => ByteStream::file(p.with_extension(format!("ntlen-tmp-{}", std::process::id())))?,
            None => ByteStream::mem(),
        };
        self.translated = Some((crate::algorithm::translation::CodonTable::new(table), st));
        Ok(())
    }

    /// Add one nucleotide sequence (ASCII) as six translated frames.
    pub fn add_nucleotide_translated(&mut self, id: &str, ascii: &[u8]) -> Result<(), IndexError> {
        let nt = crate::alphabet::encode_nucleotide(ascii);
        let (ct, _) = self.translated.as_ref().ok_or_else(|| IndexError("builder is not translated".into()))?;
        let ct = ct.clone();
        let mut fb = std::mem::take(&mut self.frame_buf);
        for f in 0..6 {
            ct.translate_codes(&nt, f, &mut fb);
            self.add_encoded(id, &fb)?;
        }
        self.frame_buf = fb;
        if let Some((_, st)) = self.translated.as_mut() {
            st.write_u64(nt.len() as u64)?;
        }
        self.nt_seqs += 1;
        self.nt_residues += nt.len() as u64;
        Ok(())
    }

    /// Add one subject given as *encoded* residues (alignment encoding).
    pub fn add_encoded(&mut self, id: &str, codes: &[u8]) -> Result<(), IndexError> {
        let nparts = self.part_off.len() - 1;
        if nparts > 0
            && (nparts >= MAX_PARTITION_SUBJECTS || self.part_codes.len() + codes.len() > self.opt.max_partition_residues)
        {
            self.flush_partition()?;
        }
        self.part_codes.extend_from_slice(codes);
        self.part_off.push(self.part_codes.len());
        self.part_max_len = self.part_max_len.max(codes.len());
        self.ids.write(id.as_bytes())?;
        self.ids_len += id.len() as u64;
        self.id_off.write_u64(self.ids_len)?;
        self.n_subjects += 1;
        self.total_residues = self
            .total_residues
            .checked_add(codes.len() as u64)
            .ok_or_else(|| IndexError("residue count overflowed u64".into()))?;
        self.max_len = self.max_len.max(codes.len() as u64);
        self.enc_off.write_u64(self.total_residues)?;
        Ok(())
    }

    /// Add one subject from ASCII.
    pub fn add_ascii(&mut self, id: &str, ascii: &[u8], scratch: &mut Vec<u8>) -> Result<(), IndexError> {
        scratch.clear();
        match self.opt.cfg.kind {
            SeqKind::Protein => scratch.extend(ascii.iter().map(|&c| crate::alphabet::encode_protein_byte(c))),
            SeqKind::Nucleotide => scratch.extend(ascii.iter().map(|&c| crate::alphabet::encode_nucleotide_byte(c))),
        }
        let codes = std::mem::take(scratch);
        let r = self.add_encoded(id, &codes);
        *scratch = codes;
        r
    }

    fn flush_partition(&mut self) -> Result<(), IndexError> {
        let nsub = self.part_off.len() - 1;
        if nsub == 0 {
            return Ok(());
        }
        let pidx = self.partitions.len();
        let codes_buf = std::mem::take(&mut self.part_codes);
        let off = std::mem::replace(&mut self.part_off, vec![0]);
        let compact = self.part_max_len < 65_536;
        let codes = codes_buf.bytes();
        // Masked runs per subject: computed once, in parallel across subjects,
        // kept as intervals and expanded per subject while seeding.
        let runs: Vec<Vec<(u32, u32)>> = match &self.masker {
            Some(t) => (0..nsub)
                .into_par_iter()
                .map_init(
                    || (TantanScratch::default(), Vec::new()),
                    |(sc, m), i| {
                        t.mask(&codes[off[i]..off[i + 1]], sc, m);
                        bool_runs(m)
                    },
                )
                .collect(),
            None => vec![Vec::new(); nsub],
        };
        self.masked_residues += runs.iter().flatten().map(|&(a, b)| (b - a) as u64).sum::<u64>();
        // Global coordinates: every subject added so far belongs to a finished
        // partition or to this one.
        let part_base = self.total_residues - codes.len() as u64;
        for (i, r) in runs.iter().enumerate() {
            for &(a, b) in r {
                self.mask_runs.push(part_base + off[i] as u64 + a as u64);
                self.mask_runs.push(part_base + off[i] as u64 + b as u64);
            }
        }
        let cfg = self.opt.cfg.clone();
        let shapes_to_index = if self.opt.postings { cfg.shapes.len() } else { 0 };
        let mut mbuf: Vec<bool> = Vec::new();
        let mut sbuf: Vec<(u32, u32)> = Vec::new();
        let mut post_counts = Vec::new();
        let mut sparse_flags = Vec::new();
        for (si, shape) in cfg.shapes.iter().enumerate().take(shapes_to_index) {
            let ncodes = cfg.num_codes(si) as usize;
            let consider_sparse = compact
                && match self.opt.offsets_layout {
                    OffsetsLayout::Auto => ncodes > SPARSE_MIN_CODES,
                    OffsetsLayout::Sparse => true,
                    OffsetsLayout::Dense => false,
                };
            // Offsets: `ends[r]` first holds the cumulative count through
            // bucket r (its end); filling postings back to front (subjects
            // and positions in descending order) decrements it to the
            // bucket's start, so postings end up ordered by (subject,
            // position) without a cursor array. `rank(code)` maps a code to
            // its bucket (identity for a dense table).
            let (mut offsets_buf, obits, total) = if consider_sparse {
                // Count in u8 cells with an overflow map (no 4-byte-per-code
                // table even transiently).
                let mut c8 = vec![0u8; ncodes];
                let mut over: std::collections::HashMap<u32, u64> = std::collections::HashMap::new();
                let mut total = 0u64;
                let mut distinct = 0usize;
                for i in 0..nsub {
                    let s = &codes[off[i]..off[i + 1]];
                    let m = expand_runs(&runs[i], s.len(), &mut mbuf);
                    for_each_seed(&cfg, shape, s, m, |_, c| {
                        total += 1;
                        let x = &mut c8[c as usize];
                        if *x == u8::MAX {
                            *over.entry(c).or_insert(0) += 1;
                        } else {
                            distinct += (*x == 0) as usize;
                            *x += 1;
                        }
                    });
                }
                if total > u32::MAX as u64 {
                    return ierr("partition postings exceed u32; reduce --index-partition-residues");
                }
                let count = |c: usize| -> u64 {
                    let x = c8[c] as u64;
                    if x == u8::MAX as u64 {
                        x + over.get(&(c as u32)).copied().unwrap_or(0)
                    } else {
                        x
                    }
                };
                let nblk = (ncodes + 63) / 64;
                let sparse_bytes = nblk * 16 + (distinct + 1) * 4;
                let dense_bytes = (ncodes + 1) * 4;
                if self.opt.offsets_layout == OffsetsLayout::Sparse || sparse_bytes * 5 < dense_bytes * 4 {
                    let mut ob_buf = AlignedBuf::zeroed_u32(nblk * 4);
                    let mut st_buf = AlignedBuf::zeroed_u32(distinct + 1);
                    {
                        let ob = ob_buf.as_u32_mut();
                        let st = st_buf.as_u32_mut();
                        let (mut r, mut acc) = (0usize, 0u64);
                        for blk in 0..nblk {
                            ob[blk * 4 + 2] = r as u32;
                            let base = blk * 64;
                            let mut word = 0u64;
                            for c in base..(base + 64).min(ncodes) {
                                if c8[c] > 0 {
                                    word |= 1u64 << (c - base);
                                    acc += count(c);
                                    st[r] = acc as u32;
                                    r += 1;
                                }
                            }
                            ob[blk * 4] = word as u32;
                            ob[blk * 4 + 1] = (word >> 32) as u32;
                        }
                        debug_assert_eq!(r, distinct);
                        st[distinct] = acc as u32;
                    }
                    (st_buf, Some(ob_buf), total)
                } else {
                    let mut ob_buf = AlignedBuf::zeroed_u32(ncodes + 1);
                    {
                        let ends = ob_buf.as_u32_mut();
                        let mut acc = 0u64;
                        for (c, e) in ends.iter_mut().enumerate().take(ncodes) {
                            acc += count(c);
                            *e = acc as u32;
                        }
                        ends[ncodes] = acc as u32;
                    }
                    (ob_buf, None, total)
                }
            } else {
                // Dense: count in place in the final table.
                let mut ob_buf = AlignedBuf::zeroed_u32(ncodes + 1);
                let mut total = 0u64;
                {
                    let ends = ob_buf.as_u32_mut();
                    for i in 0..nsub {
                        let s = &codes[off[i]..off[i + 1]];
                        let m = expand_runs(&runs[i], s.len(), &mut mbuf);
                        for_each_seed(&cfg, shape, s, m, |_, c| {
                            total += 1;
                            ends[c as usize] = ends[c as usize].wrapping_add(1);
                        });
                    }
                    if total > u32::MAX as u64 {
                        return ierr("partition postings exceed u32; reduce --index-partition-residues");
                    }
                    let mut acc = 0u32;
                    for e in ends.iter_mut().take(ncodes) {
                        acc += *e;
                        *e = acc;
                    }
                    ends[ncodes] = acc;
                }
                (ob_buf, None, total)
            };
            let total_us = total as usize;
            {
                let ends = offsets_buf.as_u32_mut();
                let ob = obits.as_ref().map(|b| bytemuck::cast_slice::<u8, u32>(b.bytes()));
                let rank = |c: u32| -> usize {
                    match ob {
                        Some(ob) => sparse_rank(ob, c).expect("seed counted in pass 1"),
                        None => c as usize,
                    }
                };
                let mut fill = |put: &mut dyn FnMut(usize, u32, u32)| {
                    for i in (0..nsub).rev() {
                        let s = &codes[off[i]..off[i + 1]];
                        let m = expand_runs(&runs[i], s.len(), &mut mbuf);
                        sbuf.clear();
                        for_each_seed(&cfg, shape, s, m, |p, c| sbuf.push((p as u32, c)));
                        for &(p, c) in sbuf.iter().rev() {
                            let r = rank(c);
                            ends[r] -= 1;
                            put(ends[r] as usize, i as u32, p);
                        }
                    }
                };
                if compact {
                    let mut post_buf = AlignedBuf::zeroed_u32(total_us);
                    {
                        let post = post_buf.as_u32_mut();
                        fill(&mut |slot, li, p| post[slot] = li | (p << 16));
                    }
                    self.sink.push(format!("post.{pidx}.{si}"), post_buf)?;
                } else {
                    let mut post = vec![0u64; total_us];
                    fill(&mut |slot, li, p| post[slot] = li as u64 | ((p as u64) << 32));
                    self.sink.push(format!("post.{pidx}.{si}"), AlignedBuf::from_u64(post))?;
                }
            }
            self.sink.push(format!("off.{pidx}.{si}"), offsets_buf)?;
            let is_sparse = obits.is_some();
            if let Some(ob) = obits {
                self.sink.push(format!("obits.{pidx}.{si}"), ob)?;
            }
            sparse_flags.push(is_sparse);
            post_counts.push(total);
            self.total_postings[si] += total;
        }
        let residues = codes.len() as u64;
        self.partitions.push(PartitionMeta {
            first_subject: self.part_first,
            num_subjects: nsub as u32,
            residues,
            max_len: self.part_max_len as u64,
            compact,
            postings: post_counts,
            sparse: sparse_flags,
        });
        // Residues of a finished partition go to the residue stream (a side
        // file for external-memory builds, so they leave RAM immediately; the
        // first in-memory partition is moved, not copied).
        drop(runs);
        self.enc_all.append_owned(codes_buf)?;
        self.part_first += nsub as u64;
        self.part_max_len = 0;
        Ok(())
    }

    fn meta(&self, sections: Vec<SectionMeta>) -> IndexMeta {
        let cfg = &self.opt.cfg;
        IndexMeta {
            format_version: FORMAT_VERSION,
            saber_version: crate::VERSION.to_string(),
            kind: cfg.kind,
            alphabet_name: cfg.alphabet.name.clone(),
            alphabet_size: cfg.alphabet.size,
            alphabet_map: cfg.alphabet.map.to_vec(),
            shapes: cfg.shapes.iter().map(|s| s.pattern.clone()).collect(),
            minimizer_window: cfg.minimizer_window as u32,
            masking: cfg.masking.to_string(),
            masking_matrix: self.opt.masking_matrix.to_string(),
            num_subjects: self.n_subjects,
            total_residues: self.total_residues,
            max_subject_len: self.max_len,
            partitions: self.partitions.clone(),
            sections,
            source: self.opt.source.clone(),
            total_postings: self.total_postings.clone(),
            translated: self.translated.is_some(),
            genetic_code: self.translated.as_ref().map(|(t, _)| format!("{:?}", t.genetic_code())).unwrap_or_default(),
            nt_sequences: self.nt_seqs,
            nt_residues: self.nt_residues,
            postings: self.opt.postings,
            scan_totals: if self.opt.scan_stats { self.scan_totals.clone() } else { Vec::new() },
            scan_freq_min: SCAN_FREQ_MIN,
        }
    }

    /// Finish an in-memory build.
    pub fn finish_in_memory(mut self) -> Result<SeedIndex, IndexError> {
        self.flush_partition()?;
        self.push_streams()?;
        let bufs = match std::mem::replace(&mut self.sink, Sink::Memory(Vec::new())) {
            Sink::Memory(v) => v,
            Sink::File { .. } => return ierr("finish_in_memory called on a file builder"),
        };
        let mut sections = Vec::new();
        let mut storage = Vec::new();
        let mut ranges = Vec::new();
        for (name, b) in bufs {
            sections.push(SectionMeta { name, offset: 0, len: b.bytes().len() as u64, crc32: 0 });
            ranges.push((0usize, b.bytes().len()));
            storage.push(b);
        }
        let meta = self.meta(sections);
        SeedIndex::from_parts(meta, Storage::Owned(storage), ranges, None)
    }

    /// Finish a file build: append the residue/id sections and metadata,
    /// write the header, fsync and atomically rename into place.
    pub fn finish_to_file(mut self) -> Result<PathBuf, IndexError> {
        self.flush_partition()?;
        self.push_streams()?;
        let sink = std::mem::replace(&mut self.sink, Sink::Memory(Vec::new()));
        let (mut f, tmp, final_path, mut pos, sections) = match sink {
            Sink::File { f, tmp, final_path, pos, sections } => (f, tmp, final_path, pos, sections),
            Sink::Memory(_) => return ierr("finish_to_file called on an in-memory builder"),
        };
        let meta = self.meta(sections);
        let json = serde_json::to_vec(&meta).map_err(|e| IndexError(e.to_string()))?;
        let pad = ((pos + 7) & !7) - pos;
        f.write_all(&[0u8; 8][..pad as usize])?;
        pos += pad;
        let meta_off = pos;
        f.write_all(&json)?;
        pos += json.len() as u64;
        let mut file = f.into_inner().map_err(|e| IndexError(e.to_string()))?;
        let mut hdr = vec![0u8; HEADER_BYTES];
        hdr[0..8].copy_from_slice(MAGIC);
        hdr[8..12].copy_from_slice(&FORMAT_VERSION.to_le_bytes());
        hdr[12..16].copy_from_slice(&ENDIAN_MARK.to_le_bytes());
        hdr[16..24].copy_from_slice(&pos.to_le_bytes());
        hdr[24..32].copy_from_slice(&meta_off.to_le_bytes());
        hdr[32..40].copy_from_slice(&(json.len() as u64).to_le_bytes());
        hdr[40..44].copy_from_slice(&crc32(&json).to_le_bytes());
        let mut c = Crc32::new();
        c.update(&hdr[0..252]);
        hdr[252..256].copy_from_slice(&c.finish().to_le_bytes());
        file.seek(SeekFrom::Start(0))?;
        file.write_all(&hdr)?;
        file.sync_all()?;
        drop(file);
        std::fs::rename(&tmp, &final_path)
            .map_err(|e| IndexError(format!("renaming {} -> {}: {e}", tmp.display(), final_path.display())))?;
        Ok(final_path)
    }

    /// Scan statistics over all residues, after the last partition: per
    /// shape (shapes in parallel, each with its own counters only while it is
    /// counted) the seedable positions and the codes above `SCAN_FREQ_MIN`.
    /// The residue and offset streams of an external-memory build are read
    /// back through a memory map.
    fn scan_stats_post_pass(&mut self) -> Result<(), IndexError> {
        for st in [&mut self.enc_all, &mut self.enc_off] {
            if let ByteStream::File { w, .. } = st {
                w.flush()?;
            }
        }
        let map = |st: &ByteStream| -> Result<Option<memmap2::Mmap>, IndexError> {
            match st {
                ByteStream::File { path, len, .. } if *len > 0 => {
                    let f = File::open(path).map_err(|e| IndexError(format!("{}: {e}", path.display())))?;
                    // SAFETY: a private side file of this build, not modified while mapped.
                    Ok(Some(unsafe { memmap2::Mmap::map(&f) }.map_err(|e| IndexError(format!("{}: {e}", path.display())))?))
                }
                _ => Ok(None),
            }
        };
        let (res_map, off_map) = (map(&self.enc_all)?, map(&self.enc_off)?);
        let res: &[u8] = match (&self.enc_all, &res_map) {
            (ByteStream::Mem(v), _) => v.bytes(),
            (_, Some(m)) => &m[..],
            _ => &[],
        };
        let off_bytes: &[u8] = match (&self.enc_off, &off_map) {
            (ByteStream::Mem(v), _) => v.bytes(),
            (_, Some(m)) => &m[..],
            _ => &[],
        };
        let offs: Vec<usize> = off_bytes.chunks_exact(8).map(|c| u64::from_le_bytes(c.try_into().expect("8 bytes")) as usize).collect();
        let cfg = self.opt.cfg.clone();
        let amap = cfg.alphabet.map;
        let asize = cfg.alphabet.size;
        let nsub = offs.len().saturating_sub(1);
        let runs = &self.mask_runs;
        let offs_ref = &offs;
        let per_shape: Vec<(u64, Vec<u64>)> = cfg
            .shapes
            .par_iter()
            .enumerate()
            .map(|(si, shape)| {
                let ncodes = cfg.num_codes(si) as usize;
                let mut cnt = vec![0u8; ncodes];
                let mut ov: std::collections::HashMap<u32, u64> = std::collections::HashMap::new();
                let inv = u32::MAX;
                let (mut cbuf, mut sb): (Vec<u32>, Vec<u8>) = (Vec::new(), Vec::new());
                let mut tot = 0u64;
                let mut ri = 0usize; // first mask run that may overlap the subject
                for i in 0..nsub {
                    let (b, e) = (offs_ref[i], offs_ref[i + 1]);
                    while ri + 1 < runs.len() && (runs[ri + 1] as usize) <= b {
                        ri += 2;
                    }
                    if e - b < shape.span {
                        continue;
                    }
                    // seed bytes of the subject: seed-alphabet code, or 0xFF
                    // (masked, or not in the seed alphabet), as code_at reads them
                    sb.clear();
                    sb.extend(res[b..e].iter().map(|&c| amap[(c as usize).min(23)]));
                    let mut k = ri;
                    while k + 1 < runs.len() && (runs[k] as usize) < e {
                        let (x, y) = ((runs[k] as usize).max(b), (runs[k + 1] as usize).min(e));
                        if y > x {
                            sb[x - b..y - b].iter_mut().for_each(|d| *d = 0xFF);
                        }
                        k += 2;
                    }
                    sb.resize(e - b + 64, 0xFF);
                    let np = e - b - shape.span + 1;
                    super::engine::seed_codes_into(&sb, np, &shape.offsets, asize, inv, &mut cbuf);
                    for k in 0..np {
                        #[cfg(target_arch = "x86_64")]
                        if let Some(&nc) = cbuf.get(k + 16) {
                            if nc != inv {
                                // SAFETY: prefetch hint inside the counters.
                                unsafe { std::arch::x86_64::_mm_prefetch(cnt.as_ptr().add(nc as usize) as *const i8, std::arch::x86_64::_MM_HINT_T0) };
                            }
                        }
                        let c = cbuf[k];
                        if c == inv {
                            continue;
                        }
                        let x = &mut cnt[c as usize];
                        if *x == u8::MAX {
                            *ov.entry(c).or_insert(0) += 1;
                        } else {
                            *x += 1;
                        }
                        tot += 1;
                    }
                }
                let mut list = Vec::new();
                for (c, &x) in cnt.iter().enumerate() {
                    if (x as u64) <= SCAN_FREQ_MIN && x != u8::MAX {
                        continue;
                    }
                    let n = x as u64 + if x == u8::MAX { ov.get(&(c as u32)).copied().unwrap_or(0) } else { 0 };
                    if n > SCAN_FREQ_MIN {
                        list.push((c as u64) << 32 | n.min(u32::MAX as u64));
                    }
                }
                (tot, list)
            })
            .collect();
        drop((res_map, off_map));
        self.scan_totals = per_shape.iter().map(|x| x.0).collect();
        self.scan_lists = Some(per_shape.into_iter().map(|x| x.1).collect());
        Ok(())
    }

    fn push_streams(&mut self) -> Result<(), IndexError> {
        if self.opt.scan_stats {
            self.scan_stats_post_pass()?;
        }
        let take = |x: &mut ByteStream| std::mem::replace(x, ByteStream::mem());
        let mut enc = take(&mut self.enc_all);
        if let ByteStream::Mem(v) = &mut enc {
            v.ensure_tail(RESIDUE_SLACK);
        }
        enc.into_section("enc", &mut self.sink)?;
        let eo = take(&mut self.enc_off);
        eo.into_section("enc_off", &mut self.sink)?;
        let ids = take(&mut self.ids);
        ids.into_section("ids", &mut self.sink)?;
        let io = take(&mut self.id_off);
        io.into_section("id_off", &mut self.sink)?;
        if self.masker.is_some() {
            let runs = std::mem::take(&mut self.mask_runs);
            self.sink.push("maskruns".to_string(), AlignedBuf::from_u64(runs))?;
        }
        if let Some(lists) = self.scan_lists.take() {
            for (si, list) in lists.into_iter().enumerate() {
                self.sink.push(format!("scanfreq.{si}"), AlignedBuf::from_u64(list))?;
            }
        } else if self.opt.scan_stats {
            if self.scan_counts.is_empty() {
                // empty database
                self.scan_totals = vec![0; self.opt.cfg.shapes.len()];
                self.scan_counts = vec![Vec::new(); self.opt.cfg.shapes.len()];
            }
            let counts = std::mem::take(&mut self.scan_counts);
            let over = std::mem::take(&mut self.scan_over);
            for (si, cnt) in counts.into_iter().enumerate() {
                // codes above the threshold: exact u8 counts, or 255 + the
                // overflow count
                let ov = over.get(si);
                let mut list: Vec<u64> = Vec::new();
                for (c, &x) in cnt.iter().enumerate() {
                    if (x as u64) <= SCAN_FREQ_MIN && x != u8::MAX {
                        continue;
                    }
                    let n = x as u64 + if x == u8::MAX { ov.and_then(|m| m.get(&(c as u32))).copied().unwrap_or(0) } else { 0 };
                    if n > SCAN_FREQ_MIN {
                        list.push((c as u64) << 32 | n.min(u32::MAX as u64));
                    }
                }
                drop(cnt);
                self.sink.push(format!("scanfreq.{si}"), AlignedBuf::from_u64(list))?;
            }
        }
        if let Some((ct, st)) = self.translated.take() {
            st.into_section("nt_len", &mut self.sink)?;
            // keep the flag for metadata
            self.translated = Some((ct, ByteStream::mem()));
        }
        Ok(())
    }

    pub fn subjects_added(&self) -> u64 {
        self.n_subjects
    }
}

impl Drop for IndexBuilder {
    fn drop(&mut self) {
        if let Sink::File { tmp, .. } = &self.sink {
            let _ = std::fs::remove_file(tmp);
        }
        for st in [&self.enc_all, &self.enc_off, &self.ids, &self.id_off] {
            st.cleanup();
        }
    }
}

/// Build an in-memory index from encoded subjects (tests, library use).
pub fn build_in_memory(opt: BuildOptions, subjects: &[(&str, &[u8])]) -> Result<SeedIndex, IndexError> {
    let mut b = IndexBuilder::in_memory(opt)?;
    for (id, s) in subjects {
        b.add_encoded(id, s)?;
    }
    b.finish_in_memory()
}

/// Helper used by tests: replace ambiguous protein codes by X.
pub fn sanitize_protein(codes: &mut [u8]) {
    for c in codes {
        if *c > 23 {
            *c = AA_X;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indexes_without_a_listing_threshold_read_as_1_2_indexes() {
        // 1.2 indexes listed codes above 4096 occurrences and had no field
        let v: serde_json::Value = serde_json::json!({ "scan_freq_min": 1 });
        assert_eq!(v["scan_freq_min"], 1);
        assert_eq!(default_scan_freq_min(), 4096);
        assert!(OCC_CAP_FLOOR > SCAN_FREQ_MIN, "automatic caps must be able to use the stored statistics");
        assert_eq!(auto_occ_cap(38_000_000, 10_000_000), 256);
        assert_eq!(auto_occ_cap(38_000_000, 1_000_000), (38.0f64 * 64.0).ceil() as u64);
    }
    use crate::alphabet::EncodedProtein;

    fn subjects() -> Vec<(String, Vec<u8>)> {
        let seqs: [&[u8]; 5] = [
            b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF",
            b"GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHABZX*",
            b"MVLSPADKTNVKAAWGKVGAHAGEY",
            b"QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ",
            b"WW",
        ];
        seqs.iter().enumerate().map(|(i, s)| (format!("s{i} desc"), EncodedProtein::from_ascii(s).as_slice().to_vec())).collect()
    }

    fn build(cfg: SeedConfig, max_part_res: usize) -> SeedIndex {
        let s = subjects();
        let mut opt = BuildOptions::new(cfg);
        opt.max_partition_residues = max_part_res;
        let v: Vec<(&str, &[u8])> = s.iter().map(|(i, c)| (i.split_whitespace().next().unwrap(), c.as_slice())).collect();
        build_in_memory(opt, &v).unwrap()
    }

    #[test]
    fn residues_round_trip_losslessly_including_bzx_stop() {
        let idx = build(SeedConfig::protein_default(), 1 << 20);
        let s = subjects();
        for (i, (_, c)) in s.iter().enumerate() {
            assert_eq!(idx.subject_codes(i), c.as_slice());
        }
        assert_eq!(idx.subject_id(1), "s1");
        idx.verify_alphabet().unwrap();
    }

    #[test]
    fn masked_homopolymer_is_not_indexed() {
        let idx = build(SeedConfig::protein_default(), 1 << 20);
        let q = EncodedProtein::from_ascii(b"QQQQ");
        let cfg = idx.seed_config().clone();
        let code = cfg.code_at(&cfg.shapes[0], q.as_slice(), None, 0).unwrap();
        let b = idx.bucket(0, 0, code);
        assert!(b.is_empty(), "poly-Q seeds must be masked in the index");
        let mut cfg_nomask = SeedConfig::protein_default();
        cfg_nomask.masking = MaskingMode::None;
        let idx2 = build(cfg_nomask, 1 << 20);
        assert!(idx2.bucket(0, 0, code).len() > 30);
    }

    #[test]
    fn partitioning_preserves_postings_and_order() {
        let mut cfg = SeedConfig::protein_default();
        cfg.masking = MaskingMode::None;
        let one = build(cfg.clone(), 1 << 20);
        let many = build(cfg, 40); // forces several partitions
        assert!(many.num_partitions() > 1);
        let ncodes = one.seed_config().num_codes(0) as u32;
        for code in 0..ncodes {
            let mut a = Vec::new();
            let b0 = one.bucket(0, 0, code);
            for i in 0..b0.len() {
                a.push(b0.get(i));
            }
            let mut b = Vec::new();
            for p in 0..many.num_partitions() {
                let first = many.partition(p).first_subject as u32;
                let bk = many.bucket(p, 0, code);
                for i in 0..bk.len() {
                    let (s, pos) = bk.get(i);
                    b.push((s + first, pos));
                }
            }
            assert_eq!(a, b, "code {code}");
            // sorted by subject within the bucket
            assert!(a.windows(2).all(|w| w[0] <= w[1]));
        }
    }

    fn all_buckets(idx: &SeedIndex, shape: usize) -> Vec<Vec<(u32, u32)>> {
        let ncodes = idx.seed_config().num_codes(shape) as u32;
        (0..ncodes)
            .map(|code| {
                let mut v = Vec::new();
                for p in 0..idx.num_partitions() {
                    let first = idx.partition(p).first_subject as u32;
                    let b = idx.bucket(p, shape, code);
                    assert_eq!(b.len(), idx.bucket_len(p, shape, code));
                    for i in 0..b.len() {
                        let (s, pos) = b.get(i);
                        v.push((s + first, pos));
                    }
                }
                v
            })
            .collect()
    }

    #[test]
    fn sparse_and_dense_offsets_give_identical_buckets() {
        let mut cfg = SeedConfig::protein_default();
        cfg.alphabet = SeedAlphabet::murphy10();
        cfg.shapes = vec![SeedShape::parse("110111").unwrap(), SeedShape::parse("1111").unwrap()];
        let s = subjects();
        let v: Vec<(&str, &[u8])> = s.iter().map(|(i, c)| (i.split_whitespace().next().unwrap(), c.as_slice())).collect();
        for part_res in [1 << 20, 40] {
            let mut dense_opt = BuildOptions::new(cfg.clone());
            dense_opt.max_partition_residues = part_res;
            dense_opt.offsets_layout = OffsetsLayout::Dense;
            let mut sparse_opt = dense_opt.clone();
            sparse_opt.offsets_layout = OffsetsLayout::Sparse;
            let dense = build_in_memory(dense_opt, &v).unwrap();
            let sparse = build_in_memory(sparse_opt.clone(), &v).unwrap();
            assert!(!dense.has_sparse_offsets());
            assert!(sparse.has_sparse_offsets());
            for shape in 0..2 {
                assert_eq!(all_buckets(&dense, shape), all_buckets(&sparse, shape));
            }
            // Sparse layout survives a file round trip.
            let path = std::env::temp_dir().join(format!("saber_sparse_{}_{part_res}.sdx", std::process::id()));
            let mut b = IndexBuilder::to_file(sparse_opt, &path).unwrap();
            for (id, c) in &v {
                b.add_encoded(id, c).unwrap();
            }
            b.finish_to_file().unwrap();
            let f = SeedIndex::open(&path).unwrap();
            f.verify_checksums().unwrap();
            for shape in 0..2 {
                assert_eq!(all_buckets(&dense, shape), all_buckets(&f, shape));
            }
            let _ = std::fs::remove_file(&path);
        }
    }

    #[test]
    fn sparse_offsets_are_validated_on_open() {
        let s = subjects();
        let mut opt = BuildOptions::new(SeedConfig::protein_default());
        opt.offsets_layout = OffsetsLayout::Sparse;
        let path = std::env::temp_dir().join(format!("saber_sparse_bad_{}.sdx", std::process::id()));
        let mut b = IndexBuilder::to_file(opt, &path).unwrap();
        for (id, c) in &s {
            b.add_encoded(id.split_whitespace().next().unwrap(), c).unwrap();
        }
        b.finish_to_file().unwrap();
        let sec = SeedIndex::open(&path).unwrap().meta().sections.iter().find(|x| x.name == "obits.0.0").unwrap().clone();
        let mut bytes = std::fs::read(&path).unwrap();
        // rank of block 1 (must equal the popcount of block 0)
        let o = sec.offset as usize + 16 + 8;
        bytes[o] = bytes[o].wrapping_add(1);
        std::fs::write(&path, &bytes).unwrap();
        let e = SeedIndex::open(&path).err().unwrap().to_string();
        assert!(e.contains("sparse ranks"), "{e}");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn large_code_spaces_use_sparse_offsets_automatically() {
        let mut cfg = SeedConfig::protein_default();
        cfg.alphabet = SeedAlphabet::murphy10();
        cfg.shapes = vec![SeedShape::parse("11101111").unwrap()]; // 10^7 codes
        cfg.masking = MaskingMode::None;
        let s = subjects();
        let v: Vec<(&str, &[u8])> = s.iter().map(|(i, c)| (i.split_whitespace().next().unwrap(), c.as_slice())).collect();
        let auto = build_in_memory(BuildOptions::new(cfg.clone()), &v).unwrap();
        assert!(auto.has_sparse_offsets());
        assert!(auto.owned_bytes() < 4 << 20, "sparse index should be far below the 40 MB dense table");
        let mut d = BuildOptions::new(cfg);
        d.offsets_layout = OffsetsLayout::Dense;
        let dense = build_in_memory(d, &v).unwrap();
        // spot-check the codes present in the subjects plus a few absent ones
        let c = auto.seed_config().clone();
        for (_, codes) in &v {
            for p in 0..codes.len() {
                if let Some(code) = c.code_at(&c.shapes[0], codes, None, p) {
                    for q in [code, code.saturating_sub(1), (code + 1).min(9_999_999)] {
                        let a: Vec<_> = (0..auto.bucket(0, 0, q).len()).map(|i| auto.bucket(0, 0, q).get(i)).collect();
                        let b: Vec<_> = (0..dense.bucket(0, 0, q).len()).map(|i| dense.bucket(0, 0, q).get(i)).collect();
                        assert_eq!(a, b);
                    }
                }
            }
        }
        assert_eq!(auto.bucket_len(0, 0, 9_999_999), dense.bucket_len(0, 0, 9_999_999));
    }

    #[test]
    fn residues_are_readable_with_slack_in_memory_and_mapped() {
        // AlignedBuf slack helpers
        let mut a = AlignedBuf::with_head(RESIDUE_SLACK);
        a.extend_from_slice(b"abc");
        a.extend_from_slice(b"defgh");
        assert_eq!(a.bytes(), b"abcdefgh");
        assert!(a.padded(RESIDUE_SLACK).is_none(), "no tail slack yet");
        a.ensure_tail(RESIDUE_SLACK);
        let p = a.padded(RESIDUE_SLACK).unwrap();
        assert_eq!(p.len(), 8 + 2 * RESIDUE_SLACK);
        assert_eq!(&p[RESIDUE_SLACK..RESIDUE_SLACK + 8], b"abcdefgh");
        assert!(AlignedBuf::from_bytes(b"xy").padded(1).is_none());
        // in-memory indexes (one and several partitions)
        for part in [1 << 20, 40] {
            let idx = build(SeedConfig::protein_default(), part);
            let r = idx.residues();
            let p = idx.residues_padded(32).expect("in-memory residues have slack");
            assert_eq!(&p[32..32 + r.len()], r);
            assert_eq!(p.len(), r.len() + 64);
        }
        // mapped index
        let s = subjects();
        let path = std::env::temp_dir().join(format!("saber_slack_{}.sdx", std::process::id()));
        let mut b = IndexBuilder::to_file(BuildOptions::new(SeedConfig::protein_default()), &path).unwrap();
        for (id, c) in &s {
            b.add_encoded(id.split_whitespace().next().unwrap(), c).unwrap();
        }
        b.finish_to_file().unwrap();
        let m = SeedIndex::open(&path).unwrap();
        let r = m.residues();
        let p = m.residues_padded(32).expect("mapped residues lie between header and metadata");
        assert_eq!(&p[32..32 + r.len()], r);
        drop(m);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn file_round_trip_and_rejections() {
        let s = subjects();
        let dir = std::env::temp_dir();
        let path = dir.join(format!("saber_v4_{}.sdx", std::process::id()));
        let mut opt = BuildOptions::new(SeedConfig::protein_default());
        opt.max_partition_residues = 50;
        let mut b = IndexBuilder::to_file(opt.clone(), &path).unwrap();
        for (id, c) in &s {
            b.add_encoded(id.split_whitespace().next().unwrap(), c).unwrap();
        }
        b.finish_to_file().unwrap();
        let m = SeedIndex::open(&path).unwrap();
        m.verify_checksums().unwrap();
        let v: Vec<(&str, &[u8])> = s.iter().map(|(i, c)| (i.split_whitespace().next().unwrap(), c.as_slice())).collect();
        let mem = build_in_memory(opt, &v).unwrap();
        assert_eq!(m.num_partitions(), mem.num_partitions());
        for i in 0..s.len() {
            assert_eq!(m.subject_codes(i), mem.subject_codes(i));
            assert_eq!(m.subject_id(i), mem.subject_id(i));
        }
        // Truncation is detected.
        let bytes = std::fs::read(&path).unwrap();
        let tpath = dir.join(format!("saber_v4_trunc_{}.sdx", std::process::id()));
        std::fs::write(&tpath, &bytes[..bytes.len() - 10]).unwrap();
        let e = SeedIndex::open(&tpath).err().unwrap().to_string();
        assert!(e.contains("truncated"), "{e}");
        // Old format versions are rejected with migration advice.
        let mut v3 = bytes.clone();
        v3[8..12].copy_from_slice(&3u32.to_le_bytes());
        std::fs::write(&tpath, &v3).unwrap();
        let e = SeedIndex::open(&tpath).err().unwrap().to_string();
        assert!(e.contains("Rebuild") || e.contains("rebuild"), "{e}");
        // Flipped byte in metadata is detected.
        let mut bad = bytes.clone();
        let n = bad.len();
        bad[n - 5] ^= 0x55;
        std::fs::write(&tpath, &bad).unwrap();
        assert!(SeedIndex::open(&tpath).is_err());
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&tpath);
    }

    #[test]
    fn spaced_and_reduced_seeds() {
        let mut cfg = SeedConfig::protein_default();
        cfg.alphabet = SeedAlphabet::murphy10();
        cfg.shapes = vec![SeedShape::parse("11011").unwrap(), SeedShape::parse("1111").unwrap()];
        cfg.masking = MaskingMode::None;
        let idx = build(cfg, 1 << 20);
        // K and R are in the same murphy10 group: "AKVL" and "ARVL" share codes.
        let a = EncodedProtein::from_ascii(b"AKVL");
        let b = EncodedProtein::from_ascii(b"ARVL");
        let c = idx.seed_config();
        assert_eq!(c.code_at(&c.shapes[1], a.as_slice(), None, 0), c.code_at(&c.shapes[1], b.as_slice(), None, 0));
        assert_eq!(idx.meta().shapes.len(), 2);
    }

    #[test]
    fn minimizers_are_subset_of_dense_positions() {
        let cfg = SeedConfig { minimizer_window: 4, masking: MaskingMode::None, ..SeedConfig::protein_default() };
        let dense = SeedConfig { minimizer_window: 1, ..cfg.clone() };
        let s = EncodedProtein::from_ascii(b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKK");
        let mut all = Vec::new();
        for_each_seed(&dense, &dense.shapes[0], s.as_slice(), None, |p, c| all.push((p, c)));
        let mut sel = Vec::new();
        for_each_seed(&cfg, &cfg.shapes[0], s.as_slice(), None, |p, c| sel.push((p, c)));
        assert!(sel.len() < all.len() && sel.len() * 3 > all.len() / 2);
        for x in &sel {
            assert!(all.contains(x));
        }
        // every window of 4 consecutive positions contains a selected one
        let selpos: Vec<usize> = sel.iter().map(|x| x.0).collect();
        for start in 0..all.len() - 3 {
            assert!((start..start + 4).any(|p| selpos.contains(&p)), "window {start}");
        }
    }
}
