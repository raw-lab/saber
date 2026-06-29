//! Inverted k-mer index over a protein or nucleotide database.
//!
//! This is the "what BLAST does" piece — built once per database, queried
//! per query in O(query_length × avg_postings_per_kmer). Replaces the
//! linear-scan prefilter with O(1) k-mer lookup.
//!
//! ## Data layout
//!
//! Flat-array postings list (BLAST/MMseqs2-style):
//!
//! ```text
//! offsets:        [hash → start in postings]   len = α^k + 1
//! postings_subj:  [subject_idx of each hit]    len = total_kmer_occurrences
//! postings_pos:   [position in subject]        len = total_kmer_occurrences
//! ```
//!
//! Looking up k-mer with hash `h`:
//! ```text
//! for i in offsets[h] .. offsets[h+1] {
//!     postings_subj[i], postings_pos[i]   ← this hit
//! }
//! ```
//!
//! Memory: 2 × 4 bytes per k-mer occurrence + 4 bytes × α^k for the
//! offsets table. For a 5003-seq protein DB with k=4, that's ~1.2 MB of
//! postings + 1.3 MB of offsets = ~2.5 MB total. For 100M residues with
//! k=4, ~800 MB — fits comfortably in RAM on any workstation.
//!
//! ## Why exact k-mers, not neighborhood words?
//!
//! BLAST uses "neighborhood words" — for each query k-mer, generate all
//! k-mers within score T of it and look up *all* of them. This catches
//! near-matches at the seed stage and is what gives BLAST its sensitivity.
//! We don't do that yet — exact k-mer match only. This means SABER's
//! seed-stage recall will be lower than BLAST's at very low identity
//! (~25%); the SIMD SW that follows will recover most of that, but if
//! recall ends up being a real issue, the next step is to add
//! neighborhood-word generation (~1 day of work).

use rayon::prelude::*;

/// Default k for protein indexing.
/// k=4 catches every significant homolog (E≤1) on the bundled benchmark
/// while keeping the candidate set small (~1000 of 5000 subjects on a
/// 5003-seq DB). k=3 has BLAST-like sensitivity but generates far more
/// candidates (≥95% of subjects) which negates the speedup. Use --index-k 3
/// for extra-sensitive remote-homology search.
pub const DEFAULT_K_PROTEIN: usize = 4;
/// Default k for nucleotide indexing — matches BLAST's blastn default.
pub const DEFAULT_K_NUCLEOTIDE: usize = 11;

/// Default minimizer window for protein subject indexing (large DBs only;
/// small DBs use dense indexing — see the adaptive logic in main). 1 = dense.
/// Tuned empirically on the UniProt bio benchmark (see CHANGELOG 2.10.0):
/// window 4 brings total peak RSS to 51 MB (below DIAMOND's 52 MB) at the
/// default preset while keeping recall at 691 hits — above BLAST (687) and
/// well above DIAMOND --sensitive (667). The query side stays dense, so the
/// sensitivity cost of sampling the subject side is minimal.
pub const DEFAULT_MINIMIZER_WINDOW_PROTEIN: usize = 4;

/// Protein alphabet size used here. Residues 0..19 are the 20 standard amino
/// acids; we exclude B/Z/X/* (indices 20..23) from k-mer generation since
/// ambiguous residues introduce spurious seed hits.
const PROTEIN_ALPHABET_FOR_INDEX: usize = 20;

/// Nucleotide alphabet size for indexing: ACGT only. N is excluded.
const NUCLEOTIDE_ALPHABET_FOR_INDEX: usize = 4;

/// One seed hit: a k-mer in the query that also appears in some subject.
#[derive(Debug, Clone, Copy)]
pub struct SeedHit {
    pub subject_idx: u32,
    /// Diagonal in the DP matrix = subject_pos − query_pos. Hits on the
    /// same diagonal are very likely part of the same alignment.
    pub diagonal: i32,
    pub q_pos: u32,
    pub s_pos: u32,
}

/// A scored candidate subject (output of [`KmerIndex::candidates`]).
#[derive(Debug, Clone, Copy)]
pub struct Candidate {
    pub subject_idx: u32,
    /// Number of distinct seed hits. Used as the prefilter score.
    pub hit_count: u32,
    /// Best diagonal estimate — the diagonal that had the most hits.
    /// Used by banded SW to restrict the DP to a strip around this diagonal.
    pub best_diagonal: i32,
    /// A representative `q_pos` of a seed on the best diagonal. Used as
    /// the anchor for ungapped X-drop extension. Picked as the median
    /// query position among the seeds on the best diagonal.
    pub seed_q_pos: u32,
}

/// Posting storage. Each posting records (subject_idx, position_in_subject).
///
/// In **Compact** mode (the common case: ≤ 65535 subjects AND every subject
/// shorter than 65535 residues) we pack both into a single u32 — 4 bytes per
/// posting. For the bundled 20K-subject UniProt benchmark this cuts the
/// postings buffer from 72 MB to 36 MB.
///
/// In **Wide** mode (huge DBs or very long subjects) we use u64s, exactly
/// the same 8 bytes per posting we used to have with two parallel u32 vecs.
///
/// The mode is chosen at build time. Iteration and lookup both dispatch on
/// the variant; the match is at a coarse granularity so branch prediction
/// gets it right consistently.
enum PostingsBuf {
    /// Bits 0..16 = subject_idx, bits 16..32 = position. 4 bytes per posting.
    Compact(Vec<u32>),
    /// Bits 0..32 = subject_idx, bits 32..64 = position. 8 bytes per posting.
    Wide(Vec<u64>),
}

impl PostingsBuf {
    /// Decode the posting at index i into (subject_idx, position).
    ///
    /// Reserved for callers that want a backing-agnostic single-element
    /// decode. The query hot path does NOT use this — it fetches the whole
    /// postings slice once via [`KmerIndex::postings_u32`] /
    /// [`KmerIndex::postings_u64`] and inlines the variant-specific bit
    /// twiddling for cache locality (see `KmerIndex::query_encoded`).
    /// Length / byte-size / compactness queries also live on `KmerIndex`,
    /// since they must work uniformly over both owned and memory-mapped
    /// backings.
    #[allow(dead_code)]
    #[inline]
    fn get(&self, i: usize) -> (u32, u32) {
        match self {
            PostingsBuf::Compact(v) => {
                let p = v[i];
                (p & 0xFFFF, p >> 16)
            }
            PostingsBuf::Wide(v) => {
                let p = v[i];
                (p as u32, (p >> 32) as u32)
            }
        }
    }
}

/// Backing storage for the index's large arrays. Either fully **owned**
/// (the in-memory build path) or **memory-mapped** from a persistent index
/// file (the `--db` path). The mmap variant lets the OS page index data in
/// on demand, so peak RSS includes only touched pages — this is what makes
/// UniRef90-scale databases workable on a workstation, and it skips the
/// rebuild entirely on repeated runs.
enum IndexData {
    Owned {
        /// Cumulative postings offsets, one per k-mer bucket + 1. Values reach
        /// `num_postings` (≈ total residues), so this is u64: a >4 GiB database
        /// would overflow u32 here.
        offsets: Vec<u64>,
        postings: PostingsBuf,
        enc_buf: Vec<u8>,
        /// Cumulative byte offsets into `enc_buf`, one per subject + 1. Values
        /// reach `total_residues`, so u64 — the field that used to overflow on
        /// multi-gigabyte databases.
        enc_off: Vec<u64>,
        /// Per-subject length (residues). Length = num_subjects. A single
        /// sequence is always < 4 GiB, so u32 is sufficient here.
        subj_lens: Vec<u32>,
        /// Concatenated subject ID strings.
        id_blob: Vec<u8>,
        /// id_offsets[i]..id_offsets[i+1] is subject i's ID in id_blob. u64
        /// because the concatenated ID blob can exceed 4 GiB on huge databases.
        id_offsets: Vec<u64>,
    },
    Mapped {
        mmap: memmap2::Mmap,
        postings_compact: bool,
        offsets: ByteRange,
        postings: ByteRange,
        enc_buf: ByteRange,
        enc_off: ByteRange,
        subj_lens: ByteRange,
        id_blob: ByteRange,
        id_offsets: ByteRange,
    },
}

/// (start, len) byte range within the mmap.
#[derive(Clone, Copy)]
struct ByteRange {
    start: usize,
    len: usize,
}

/// Inverted k-mer index. Build with [`KmerIndex::build_protein`] /
/// [`KmerIndex::build_protein_windowed`] (in-memory), or load a persistent
/// one with [`KmerIndex::open`] (memory-mapped). Query with
/// [`KmerIndex::query_protein`] / [`KmerIndex::query_nucleotide`].
pub struct KmerIndex {
    k: usize,
    alphabet_size: usize,
    num_subjects: usize,
    total_residues: usize,
    data: IndexData,
}

impl KmerIndex {
    /// Build an inverted k-mer index over `subjects` (raw ASCII sequences).
    /// Use [`DEFAULT_K_PROTEIN`] (= 4) for typical protein databases.
    ///
    /// Dense indexing (every k-mer position). Equivalent to
    /// `build_protein_windowed(subjects, k, 1)`.
    pub fn build_protein(subjects: &[&[u8]], k: usize) -> Self {
        Self::build(subjects, k, PROTEIN_ALPHABET_FOR_INDEX, encode_protein_byte, 1)
    }

    pub fn build_nucleotide(subjects: &[&[u8]], k: usize) -> Self {
        Self::build(subjects, k, NUCLEOTIDE_ALPHABET_FOR_INDEX, encode_nucleotide_byte, 1)
    }

    /// Build a **sparse minimizer index** over protein `subjects`. Only the
    /// minimizer k-mer of each `window`-length window is indexed, cutting the
    /// postings array (and thus peak memory) by roughly `2/(window+1)`.
    /// `window = 1` is identical to [`Self::build_protein`] (dense).
    ///
    /// Sensitivity is preserved on the query side: the search still looks up
    /// *every* query k-mer (and optional neighborhood), so seeds are missed
    /// only when a conserved k-mer was never selected as a subject minimizer.
    /// The two-hit / diagonal-clustering stage tolerates this because
    /// conserved regions span several windows and contribute several
    /// minimizers.
    pub fn build_protein_windowed(subjects: &[&[u8]], k: usize, window: usize) -> Self {
        Self::build(subjects, k, PROTEIN_ALPHABET_FOR_INDEX, encode_protein_byte, window.max(1))
    }

    pub fn build_nucleotide_windowed(subjects: &[&[u8]], k: usize, window: usize) -> Self {
        Self::build(subjects, k, NUCLEOTIDE_ALPHABET_FOR_INDEX, encode_nucleotide_byte, window.max(1))
    }

    fn build<F>(subjects: &[&[u8]], k: usize, alphabet_size: usize, encode_byte: F, window: usize) -> Self
    where
        F: Fn(u8) -> Option<u8> + Sync + Send,
    {
        // For protein (alphabet=20) k>6 already exhausts ~64M offsets; for
        // nucleotide (alphabet=4) k up to 13 is fine (~64M offsets).
        assert!(k >= 2 && k <= 13, "k must be in [2, 13]; got {}", k);
        assert!(
            (alphabet_size as u64).checked_pow(k as u32).map_or(false, |n| n < (u32::MAX as u64) / 2),
            "k={} too large for alphabet_size={} (offsets would exceed u32)",
            k, alphabet_size
        );
        let num_kmers = (alphabet_size as u64).pow(k as u32) as usize;

        // ---- Build the flat encoded buffer FIRST. ----
        // We used to do this:
        //   1. encode every subject into a Vec<Vec<u8>> `encoded` (9 MB)
        //   2. count k-mers from `encoded`
        //   3. allocate postings
        //   4. fill postings from `encoded`
        //   5. flatten `encoded` into enc_buf (another 9 MB)
        //   6. drop encoded
        // — that peaked with `encoded` + `enc_buf` simultaneously alive
        // (18 MB of redundant sequence storage). The new flow encodes
        // directly into `enc_buf` and uses `enc_off` boundaries for all
        // later passes, eliminating the intermediate Vec<Vec<u8>>.
        let total_residues: usize = subjects.iter().map(|s| s.len()).sum();
        let mut enc_buf: Vec<u8> = Vec::with_capacity(total_residues);
        let mut enc_off: Vec<u64> = Vec::with_capacity(subjects.len() + 1);
        enc_off.push(0u64);
        let mut acc: u64 = 0;
        let mut max_subj_len: usize = 0;
        for s in subjects {
            // Encode in-place: each ASCII byte → encoded byte (0..alphabet-1)
            // or 255 for ambiguous. Direct push into enc_buf avoids the
            // intermediate Vec<u8>.
            for &c in *s {
                enc_buf.push(encode_byte(c).unwrap_or(255));
            }
            // u64 accumulator: a multi-gigabyte database used to overflow the
            // old u32 here ("encoded subject offsets overflowed u32"). u64 holds
            // any realistic total residue count.
            acc = acc.checked_add(s.len() as u64)
                .expect("encoded subject offsets overflowed u64 (impossibly large database)");
            enc_off.push(acc);
            if s.len() > max_subj_len { max_subj_len = s.len(); }
        }

        // Helper: subject i lives at enc_buf[enc_off[i]..enc_off[i+1]].
        let subj_slice = |i: usize| -> &[u8] {
            let lo = enc_off[i] as usize;
            let hi = enc_off[i + 1] as usize;
            &enc_buf[lo..hi]
        };

        // ---- Pass 1: count k-mers per hash bucket. ----
        // Parallel reduction over subject indices, accumulating into per-thread
        // count arrays then summing. For small DBs this is overkill but it
        // costs nothing when the DB is big.
        let counts: Vec<u32> = (0..subjects.len())
            .into_par_iter()
            .fold(
                || vec![0u32; num_kmers],
                |mut acc, i| {
                    count_kmers_into(&mut acc, subj_slice(i), k, alphabet_size, window);
                    acc
                },
            )
            .reduce(
                || vec![0u32; num_kmers],
                |mut a, b| {
                    for i in 0..num_kmers {
                        a[i] = a[i].saturating_add(b[i]);
                    }
                    a
                },
            );

        // ---- Build offsets via prefix sum. ----
        // u64 offsets: `sum` reaches the total number of k-mer occurrences,
        // which for a multi-gigabyte database exceeds u32::MAX. The old u32
        // offsets asserted-and-aborted here ("index too large for u32
        // postings"); u64 scales to arbitrarily large databases.
        let mut offsets: Vec<u64> = Vec::with_capacity(num_kmers + 1);
        offsets.push(0u64);
        let mut sum: u64 = 0;
        for &c in &counts {
            sum += c as u64;
            offsets.push(sum);
        }
        drop(counts);

        let total: usize = sum as usize;

        // ---- Choose compact vs wide postings encoding. ----
        // Compact: subject_idx + position both fit in u16 → 4 bytes per posting.
        // Wide: anything bigger → 8 bytes per posting.
        let use_compact = subjects.len() <= u16::MAX as usize
                       && max_subj_len <= u16::MAX as usize;

        let mut postings = if use_compact {
            PostingsBuf::Compact(vec![0u32; total])
        } else {
            PostingsBuf::Wide(vec![0u64; total])
        };

        // Write cursor: starts at offsets, advances as we fill.
        let mut cursor = offsets.clone();

        // ---- Pass 2: fill postings. ----
        // Iterate subject by subject, slicing into enc_buf via enc_off.
        // Uses the SAME for_each_selected_kmer helper as the counting pass,
        // so the set of (pos, code) emitted is identical — critical, since
        // the cursor was sized from the counts.
        let alpha = alphabet_size as u32;
        for sidx in 0..subjects.len() {
            let s = subj_slice(sidx);
            if s.len() < k { continue; }
            for_each_selected_kmer(s, k, alpha, window, |pos, h| {
                let h = h as usize;
                let wp = cursor[h] as usize;
                match &mut postings {
                    PostingsBuf::Compact(v) => {
                        v[wp] = (sidx as u32) | ((pos as u32) << 16);
                    }
                    PostingsBuf::Wide(v) => {
                        v[wp] = (sidx as u64) | ((pos as u64) << 32);
                    }
                }
                cursor[h] += 1;
            });
        }
        drop(cursor);

        Self {
            k,
            alphabet_size,
            num_subjects: subjects.len(),
            total_residues,
            data: IndexData::Owned {
                offsets,
                postings,
                enc_buf,
                enc_off,
                // Subject metadata is attached separately via
                // set_subject_meta() by the caller (it has the FASTA records);
                // empty here keeps the build API sequence-only.
                subj_lens: Vec::new(),
                id_blob: Vec::new(),
                id_offsets: Vec::new(),
            },
        }
    }

    // ---- Data accessors (uniform over owned / mmap backing) ------------

    #[inline]
    fn offsets(&self) -> &[u64] {
        match &self.data {
            IndexData::Owned { offsets, .. } => offsets,
            IndexData::Mapped { mmap, offsets, .. } => {
                bytemuck::cast_slice(&mmap[offsets.start..offsets.start + offsets.len])
            }
        }
    }

    #[inline]
    fn enc_buf(&self) -> &[u8] {
        match &self.data {
            IndexData::Owned { enc_buf, .. } => enc_buf,
            IndexData::Mapped { mmap, enc_buf, .. } => &mmap[enc_buf.start..enc_buf.start + enc_buf.len],
        }
    }

    #[inline]
    fn enc_off(&self) -> &[u64] {
        match &self.data {
            IndexData::Owned { enc_off, .. } => enc_off,
            IndexData::Mapped { mmap, enc_off, .. } => {
                bytemuck::cast_slice(&mmap[enc_off.start..enc_off.start + enc_off.len])
            }
        }
    }

    /// Postings as a u32 slice (compact mode) or None (wide mode).
    #[inline]
    fn postings_u32(&self) -> Option<&[u32]> {
        match &self.data {
            IndexData::Owned { postings: PostingsBuf::Compact(v), .. } => Some(v),
            IndexData::Owned { .. } => None,
            IndexData::Mapped { postings_compact: false, .. } => None,
            IndexData::Mapped { mmap, postings, .. } => {
                Some(bytemuck::cast_slice(&mmap[postings.start..postings.start + postings.len]))
            }
        }
    }

    /// Postings as a u64 slice (wide mode) or None (compact mode).
    #[inline]
    fn postings_u64(&self) -> Option<&[u64]> {
        match &self.data {
            IndexData::Owned { postings: PostingsBuf::Wide(v), .. } => Some(v),
            IndexData::Owned { .. } => None,
            IndexData::Mapped { postings_compact: true, .. } => None,
            IndexData::Mapped { mmap, postings, .. } => {
                Some(bytemuck::cast_slice(&mmap[postings.start..postings.start + postings.len]))
            }
        }
    }

    fn subj_lens(&self) -> &[u32] {
        match &self.data {
            IndexData::Owned { subj_lens, .. } => subj_lens,
            IndexData::Mapped { mmap, subj_lens, .. } => {
                bytemuck::cast_slice(&mmap[subj_lens.start..subj_lens.start + subj_lens.len])
            }
        }
    }

    fn id_blob(&self) -> &[u8] {
        match &self.data {
            IndexData::Owned { id_blob, .. } => id_blob,
            IndexData::Mapped { mmap, id_blob, .. } => &mmap[id_blob.start..id_blob.start + id_blob.len],
        }
    }

    fn id_offsets(&self) -> &[u64] {
        match &self.data {
            IndexData::Owned { id_offsets, .. } => id_offsets,
            IndexData::Mapped { mmap, id_offsets, .. } => {
                bytemuck::cast_slice(&mmap[id_offsets.start..id_offsets.start + id_offsets.len])
            }
        }
    }

    /// Subject length in residues. Falls back to the decoded length when no
    /// metadata is attached (in-memory build before set_subject_meta()).
    pub fn subject_len(&self, i: usize) -> usize {
        let lens = self.subj_lens();
        if i < lens.len() {
            lens[i] as usize
        } else {
            // No metadata attached — derive from encoded buffer.
            let eo = self.enc_off();
            (eo[i + 1] - eo[i]) as usize
        }
    }

    /// Subject ID string. Returns "" when no metadata is attached.
    pub fn subject_id(&self, i: usize) -> &str {
        let off = self.id_offsets();
        if i + 1 < off.len() {
            let blob = self.id_blob();
            let s = off[i] as usize;
            let e = off[i + 1] as usize;
            std::str::from_utf8(&blob[s..e]).unwrap_or("")
        } else {
            ""
        }
    }

    /// True once subject metadata (ids + lengths) is available.
    pub fn has_subject_meta(&self) -> bool {
        !self.subj_lens().is_empty()
    }

    /// Attach subject metadata (IDs + lengths) to an in-memory index. The
    /// caller passes them straight from the FASTA records; this lets the
    /// search read all per-subject output data from the index and drop the
    /// FASTA `Database` entirely.
    pub fn set_subject_meta(&mut self, ids: &[String], lens: &[u32]) {
        if let IndexData::Owned { subj_lens, id_blob, id_offsets, .. } = &mut self.data {
            subj_lens.clear();
            subj_lens.extend_from_slice(lens);
            id_blob.clear();
            id_offsets.clear();
            id_offsets.push(0u64);
            let mut acc: u64 = 0;
            for id in ids {
                id_blob.extend_from_slice(id.as_bytes());
                acc += id.len() as u64;
                id_offsets.push(acc);
            }
        }
    }

    // ---- Persistent on-disk index (memory-mapped) ----------------------
    //
    // File layout (all little-endian; an arch/endian marker triggers a clean
    // rebuild on mismatch rather than silent corruption):
    //
    //   [Header: fixed 192 bytes]
    //   [offsets: (num_kmers+1) * u64]      8-byte aligned
    //   [postings: num_postings * (4|8)]    8-byte aligned
    //   [enc_buf: total_residues * u8]      8-byte aligned
    //   [enc_off: (num_subjects+1) * u64]   8-byte aligned
    //   [subj_lens: num_subjects * u32]     8-byte aligned
    //   [id_offsets: (num_subjects+1) * u64]8-byte aligned
    //   [id_blob: bytes]                    8-byte aligned
    //
    // offsets / enc_off / id_offsets are u64 (format v3) so the index scales
    // past 4 GiB of residues; subj_lens stays u32 (one sequence < 4 GiB).
    // Sections are padded to 8-byte boundaries so bytemuck can reinterpret
    // mapped bytes as &[u32]/&[u64] without copying or misalignment.

    /// Serialize this (owned) index to `path` for later memory-mapped reuse.
    pub fn save(&self, path: &std::path::Path) -> std::io::Result<()> {
        use std::io::Write;

        let offsets = self.offsets();
        let enc_buf = self.enc_buf();
        let enc_off = self.enc_off();
        let subj_lens = self.subj_lens();
        let id_offsets = self.id_offsets();
        let id_blob = self.id_blob();
        let postings_compact = self.postings_are_compact();

        // Serialize postings to raw bytes (u32 or u64 little-endian).
        let postings_bytes: &[u8] = if let Some(p) = self.postings_u32() {
            bytemuck::cast_slice(p)
        } else if let Some(p) = self.postings_u64() {
            bytemuck::cast_slice(p)
        } else {
            &[]
        };

        let f = std::fs::File::create(path)?;
        let mut w = std::io::BufWriter::new(f);

        // Compute section offsets (header is 192 bytes; pad each section to 8).
        let align8 = |x: usize| (x + 7) & !7usize;
        let mut cursor = 192usize;
        let sec = |cursor: &mut usize, len: usize| -> (u64, u64) {
            let start = align8(*cursor);
            *cursor = start + len;
            (start as u64, len as u64)
        };
        let (off_offsets, len_offsets) = sec(&mut cursor, offsets.len() * 8);
        let (off_post, len_post) = sec(&mut cursor, postings_bytes.len());
        let (off_enc, len_enc) = sec(&mut cursor, enc_buf.len());
        let (off_encoff, len_encoff) = sec(&mut cursor, enc_off.len() * 8);
        let (off_lens, len_lens) = sec(&mut cursor, subj_lens.len() * 4);
        let (off_idoff, len_idoff) = sec(&mut cursor, id_offsets.len() * 8);
        let (off_idblob, len_idblob) = sec(&mut cursor, id_blob.len());

        // ---- Header (192 bytes) ----
        let mut hdr = Vec::with_capacity(128);
        hdr.extend_from_slice(b"SABERIDX");                 // 8  magic
        hdr.extend_from_slice(&3u32.to_le_bytes());          // 4  format version (3 = u64 offsets)
        hdr.extend_from_slice(&0x0102_0304u32.to_le_bytes()); // 4  endian marker
        hdr.extend_from_slice(&(self.k as u32).to_le_bytes());
        hdr.extend_from_slice(&(self.alphabet_size as u32).to_le_bytes());
        hdr.extend_from_slice(&(self.num_subjects as u64).to_le_bytes());
        hdr.extend_from_slice(&(self.total_residues as u64).to_le_bytes());
        hdr.extend_from_slice(&(postings_compact as u32).to_le_bytes());
        hdr.extend_from_slice(&0u32.to_le_bytes());          // reserved/pad
        for (o, l) in [
            (off_offsets, len_offsets), (off_post, len_post), (off_enc, len_enc),
            (off_encoff, len_encoff), (off_lens, len_lens), (off_idoff, len_idoff),
            (off_idblob, len_idblob),
        ] {
            hdr.extend_from_slice(&o.to_le_bytes());
            hdr.extend_from_slice(&l.to_le_bytes());
        }
        hdr.resize(192, 0);
        w.write_all(&hdr)?;

        // ---- Sections, each padded up to its 8-byte start ----
        let mut pos = 192usize;
        let write_section = |w: &mut std::io::BufWriter<std::fs::File>, start: u64, bytes: &[u8], pos: &mut usize| -> std::io::Result<()> {
            while *pos < start as usize {
                w.write_all(&[0u8])?;
                *pos += 1;
            }
            w.write_all(bytes)?;
            *pos += bytes.len();
            Ok(())
        };
        write_section(&mut w, off_offsets, bytemuck::cast_slice(offsets), &mut pos)?;
        write_section(&mut w, off_post, postings_bytes, &mut pos)?;
        write_section(&mut w, off_enc, enc_buf, &mut pos)?;
        write_section(&mut w, off_encoff, bytemuck::cast_slice(enc_off), &mut pos)?;
        write_section(&mut w, off_lens, bytemuck::cast_slice(subj_lens), &mut pos)?;
        write_section(&mut w, off_idoff, bytemuck::cast_slice(id_offsets), &mut pos)?;
        write_section(&mut w, off_idblob, id_blob, &mut pos)?;
        w.flush()?;
        Ok(())
    }

    /// Open a persistent index by memory-mapping `path`. Returns an error if
    /// the file is not a SABER index, the format version differs, or the
    /// byte order doesn't match this machine (caller should rebuild).
    pub fn open(path: &std::path::Path) -> std::io::Result<Self> {
        let file = std::fs::File::open(path)?;
        // SAFETY: standard mmap of a read-only file we treat as immutable.
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        if mmap.len() < 192 || &mmap[0..8] != b"SABERIDX" {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "not a SABER index file"));
        }
        let rd_u32 = |o: usize| u32::from_le_bytes(mmap[o..o + 4].try_into().unwrap());
        let rd_u64 = |o: usize| u64::from_le_bytes(mmap[o..o + 8].try_into().unwrap());

        let version = rd_u32(8);
        if version != 3 {
            let hint = if version == 2 {
                " (v2 used 32-bit offsets; rebuild with --makedb to get 64-bit offsets)"
            } else {
                "; rebuild"
            };
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData,
                format!("SABER index format version {version} != 3{hint}")));
        }
        if rd_u32(12) != 0x0102_0304 {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData,
                "SABER index byte order mismatch; rebuild on this machine"));
        }
        let k = rd_u32(16) as usize;
        let alphabet_size = rd_u32(20) as usize;
        let num_subjects = rd_u64(24) as usize;
        let total_residues = rd_u64(32) as usize;
        let postings_compact = rd_u32(40) != 0;

        // 7 (offset,len) u64 pairs start at byte 48.
        let mut base = 48usize;
        let mut next_range = || {
            let start = rd_u64(base) as usize;
            let len = rd_u64(base + 8) as usize;
            base += 16;
            ByteRange { start, len }
        };
        let offsets = next_range();
        let postings = next_range();
        let enc_buf = next_range();
        let enc_off = next_range();
        let subj_lens = next_range();
        let id_offsets = next_range();
        let id_blob = next_range();

        // Validate ranges are inside the file and aligned for casting.
        for r in [offsets, postings, enc_buf, enc_off, subj_lens, id_offsets, id_blob] {
            if r.start + r.len > mmap.len() {
                return Err(std::io::Error::new(std::io::ErrorKind::InvalidData,
                    "SABER index section out of bounds"));
            }
        }
        // offsets / enc_off / id_offsets are u64 → need 8-byte alignment to
        // reinterpret as &[u64]; subj_lens is u32 → 4-byte. All sections are
        // 8-padded on write, so both hold.
        for r in [offsets, enc_off, id_offsets] {
            if r.start % 8 != 0 {
                return Err(std::io::Error::new(std::io::ErrorKind::InvalidData,
                    "SABER index u64 section misaligned"));
            }
        }
        if subj_lens.start % 4 != 0 {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData,
                "SABER index u32 section misaligned"));
        }
        if !postings_compact && postings.start % 8 != 0 {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData,
                "SABER index u64 postings misaligned"));
        }

        Ok(Self {
            k,
            alphabet_size,
            num_subjects,
            total_residues,
            data: IndexData::Mapped {
                mmap, postings_compact,
                offsets, postings, enc_buf, enc_off, subj_lens, id_offsets, id_blob,
            },
        })
    }

    pub fn k(&self) -> usize { self.k }
    pub fn num_subjects(&self) -> usize { self.num_subjects }
    pub fn num_postings(&self) -> usize {
        self.postings_u32().map(|s| s.len())
            .or_else(|| self.postings_u64().map(|s| s.len()))
            .unwrap_or(0)
    }
    pub fn total_residues(&self) -> usize { self.total_residues }

    /// Whether this index uses the compact (4-byte) posting encoding.
    /// Returns false on huge DBs where wide (8-byte) postings were needed.
    pub fn postings_are_compact(&self) -> bool {
        self.postings_u32().is_some()
    }

    /// Access the i-th subject as encoded bytes (0..alphabet_size or 255 for
    /// ambiguous residues). Caller is expected to translate 255 to a
    /// substitution-matrix-friendly value before passing to SW (the SIMD
    /// kernel clamps anything > 23 to X = 22).
    pub fn encoded_subject(&self, i: usize) -> &[u8] {
        let eo = self.enc_off();
        let lo = eo[i] as usize;
        let hi = eo[i + 1] as usize;
        &self.enc_buf()[lo..hi]
    }

    /// Reconstruct the ASCII residue sequence for subject `i` from the
    /// encoded buffer. Used by the traceback step after the database's
    /// ASCII sequences have been dropped to save memory. Ambiguous
    /// residues (encoded as 255) reconstruct to ASCII `X`.
    ///
    /// Allocates a new Vec<u8>; only called for the small number of
    /// candidates that survive scoring and need a traceback (≤ max_aligns
    /// per query, typically ≤ 10–100).
    pub fn ascii_subject(&self, i: usize) -> Vec<u8> {
        let enc = self.encoded_subject(i);
        let mut out = Vec::with_capacity(enc.len());
        for &b in enc {
            out.push(decode_residue(b, self.alphabet_size));
        }
        out
    }

    /// Memory footprint in bytes of the index's large arrays. For an mmap'd
    /// index this is the on-disk/virtual size, not resident memory (the OS
    /// pages it in on demand).
    pub fn memory_bytes(&self) -> usize {
        let postings_bytes = self.postings_u32().map(|s| s.len() * 4)
            .or_else(|| self.postings_u64().map(|s| s.len() * 8))
            .unwrap_or(0);
        self.offsets().len() * 8
            + postings_bytes
            + self.enc_buf().len()
            + self.enc_off().len() * 8
    }

    /// Query the index with a raw ASCII protein query. Returns one
    /// [`SeedHit`] per (query_kmer_position, matching_db_occurrence).
    pub fn query_protein(&self, query: &[u8]) -> Vec<SeedHit> {
        let encoded: Vec<u8> = query
            .iter()
            .map(|&c| encode_protein_byte(c).unwrap_or(255))
            .collect();
        self.query_encoded(&encoded)
    }

    pub fn query_nucleotide(&self, query: &[u8]) -> Vec<SeedHit> {
        let encoded: Vec<u8> = query
            .iter()
            .map(|&c| encode_nucleotide_byte(c).unwrap_or(255))
            .collect();
        self.query_encoded(&encoded)
    }

    fn query_encoded(&self, query: &[u8]) -> Vec<SeedHit> {
        if query.len() < self.k {
            return Vec::new();
        }
        let alpha = self.alphabet_size as u32;
        let mut hits = Vec::new();
        let stop = query.len() - self.k + 1;

        // Fetch the backing slices ONCE (owned Vec or mmap region), then run
        // a branch-free hot loop over them. The postings variant is dispatched
        // a single time, outside the loop.
        let offsets = self.offsets();
        if let Some(buf) = self.postings_u32() {
            for qp in 0..stop {
                let Some(h) = hash_kmer(&query[qp..qp + self.k], alpha) else { continue };
                let h = h as usize;
                let start = offsets[h] as usize;
                let end = offsets[h + 1] as usize;
                for &p in &buf[start..end] {
                    let sidx = p & 0xFFFF;
                    let sp = p >> 16;
                    hits.push(SeedHit {
                        subject_idx: sidx,
                        diagonal: sp as i32 - qp as i32,
                        q_pos: qp as u32,
                        s_pos: sp,
                    });
                }
            }
        } else if let Some(buf) = self.postings_u64() {
            for qp in 0..stop {
                let Some(h) = hash_kmer(&query[qp..qp + self.k], alpha) else { continue };
                let h = h as usize;
                let start = offsets[h] as usize;
                let end = offsets[h + 1] as usize;
                for &p in &buf[start..end] {
                    let sidx = p as u32;
                    let sp = (p >> 32) as u32;
                    hits.push(SeedHit {
                        subject_idx: sidx,
                        diagonal: sp as i32 - qp as i32,
                        q_pos: qp as u32,
                        s_pos: sp,
                    });
                }
            }
        }
        hits
    }

    /// Query the index using **neighborhood word expansion** — for each
    /// query k-mer, look up not just the exact match but every k-mer
    /// scoring ≥ T against it (per [`NeighborhoodGenerator`]). This is
    /// BLAST's seed-stage approach; the key effect is that real homologs
    /// produce dense seed clusters on the right diagonal even at low
    /// sequence identity, while random subjects don't.
    ///
    /// The caller is responsible for matching `nbgen.k()` to `self.k()`
    /// and using a `nbgen` built against the same scoring matrix the
    /// downstream SW will use.
    pub fn query_protein_with_neighborhood(
        &self,
        query: &[u8],
        nbgen: &super::neighborhood::NeighborhoodGenerator,
    ) -> Vec<SeedHit> {
        debug_assert_eq!(self.k, nbgen.k());
        let encoded: Vec<u8> = query
            .iter()
            .map(|&c| encode_protein_byte(c).unwrap_or(255))
            .collect();
        if encoded.len() < self.k {
            return Vec::new();
        }

        let alpha = self.alphabet_size as u32;
        let mut hits = Vec::new();
        let mut neighbors: Vec<u32> = Vec::new();
        let stop = encoded.len() - self.k + 1;

        // Fetch backing slices once; dispatch postings variant once.
        let offsets = self.offsets();
        if let Some(buf) = self.postings_u32() {
            for qp in 0..stop {
                let kmer = &encoded[qp..qp + self.k];
                if kmer.iter().any(|&b| (b as usize) >= 20) { continue; }
                neighbors.clear();
                nbgen.expand(kmer, alpha, &mut neighbors);
                for &nh in &neighbors {
                    let h = nh as usize;
                    let start = offsets[h] as usize;
                    let end = offsets[h + 1] as usize;
                    for &p in &buf[start..end] {
                        let sidx = p & 0xFFFF;
                        let sp = p >> 16;
                        hits.push(SeedHit {
                            subject_idx: sidx,
                            diagonal: sp as i32 - qp as i32,
                            q_pos: qp as u32,
                            s_pos: sp,
                        });
                    }
                }
            }
        } else if let Some(buf) = self.postings_u64() {
            for qp in 0..stop {
                let kmer = &encoded[qp..qp + self.k];
                if kmer.iter().any(|&b| (b as usize) >= 20) { continue; }
                neighbors.clear();
                nbgen.expand(kmer, alpha, &mut neighbors);
                for &nh in &neighbors {
                    let h = nh as usize;
                    let start = offsets[h] as usize;
                    let end = offsets[h + 1] as usize;
                    for &p in &buf[start..end] {
                        let sidx = p as u32;
                        let sp = (p >> 32) as u32;
                        hits.push(SeedHit {
                            subject_idx: sidx,
                            diagonal: sp as i32 - qp as i32,
                            q_pos: qp as u32,
                            s_pos: sp,
                        });
                    }
                }
            }
        }
        hits
    }

    /// Score candidates using **BLAST's two-hit heuristic**.
    ///
    /// For each subject, group seed hits by (exact) diagonal, then within
    /// each diagonal look for *pairs of hits whose query positions are
    /// within `window` residues of each other*. A subject becomes a
    /// candidate iff it has ≥ 1 such pair (i.e., at least one "double hit"
    /// on a diagonal). This is much more selective than my old
    /// "diagonal bucket count ≥ min_hits" approach — random subjects
    /// rarely produce close-on-diagonal pairs even with neighborhood
    /// expansion, while real homologs produce many.
    ///
    /// Candidates are then ranked by their pair count; the top `top_n`
    /// are returned. `best_diagonal` is the diagonal of the best-scoring
    /// pair cluster, used by downstream banded SW.
    ///
    /// Parameters per BLAST defaults: `window = 40` (matches the BLAST
    /// `-window_size` default for blastp).
    pub fn candidates_two_hit(
        &self,
        hits: &[SeedHit],
        window: u32,
        top_n: usize,
    ) -> Vec<Candidate> {
        if hits.is_empty() {
            return Vec::new();
        }

        // Sort by (subject, diagonal, q_pos).
        let mut sorted = hits.to_vec();
        sorted.sort_unstable_by_key(|h| (h.subject_idx, h.diagonal, h.q_pos));

        let mut out: Vec<Candidate> = Vec::new();
        let mut i = 0;
        while i < sorted.len() {
            let sidx = sorted[i].subject_idx;
            // Find this subject's range in the sorted array.
            let mut j = i;
            while j < sorted.len() && sorted[j].subject_idx == sidx {
                j += 1;
            }

            // Within this subject, walk through diagonals; for each, count
            // pairs (k, k+1) with q_pos distance ≤ window.
            let mut best_pairs: u32 = 0;
            let mut best_diag: i32 = 0;
            // Track the median q_pos on the best-pair diagonal as the
            // anchor for ungapped extension downstream.
            let mut best_anchor_qpos: u32 = 0;
            let mut diag_start = i;
            while diag_start < j {
                let diag = sorted[diag_start].diagonal;
                let mut diag_end = diag_start + 1;
                while diag_end < j && sorted[diag_end].diagonal == diag {
                    diag_end += 1;
                }
                // Count close-pair hits within this diagonal.
                let mut pairs: u32 = 0;
                for k in diag_start..diag_end.saturating_sub(1) {
                    if sorted[k + 1].q_pos.saturating_sub(sorted[k].q_pos) <= window {
                        pairs += 1;
                    }
                }
                if pairs > best_pairs {
                    best_pairs = pairs;
                    best_diag = diag;
                    // Anchor: median q_pos on this diagonal.
                    let mid = diag_start + (diag_end - diag_start) / 2;
                    best_anchor_qpos = sorted[mid].q_pos;
                }
                diag_start = diag_end;
            }

            if best_pairs >= 1 {
                out.push(Candidate {
                    subject_idx: sidx,
                    hit_count: best_pairs,
                    best_diagonal: best_diag,
                    seed_q_pos: best_anchor_qpos,
                });
            }
            i = j;
        }
        out.sort_unstable_by(|a, b| b.hit_count.cmp(&a.hit_count));
        out.truncate(top_n);
        out
    }

    /// Score candidates from raw seed hits using diagonal-bucket clustering.
    /// Used when neighborhood expansion is OFF (T=0); for T>0 use
    /// [`Self::candidates_two_hit`] which is much more selective.
    pub fn candidates(
        &self,
        hits: &[SeedHit],
        min_bucket_hits: u32,
        diag_bucket_size: i32,
        top_n: usize,
    ) -> Vec<Candidate> {
        if hits.is_empty() {
            return Vec::new();
        }
        // Sort by (subject_idx, diagonal) to enable a single linear pass.
        let mut sorted = hits.to_vec();
        sorted.sort_unstable_by_key(|h| (h.subject_idx, h.diagonal));

        let mut out: Vec<Candidate> = Vec::new();
        let mut i = 0;
        while i < sorted.len() {
            let sidx = sorted[i].subject_idx;
            // Scan this subject's contiguous run.
            let mut j = i;
            while j < sorted.len() && sorted[j].subject_idx == sidx {
                j += 1;
            }
            // Within this subject's hits, find the best (highest-count)
            // diagonal bucket.
            let mut best_count: u32 = 0;
            let mut best_diag: i32 = 0;
            let mut best_anchor_qpos: u32 = 0;
            let mut bucket_start = i;
            while bucket_start < j {
                let bucket_id = sorted[bucket_start].diagonal / diag_bucket_size;
                let mut bucket_end = bucket_start + 1;
                let mut count: u32 = 1;
                let mut diag_sum: i64 = sorted[bucket_start].diagonal as i64;
                while bucket_end < j
                    && sorted[bucket_end].diagonal / diag_bucket_size == bucket_id
                {
                    count += 1;
                    diag_sum += sorted[bucket_end].diagonal as i64;
                    bucket_end += 1;
                }
                if count > best_count {
                    best_count = count;
                    best_diag = (diag_sum / count as i64) as i32;
                    // Use the median hit's q_pos as the anchor.
                    let mid = bucket_start + (bucket_end - bucket_start) / 2;
                    best_anchor_qpos = sorted[mid].q_pos;
                }
                bucket_start = bucket_end;
            }
            if best_count >= min_bucket_hits {
                out.push(Candidate {
                    subject_idx: sidx,
                    hit_count: best_count,
                    best_diagonal: best_diag,
                    seed_q_pos: best_anchor_qpos,
                });
            }
            i = j;
        }
        // Rank by hit count descending, take top_n.
        out.sort_unstable_by(|a, b| b.hit_count.cmp(&a.hit_count));
        out.truncate(top_n);
        out
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn count_kmers_into(acc: &mut [u32], s: &[u8], k: usize, alphabet_size: usize, window: usize) {
    let alpha = alphabet_size as u32;
    for_each_selected_kmer(s, k, alpha, window, |_pos, h| {
        acc[h as usize] = acc[h as usize].saturating_add(1);
    });
}

#[inline]
fn hash_kmer(kmer: &[u8], alphabet_size: u32) -> Option<u32> {
    let mut h: u32 = 0;
    for &b in kmer {
        if (b as u32) >= alphabet_size {
            return None; // ambiguous residue → skip this k-mer
        }
        h = h * alphabet_size + b as u32;
    }
    Some(h)
}

/// Integer finalizer (murmur3 fmix32) used to give every k-mer code a
/// well-distributed *ranking* hash for minimizer selection. Without mixing,
/// the lexicographic k-mer code biases minimizer selection toward k-mers
/// starting with low-index residues (A, R, N…), which clusters seeds and
/// hurts sensitivity. Mixing makes minimizer selection effectively uniform.
#[inline]
fn mix32(mut x: u32) -> u32 {
    x ^= x >> 16;
    x = x.wrapping_mul(0x7feb_352d);
    x ^= x >> 15;
    x = x.wrapping_mul(0x846c_a68b);
    x ^= x >> 16;
    x
}

/// Invoke `f(pos, code)` for each *selected* k-mer of `s`, where selection
/// depends on the minimizer window `w`:
///
/// * `w <= 1` — dense: every valid k-mer position is emitted (this is the
///   original, no-sampling behavior).
/// * `w >= 2` — sparse: within each sliding window of `w` consecutive
///   k-mers, only the k-mer with the smallest *mixed* hash is emitted. A
///   minimizer selected by overlapping windows is emitted once. This is the
///   standard (w, k)-minimizer scheme. It cuts the number of indexed
///   positions to roughly `2 / (w + 1)` of the dense count.
///
/// `code` is always the canonical [`hash_kmer`] code (the bucket key), so a
/// dense query lookup by canonical code still finds these minimizers. The
/// mixed hash is used only for *ranking* within the window.
///
/// k-mers containing ambiguous residues (where [`hash_kmer`] returns `None`)
/// are skipped and act as hard window boundaries.
#[inline]
fn for_each_selected_kmer<F: FnMut(usize, u32)>(
    s: &[u8],
    k: usize,
    alphabet_size: u32,
    w: usize,
    mut f: F,
) {
    if s.len() < k {
        return;
    }
    let stop = s.len() - k + 1;

    if w <= 1 {
        for pos in 0..stop {
            if let Some(h) = hash_kmer(&s[pos..pos + k], alphabet_size) {
                f(pos, h);
            }
        }
        return;
    }

    // Sliding-window minimizer. We keep it simple and robust against
    // ambiguous-residue breaks: scan each length-`w` window, pick the min
    // mixed-hash valid k-mer in it, emit if it differs from the last emit.
    // O(stop * w); w is small (≤ ~16) and this runs only at build time.
    let mut last_emitted: i64 = -1;
    for win_start in 0..stop {
        let win_end = (win_start + w).min(stop); // exclusive
        let mut best_pos: Option<usize> = None;
        let mut best_rank: u32 = u32::MAX;
        let mut best_code: u32 = 0;
        for pos in win_start..win_end {
            if let Some(code) = hash_kmer(&s[pos..pos + k], alphabet_size) {
                let r = mix32(code);
                // Tie-break on position (smaller pos wins) for determinism.
                if r < best_rank {
                    best_rank = r;
                    best_pos = Some(pos);
                    best_code = code;
                }
            }
        }
        if let Some(p) = best_pos {
            if p as i64 != last_emitted {
                f(p, best_code);
                last_emitted = p as i64;
            }
        }
    }
}

#[inline]
fn encode_protein_byte(c: u8) -> Option<u8> {
    // The 20 standard amino acids ordered as in the BLOSUM 24-letter alphabet's
    // first 20 entries: ARNDCQEGHILKMFPSTWYV. Lower-case folded.
    match c {
        b'A' | b'a' => Some(0),
        b'R' | b'r' => Some(1),
        b'N' | b'n' => Some(2),
        b'D' | b'd' => Some(3),
        b'C' | b'c' => Some(4),
        b'Q' | b'q' => Some(5),
        b'E' | b'e' => Some(6),
        b'G' | b'g' => Some(7),
        b'H' | b'h' => Some(8),
        b'I' | b'i' => Some(9),
        b'L' | b'l' => Some(10),
        b'K' | b'k' => Some(11),
        b'M' | b'm' => Some(12),
        b'F' | b'f' => Some(13),
        b'P' | b'p' => Some(14),
        b'S' | b's' => Some(15),
        b'T' | b't' => Some(16),
        b'W' | b'w' => Some(17),
        b'Y' | b'y' => Some(18),
        b'V' | b'v' => Some(19),
        // B (Asx), Z (Glx), X (any), * (stop), J (Leu/Ile) → ambiguous, excluded
        _ => None,
    }
}

#[inline]
fn encode_nucleotide_byte(c: u8) -> Option<u8> {
    match c.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// Inverse of `encode_protein_byte` / `encode_nucleotide_byte`.
/// Maps an encoded residue back to its ASCII character. Uppercase only.
/// Ambiguous residues (encoded as 255 by the sentinel-encoding) decode to `X`
/// for protein and `N` for nucleotide.
#[inline]
fn decode_residue(b: u8, alphabet_size: usize) -> u8 {
    if alphabet_size == PROTEIN_ALPHABET_FOR_INDEX {
        // First 20 of BLOSUM 24-letter alphabet.
        const PROTEIN_DECODE: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";
        if (b as usize) < PROTEIN_DECODE.len() {
            PROTEIN_DECODE[b as usize]
        } else {
            b'X'
        }
    } else {
        // Nucleotide.
        const NUC_DECODE: &[u8; 4] = b"ACGT";
        if (b as usize) < NUC_DECODE.len() {
            NUC_DECODE[b as usize]
        } else {
            b'N'
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn persistent_index_roundtrip_matches_in_memory() {
        let subjects: Vec<&[u8]> = vec![
            b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF",
            b"GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHA",
            b"MVLSPADKTNVKAAWGKVGAHAGEY",
            b"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW",
        ];
        let ids = vec![
            "sp|P1".to_string(), "sp|P2".to_string(),
            "tr|P3".to_string(), "junk|P4".to_string(),
        ];
        let lens: Vec<u32> = subjects.iter().map(|s| s.len() as u32).collect();

        let mut mem = KmerIndex::build_protein_windowed(&subjects, 4, 2);
        mem.set_subject_meta(&ids, &lens);

        let mut path = std::env::temp_dir();
        path.push(format!("saber_test_idx_{}.sdx", std::process::id()));
        mem.save(&path).expect("save");
        let mapped = KmerIndex::open(&path).expect("open");

        assert_eq!(mem.num_subjects(), mapped.num_subjects());
        assert_eq!(mem.num_postings(), mapped.num_postings());
        assert_eq!(mem.k(), mapped.k());
        assert_eq!(mem.total_residues(), mapped.total_residues());
        assert!(mapped.has_subject_meta());

        for i in 0..subjects.len() {
            assert_eq!(mem.subject_id(i), mapped.subject_id(i), "id {i}");
            assert_eq!(mem.subject_len(i), mapped.subject_len(i), "len {i}");
            assert_eq!(mem.encoded_subject(i), mapped.encoded_subject(i), "enc {i}");
            assert_eq!(mem.ascii_subject(i), mapped.ascii_subject(i), "ascii {i}");
        }

        let q = b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF";
        let h_mem = mem.query_protein(q);
        let h_map = mapped.query_protein(q);
        assert_eq!(h_mem.len(), h_map.len(), "hit count");
        for (a, b) in h_mem.iter().zip(h_map.iter()) {
            assert_eq!(a.subject_idx, b.subject_idx);
            assert_eq!(a.q_pos, b.q_pos);
            assert_eq!(a.s_pos, b.s_pos);
            assert_eq!(a.diagonal, b.diagonal);
        }

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn open_rejects_non_index_file() {
        let mut path = std::env::temp_dir();
        path.push(format!("saber_test_bad_{}.sdx", std::process::id()));
        std::fs::write(&path, b"this is not a SABER index file at all, padding..............").unwrap();
        assert!(KmerIndex::open(&path).is_err());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn index_uses_u64_offsets_format_v3() {
        // Regression guard for the "encoded subject offsets overflowed u32"
        // bug: a >4 GiB database overflowed the old 32-bit offset arrays.
        // The fix widens offsets / enc_off / id_offsets to u64 and bumps the
        // on-disk format to v3. We can't allocate 4 GiB in a unit test, so we
        // assert the *format* instead: the persisted offset sections must use
        // 8-byte (u64) elements, and old v2 files must be rejected.
        let subjects: Vec<&[u8]> = vec![
            b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF",
            b"GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHA",
            b"MVLSPADKTNVKAAWGKVGAHAGEY",
        ];
        let ids = vec!["sp|P1".to_string(), "sp|P2".to_string(), "tr|P3".to_string()];
        let lens: Vec<u32> = subjects.iter().map(|s| s.len() as u32).collect();
        let mut mem = KmerIndex::build_protein(&subjects, 4);
        mem.set_subject_meta(&ids, &lens);

        let mut path = std::env::temp_dir();
        path.push(format!("saber_u64fmt_{}.sdx", std::process::id()));
        mem.save(&path).expect("save");
        let bytes = std::fs::read(&path).expect("read");

        // Header: magic[0..8], format version u32 @8 (must be 3).
        assert_eq!(&bytes[0..8], b"SABERIDX");
        let version = u32::from_le_bytes(bytes[8..12].try_into().unwrap());
        assert_eq!(version, 3, "persistent format must be v3 (u64 offsets)");

        // Section table: 7 (off u64, len u64) pairs starting at byte 48.
        // enc_off is pair index 3 → its length lives at byte 48 + 3*16 + 8.
        let enc_off_len = u64::from_le_bytes(
            bytes[48 + 3 * 16 + 8..48 + 3 * 16 + 16].try_into().unwrap());
        assert_eq!(enc_off_len as usize, (subjects.len() + 1) * 8,
            "enc_off must be stored as u64 (8 bytes/elem), not u32");
        // id_offsets is pair index 5 → 8 bytes/elem too.
        let id_off_len = u64::from_le_bytes(
            bytes[48 + 5 * 16 + 8..48 + 5 * 16 + 16].try_into().unwrap());
        assert_eq!(id_off_len as usize, (subjects.len() + 1) * 8,
            "id_offsets must be stored as u64");

        // A v2-stamped file (32-bit offsets) must be cleanly rejected.
        let mut v2 = bytes.clone();
        v2[8..12].copy_from_slice(&2u32.to_le_bytes());
        let mut p2 = std::env::temp_dir();
        p2.push(format!("saber_u64fmt_v2_{}.sdx", std::process::id()));
        std::fs::write(&p2, &v2).unwrap();
        let err = match KmerIndex::open(&p2) {
            Ok(_) => panic!("a v2 (32-bit offset) index must be rejected"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("32-bit") || err.contains("rebuild"),
            "v2 rejection should hint at rebuilding for 64-bit offsets: {err}");

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&p2);
    }

    #[test]
    fn protein_index_finds_perfect_match() {
        let subjects: Vec<&[u8]> = vec![
            b"MVLSPADKTNVKAAWGKVGAHAGEY",
            b"PPPPPPPPPPPPPPPPPPPPPPPPP",
            b"MVLSPADKTNVKAAWGKVGAHAGEY", // duplicate of #0
        ];
        let idx = KmerIndex::build_protein(&subjects, 4);
        assert_eq!(idx.num_subjects(), 3);
        let query = b"MVLSPADKTNVKAAWGKVGAHAGEY";
        let hits = idx.query_protein(query);

        // Subjects 0 and 2 are identical to the query → should generate many hits each.
        let candidates = idx.candidates(&hits, 2, 16, 100);
        assert!(!candidates.is_empty());
        // Top two candidates should be subjects 0 and 2 (in some order).
        let top_subj: Vec<u32> = candidates.iter().take(2).map(|c| c.subject_idx).collect();
        assert!(top_subj.contains(&0));
        assert!(top_subj.contains(&2));
        // Junk subject should not be in candidates (with min_bucket_hits=2).
        assert!(!candidates.iter().any(|c| c.subject_idx == 1));
    }

    #[test]
    fn protein_index_orders_by_hit_count() {
        let subjects: Vec<&[u8]> = vec![
            b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFL",          // strong overlap
            b"MVLSPAD",                                       // weak overlap
            b"YYYYYYYYYYYYYY",                                // no overlap
        ];
        let idx = KmerIndex::build_protein(&subjects, 4);
        let hits = idx.query_protein(b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFL");
        let candidates = idx.candidates(&hits, 1, 16, 10);
        assert!(candidates.len() >= 1);
        assert_eq!(candidates[0].subject_idx, 0); // strongest first
    }

    #[test]
    fn nucleotide_index_basic() {
        let subjects: Vec<&[u8]> = vec![
            b"ACGTACGTACGTACGTACGTACGT",
            b"TTTTTTTTTTTTTTTTTTTT",
        ];
        let idx = KmerIndex::build_nucleotide(&subjects, 5);
        let hits = idx.query_nucleotide(b"ACGTACGTACGTACGTACGT");
        let candidates = idx.candidates(&hits, 1, 16, 10);
        assert!(candidates.iter().any(|c| c.subject_idx == 0));
        assert!(!candidates.iter().any(|c| c.subject_idx == 1));
    }

    #[test]
    fn empty_query_returns_no_hits() {
        let subjects: Vec<&[u8]> = vec![b"MVLSPADKTNVKAAWGK"];
        let idx = KmerIndex::build_protein(&subjects, 4);
        assert!(idx.query_protein(b"").is_empty());
        assert!(idx.query_protein(b"MV").is_empty()); // shorter than k
    }

    #[test]
    fn ambiguous_residues_skip_kmer() {
        // Subject contains X — k-mers spanning the X should not be indexed.
        let subjects: Vec<&[u8]> = vec![b"MVLSPAXXXXXXKVGAHAGEY"];
        let idx = KmerIndex::build_protein(&subjects, 4);
        // Querying MVLS should still match (it doesn't span an X).
        let hits = idx.query_protein(b"MVLS");
        assert!(!hits.is_empty(), "MVLS should match subject prefix");
        // Querying SPAX would have been a near-match before — but k-mers
        // with X are excluded, so SPAX produces no hit.
        let hits_x = idx.query_protein(b"SPAX");
        assert!(hits_x.is_empty(), "X-containing k-mer should produce no hits");
    }

    #[test]
    fn diagonal_clustering_separates_real_hits_from_random() {
        // Two subjects, one with hits clustered on a diagonal, one with scattered hits.
        // The clustered one should be a candidate; the scattered one should not.
        let subjects: Vec<&[u8]> = vec![
            // Subject 0: contains MVLSPADK at position 5 → same diagonal across multiple k-mers
            b"AAAAAMVLSPADKAAAAAA",
            // Subject 1: contains MVLS at pos 5 and PADK at pos 30 → different diagonals
            b"AAAAAMVLSAAAAAAAAAAAAAAAAAAAAAAAAAPADKAAAA",
        ];
        let idx = KmerIndex::build_protein(&subjects, 4);
        let hits = idx.query_protein(b"MVLSPADK");
        // Require at least 2 hits in the same diagonal bucket of size 16.
        let cands = idx.candidates(&hits, 2, 16, 10);
        // Subject 0 has MVLS at qpos=0,subjpos=5 (diag=5) AND PADK at qpos=4,subjpos=9 (diag=5)
        // → both on diagonal 5, same bucket. Should pass the 2-hit filter.
        assert!(cands.iter().any(|c| c.subject_idx == 0), "clustered hits should pass: {:?}", cands);
        // Subject 1's hits are on different diagonals, no bucket has 2 → fails filter.
        assert!(!cands.iter().any(|c| c.subject_idx == 1));
    }
}
