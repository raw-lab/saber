//! The SABER search engine (v1.2; scan schedule and per-subject E-values since v1.3).
//!
//! One implementation serves the CLI (`pp`, `nn`, `nx`, `pn` modes) and the
//! library API. Per query:
//!
//! 1. **Preparation** — encode, compute the composition-bias correction and
//!    the tantan seed mask, derive per-query statistics and the exact score
//!    floor implied by `-E` / `--min-bitscore`.
//! 2. **Seeding, bounded** — for every partition of the index, look up the
//!    query's seeds (skipping masked, ambiguous and over-frequent seeds), count
//!    postings per subject in one pass, then scatter them into a worker-local
//!    buffer of bounded capacity grouped by subject. A partition whose hits
//!    exceed the buffer is processed in consecutive subject ranges (postings
//!    are subject-sorted within each bucket, so a range is found by binary
//!    search); a single subject that alone exceeds the buffer is streamed
//!    without materialising its hits. No evidence is discarded to bound
//!    memory.
//! 3. **Diagonal evidence and ungapped extension** — per subject, hits are
//!    visited in query order through a diagonal table. A hit that falls inside
//!    a segment already extended on its diagonal is skipped (counted as an
//!    avoided extension); otherwise an X-drop ungapped extension is run, with
//!    the bias correction applied. Optional two-hit triggering requires a
//!    second, non-overlapping hit on the same diagonal within a window.
//! 4. **Candidate regions** — HSPs of a subject are clustered by diagonal into
//!    regions (several per subject, so multi-domain matches survive); regions
//!    enter a bounded top-N heap with deterministic tie-breaking. When the
//!    heap drops regions, the drop is counted and reported — a capped search
//!    is never described as exhaustive.
//! 5. **Gapped scoring** — every region's subject window (the region's
//!    diagonal band projected over the full query, plus slack) is scored by
//!    the SIMD kernel with the bias correction; saturated lanes are recomputed
//!    exactly. Candidates below the exact E-value floor are dropped here.
//! 6. **Progressive traceback** — candidates are traced back in score order
//!    with a genuine banded DP that is widened adaptively until it reproduces
//!    the window score (exactness w.r.t. the window); filters on identity and
//!    coverage are applied afterwards and the walk continues until
//!    `max_target_seqs` subjects pass or candidates run out — no fixed
//!    overshoot.

use rayon::prelude::*;
use std::time::Instant;

use crate::algorithm::banded::{banded_align_until, banded_score, Band, BandScratch, LocalAlignment, Scheme};
use crate::algorithm::sw_simd::{score_many, KernelCounters, SimdScratch};
use crate::alphabet::{code_index, decode_protein, ProteinCodes, NUCLEOTIDE_ALPHABET};
use crate::masking::{MaskingMode, Tantan, TantanScratch};
use crate::scoring::bias::composition_bias;
use crate::scoring::matrices::ScoreMatrix;
use crate::scoring::stats::{EvalueModel, GapCosts, QueryStatistics, SearchStatistics};

use super::index::{Postings, SeedIndex, SeqKind};

mod scan;
pub(crate) use scan::seed_codes_into;
pub use scan::ScanSeeds;

/// Search parameters (CLI flags map 1:1 onto these).
#[derive(Debug, Clone)]
pub struct SearchParams {
    pub evalue: f64,
    pub min_bitscore: Option<f64>,
    /// Subjects reported per query (`-k`).
    pub max_target_seqs: usize,
    /// Alignments reported per subject.
    pub max_hsps: usize,
    pub min_identity: f64,
    pub min_query_cover: f64,
    pub min_subject_cover: f64,
    /// Minimum (bias-corrected) ungapped HSP raw score for a region.
    pub hsp_threshold: i32,
    /// Region threshold: a candidate region is kept only if one of its
    /// ungapped segments, or the sum of two or more, reaches this fraction
    /// of the query's exact raw-score floor at `-E` or at
    /// [`REGION_FLOOR_EVALUE`], whichever is more permissive (a lone segment
    /// far below the reporting floor rarely grows into a reportable gapped
    /// alignment). 0 = use `hsp_threshold` alone.
    pub hsp_floor_fraction: f64,
    pub x_drop_ungapped: i32,
    /// Two-hit window (0 = single-hit triggering).
    pub two_hit_window: u32,
    /// Hits whose own seed score (Σ S over the seed's care positions, with
    /// bias) is below this are ignored (0 = keep all).
    pub seed_min_score: i32,
    /// With two-hit triggering, a single hit whose seed score reaches this
    /// value is extended immediately (0 = never).
    pub single_hit_score: i32,
    /// Hit filter (DIAMOND-style): before extending a hit, count identities
    /// in the seed alphabet over a 32-residue window centred on the seed;
    /// hits with fewer than this many are skipped (0 = off).
    pub hit_filter_min_ids: u32,
    /// Weight of exact (20-letter) identities in the window filter: the
    /// filter value is `seed-alphabet identities + w × exact identities`
    /// (0 = seed-alphabet identities only).
    pub hit_filter_exact_weight: u32,
    /// Initial traceback band half-width beyond the region's diagonal span.
    pub band_pad: i64,
    /// Subject-window slack around the projected diagonal band.
    pub window_slack: i64,
    /// Score every candidate region inside the diagonal band spanned by its
    /// HSPs ± this many diagonals instead of the whole window rectangle
    /// (0 = window scoring).
    pub score_band: i64,
    /// Two-stage banded scoring: every banded region is first scored inside
    /// its HSP diagonals ± this many diagonals (0 = off), and only regions
    /// whose narrow-band score reaches `score_prefilter_frac` of the
    /// reporting floor are scored in the full band (heuristic: an alignment
    /// that needs the wide band for most of its score is dropped).
    pub score_prefilter_band: i64,
    pub score_prefilter_frac: f64,
    /// Candidate regions kept per query (heuristic limit; drops are counted).
    pub max_candidates: usize,
    pub max_regions_per_subject: usize,
    /// Consecutive engine subjects that form one reported sequence (6 for a
    /// translated index: the six frames of one nucleotide sequence). `-k`
    /// and `--max-hsps` count sequences, not frames.
    pub subject_group: u64,
    /// Skip seeds occurring more often than this in the database (0 = auto).
    pub seed_max_occ: usize,
    pub query_masking: MaskingMode,
    pub comp_bias: bool,
    /// Hit-buffer capacity per worker (hits; 8 bytes each).
    pub hit_buffer: usize,
    /// Direction-matrix cap per traceback (bytes).
    pub traceback_max_bytes: usize,
    /// Reference mode: score every subject, no seeding.
    pub exhaustive: bool,
    /// Build the gapped alignment strings (query, subject, midline) of each
    /// hit. Only pairwise (m0) and SAM output need them; statistics, CIGAR
    /// and coordinates are always computed.
    pub alignment_strings: bool,
    /// Neighborhood-word threshold T (full 20-letter seed alphabet only):
    /// every query seed is expanded to all words scoring >= T against it
    /// (BLAST-style). 0 = exact seed matches only.
    pub neighborhood_t: i32,
    /// At most this many neighbor words per query seed position (the
    /// highest-scoring ones are kept); bounds neighborhood expansion.
    pub neighborhood_max_words: usize,
    /// Memory budget of the shared neighborhood cache (bytes; 0 = no cache).
    pub neighborhood_cache_bytes: usize,
    /// Scan schedule: a (query, subject) pair whose window-filter survivors
    /// would exceed this many is scanned and visited on its own, in slices of
    /// query positions of about this many survivors each, so that its
    /// survivors are never held at once (repeat-rich pairs with masking and
    /// the seed-frequency cap disabled). Same results; 0 = default (2^22).
    pub scan_pair_survivors: usize,
}

impl Default for SearchParams {
    fn default() -> Self {
        Self {
            evalue: 10.0,
            min_bitscore: None,
            max_target_seqs: 10,
            max_hsps: 1,
            min_identity: 0.0,
            min_query_cover: 0.0,
            min_subject_cover: 0.0,
            hsp_threshold: 40,
            hsp_floor_fraction: 0.0,
            x_drop_ungapped: 20,
            two_hit_window: 0,
            seed_min_score: 0,
            single_hit_score: 0,
            hit_filter_min_ids: 0,
            hit_filter_exact_weight: 0,
            band_pad: 16,
            window_slack: 64,
            score_band: 0,
            score_prefilter_band: 0,
            score_prefilter_frac: 0.75,
            max_candidates: 5000,
            max_regions_per_subject: 4,
            subject_group: 1,
            seed_max_occ: 0,
            query_masking: MaskingMode::Tantan,
            comp_bias: true,
            hit_buffer: 1 << 18,
            traceback_max_bytes: 32 << 20,
            exhaustive: false,
            alignment_strings: true,
            neighborhood_t: 0,
            neighborhood_max_words: 256,
            neighborhood_cache_bytes: 64 << 20,
            scan_pair_survivors: 0,
        }
    }
}

/// Work and timing counters, summed over queries and threads.
#[derive(Debug, Default, Clone, Copy, serde::Serialize)]
pub struct Counters {
    pub queries: u64,
    pub query_residues: u64,
    pub query_residues_masked: u64,
    pub seed_positions: u64,
    pub seeds_frequent_skipped: u64,
    pub postings_visited: u64,
    pub hits: u64,
    pub subject_range_passes: u64,
    pub streamed_subjects: u64,
    pub ungapped_extensions: u64,
    pub extensions_avoided: u64,
    pub two_hit_waits: u64,
    pub seeds_low_score: u64,
    pub hits_filtered: u64,
    pub hsps: u64,
    pub regions: u64,
    pub regions_dropped_by_cap: u64,
    pub regions_below_evalue: u64,
    pub simd_cells: u64,
    pub scalar_cells: u64,
    pub saturation_retries: u64,
    pub banded_cells: u64,
    pub band_widenings: u64,
    pub band_capped: u64,
    pub tracebacks: u64,
    pub filtered_identity_coverage: u64,
    pub alignments: u64,
    pub t_prep_ns: u64,
    pub t_seed_ns: u64,
    pub t_extend_ns: u64,
    pub t_score_ns: u64,
    pub t_traceback_ns: u64,
    /// Extra seeds generated by neighborhood expansion (--neighborhood-t).
    pub neighbor_seeds: u64,
    /// Windows scored by an external score backend (GPU).
    pub backend_windows: u64,
    /// Backend calls that failed and were re-scored on the CPU.
    pub backend_fallbacks: u64,
    /// Neighborhood lists served from / added to the shared cache.
    pub neighbor_cache_hits: u64,
    pub neighbor_cache_misses: u64,
    /// Neighbor words dropped by `neighborhood_max_words`.
    pub neighbor_words_capped: u64,
    /// Candidate regions scored with banded DP (long sequence pairs).
    pub band_scored_regions: u64,
    /// Regions with a single weak HSP below the floor-relative threshold.
    pub regions_weak_dropped: u64,
    /// Tracebacks that hit the memory cap over the region's diagonal range and
    /// were redone in a band around the optimum's end cell.
    pub end_anchored_tracebacks: u64,
    /// Windows/bands whose 8-bit first-pass score reached the 8-bit limit and
    /// were re-scored exactly with wider lanes.
    pub u8_rescored: u64,
    /// Scan schedule, CPU cycles (TSC) per stage: seed codes and table
    /// lookups, table ranges and prefetches, window filter, survivor visits.
    pub c_scan_cand: u64,
    pub c_scan_range: u64,
    pub c_scan_filter: u64,
    pub c_scan_visit: u64,
    /// Two-stage scoring: regions scored in the full band after the narrow pass.
    pub band_rescored: u64,
    /// Within `c_scan_visit`: the per-subject sort and grouping.
    pub c_scan_sort: u64,
    /// Expanding the entry runs into hits (not part of `c_scan_range`).
    pub c_scan_expand: u64,
    /// Scan schedule: database positions whose seed code is in the batch.
    pub scan_candidates: u64,
    /// Scan schedule: X-drop steps (residue pairs scored) and (subject, query) groups visited.
    pub xdrop_steps: u64,
    pub scan_groups: u64,
    /// Scan schedule: work items cut (by chunks, subjects or queries) to
    /// bound the window-filter survivors held in memory.
    pub scan_splits: u64,
    /// Scan schedule: (query, subject) pairs scanned and visited on their own
    /// in query-position slices (see `SearchParams::scan_pair_survivors`).
    pub scan_pairs_streamed: u64,
}

impl Counters {
    pub fn sub(&mut self, o: &Counters) {
        macro_rules! d { ($($f:ident),*) => { $( self.$f -= o.$f; )* } }
        d!(queries, query_residues, query_residues_masked, seed_positions, seeds_frequent_skipped,
           postings_visited, hits, subject_range_passes, streamed_subjects, ungapped_extensions,
           extensions_avoided, two_hit_waits, seeds_low_score, hits_filtered, hsps, regions, regions_dropped_by_cap, regions_below_evalue,
           simd_cells, scalar_cells, saturation_retries, banded_cells, band_widenings, band_capped,
           tracebacks, filtered_identity_coverage, alignments, t_prep_ns, t_seed_ns, t_extend_ns,
           t_score_ns, t_traceback_ns, neighbor_seeds, backend_windows, backend_fallbacks, neighbor_cache_hits,
           neighbor_cache_misses, neighbor_words_capped, band_scored_regions, regions_weak_dropped,
           end_anchored_tracebacks, u8_rescored, c_scan_cand, c_scan_range, c_scan_filter, c_scan_visit, c_scan_sort, c_scan_expand, scan_candidates, xdrop_steps, scan_groups, band_rescored, scan_splits, scan_pairs_streamed);
    }

    pub fn add(&mut self, o: &Counters) {
        macro_rules! a { ($($f:ident),*) => { $( self.$f += o.$f; )* } }
        a!(queries, query_residues, query_residues_masked, seed_positions, seeds_frequent_skipped,
           postings_visited, hits, subject_range_passes, streamed_subjects, ungapped_extensions,
           extensions_avoided, two_hit_waits, seeds_low_score, hits_filtered, hsps, regions, regions_dropped_by_cap, regions_below_evalue,
           simd_cells, scalar_cells, saturation_retries, banded_cells, band_widenings, band_capped,
           tracebacks, filtered_identity_coverage, alignments, t_prep_ns, t_seed_ns, t_extend_ns,
           t_score_ns, t_traceback_ns, neighbor_seeds, backend_windows, backend_fallbacks, neighbor_cache_hits,
           neighbor_cache_misses, neighbor_words_capped, band_scored_regions, regions_weak_dropped,
           end_anchored_tracebacks, u8_rescored, c_scan_cand, c_scan_range, c_scan_filter, c_scan_visit, c_scan_sort, c_scan_expand, scan_candidates, xdrop_steps, scan_groups, band_rescored, scan_splits, scan_pairs_streamed);
    }
}

/// One reported alignment (coordinates 1-based, inclusive).
#[derive(Debug, Clone, serde::Serialize)]
pub struct Hit {
    pub subject: u64,
    pub score: i32,
    pub bitscore: f64,
    pub evalue: f64,
    pub q_start: usize,
    pub q_end: usize,
    pub s_start: usize,
    pub s_end: usize,
    pub columns: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub positives: usize,
    pub gap_opens: usize,
    pub gap_columns: usize,
    pub q_len: usize,
    pub s_len: usize,
    /// SAM-convention ops relative to the subject.
    pub ops: Vec<(u8, u32)>,
    pub q_aln: String,
    pub s_aln: String,
    pub midline: String,
    /// False when the traceback band hit the memory cap before reproducing
    /// the window score (explicit fallback; counted in `band_capped`).
    pub exact_in_window: bool,
}

impl Hit {
    pub fn identity(&self) -> f64 {
        if self.columns == 0 { 0.0 } else { 100.0 * self.matches as f64 / self.columns as f64 }
    }
    pub fn query_cover(&self) -> f64 {
        if self.q_len == 0 { 0.0 } else { 100.0 * (self.q_end + 1 - self.q_start) as f64 / self.q_len as f64 }
    }
    pub fn subject_cover(&self) -> f64 {
        if self.s_len == 0 { 0.0 } else { 100.0 * (self.s_end + 1 - self.s_start) as f64 / self.s_len as f64 }
    }
    pub fn cigar(&self) -> String {
        let mut s = String::new();
        for &(op, l) in &self.ops {
            s.push_str(&l.to_string());
            s.push(op as char);
        }
        s
    }
}

/// Worker-local scratch buffers (one per thread; never shared).
pub struct Worker {
    pub counters: Counters,
    simd: SimdScratch,
    band: BandScratch,
    band8: crate::algorithm::band_simd8::Band8Scratch,
    tantan: TantanScratch,
    mask: Vec<bool>,
    bias: Vec<i8>,
    // seeding
    seeds: Vec<(u32, u16, u32)>, // (qpos, shape, code)
    counts: Vec<u32>,
    touched: Vec<u16>,
    starts: Vec<u32>,
    hits: Vec<(u32, u32)>, // (qpos, spos) grouped by subject
    /// Hits that passed the window filter, in seed order: (qpos, posting).
    surv: Vec<(u32, u32)>,
    diag_epoch: Vec<u32>,
    diag_cover: Vec<u32>,
    diag_last: Vec<u32>,
    epoch: u32,
    hsps: Vec<Hsp>,
    regions: Vec<Region>,
    /// Query in the seed alphabet, padded by HF_PAD sentinels each side.
    qred: Vec<u8>,
    /// Query in the alignment encoding with the same padding/sentinels as
    /// `qred` (exact-identity part of the window filter).
    qraw: Vec<u8>,
    /// Ungapped scoring profile: prof[i*32 + c] = S(q_i, c) + bias_i.
    prof: Vec<i16>,
    sbuf: [u8; 32],
    /// HSP threshold of the current query (ungapped segments below it are
    /// not kept).
    hsp_thr: i32,
    /// Region threshold of the current query: a region is kept if its best
    /// HSP reaches it, or if it has ≥ 2 HSPs whose scores sum to it.
    region_thr: i32,
}

/// Candidate windows whose rectangle (query length × subject window) exceeds
/// this many cells are scored with banded DP around the seed diagonals when
/// that is at least `BAND_COST_RATIO` times cheaper (long sequence pairs).
const BAND_SCORE_MIN_CELLS: u64 = 1 << 20;
const BAND_COST_RATIO: u64 = 12;
/// Regions of one subject closer than this many diagonals share a band.
const BAND_BRIDGE: i64 = 1024;
/// Band-scratch buffers above this size are released after a traceback.
const BAND_SCRATCH_KEEP: usize = 8 << 20;

/// Subjects ahead to prefetch while visiting hits.
const PREFETCH_AHEAD: usize = 4;

/// Padding around the reduced query (hit filter window half-width).
const HF_PAD: usize = 32;

impl Default for Worker {
    fn default() -> Self {
        Self::new()
    }
}

impl Worker {
    pub fn new() -> Self {
        Self {
            counters: Counters::default(),
            simd: SimdScratch::new(),
            band: BandScratch::new(),
            band8: crate::algorithm::band_simd8::Band8Scratch::new(),
            tantan: TantanScratch::default(),
            mask: Vec::new(),
            bias: Vec::new(),
            seeds: Vec::new(),
            counts: vec![0; 65_536],
            touched: Vec::new(),
            starts: vec![0; 65_537],
            hits: Vec::new(),
            surv: Vec::new(),
            diag_epoch: Vec::new(),
            diag_cover: Vec::new(),
            diag_last: Vec::new(),
            epoch: 0,
            hsps: Vec::new(),
            regions: Vec::new(),
            qred: Vec::new(),
            qraw: Vec::new(),
            prof: Vec::new(),
            sbuf: [0xFF; 32],
            hsp_thr: 0,
            region_thr: 0,
        }
    }

    /// Approximate bytes held by this worker's scratch space.
    pub fn scratch_bytes(&self) -> usize {
        self.hits.capacity() * 8 + self.surv.capacity() * 8 + self.counts.capacity() * 4 + self.starts.capacity() * 4
            + self.diag_epoch.capacity() * 12 + self.seeds.capacity() * 12 + self.mask.capacity()
    }
}

#[derive(Debug, Clone, Copy)]
#[allow(dead_code)] // q_lo/q_hi are kept for debugging output
struct Hsp {
    diag: i64,
    q_lo: u32,
    q_hi: u32,
    score: i32,
}

#[derive(Debug, Clone, Copy)]
struct Region {
    subject: u64,
    d_min: i64,
    d_max: i64,
    best: i32,
}

/// Where a candidate region is scored and traced back: the subject window
/// `[lo, hi)` and, for long sequence pairs, the diagonal band inside it
/// (relative to the window). `band == None` = the whole window rectangle.
#[derive(Debug, Clone, Copy)]
struct Plan {
    lo: usize,
    hi: usize,
    band: Option<Band>,
}

/// Bounded min-heap of candidate regions. Elements are ordered so that the
/// heap minimum is the *weakest* region under the total order used for
/// ranking everywhere: score desc, subject asc, d_min asc, d_max asc.
/// Top-N under a total order does not depend on how the database was split,
/// so sequential and parallel searches keep exactly the same regions.
type RegionKey = std::cmp::Reverse<(i32, std::cmp::Reverse<u64>, std::cmp::Reverse<i64>, std::cmp::Reverse<i64>)>;
type RegionHeap = std::collections::BinaryHeap<RegionKey>;

fn region_order(a: &Region, b: &Region) -> std::cmp::Ordering {
    b.best.cmp(&a.best).then(a.subject.cmp(&b.subject)).then(a.d_min.cmp(&b.d_min)).then(a.d_max.cmp(&b.d_max))
}

/// Per-query state shared by all stages (and all threads in
/// [`Engine::search_parallel`]).
pub struct QueryPrep {
    pub bias: Vec<i8>,
    pub mask: Vec<bool>,
    pub has_bias: bool,
    pub qstats: QueryStatistics,
    /// Exact raw-score floor implied by -E / --min-bitscore.
    pub floor: i64,
    /// The floor the region threshold is relative to: at
    /// [`REGION_FLOOR_EVALUE`] or the search's -E, whichever is larger.
    pub region_floor: i64,
}

/// The region threshold (`hsp_floor_fraction`) and the two-stage scoring
/// threshold (`score_prefilter_frac`) are fractions of the query's reporting
/// floor at this E-value or at the search's own, whichever is more
/// permissive: a more stringent `-E` (or `--min-bitscore`) then drops no more
/// candidate regions than the setting the presets were tuned at (`-E 10`), so
/// it only filters what an `-E 10` search reports.
pub const REGION_FLOOR_EVALUE: f64 = 10.0;

/// Search engine bound to one index and one scoring configuration.
pub struct Engine<'a> {
    pub index: &'a SeedIndex,
    pub matrix: ScoreMatrix,
    pub gaps: GapCosts,
    pub stats: SearchStatistics,
    pub params: SearchParams,
    tantan: Option<Tantan>,
    occ_cap: Vec<u64>,
    is_protein: bool,
    /// Two-hit: minimum query distance between the two hits (non-overlap).
    two_hit_min: u32,
    /// Span used for seed scoring (shortest shape span).
    seed_span: usize,
    avx2: bool,
    /// Seed-alphabet map as two 16-byte shuffle tables (codes 0..15, 16..23).
    red_lo: [u8; 16],
    red_hi: [u8; 16],
    /// Per-shape neighborhood generators (empty when T = 0).
    nbr: Vec<Option<crate::algorithm::neighborhood::NeighborhoodGenerator>>,
    /// Shared cache of neighbor lists keyed by (shape, word code), bounded by
    /// `neighborhood_cache_bytes` (entries are no longer added once full).
    nbr_cache: std::sync::RwLock<(std::collections::HashMap<(u16, u32), std::sync::Arc<[u32]>>, usize)>,
    /// Optional external backend for the region-scoring stage (GPU).
    backend: Option<std::sync::Arc<dyn crate::gpu::ScoreBackend>>,
    /// Seeding scheme of the scan schedule ([`Self::search_scan`]).
    scan_seeds: Option<ScanSeeds>,
    /// Database-side scan data, built on first use.
    scan_db: std::sync::OnceLock<scan::ScanDb<'a>>,
    /// Substitution scores `mflat[a * 32 + c]` (c >= 24 scored as X).
    mflat: [i32; 768],
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt,bmi1")]
unsafe fn ids_avx2(q: &[u8], s: &[u8], lo: &[u8; 16], hi: &[u8; 16]) -> u32 {
    use std::arch::x86_64::*;
    let sv = _mm256_loadu_si256(s.as_ptr() as *const __m256i);
    let qv = _mm256_loadu_si256(q.as_ptr() as *const __m256i);
    let tlo = _mm256_broadcastsi128_si256(_mm_loadu_si128(lo.as_ptr() as *const __m128i));
    let thi = _mm256_broadcastsi128_si256(_mm_loadu_si128(hi.as_ptr() as *const __m128i));
    let idx = _mm256_and_si256(sv, _mm256_set1_epi8(0x0f));
    let ishi = _mm256_cmpgt_epi8(sv, _mm256_set1_epi8(15));
    let red = _mm256_blendv_epi8(_mm256_shuffle_epi8(tlo, idx), _mm256_shuffle_epi8(thi, idx), ishi);
    let eq = _mm256_cmpeq_epi8(red, qv);
    (_mm256_movemask_epi8(eq) as u32).count_ones()
}

impl<'a> Engine<'a> {
    pub fn new(index: &'a SeedIndex, matrix: ScoreMatrix, gaps: GapCosts, stats: SearchStatistics, params: SearchParams) -> Self {
        let is_protein = index.kind() == SeqKind::Protein;
        let tantan = if is_protein && params.query_masking == MaskingMode::Tantan {
            Some(Tantan::protein(&ScoreMatrix::new(crate::scoring::ScoringSystem::Blosum62)))
        } else {
            None
        };
        let cfg = index.seed_config();
        let occ_cap = (0..cfg.shapes.len())
            .map(|s| {
                if params.seed_max_occ > 0 {
                    params.seed_max_occ as u64
                } else {
                    // Auto: 64× the mean bucket size, at least OCC_CAP_FLOOR.
                    super::index::auto_occ_cap(index.meta().total_postings.get(s).copied().unwrap_or(0), cfg.num_codes(s))
                }
            })
            .collect();
        let two_hit_min = cfg.shapes.iter().map(|s| s.span as u32).min().unwrap_or(1);
        let seed_span = cfg.shapes.iter().map(|s| s.span).min().unwrap_or(4);
        let mut red_lo = [0xFFu8; 16];
        let mut red_hi = [0xFFu8; 16];
        for c in 0..24usize {
            let r = cfg.alphabet.map[c];
            if c < 16 { red_lo[c] = r } else { red_hi[c - 16] = r }
        }
        #[cfg(target_arch = "x86_64")]
        let avx2 = is_x86_feature_detected!("avx2") && is_x86_feature_detected!("popcnt") && is_x86_feature_detected!("bmi1");
        #[cfg(not(target_arch = "x86_64"))]
        let avx2 = false;
        let nbr = if params.neighborhood_t > 0 && is_protein && cfg.alphabet.size == 20 && cfg.alphabet.name == "full" {
            cfg.shapes
                .iter()
                .map(|sh| {
                    let w = sh.weight();
                    if (2..=8).contains(&w) {
                        Some(crate::algorithm::neighborhood::NeighborhoodGenerator::new(&matrix, w, params.neighborhood_t))
                    } else {
                        None
                    }
                })
                .collect()
        } else {
            Vec::new()
        };
        let mut mflat = [0i32; 768];
        for a in 0..24 {
            for c in 0..32 {
                mflat[a * 32 + c] = matrix.matrix[a][if c < 24 { c } else { crate::alphabet::AA_X as usize }];
            }
        }
        Self {
            index,
            matrix,
            gaps,
            stats,
            params,
            tantan,
            occ_cap,
            is_protein,
            two_hit_min,
            seed_span,
            avx2,
            red_lo,
            red_hi,
            nbr,
            nbr_cache: std::sync::RwLock::new((std::collections::HashMap::new(), 0)),
            backend: None,
            scan_seeds: None,
            scan_db: std::sync::OnceLock::new(),
            mflat,
        }
    }

    /// Use `backend` (e.g. a GPU) for the region-scoring stage. Scores must
    /// be exact; on a backend error the CPU kernel re-scores the batch.
    pub fn with_backend(mut self, backend: std::sync::Arc<dyn crate::gpu::ScoreBackend>) -> Self {
        self.backend = Some(backend);
        self
    }

    /// Neighbor words of the seed `code` at query position `p` (bounded by
    /// `neighborhood_max_words`, highest-scoring first), via the shared cache.
    #[allow(clippy::too_many_arguments)]
    fn neighbors(
        &self,
        g: &crate::algorithm::neighborhood::NeighborhoodGenerator,
        shape_idx: u16,
        code: u32,
        shape: &crate::alphabet::SeedShape,
        query: &[u8],
        p: usize,
        word: &mut Vec<u8>,
        nb: &mut Vec<u32>,
        ctr: &mut Counters,
    ) -> std::sync::Arc<[u32]> {
        if self.params.neighborhood_cache_bytes > 0 {
            if let Ok(c) = self.nbr_cache.read() {
                if let Some(v) = c.0.get(&(shape_idx, code)) {
                    ctr.neighbor_cache_hits += 1;
                    return v.clone();
                }
            }
        }
        ctr.neighbor_cache_misses += 1;
        word.clear();
        word.extend(shape.offsets.iter().map(|&o| query[p + o as usize]));
        nb.clear();
        let alpha = self.index.seed_config().alphabet.size as u32;
        g.expand(word, alpha, nb);
        let cap = self.params.neighborhood_max_words.max(1);
        if nb.len() > cap {
            // keep the highest-scoring neighbors (score each word against the query word)
            let score = |code: u32| -> i32 {
                let mut c = code;
                let mut t = 0;
                for &qc in word.iter().rev() {
                    let r = (c % alpha) as usize;
                    c /= alpha;
                    t += self.matrix.matrix[code_index(qc)][r];
                }
                t
            };
            nb.sort_by_cached_key(|&n| (std::cmp::Reverse(score(n)), n));
            ctr.neighbor_words_capped += (nb.len() - cap) as u64;
            nb.truncate(cap);
        }
        let v: std::sync::Arc<[u32]> = std::sync::Arc::from(nb.as_slice());
        if self.params.neighborhood_cache_bytes > 0 {
            if let Ok(mut c) = self.nbr_cache.write() {
                let bytes = v.len() * 4 + 48;
                if c.1 + bytes <= self.params.neighborhood_cache_bytes {
                    c.1 += bytes;
                    c.0.insert((shape_idx, code), v.clone());
                }
            }
        }
        v
    }

    /// Effective ungapped HSP threshold for a query (see `hsp_floor_fraction`).
    fn hsp_threshold_for(&self, prep: &QueryPrep) -> i32 {
        let f = self.params.hsp_floor_fraction;
        if f > 0.0 && prep.region_floor > 0 {
            self.params.hsp_threshold.max((prep.region_floor as f64 * f).ceil().min(i32::MAX as f64) as i32)
        } else {
            self.params.hsp_threshold
        }
    }

    /// Disable the AVX2 hit-visiting path (testing: the scalar path must give
    /// identical results).
    #[doc(hidden)]
    pub fn force_scalar_hit_path(&mut self) {
        self.avx2 = false;
    }

    /// True when neighborhood expansion is active.
    pub fn neighborhood_active(&self) -> bool {
        self.nbr.iter().any(|g| g.is_some())
    }

    pub fn scheme(&self) -> Scheme<'_> {
        Scheme { matrix: &self.matrix.matrix, gap_open: self.gaps.open, gap_extend: self.gaps.extend }
    }

    pub fn occ_caps(&self) -> &[u64] {
        &self.occ_cap
    }

    /// Per-query preparation: composition bias, tantan seed mask, statistics.
    pub fn prepare(&self, query: &[u8], w: &mut Worker) -> QueryPrep {
        let t0 = Instant::now();
        let qlen = query.len();
        let mut bias = Vec::new();
        if self.params.comp_bias && self.is_protein {
            composition_bias(query, &self.matrix, &mut bias);
        } else {
            bias.resize(qlen, 0);
        }
        let has_bias = bias.iter().any(|&b| b != 0);
        let mut mask = Vec::new();
        if let Some(t) = &self.tantan {
            let n = t.mask(query, &mut w.tantan, &mut mask);
            w.counters.query_residues_masked += n as u64;
        } else {
            mask.resize(qlen, false);
        }
        let qstats = self.stats.for_query(qlen.max(1));
        let floor = qstats.min_score(self.params.evalue, self.params.min_bitscore);
        let region_floor = if self.params.evalue >= REGION_FLOOR_EVALUE && self.params.min_bitscore.is_none() {
            floor
        } else {
            qstats.min_score(self.params.evalue.max(REGION_FLOOR_EVALUE), None)
        };
        w.counters.t_prep_ns += t0.elapsed().as_nanos() as u64;
        QueryPrep { bias, mask, has_bias, qstats, floor, region_floor }
    }

    /// Search one encoded query on the calling thread.
    /// `query` must be in the index's encoding.
    pub fn search(&self, query: &[u8], w: &mut Worker) -> Vec<Hit> {
        w.counters.queries += 1;
        w.counters.query_residues += query.len() as u64;
        if query.is_empty() || self.index.num_subjects() == 0 {
            return Vec::new();
        }
        let prep = self.prepare(query, w);
        let mut regions = std::mem::take(&mut w.regions);
        regions.clear();
        if self.params.exhaustive {
            self.all_regions(query.len(), &mut regions);
        } else {
            let ranges: Vec<(usize, u32, u32)> =
                (0..self.index.num_partitions()).map(|p| (p, 0u32, self.index.partition(p).num_subjects)).collect();
            self.collect_regions(query, &prep, w, &ranges, &mut regions);
        }
        w.counters.regions += regions.len() as u64;
        let t2 = Instant::now();
        let windows = self.plan_windows(&regions, query.len());
        let scores = self.score_windows(query, &prep, &regions, &windows, w);
        w.counters.t_score_ns += t2.elapsed().as_nanos() as u64;
        let out = self.finish(query, &prep, &regions, &windows, &scores, w);
        w.regions = regions;
        out
    }

    /// Search one query using all threads of the current rayon pool: the
    /// database is split into subject ranges searched in parallel, region
    /// scoring is chunked across threads, traceback runs on the caller.
    /// Results are identical to [`Self::search`].
    pub fn search_parallel(&self, query: &[u8]) -> (Vec<Hit>, Counters) {
        let mut main = Worker::new();
        main.counters.queries += 1;
        main.counters.query_residues += query.len() as u64;
        if query.is_empty() || self.index.num_subjects() == 0 {
            return (Vec::new(), main.counters);
        }
        let prep = self.prepare(query, &mut main);
        let mut regions: Vec<Region> = Vec::new();
        let mut ctr = Counters::default();
        if self.params.exhaustive {
            self.all_regions(query.len(), &mut regions);
        } else {
            let threads = rayon::current_num_threads().max(1);
            self.build_seeds(query, &prep, &mut main);
            let tasks = self.plan_tasks(threads);
            let seed_src = &main;
            let parts: Vec<(Vec<Region>, Counters)> = tasks
                .par_iter()
                .map_init(Worker::new, |w, ranges| {
                    let before = w.counters;
                    Self::copy_seed_state(seed_src, w);
                    let mut out = Vec::new();
                    self.scan_ranges(query, w, ranges, &mut out);
                    let mut d = w.counters;
                    d.sub(&before);
                    (out, d)
                })
                .collect();
            for (r, c) in parts {
                regions.extend(r);
                ctr.add(&c);
            }
            regions.sort_unstable_by(region_order);
            let cap = self.params.max_candidates.max(1);
            if regions.len() > cap {
                ctr.regions_dropped_by_cap += (regions.len() - cap) as u64;
                regions.truncate(cap);
            }
        }
        main.counters.regions += regions.len() as u64;
        let t2 = Instant::now();
        let windows = self.plan_windows(&regions, query.len());
        let chunk = 64usize;
        let idxs: Vec<usize> = (0..regions.len()).step_by(chunk).collect();
        let scored: Vec<(Vec<i32>, Counters)> = idxs
            .par_iter()
            .map(|&st| {
                let en = (st + chunk).min(regions.len());
                let mut w = Worker::new();
                let sc = self.score_windows(query, &prep, &regions[st..en], &windows[st..en], &mut w);
                (sc, w.counters)
            })
            .collect();
        let mut scores = Vec::with_capacity(regions.len());
        for (sc, c) in scored {
            scores.extend(sc);
            ctr.add(&c);
        }
        main.counters.t_score_ns += t2.elapsed().as_nanos() as u64;
        let out = self.finish(query, &prep, &regions, &windows, &scores, &mut main);
        ctr.add(&main.counters);
        (out, ctr)
    }

    /// Search a batch of queries *partition-major*: every index partition is
    /// scanned once for the whole batch (all queries in parallel) before the
    /// next one is touched, so each partition's postings and residues are
    /// paged in / cached once per batch rather than once per query. This is
    /// the blocked "batch seed-join" schedule for large query sets against
    /// large (e.g. memory-mapped, larger-than-RAM) indexes. Per-query
    /// candidate lists stay bounded by `max_candidates`; results are
    /// identical to [`Self::search`] for every query.
    pub fn search_batch(&self, queries: &[&[u8]]) -> Vec<(Vec<Hit>, Counters)> {
        struct State {
            prep: QueryPrep,
            seeds: Vec<(u32, u16, u32)>,
            qred: Vec<u8>,
            qraw: Vec<u8>,
            regions: Vec<Region>,
            ctr: Counters,
        }
        let empty_db = self.index.num_subjects() == 0;
        let mut states: Vec<State> = queries
            .par_iter()
            .map_init(Worker::new, |w, q| {
                let before = w.counters;
                w.counters.queries += 1;
                w.counters.query_residues += q.len() as u64;
                let prep = self.prepare(q, w);
                let (seeds, qred, qraw) = if q.is_empty() || empty_db || self.params.exhaustive {
                    (Vec::new(), Vec::new(), Vec::new())
                } else {
                    self.load_prep(q, &prep, w);
                    self.compute_seeds(q, w);
                    (w.seeds.clone(), w.qred.clone(), w.qraw.clone())
                };
                let mut d = w.counters;
                d.sub(&before);
                State { prep, seeds, qred, qraw, regions: Vec::new(), ctr: d }
            })
            .collect();
        let cap = self.params.max_candidates.max(1);
        if !self.params.exhaustive && !empty_db {
            for p in 0..self.index.num_partitions() {
                let nsub = self.index.partition(p).num_subjects;
                states.par_iter_mut().zip(queries.par_iter()).for_each_init(Worker::new, |w, (st, q)| {
                    if st.seeds.is_empty() {
                        return;
                    }
                    let before = w.counters;
                    self.load_prep(q, &st.prep, w);
                    w.seeds.clone_from(&st.seeds);
                    w.qred.clone_from(&st.qred);
                    w.qraw.clone_from(&st.qraw);
                    let mut out = Vec::new();
                    self.scan_ranges(q, w, &[(p, 0, nsub)], &mut out);
                    st.regions.extend(out);
                    if st.regions.len() > 2 * cap {
                        st.regions.sort_unstable_by(region_order);
                        st.ctr.regions_dropped_by_cap += (st.regions.len() - cap) as u64;
                        st.regions.truncate(cap);
                    }
                    let mut d = w.counters;
                    d.sub(&before);
                    st.ctr.add(&d);
                });
            }
        }
        states
            .into_par_iter()
            .zip(queries.par_iter())
            .map_init(Worker::new, |w, (mut st, q)| {
                if q.is_empty() || empty_db {
                    return (Vec::new(), st.ctr);
                }
                let before = w.counters;
                let mut regions = std::mem::take(&mut st.regions);
                if self.params.exhaustive {
                    self.all_regions(q.len(), &mut regions);
                } else {
                    regions.sort_unstable_by(region_order);
                    if regions.len() > cap {
                        st.ctr.regions_dropped_by_cap += (regions.len() - cap) as u64;
                        regions.truncate(cap);
                    }
                }
                w.counters.regions += regions.len() as u64;
                let t2 = Instant::now();
                let windows = self.plan_windows(&regions, q.len());
                let scores = self.score_windows(q, &st.prep, &regions, &windows, w);
                w.counters.t_score_ns += t2.elapsed().as_nanos() as u64;
                let hits = self.finish(q, &st.prep, &regions, &windows, &scores, w);
                let mut d = w.counters;
                d.sub(&before);
                st.ctr.add(&d);
                (hits, st.ctr)
            })
            .collect()
    }

    /// Split the database into about `4 × threads` tasks: groups of whole
    /// partitions when there are many, subject sub-ranges otherwise.
    fn plan_tasks(&self, threads: usize) -> Vec<Vec<(usize, u32, u32)>> {
        let np = self.index.num_partitions();
        let target = (threads * 4).max(1);
        let mut tasks = Vec::new();
        if np == 0 {
            return tasks;
        }
        if np >= target {
            let per = (np + target - 1) / target;
            let mut p = 0;
            while p < np {
                let e = (p + per).min(np);
                tasks.push((p..e).map(|q| (q, 0u32, self.index.partition(q).num_subjects)).collect());
                p = e;
            }
        } else {
            let per_part = ((target + np - 1) / np) as u32;
            for p in 0..np {
                let n = self.index.partition(p).num_subjects;
                let step = ((n + per_part - 1) / per_part).max(1);
                let mut lo = 0u32;
                while lo < n {
                    let hi = (lo + step).min(n);
                    tasks.push(vec![(p, lo, hi)]);
                    lo = hi;
                }
            }
        }
        tasks
    }

    fn all_regions(&self, qlen: usize, regions: &mut Vec<Region>) {
        for s in 0..self.index.num_subjects() {
            let sl = self.index.subject_len(s) as i64;
            regions.push(Region { subject: s as u64, d_min: -(qlen as i64), d_max: sl, best: 0 });
        }
    }

    fn score_windows(&self, query: &[u8], prep: &QueryPrep, regions: &[Region], windows: &[Plan], w: &mut Worker) -> Vec<i32> {
        let bias: Option<&[i8]> = if prep.has_bias { Some(&prep.bias) } else { None };
        // Banded plans (narrow-band scoring, long sequence pairs): one
        // inter-sequence 8-bit pass over all of them, exact fallback for wide
        // bands and saturated lanes (same scores as `banded_score`).
        let mut out = vec![0i32; regions.len()];
        let mut simd_idx: Vec<usize> = Vec::with_capacity(regions.len());
        let mut band_idx: Vec<usize> = Vec::new();
        for (i, pl) in windows.iter().enumerate() {
            match pl.band {
                Some(_) => band_idx.push(i),
                None => simd_idx.push(i),
            }
        }
        if !band_idx.is_empty() {
            let mut kc = KernelCounters::default();
            let mut bc = crate::algorithm::band_simd8::Band8Counters::default();
            let wn = self.params.score_prefilter_band;
            let full_band = self.params.score_band;
            // Two-stage: the narrow band first, the full band only for
            // regions that reach a fraction of their reporting floor.
            let mut full_idx: Vec<usize> = band_idx.clone();
            if wn > 0 && wn < full_band {
                let qs = &prep.qstats;
                let jobs: Vec<crate::algorithm::band_simd8::BandJob<'_>> = band_idx
                    .iter()
                    .map(|&i| {
                        let pl = windows[i];
                        let r = &regions[i];
                        let subj = &self.index.subject_codes(r.subject as usize)[pl.lo..pl.hi];
                        let band = Band { lo: r.d_min - wn - pl.lo as i64, hi: r.d_max + wn - pl.lo as i64 }.clip(query.len(), subj.len());
                        crate::algorithm::band_simd8::BandJob { subject: subj, band }
                    })
                    .collect();
                let sc = crate::algorithm::band_simd8::score_banded_many(query, bias, &jobs, self.scheme(), &mut w.band8, &mut kc, &mut bc);
                full_idx.clear();
                for (k, &i) in band_idx.iter().enumerate() {
                    let n = self.index.subject_len(regions[i].subject as usize) as usize;
                    // relative to the floor at -E or E = 10, whichever is
                    // lower (as the region threshold): a stringent -E only
                    // filters what an -E 10 search reports
                    let thr = self.params.score_prefilter_frac
                        * qs.floor_lower_bound(self.params.evalue.max(REGION_FLOOR_EVALUE), None, n, prep.region_floor);
                    if sc[k] as f64 >= thr {
                        full_idx.push(i);
                    } else {
                        out[i] = sc[k];
                    }
                }
                w.counters.band_rescored += full_idx.len() as u64;
            }
            let jobs: Vec<crate::algorithm::band_simd8::BandJob<'_>> = full_idx
                .iter()
                .map(|&i| {
                    let pl = windows[i];
                    crate::algorithm::band_simd8::BandJob {
                        subject: &self.index.subject_codes(regions[i].subject as usize)[pl.lo..pl.hi],
                        band: pl.band.expect("banded plan"),
                    }
                })
                .collect();
            let sc = crate::algorithm::band_simd8::score_banded_many(query, bias, &jobs, self.scheme(), &mut w.band8, &mut kc, &mut bc);
            w.counters.banded_cells += kc.simd_cells + kc.scalar_cells;
            w.counters.band_scored_regions += band_idx.len() as u64;
            w.counters.u8_rescored += bc.u8_overflows;
            for (k, &i) in full_idx.iter().enumerate() {
                out[i] = sc[k];
            }
        }
        if simd_idx.is_empty() {
            return out;
        }
        let sub_regions: Vec<Region> = simd_idx.iter().map(|&i| regions[i]).collect();
        let sub_windows: Vec<(usize, usize)> = simd_idx.iter().map(|&i| (windows[i].lo, windows[i].hi)).collect();
        let sc = self.score_windows_simd(query, bias, &sub_regions, &sub_windows, w);
        for (k, &i) in simd_idx.iter().enumerate() {
            out[i] = sc[k];
        }
        out
    }

    fn score_windows_simd(&self, query: &[u8], bias: Option<&[i8]>, regions: &[Region], windows: &[(usize, usize)], w: &mut Worker) -> Vec<i32> {
        let views: Vec<ProteinCodes<'_>> = regions
            .iter()
            .zip(windows.iter())
            .map(|(r, &(lo, hi))| ProteinCodes::new_unchecked(&self.index.subject_codes(r.subject as usize)[lo..hi]))
            .collect();
        if let Some(be) = &self.backend {
            let slices: Vec<&[u8]> = views.iter().map(|v| v.as_slice()).collect();
            match be.score_windows(query, bias, &slices, &self.matrix.matrix, self.gaps) {
                Ok(sc) if sc.len() == slices.len() => {
                    w.counters.backend_windows += slices.len() as u64;
                    return sc;
                }
                _ => w.counters.backend_fallbacks += 1,
            }
        }
        let mut kc = KernelCounters::default();
        let scores = score_many(ProteinCodes::new_unchecked(query), bias, &views, &self.matrix.matrix, self.gaps, &mut w.simd, &mut kc);
        w.counters.simd_cells += kc.simd_cells;
        w.counters.scalar_cells += kc.scalar_cells;
        w.counters.saturation_retries += kc.saturation_retries;
        w.counters.u8_rescored += kc.u8_rescored;
        scores
    }

    /// E-value floor, ranking and progressive traceback.
    fn finish(&self, query: &[u8], prep: &QueryPrep, regions: &[Region], windows: &[Plan], scores: &[i32], w: &mut Worker) -> Vec<Hit> {
        let floor = prep.floor;
        let qs = &prep.qstats;
        let ecut = if self.params.evalue > 0.0 && self.params.evalue.is_finite() { self.params.evalue } else { f64::INFINITY };
        let min_bits = self.params.min_bitscore;
        let bias: Option<&[i8]> = if prep.has_bias { Some(&prep.bias) } else { None };
        // E-value of every region that passes the cut-offs (the traceback
        // reproduces the region's score, so this is the hit's E-value), in
        // E-value order. Per-subject E-values depend on the subject length: a
        // one-logarithm lower bound rejects most regions before the exact
        // length adjustment of the pair is computed.
        let mut order: Vec<(f64, usize)> = Vec::new();
        for (i, r) in regions.iter().enumerate() {
            let s = scores[i] as i64;
            match qs.model {
                EvalueModel::Database => {
                    if s >= floor {
                        order.push((qs.evalue(s), i));
                    }
                }
                EvalueModel::PerSubject => {
                    if s <= 0 || min_bits.map_or(false, |b| qs.bitscore(s) < b) {
                        continue;
                    }
                    let n = self.index.subject_len(r.subject as usize) as usize;
                    if qs.evalue_in(s, qs.space_lower_bound(n)) > ecut {
                        continue;
                    }
                    let e = qs.evalue_for(s, n);
                    if e <= ecut {
                        order.push((e, i));
                    }
                }
            }
        }
        w.counters.regions_below_evalue += (regions.len() - order.len()) as u64;
        order.sort_unstable_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(scores[b.1].cmp(&scores[a.1]))
                .then(region_order(&regions[a.1], &regions[b.1]))
        });
        let t3 = Instant::now();
        let mut out: Vec<Hit> = Vec::new();
        let mut per_subject: std::collections::HashMap<u64, usize> = std::collections::HashMap::new();
        let mut reported_subjects = 0usize;
        for &(_, ri) in &order {
            let group = self.params.subject_group.max(1);
            if reported_subjects >= self.params.max_target_seqs && !per_subject.contains_key(&(regions[ri].subject / group)) {
                break;
            }
            let r = regions[ri];
            let cnt = per_subject.get(&(r.subject / group)).copied().unwrap_or(0);
            if cnt >= self.params.max_hsps {
                continue;
            }
            let Plan { lo, hi, band: limit } = windows[ri];
            let subj = &self.index.subject_codes(r.subject as usize)[lo..hi];
            let tb = self.traceback(query, subj, bias, &r, lo, scores[ri], limit, w);
            w.band.trim(BAND_SCRATCH_KEEP);
            let Some((al, exact)) = tb else { continue };
            let hit = self.make_hit(query, self.index.subject_codes(r.subject as usize), r.subject, lo, &al, exact, &prep.qstats);
            if cnt > 0 && out.iter().any(|h| h.subject == hit.subject && overlaps(h, &hit)) {
                continue;
            }
            let passes = match qs.model {
                EvalueModel::Database => hit.score as i64 >= floor,
                EvalueModel::PerSubject => hit.evalue <= ecut && min_bits.map_or(true, |b| hit.bitscore >= b),
            };
            if hit.identity() < self.params.min_identity
                || hit.query_cover() < self.params.min_query_cover
                || hit.subject_cover() < self.params.min_subject_cover
                || !passes
            {
                w.counters.filtered_identity_coverage += 1;
                continue;
            }
            if cnt == 0 {
                reported_subjects += 1;
            }
            per_subject.insert(r.subject / group, cnt + 1);
            out.push(hit);
        }
        w.counters.t_traceback_ns += t3.elapsed().as_nanos() as u64;
        out.sort_by(|a, b| {
            a.evalue
                .partial_cmp(&b.evalue)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(b.score.cmp(&a.score))
                .then(a.subject.cmp(&b.subject))
                .then(a.s_start.cmp(&b.s_start))
        });
        w.counters.alignments += out.len() as u64;
        out
    }

    fn window(&self, r: &Region, qlen: usize) -> (usize, usize) {
        let sl = self.index.subject_len(r.subject as usize) as i64;
        let slack = self.params.window_slack;
        let lo = (r.d_min - slack).max(0).min(sl);
        let hi = (qlen as i64 + r.d_max + slack).min(sl).max(lo);
        // A window covering most of the subject costs no more than the subject.
        if (hi - lo) * 4 >= sl * 3 {
            (0, sl as usize)
        } else {
            (lo as usize, hi as usize)
        }
    }

    /// Scoring/traceback plan for every candidate region of one query.
    ///
    /// Normally a region is scored against the subject window spanned by its
    /// diagonals ± `window_slack` with the inter-sequence SIMD kernel (the
    /// whole rectangle: query × window). For long sequence pairs, where that
    /// rectangle exceeds `BAND_SCORE_MIN_CELLS` and banded DP is at least
    /// `BAND_COST_RATIO`× cheaper, the region is scored and traced back inside
    /// the diagonal band spanned by the chain of that subject's regions within
    /// `BAND_BRIDGE` diagonals of each other, ± `window_slack` (so an indel of
    /// up to ~1,000 residues between two seeded parts stays inside the band).
    /// Depends only on the full region list (deterministic for every schedule).
    fn plan_windows(&self, regions: &[Region], qlen: usize) -> Vec<Plan> {
        // Chain each subject's regions whose diagonal ranges lie within
        // `BAND_BRIDGE` of each other (single linkage): a long indel between
        // two seeded parts of one alignment then falls inside one band, while
        // a distant spurious region does not widen it.
        let mut by_subject: std::collections::HashMap<u64, Vec<usize>> = std::collections::HashMap::new();
        for (i, r) in regions.iter().enumerate() {
            by_subject.entry(r.subject).or_default().push(i);
        }
        let mut span: Vec<(i64, i64)> = regions.iter().map(|r| (r.d_min, r.d_max)).collect();
        for idx in by_subject.values_mut() {
            // Only strong regions are bridged (a weak, typically spurious HSP
            // must not widen the band of a real alignment).
            let top = idx.iter().map(|&i| regions[i].best).max().unwrap_or(0);
            let strong = (2 * self.params.hsp_threshold).max(top / 4);
            idx.retain(|&i| regions[i].best >= strong);
            idx.sort_unstable_by_key(|&i| (regions[i].d_min, regions[i].d_max, i));
            let mut start = 0;
            while start < idx.len() {
                let (mut lo, mut hi) = (regions[idx[start]].d_min, regions[idx[start]].d_max);
                let mut end = start + 1;
                while end < idx.len() && regions[idx[end]].d_min - hi <= BAND_BRIDGE {
                    hi = hi.max(regions[idx[end]].d_max);
                    end += 1;
                }
                lo = lo.min(regions[idx[start]].d_min);
                for &i in &idx[start..end] {
                    span[i] = (lo, hi);
                }
                start = end;
            }
        }
        regions
            .iter()
            .enumerate()
            .map(|(ri, r)| {
                let (lo, hi) = self.window(r, qlen);
                let own = Plan { lo, hi, band: None };
                if self.params.exhaustive || r.best == 0 || self.params.window_slack >= i64::MAX / 8 {
                    return own;
                }
                if self.params.score_band > 0 {
                    // the region's diagonals, or those of its chain of strong
                    // regions (a long indel between two seeded parts of one
                    // alignment stays inside the band), ± score_band
                    let w = self.params.score_band;
                    let (dmin, dmax) = span[ri];
                    let sl = self.index.subject_len(r.subject as usize) as i64;
                    let lo = (dmin - w).max(0).min(sl);
                    let hi = (qlen as i64 + dmax + w).min(sl).max(lo);
                    let n = (hi - lo) as usize;
                    let band = Band { lo: dmin - w - lo, hi: dmax + w - lo }.clip(qlen, n);
                    if n == 0 || band.covers(qlen, n) {
                        return Plan { lo: lo as usize, hi: hi as usize, band: None };
                    }
                    return Plan { lo: lo as usize, hi: hi as usize, band: Some(band) };
                }
                if (qlen as u64) * ((hi - lo) as u64) < BAND_SCORE_MIN_CELLS {
                    return own;
                }
                let (dmin, dmax) = span[ri];
                let u = Region { subject: r.subject, d_min: dmin, d_max: dmax, best: r.best };
                let (ulo, uhi) = self.window(&u, qlen);
                let n = uhi - ulo;
                let band = Band { lo: dmin - ulo as i64, hi: dmax - ulo as i64 }.widen(self.params.window_slack).clip(qlen, n);
                if band.covers(qlen, n) {
                    return own;
                }
                let rows = qlen.min(n + band.width()) as u64;
                let rect = qlen as u64 * n as u64;
                if band.width() as u64 * rows * BAND_COST_RATIO < rect {
                    Plan { lo: ulo, hi: uhi, band: Some(band) }
                } else {
                    own
                }
            })
            .collect()
    }

    /// Banded traceback, widened until it reproduces the window score.
    #[allow(clippy::too_many_arguments)]
    fn traceback(
        &self,
        query: &[u8],
        subj: &[u8],
        bias: Option<&[i8]>,
        r: &Region,
        lo: usize,
        want: i32,
        limit: Option<Band>,
        w: &mut Worker,
    ) -> Option<(LocalAlignment, bool)> {
        let m = query.len();
        let n = subj.len();
        let mut pad = self.params.band_pad.max(4);
        let full = Band::full(m, n);
        let base = if self.params.exhaustive || r.best == 0 {
            full
        } else {
            Band { lo: r.d_min - lo as i64, hi: r.d_max - lo as i64 }
        };
        let mut best: Option<LocalAlignment> = None;
        loop {
            let mut band = if base == full { full } else { base.widen(pad) };
            let mut covers = band.covers(m, n);
            if let Some(l) = limit {
                // the scoring band bounds the search: never widen past it
                band = Band { lo: band.lo.max(l.lo), hi: band.hi.min(l.hi) };
                covers = covers || (band.lo <= l.lo && band.hi >= l.hi);
            }
            w.counters.tracebacks += 1;
            // `want` is the exact score of the window that contains the band,
            // so rows after the first cell reaching it need not be computed.
            match banded_align_until(query, subj, bias, self.scheme(), band, &mut w.band, self.params.traceback_max_bytes, want) {
                Ok(Some(al)) => {
                    w.counters.banded_cells += al.cells;
                    let done = al.score >= want || covers;
                    best = Some(al);
                    if done {
                        let exact = best.as_ref().map_or(false, |a| a.score >= want);
                        return best.map(|a| (a, exact));
                    }
                }
                Ok(None) => {
                    if covers {
                        return None;
                    }
                }
                Err(_) => {
                    // The band over the region's diagonal range needs more than
                    // the memory cap (a long pair whose region spans many
                    // diagonals, e.g. internal repeats). Locate the optimum's
                    // end cell with a score-only pass and trace back from a
                    // band around its diagonal instead.
                    if let Some(found) = self.traceback_from_end(query, subj, bias, want, limit, w) {
                        return Some(found);
                    }
                    // Still over the cap: report the best alignment found so
                    // far, flagged.
                    w.counters.band_capped += 1;
                    return best.map(|a| (a, false));
                }
            }
            w.counters.band_widenings += 1;
            pad *= 2;
        }
    }

    /// Fallback of [`Self::traceback`] when the band over the region's
    /// diagonals exceeds the memory cap: a linear-memory score pass over the
    /// scoring window (`limit`, or the whole window) gives the first optimal
    /// end cell, and the traceback band starts on that cell's diagonal and
    /// widens until it reproduces the window score. Returns `None` when this
    /// band also reaches the cap before the score is reproduced (or the
    /// window holds no positive alignment).
    #[allow(clippy::too_many_arguments)]
    fn traceback_from_end(
        &self,
        query: &[u8],
        subj: &[u8],
        bias: Option<&[i8]>,
        want: i32,
        limit: Option<Band>,
        w: &mut Worker,
    ) -> Option<(LocalAlignment, bool)> {
        let (m, n) = (query.len(), subj.len());
        let scope = limit.unwrap_or_else(|| Band::full(m, n));
        let end = banded_score(query, subj, bias, self.scheme(), scope, &mut w.band);
        w.counters.banded_cells += end.cells;
        if end.score <= 0 {
            return None;
        }
        w.counters.end_anchored_tracebacks += 1;
        let target = end.score.max(want);
        let d_end = end.s_end as i64 - end.q_end as i64;
        let mut pad = self.params.band_pad.max(4);
        loop {
            let band = Band { lo: (d_end - pad).max(scope.lo), hi: (d_end + pad).min(scope.hi) };
            let covers = band.lo <= scope.lo && band.hi >= scope.hi;
            w.counters.tracebacks += 1;
            match banded_align_until(query, subj, bias, self.scheme(), band, &mut w.band, self.params.traceback_max_bytes, target) {
                Ok(Some(al)) => {
                    w.counters.banded_cells += al.cells;
                    if al.score >= target || covers {
                        let exact = al.score >= want;
                        return Some((al, exact));
                    }
                }
                Ok(None) => {
                    if covers {
                        return None;
                    }
                }
                Err(_) => return None,
            }
            w.counters.band_widenings += 1;
            pad *= 2;
        }
    }

    fn make_hit(&self, q: &[u8], s_full: &[u8], subject: u64, s_off: usize, al: &LocalAlignment, exact: bool, qs: &QueryStatistics) -> Hit {
        let s = &s_full[s_off..];
        let strings = self.params.alignment_strings;
        let cap = if strings { al.columns() } else { 0 };
        let mut q_aln = String::with_capacity(cap);
        let mut s_aln = String::with_capacity(cap);
        let mut mid = String::with_capacity(cap);
        let (mut i, mut j) = (al.q_start, al.s_start);
        let (mut matches, mut mism, mut pos, mut gopen, mut gcols) = (0usize, 0usize, 0usize, 0usize, 0usize);
        let dec = |c: u8| -> char {
            if self.is_protein { decode_protein(c) as char } else { NUCLEOTIDE_ALPHABET[(c as usize).min(4)] as char }
        };
        for &(op, l) in &al.ops {
            match op {
                b'M' => {
                    for _ in 0..l {
                        let (a, b) = (q[i], s[j]);
                        let identical = a == b && (self.is_protein && a < 20 || !self.is_protein && a < 4);
                        let positive = identical || self.matrix.matrix[code_index(a)][code_index(b)] > 0;
                        if identical {
                            matches += 1;
                        } else {
                            mism += 1;
                        }
                        pos += positive as usize;
                        if strings {
                            q_aln.push(dec(a));
                            s_aln.push(dec(b));
                            mid.push(if identical { dec(a) } else if positive { '+' } else { ' ' });
                        }
                        i += 1;
                        j += 1;
                    }
                }
                b'I' => {
                    gopen += 1;
                    gcols += l as usize;
                    if strings {
                        for k in 0..l as usize {
                            q_aln.push(dec(q[i + k]));
                            s_aln.push('-');
                            mid.push(' ');
                        }
                    }
                    i += l as usize;
                }
                _ => {
                    gopen += 1;
                    gcols += l as usize;
                    if strings {
                        for k in 0..l as usize {
                            q_aln.push('-');
                            s_aln.push(dec(s[j + k]));
                            mid.push(' ');
                        }
                    }
                    j += l as usize;
                }
            }
        }
        Hit {
            subject,
            score: al.score,
            bitscore: qs.bitscore(al.score as i64),
            evalue: qs.evalue_for(al.score as i64, s_full.len()),
            q_start: al.q_start + 1,
            q_end: al.q_end,
            s_start: al.s_start + s_off + 1,
            s_end: al.s_end + s_off,
            columns: al.columns(),
            matches,
            mismatches: mism,
            positives: pos,
            gap_opens: gopen,
            gap_columns: gcols,
            q_len: q.len(),
            s_len: s_full.len(),
            ops: al.ops.clone(),
            q_aln,
            s_aln,
            midline: mid,
            exact_in_window: exact,
        }
    }

    // -----------------------------------------------------------------------
    // Seeding and candidate regions
    // -----------------------------------------------------------------------

    fn collect_regions(&self, query: &[u8], prep: &QueryPrep, w: &mut Worker, ranges: &[(usize, u32, u32)], out: &mut Vec<Region>) {
        self.build_seeds(query, prep, w);
        self.scan_ranges(query, w, ranges, out);
    }

    /// Copy the per-query seeding state of `from` into `to` (parallel search).
    fn copy_seed_state(from: &Worker, to: &mut Worker) {
        to.hsp_thr = from.hsp_thr;
        to.region_thr = from.region_thr;
        to.bias.clone_from(&from.bias);
        to.mask.clone_from(&from.mask);
        to.prof.clone_from(&from.prof);
        to.seeds.clone_from(&from.seeds);
        to.qred.clone_from(&from.qred);
        to.qraw.clone_from(&from.qraw);
    }

    /// Query profile, seed list (minus masked and over-frequent seeds) and
    /// reduced query for the window identity filter.
    fn build_seeds(&self, query: &[u8], prep: &QueryPrep, w: &mut Worker) {
        self.load_prep(query, prep, w);
        self.compute_seeds(query, w);
    }

    /// Copy bias and mask into the worker and build the ungapped profile.
    fn load_prep(&self, query: &[u8], prep: &QueryPrep, w: &mut Worker) {
        w.hsp_thr = self.params.hsp_threshold;
        w.region_thr = self.hsp_threshold_for(prep);
        w.bias.clear();
        w.bias.extend_from_slice(&prep.bias);
        w.mask.clear();
        w.mask.extend_from_slice(&prep.mask);
        w.prof.clear();
        w.prof.resize(query.len() * 32, 0);
        for (i, &qc) in query.iter().enumerate() {
            let row = &self.matrix.matrix[code_index(qc)];
            let b = w.bias[i] as i32;
            for c in 0..32 {
                w.prof[i * 32 + c] = (row[if c < 24 { c } else { crate::alphabet::AA_X as usize }] + b) as i16;
            }
        }
    }

    /// Seed list (minus masked and over-frequent seeds) and the reduced
    /// query for the window identity filter. Requires [`Self::load_prep`].
    fn compute_seeds(&self, query: &[u8], w: &mut Worker) {
        let t0 = Instant::now();
        let cfg = self.index.seed_config();
        let qlen = query.len();
        // Query seeds sorted by (qpos, shape): hits per subject then arrive
        // in query order, which the diagonal table relies on.
        w.seeds.clear();
        let mut word: Vec<u8> = Vec::new();
        let mut nb: Vec<u32> = Vec::new();
        let mut extra = 0u64;
        for (si, shape) in cfg.shapes.iter().enumerate() {
            if qlen < shape.span {
                continue;
            }
            let gen = self.nbr.get(si).and_then(|g| g.as_ref());
            for p in 0..=qlen - shape.span {
                if let Some(c) = cfg.code_at(shape, query, Some(&w.mask), p) {
                    w.seeds.push((p as u32, si as u16, c));
                    if let Some(g) = gen {
                        let list = self.neighbors(g, si as u16, c, shape, query, p, &mut word, &mut nb, &mut w.counters);
                        for &n in list.iter() {
                            if n != c {
                                w.seeds.push((p as u32, si as u16, n));
                                extra += 1;
                            }
                        }
                    }
                }
            }
        }
        w.counters.neighbor_seeds += extra;
        w.counters.seed_positions += w.seeds.len() as u64;
        if self.params.hit_filter_min_ids > 0 {
            // Reduced query with sentinels: masked/ambiguous positions and the
            // padding (0xFE) never match a subject byte (0xFF for non-seedable).
            w.qred.clear();
            w.qred.resize(qlen + 2 * HF_PAD, 0xFE);
            w.qraw.clear();
            if self.params.hit_filter_exact_weight > 0 {
                w.qraw.resize(qlen + 2 * HF_PAD, 0xFE);
            }
            for i in 0..qlen {
                if !w.mask[i] {
                    if let Some(r) = cfg.alphabet.reduce(query[i]) {
                        w.qred[HF_PAD + i] = r;
                        if self.params.hit_filter_exact_weight > 0 {
                            w.qraw[HF_PAD + i] = query[i];
                        }
                    }
                }
            }
        }
        // Over-frequent seeds (total occurrences across partitions).
        if self.index.num_partitions() > 0 {
            let caps = &self.occ_cap;
            let idx = self.index;
            let before = w.seeds.len();
            w.seeds.retain(|&(_, s, c)| {
                let tot: u64 = (0..idx.num_partitions()).map(|p| idx.bucket_len(p, s as usize, c) as u64).sum();
                tot <= caps[s as usize]
            });
            w.counters.seeds_frequent_skipped += (before - w.seeds.len()) as u64;
        }
        w.seeds.sort_unstable_by_key(|&(p, s, _)| (p, s));
        w.counters.t_seed_ns += t0.elapsed().as_nanos() as u64;
    }

    /// Scan subject ranges `(partition, lo, hi)` (local subject ids) and keep
    /// the best `max_candidates` regions. Requires [`Self::build_seeds`].
    fn scan_ranges(&self, query: &[u8], w: &mut Worker, ranges: &[(usize, u32, u32)], out: &mut Vec<Region>) {
        let qlen = query.len();
        // Diagonal table sized for this query and the longest subject.
        let max_s = self.index.meta().max_subject_len as usize;
        let need = qlen + max_s + 2;
        if w.diag_epoch.len() < need {
            w.diag_epoch.resize(need, 0);
            w.diag_cover.resize(need, 0);
            w.diag_last.resize(need, 0);
        }

        let mut heap: RegionHeap = std::collections::BinaryHeap::new();
        let cap = self.params.max_candidates.max(1);

        for &(part, s_lo, s_hi) in ranges {
            self.process_partition(query, part, s_lo, s_hi, w, &mut heap, cap);
        }
        out.clear();
        for std::cmp::Reverse((best, std::cmp::Reverse(subject), std::cmp::Reverse(d_min), std::cmp::Reverse(d_max))) in heap.into_iter() {
            out.push(Region { subject, d_min, d_max, best });
        }
        out.sort_unstable_by(region_order);
    }

    #[allow(clippy::too_many_arguments)]
    fn process_partition(
        &self,
        query: &[u8],
        part: usize,
        r_lo: u32,
        r_hi: u32,
        w: &mut Worker,
        heap: &mut RegionHeap,
        cap: usize,
    ) {
        let t0 = Instant::now();
        let idx = self.index;
        let pm = idx.partition(part);
        let nsub = pm.num_subjects as usize;
        if nsub == 0 || w.seeds.is_empty() {
            return;
        }
        #[cfg(target_arch = "x86_64")]
        {
            let whole = r_lo == 0 && r_hi as usize >= nsub;
            if whole
                && self.avx2
                && pm.compact
                && self.params.hit_filter_min_ids > 0
                && self.params.two_hit_window == 0
                && self.params.seed_min_score == 0
                && self.few_hits_per_subject(part, nsub, w)
            {
                // SAFETY: AVX2 (and popcnt/bmi1) detected at engine construction.
                if unsafe { self.filter_first(query, part, w, heap, cap) } {
                    return;
                }
            }
        }
        // Pass 1: count hits per local subject. `touched` collects every
        // subject whose count goes 0 -> 1 (branch-free append).
        let mut total: u64 = 0;
        let whole = r_lo == 0 && r_hi as usize >= nsub;
        if w.touched.len() < 65_536 {
            w.touched.resize(65_536, 0);
        }
        let mut nt = 0usize;
        {
            let counts = &mut w.counts[..];
            let touched = &mut w.touched[..];
            for &(_, s, c) in &w.seeds {
                match idx.bucket(part, s as usize, c) {
                    Postings::Compact(v) => {
                        let (st, en) = if whole { (0, v.len()) } else { (lower_bound_compact(v, r_lo), lower_bound_compact(v, r_hi)) };
                        total += (en - st) as u64;
                        for &p in &v[st..en] {
                            let ls = (p & 0xFFFF) as usize;
                            let cnt = counts[ls];
                            touched[nt] = ls as u16;
                            nt += (cnt == 0) as usize;
                            counts[ls] = cnt + 1;
                        }
                    }
                    b @ Postings::Wide(_) => {
                        let (st, en) = if whole { (0, b.len()) } else { (b.lower_bound(r_lo), b.lower_bound(r_hi)) };
                        total += (en - st) as u64;
                        for i in st..en {
                            let ls = b.get(i).0 as usize;
                            let cnt = counts[ls];
                            touched[nt] = ls as u16;
                            nt += (cnt == 0) as usize;
                            counts[ls] = cnt + 1;
                        }
                    }
                }
            }
        }
        w.counters.postings_visited += total;
        w.counters.hits += total;
        // Sorted list of touched subjects: a (SIMD) scan of the count array
        // is cheaper than sorting when much of the partition was touched.
        if whole && nt * 8 >= nsub {
            nt = nonzero_indices(&w.counts[..nsub], &mut w.touched[..]);
        } else {
            w.touched[..nt].sort_unstable();
        }
        let cap_hits = self.params.hit_buffer.max(1024);
        // Plan subject ranges whose hits fit the buffer.
        let mut ranges: Vec<(usize, usize, u64)> = Vec::new(); // [ti_lo, ti_hi) over touched, hits
        {
            let mut lo = 0usize;
            let mut acc = 0u64;
            for ti in 0..nt {
                let c = w.counts[w.touched[ti] as usize] as u64;
                if acc > 0 && acc + c > cap_hits as u64 {
                    ranges.push((lo, ti, acc));
                    lo = ti;
                    acc = 0;
                }
                acc += c;
            }
            if lo < nt {
                ranges.push((lo, nt, acc));
            }
        }
        w.counters.t_seed_ns += t0.elapsed().as_nanos() as u64;
        // `touched_all` keeps its full length (no re-zeroing per call).
        let touched_all = std::mem::take(&mut w.touched);
        let touched = &touched_all[..nt];
        for &(tlo, thi, nh) in &ranges {
            w.counters.subject_range_passes += 1;
            if nh > cap_hits as u64 && thi - tlo == 1 {
                // One subject alone exceeds the buffer: stream it.
                w.counters.streamed_subjects += 1;
                self.stream_subject(query, part, touched[tlo] as u32, w, heap, cap);
                continue;
            }
            let t1 = Instant::now();
            // Pass 2: scatter hits of subjects touched[tlo..thi] grouped by subject.
            let s_lo = touched[tlo] as u32;
            let s_hi = touched[thi - 1] as u32 + 1;
            let mut acc = 0u32;
            for ti in tlo..thi {
                let ls = touched[ti] as usize;
                w.starts[ls] = acc;
                acc += w.counts[ls];
            }
            if w.hits.len() < acc as usize {
                w.hits.resize(acc as usize, (0, 0));
            }
            let full_range = whole && tlo == 0 && thi == touched.len();
            {
                let starts = &mut w.starts[..];
                let hits = &mut w.hits[..acc as usize];
                for &(qp, s, c) in &w.seeds {
                    match idx.bucket(part, s as usize, c) {
                        Postings::Compact(v) => {
                            let start = if full_range { 0 } else { lower_bound_compact(v, s_lo) };
                            for &p in &v[start..] {
                                let ls = p & 0xFFFF;
                                if ls >= s_hi {
                                    break;
                                }
                                let slot = starts[ls as usize];
                                hits[slot as usize] = (qp, p >> 16);
                                starts[ls as usize] = slot + 1;
                            }
                        }
                        b @ Postings::Wide(_) => {
                            let start = if full_range { 0 } else { b.lower_bound(s_lo) };
                            for i in start..b.len() {
                                let (ls, sp) = b.get(i);
                                if ls >= s_hi {
                                    break;
                                }
                                let slot = starts[ls as usize];
                                hits[slot as usize] = (qp, sp);
                                starts[ls as usize] = slot + 1;
                            }
                        }
                    }
                }
            }
            w.counters.t_seed_ns += t1.elapsed().as_nanos() as u64;
            // Process subjects.
            let t2 = Instant::now();
            let mut begin = 0usize;
            let hits = std::mem::take(&mut w.hits);
            for ti in tlo..thi {
                let ls = touched[ti] as usize;
                let n = w.counts[ls] as usize;
                let subject = pm.first_subject + ls as u64;
                #[cfg(target_arch = "x86_64")]
                if ti + PREFETCH_AHEAD < thi {
                    // Warm a later subject's residues (random access into
                    // the residue buffer is the main stall in hit filtering).
                    let nx = idx.subject_codes((pm.first_subject + touched[ti + PREFETCH_AHEAD] as u64) as usize);
                    let mut o = 0;
                    while o < nx.len().min(256) {
                        // SAFETY: prefetch is a hint; the address is inside `nx`.
                        unsafe { std::arch::x86_64::_mm_prefetch(nx.as_ptr().add(o) as *const i8, std::arch::x86_64::_MM_HINT_T0) };
                        o += 64;
                    }
                }
                self.process_subject(query, subject, &hits[begin..begin + n], w, heap, cap);
                begin += n;
            }
            w.hits = hits;
            w.counters.t_extend_ns += t2.elapsed().as_nanos() as u64;
        }
        for &ls in touched {
            w.counts[ls as usize] = 0;
        }
        w.touched = touched_all;
    }

    /// Choose between the two hit pipelines (identical results). Grouping all
    /// hits by subject first pays off when a subject receives many hits (its
    /// residues stay in cache while they are filtered); with few hits per
    /// subject the counting and scattering passes cost more than they save
    /// and [`Self::filter_first`] is faster. The hits per touched subject are
    /// estimated from the bucket sizes (Poisson occupancy of the partition's
    /// subjects).
    fn few_hits_per_subject(&self, part: usize, nsub: usize, w: &Worker) -> bool {
        const MAX_HITS_PER_SUBJECT: f64 = 6.0;
        let total: u64 = w.seeds.iter().map(|&(_, s, c)| self.index.bucket_len(part, s as usize, c) as u64).sum();
        let n = nsub as f64;
        let touched = n * (1.0 - (-(total as f64) / n).exp());
        (total as f64) <= MAX_HITS_PER_SUBJECT * touched.max(1.0)
    }

    /// Filter-first hit processing for a whole compact partition: every
    /// posting of every seed is tested with the window identity filter in
    /// posting order (the residues of a partition are read in place; the
    /// window of a later hit is prefetched), and only the survivors — about
    /// a tenth of the hits — are grouped by subject and visited (diagonal
    /// coverage check and X-drop extension). Each hit is extended exactly
    /// when the subject-grouped path would extend it: the filter is a pure
    /// function of the hit, and survivors are visited in the same (seed)
    /// order, so the coverage state evolves identically. Returns `false`
    /// without side effects on the regions when the survivors would exceed
    /// the hit buffer; the caller then uses the bounded range path.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2,popcnt,bmi1")]
    unsafe fn filter_first(&self, query: &[u8], part: usize, w: &mut Worker, heap: &mut RegionHeap, cap: usize) -> bool {
        use std::arch::x86_64::*;
        const PF: usize = 12;
        let t0 = Instant::now();
        let idx = self.index;
        let pm = idx.partition(part);
        let first = pm.first_subject as usize;
        let nsub = pm.num_subjects as usize;
        let offs = &idx.subject_offsets()[first..first + nsub + 1];
        let all = idx.residues();
        let all_len = all.len() as i64;
        let min_ids = self.params.hit_filter_min_ids;
        let c = (self.seed_span / 2) as i64;
        let cap_hits = self.params.hit_buffer.max(1024);
        let tlo = _mm256_broadcastsi128_si256(_mm_loadu_si128(self.red_lo.as_ptr() as *const __m128i));
        let thi = _mm256_broadcastsi128_si256(_mm_loadu_si128(self.red_hi.as_ptr() as *const __m128i));
        let nib = _mm256_set1_epi8(0x0f);
        let fifteen = _mm256_set1_epi8(15);
        let map = &idx.seed_config().alphabet.map;
        if w.touched.len() < 65_536 {
            w.touched.resize(65_536, 0);
        }
        let mut surv = std::mem::take(&mut w.surv);
        surv.clear();
        let mut nt = 0usize;
        let mut total = 0u64;
        let mut filtered = 0u64;
        let mut ok = true;
        {
            let counts = &mut w.counts[..];
            let touched = &mut w.touched[..];
            let qred = &w.qred[..];
            let kx = self.params.hit_filter_exact_weight;
            let qraw = &w.qraw[..];
            'seeds: for &(qp, s, code) in &w.seeds {
                let v = match idx.bucket(part, s as usize, code) {
                    Postings::Compact(v) => v,
                    Postings::Wide(_) => {
                        ok = false;
                        break 'seeds;
                    }
                };
                total += v.len() as u64;
                let qo = (qp as i64 + c - 16 + HF_PAD as i64) as usize;
                let qv = _mm256_loadu_si256(qred.as_ptr().add(qo) as *const __m256i);
                let qrv = if kx > 0 { _mm256_loadu_si256(qraw.as_ptr().add(qo) as *const __m256i) } else { qv };
                for k in 0..v.len() {
                    if k + PF < v.len() {
                        let pp = *v.get_unchecked(k + PF);
                        let o = *offs.get_unchecked((pp & 0xFFFF) as usize) as i64 + (pp >> 16) as i64 + c - 16;
                        if o >= 0 && o < all_len {
                            _mm_prefetch(all.as_ptr().add(o as usize) as *const i8, _MM_HINT_T0);
                        }
                    }
                    let p = *v.get_unchecked(k);
                    let ls = (p & 0xFFFF) as usize;
                    if ls >= nsub {
                        continue; // defensive: corrupt posting
                    }
                    let sp = (p >> 16) as i64;
                    let base = offs[ls] as i64;
                    let sl = offs[ls + 1] as i64 - base;
                    if sp >= sl {
                        continue; // defensive: corrupt posting
                    }
                    let ss = sp + c - 16;
                    let g = base + ss;
                    let ids = if g >= 0 && g + 32 <= all_len {
                        // The load may extend into neighbouring subjects; those lanes are masked.
                        let sv = _mm256_loadu_si256(all.as_ptr().add(g as usize) as *const __m256i);
                        let ix = _mm256_and_si256(sv, nib);
                        let ishi = _mm256_cmpgt_epi8(sv, fifteen);
                        let red = _mm256_blendv_epi8(_mm256_shuffle_epi8(tlo, ix), _mm256_shuffle_epi8(thi, ix), ishi);
                        let eq = _mm256_movemask_epi8(_mm256_cmpeq_epi8(red, qv)) as u32;
                        let lo = (-ss).max(0) as u32;
                        let hi = (sl - ss).min(32) as u32;
                        let lane = (if hi >= 32 { u32::MAX } else { (1u32 << hi) - 1 }) & !((1u32 << lo) - 1);
                        let r = (eq & lane).count_ones();
                        if kx > 0 {
                            let ex = _mm256_movemask_epi8(_mm256_cmpeq_epi8(sv, qrv)) as u32;
                            r + kx * (ex & eq & lane).count_ones()
                        } else {
                            r
                        }
                    } else {
                        let subj = &all[base as usize..(base + sl) as usize];
                        let mut n = 0u32;
                        for t in 0..32i64 {
                            let j = ss + t;
                            if j >= 0 && j < sl {
                                let hit = map[subj[j as usize] as usize] == qred[qo + t as usize];
                                n += hit as u32;
                                if kx > 0 && hit && subj[j as usize] == qraw[qo + t as usize] {
                                    n += kx;
                                }
                            }
                        }
                        n
                    };
                    if ids >= min_ids {
                        surv.push((qp, p));
                        let cnt = counts[ls];
                        touched[nt] = ls as u16;
                        nt += (cnt == 0) as usize;
                        counts[ls] = cnt + 1;
                        if surv.len() > cap_hits {
                            ok = false;
                            break 'seeds;
                        }
                    } else {
                        filtered += 1;
                    }
                }
            }
        }
        if !ok {
            for &ls in &w.touched[..nt] {
                w.counts[ls as usize] = 0;
            }
            w.surv = surv;
            w.counters.t_extend_ns += t0.elapsed().as_nanos() as u64;
            return false;
        }
        w.counters.postings_visited += total;
        w.counters.hits += total;
        w.counters.hits_filtered += filtered;
        // Group survivors by subject (ascending local id), keeping seed order.
        w.touched[..nt].sort_unstable();
        let mut acc = 0u32;
        for ti in 0..nt {
            let ls = w.touched[ti] as usize;
            w.starts[ls] = acc;
            acc += w.counts[ls];
        }
        if w.hits.len() < acc as usize {
            w.hits.resize(acc as usize, (0, 0));
        }
        for &(qp, p) in &surv {
            let ls = (p & 0xFFFF) as usize;
            let slot = w.starts[ls];
            w.hits[slot as usize] = (qp, p >> 16);
            w.starts[ls] = slot + 1;
        }
        let hits = std::mem::take(&mut w.hits);
        let touched = std::mem::take(&mut w.touched);
        let mut begin = 0usize;
        for &ls in &touched[..nt] {
            let n = w.counts[ls as usize] as usize;
            let subject = pm.first_subject + ls as u64;
            let subj = idx.subject_codes(subject as usize);
            self.begin_subject(w);
            self.visit_hits_fast::<false>(query, subj, &hits[begin..begin + n], w);
            self.finish_subject(subject, subj.len(), w, heap, cap);
            begin += n;
        }
        for &ls in &touched[..nt] {
            w.counts[ls as usize] = 0;
        }
        w.hits = hits;
        w.touched = touched;
        w.surv = surv;
        w.counters.t_extend_ns += t0.elapsed().as_nanos() as u64;
        true
    }

    /// Stream one subject's hits (query order) without materialising them.
    fn stream_subject(
        &self,
        query: &[u8],
        part: usize,
        ls: u32,
        w: &mut Worker,
        heap: &mut RegionHeap,
        cap: usize,
    ) {
        let t2 = Instant::now();
        let idx = self.index;
        let subject = idx.partition(part).first_subject + ls as u64;
        let subj = idx.subject_codes(subject as usize);
        self.begin_subject(w);
        let seeds = std::mem::take(&mut w.seeds);
        let mut one = Vec::with_capacity(64);
        for &(qp, s, c) in &seeds {
            let b = idx.bucket(part, s as usize, c);
            let mut i = b.lower_bound(ls);
            one.clear();
            while i < b.len() {
                let (l, sp) = b.get(i);
                if l != ls {
                    break;
                }
                one.push((qp, sp));
                i += 1;
            }
            for &h in &one {
                self.visit_hit(query, subj, h, w);
            }
        }
        w.seeds = seeds;
        self.finish_subject(subject, subj.len(), w, heap, cap);
        w.counters.t_extend_ns += t2.elapsed().as_nanos() as u64;
    }

    fn process_subject(
        &self,
        query: &[u8],
        subject: u64,
        hits: &[(u32, u32)],
        w: &mut Worker,
        heap: &mut RegionHeap,
        cap: usize,
    ) {
        let subj = self.index.subject_codes(subject as usize);
        self.begin_subject(w);
        let simple = self.params.two_hit_window == 0 && self.params.seed_min_score == 0;
        #[cfg(target_arch = "x86_64")]
        {
            if simple && self.avx2 {
                // SAFETY: AVX2 detected at engine construction.
                unsafe { self.visit_hits_fast::<true>(query, subj, hits, w) };
                self.finish_subject(subject, subj.len(), w, heap, cap);
                return;
            }
        }
        let _ = simple;
        for &h in hits {
            self.visit_hit(query, subj, h, w);
        }
        self.finish_subject(subject, subj.len(), w, heap, cap);
    }

    /// Fast path for single-hit triggering: diagonal coverage check, SIMD
    /// window-identity filter and profile-based X-drop extension, with all
    /// buffers held in local slices. Same decisions as [`Self::visit_hit`].
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2,popcnt,bmi1")]
    unsafe fn visit_hits_fast<const FILTER: bool>(&self, query: &[u8], subj: &[u8], hits: &[(u32, u32)], w: &mut Worker) {
        use std::arch::x86_64::*;
        let qlen = query.len();
        let sl = subj.len();
        let epoch = w.epoch;
        let min_ids = self.params.hit_filter_min_ids;
        let thr = w.hsp_thr;
        let x = self.params.x_drop_ungapped;
        let c = self.seed_span / 2;
        let tlo = _mm256_broadcastsi128_si256(_mm_loadu_si128(self.red_lo.as_ptr() as *const __m128i));
        let thi = _mm256_broadcastsi128_si256(_mm_loadu_si128(self.red_hi.as_ptr() as *const __m128i));
        let nib = _mm256_set1_epi8(0x0f);
        let fifteen = _mm256_set1_epi8(15);
        let map = &self.index.seed_config().alphabet.map;
        let (mut avoided, mut filtered, mut ext) = (0u64, 0u64, 0u64);
        let de = &mut w.diag_epoch[..];
        let dc = &mut w.diag_cover[..];
        let dl = &mut w.diag_last[..];
        let qred = &w.qred[..];
        let kx = self.params.hit_filter_exact_weight;
        let qraw = &w.qraw[..];
        let prof = &w.prof[..];
        let hsps = &mut w.hsps;
        let all = self.index.residues();
        let base = subj.as_ptr() as i64 - all.as_ptr() as i64;
        let all_len = all.len() as i64;
        for &(qp, sp) in hits {
            let d = sp as i64 - qp as i64;
            let di = (d + qlen as i64) as usize;
            if di >= de.len() || sp as usize >= sl {
                continue;
            }
            if de[di] != epoch {
                de[di] = epoch;
                dc[di] = 0;
                dl[di] = u32::MAX;
            }
            if FILTER && min_ids > 0 {
                let qs = (qp as usize + c) as i64 - 16;
                let ss = (sp as usize + c) as i64 - 16;
                let qo = (qs + HF_PAD as i64) as usize;
                let g = base + ss;
                let ids = if g >= 0 && g + 32 <= all_len {
                    // The 32-byte load may extend into neighbouring subjects
                    // (still inside the residue buffer); those lanes are masked.
                    let sv = _mm256_loadu_si256(all.as_ptr().add(g as usize) as *const __m256i);
                    let qv = _mm256_loadu_si256(qred.as_ptr().add(qo) as *const __m256i);
                    let idx = _mm256_and_si256(sv, nib);
                    let ishi = _mm256_cmpgt_epi8(sv, fifteen);
                    let red = _mm256_blendv_epi8(_mm256_shuffle_epi8(tlo, idx), _mm256_shuffle_epi8(thi, idx), ishi);
                    let eq = _mm256_movemask_epi8(_mm256_cmpeq_epi8(red, qv)) as u32;
                    let lo = (-ss).max(0) as u32;
                    let hi = (sl as i64 - ss).min(32) as u32;
                    let lane = (if hi >= 32 { u32::MAX } else { (1u32 << hi) - 1 }) & !((1u32 << lo) - 1);
                    let r = (eq & lane).count_ones();
                    if kx > 0 {
                        let qrv = _mm256_loadu_si256(qraw.as_ptr().add(qo) as *const __m256i);
                        let ex = _mm256_movemask_epi8(_mm256_cmpeq_epi8(sv, qrv)) as u32;
                        r + kx * (ex & eq & lane).count_ones()
                    } else {
                        r
                    }
                } else {
                    let mut n = 0u32;
                    for t in 0..32i64 {
                        let j = ss + t;
                        if j >= 0 && j < sl as i64 {
                            let hit = map[subj[j as usize] as usize] == qred[qo + t as usize];
                            n += hit as u32;
                            if kx > 0 && hit && subj[j as usize] == qraw[qo + t as usize] {
                                n += kx;
                            }
                        }
                    }
                    n
                };
                if ids < min_ids {
                    filtered += 1;
                    continue;
                }
            }
            if qp < dc[di] {
                avoided += 1;
                continue;
            }
            ext += 1;
            let (score, lo, hi) = ungapped_xdrop_prof(prof, qlen, subj, qp as usize, sp as usize, x);
            dc[di] = hi as u32;
            dl[di] = qp;
            if score >= thr {
                hsps.push(Hsp { diag: d, q_lo: lo as u32, q_hi: hi as u32, score });
            }
        }
        w.counters.extensions_avoided += avoided;
        w.counters.hits_filtered += filtered;
        w.counters.ungapped_extensions += ext;
    }

    #[inline]
    fn begin_subject(&self, w: &mut Worker) {
        w.epoch = w.epoch.wrapping_add(1);
        if w.epoch == 0 {
            w.diag_epoch.iter_mut().for_each(|e| *e = 0);
            w.epoch = 1;
        }
        w.hsps.clear();
    }

    #[inline]
    fn visit_hit(&self, query: &[u8], subj: &[u8], (qp, sp): (u32, u32), w: &mut Worker) {
        let qlen = query.len();
        let d = sp as i64 - qp as i64;
        let di = (d + qlen as i64) as usize;
        if di >= w.diag_epoch.len() || sp as usize >= subj.len() {
            return; // defensive: corrupt posting
        }
        if w.diag_epoch[di] != w.epoch {
            w.diag_epoch[di] = w.epoch;
            w.diag_cover[di] = 0;
            w.diag_last[di] = u32::MAX;
        }
        // Window filter before the coverage check: both only reject, so the
        // order does not change which hits are extended, and this order makes
        // the counters (`hits_filtered` = all hits failing the filter,
        // `extensions_avoided` = passing hits on an already extended stretch)
        // the same for every hit pipeline.
        if self.params.hit_filter_min_ids > 0 {
            let ids = self.window_identities(subj, qp as usize, sp as usize, w);
            if ids < self.params.hit_filter_min_ids {
                w.counters.hits_filtered += 1;
                return;
            }
        }
        if qp < w.diag_cover[di] {
            w.counters.extensions_avoided += 1;
            return;
        }
        let need_score = self.params.seed_min_score > 0 || (self.params.two_hit_window > 0 && self.params.single_hit_score > 0);
        let sscore = if need_score { self.seed_score(query, subj, qp as usize, sp as usize, &w.bias) } else { 0 };
        if self.params.seed_min_score > 0 && sscore < self.params.seed_min_score {
            w.counters.seeds_low_score += 1;
            return;
        }
        let win = self.params.two_hit_window;
        if win > 0 && !(self.params.single_hit_score > 0 && sscore >= self.params.single_hit_score) {
            let last = w.diag_last[di];
            let ok = last != u32::MAX && qp >= last + self.two_hit_min && qp - last <= win;
            if !ok {
                w.diag_last[di] = qp;
                w.counters.two_hit_waits += 1;
                return;
            }
        }
        // Ungapped X-drop extension (bias-corrected) through the seed start.
        w.counters.ungapped_extensions += 1;
        let (score, lo, hi) = ungapped_xdrop_prof(&w.prof, query.len(), subj, qp as usize, sp as usize, self.params.x_drop_ungapped);
        w.diag_cover[di] = hi as u32;
        w.diag_last[di] = qp;
        if score >= w.hsp_thr {
            w.hsps.push(Hsp { diag: d, q_lo: lo as u32, q_hi: hi as u32, score });
        }
    }

    /// Identities in the seed alphabet between the 32-residue windows centred
    /// on the seed in query and subject.
    #[inline]
    fn window_identities(&self, subj: &[u8], qp: usize, sp: usize, w: &mut Worker) -> u32 {
        let c = self.seed_span / 2;
        // window start offsets (may be negative near sequence starts)
        let qs = (qp + c) as i64 - 16;
        let ss = (sp + c) as i64 - 16;
        let qo = (qs + HF_PAD as i64) as usize;
        let qwin = &w.qred[qo..qo + 32];
        let kx = self.params.hit_filter_exact_weight;
        let map = &self.index.seed_config().alphabet.map;
        let sl = subj.len() as i64;
        let mut n = 0u32;
        if ss >= 0 && ss + 32 <= sl && kx == 0 {
            let sw = &subj[ss as usize..ss as usize + 32];
            #[cfg(target_arch = "x86_64")]
            {
                if self.avx2 {
                    // SAFETY: AVX2 detected at engine construction.
                    return unsafe { ids_avx2(qwin, sw, &self.red_lo, &self.red_hi) };
                }
            }
            for t in 0..32 {
                n += (map[sw[t] as usize] == qwin[t]) as u32;
            }
            return n;
        }
        for t in 0..32i64 {
            let j = ss + t;
            if j >= 0 && j < sl {
                let hit = map[subj[j as usize] as usize] == qwin[t as usize];
                n += hit as u32;
                if kx > 0 && hit && subj[j as usize] == w.qraw[qo + t as usize] {
                    n += kx;
                }
            }
        }
        let _ = &mut w.sbuf;
        n
    }

    /// Score of the seed's care positions (first shape's span, all positions
    /// of the span for contiguous seeds).
    #[inline]
    fn seed_score(&self, q: &[u8], s: &[u8], qp: usize, sp: usize, bias: &[i8]) -> i32 {
        let span = self.seed_span;
        if qp + span > q.len() || sp + span > s.len() {
            return 0;
        }
        let mut t = 0i32;
        for o in 0..span {
            t += self.matrix.matrix[code_index(q[qp + o])][code_index(s[sp + o])] + bias[qp + o] as i32;
        }
        t
    }

    fn finish_subject(
        &self,
        subject: u64,
        _slen: usize,
        w: &mut Worker,
        heap: &mut RegionHeap,
        cap: usize,
    ) {
        if w.hsps.is_empty() {
            return;
        }
        w.counters.hsps += w.hsps.len() as u64;
        // Cluster HSPs by diagonal: HSPs closer than the window slack share a
        // region (their windows would overlap anyway).
        w.hsps.sort_unstable_by_key(|h| h.diag);
        let join = 2 * self.params.window_slack;
        let mut regs: Vec<Region> = Vec::new();
        let mut evidence: Vec<(u32, i32)> = Vec::new(); // (HSP count, score sum) per region
        for h in w.hsps.iter() {
            match regs.last_mut() {
                Some(r) if h.diag - r.d_max <= join => {
                    r.d_max = r.d_max.max(h.diag);
                    r.best = r.best.max(h.score);
                    let e = evidence.last_mut().expect("parallel vectors");
                    e.0 += 1;
                    e.1 = e.1.saturating_add(h.score);
                }
                _ => {
                    regs.push(Region { subject, d_min: h.diag, d_max: h.diag, best: h.score });
                    evidence.push((1, h.score));
                }
            }
        }
        // Coherent evidence: keep a region whose best HSP reaches the region
        // threshold, or that has at least two HSPs whose scores add up to it
        // (a gapped alignment joining several segments).
        let thr = w.region_thr;
        if thr > self.params.hsp_threshold {
            let before = regs.len();
            let mut k = 0;
            regs.retain(|r| {
                let (n, sum) = evidence[k];
                k += 1;
                r.best >= thr || (n >= 2 && sum >= thr)
            });
            w.counters.regions_weak_dropped += (before - regs.len()) as u64;
        }
        regs.sort_unstable_by(|a, b| b.best.cmp(&a.best).then(a.d_min.cmp(&b.d_min)));
        regs.truncate(self.params.max_regions_per_subject.max(1));
        for r in regs {
            let key = std::cmp::Reverse((r.best, std::cmp::Reverse(r.subject), std::cmp::Reverse(r.d_min), std::cmp::Reverse(r.d_max)));
            if heap.len() < cap {
                heap.push(key);
            } else if let Some(top) = heap.peek() {
                // `top` is the weakest region kept so far.
                if key < *top {
                    heap.pop();
                    heap.push(key);
                }
                w.counters.regions_dropped_by_cap += 1;
            }
        }
    }
}

fn overlaps(a: &Hit, b: &Hit) -> bool {
    a.q_start <= b.q_end && b.q_start <= a.q_end && a.s_start <= b.s_end && b.s_start <= a.s_end
}

/// Ungapped X-drop extension using a per-query profile (`prof[i*32 + c]` =
/// S(q_i, c) + bias_i). Identical results to [`ungapped_xdrop`]. The best
/// score and its position are updated with selects rather than branches
/// (the improvement test is data-dependent and mispredicts often); the only
/// branch left is the X-drop exit.
#[inline]
pub fn ungapped_xdrop_prof(prof: &[i16], qlen: usize, s: &[u8], qp: usize, sp: usize, x: i32) -> (i32, usize, usize) {
    let n = (qlen - qp).min(s.len() - sp);
    let pr = &prof[qp * 32..(qp + n) * 32];
    let sr = &s[sp..sp + n];
    let (mut best_r, mut cur, mut end_t) = (0i32, 0i32, 0usize);
    for t in 0..n {
        // SAFETY: t < n, so t*32 + (c & 31) < n*32 = pr.len(), and t < sr.len().
        cur += unsafe { *pr.get_unchecked(t * 32 + (*sr.get_unchecked(t) as usize & 31)) } as i32;
        let better = cur > best_r;
        end_t = if better { t + 1 } else { end_t };
        best_r = best_r.max(cur);
        if cur < best_r - x {
            break;
        }
    }
    let nl = qp.min(sp);
    let pl = &prof[(qp - nl) * 32..qp * 32];
    let sl = &s[sp - nl..sp];
    let (mut best_l, mut cur, mut st_t) = (0i32, 0i32, 0usize);
    for t in 1..=nl {
        // SAFETY: 1 <= t <= nl: (nl-t)*32 + (c & 31) < nl*32 = pl.len(), nl-t < sl.len().
        cur += unsafe { *pl.get_unchecked((nl - t) * 32 + (*sl.get_unchecked(nl - t) as usize & 31)) } as i32;
        let better = cur > best_l;
        st_t = if better { t } else { st_t };
        best_l = best_l.max(cur);
        if cur < best_l - x {
            break;
        }
    }
    (best_r + best_l, qp - st_t, qp + end_t)
}

/// Ungapped X-drop extension through (qp, sp) with per-query-position bias.
/// Returns (score, q_lo, q_hi) with q range half-open.
#[inline]
pub fn ungapped_xdrop(q: &[u8], s: &[u8], qp: usize, sp: usize, m: &[[i32; 24]; 24], bias: &[i8], x: i32) -> (i32, usize, usize) {
    // right (including the seed start)
    let mut best_r = 0i32;
    let mut cur = 0i32;
    let mut end = qp;
    let (mut i, mut j) = (qp, sp);
    while i < q.len() && j < s.len() {
        cur += m[code_index(q[i])][code_index(s[j])] + bias[i] as i32;
        i += 1;
        j += 1;
        if cur > best_r {
            best_r = cur;
            end = i;
        } else if cur < best_r - x {
            break;
        }
    }
    // left
    let mut best_l = 0i32;
    let mut cur = 0i32;
    let mut start = qp;
    let (mut i, mut j) = (qp as isize - 1, sp as isize - 1);
    while i >= 0 && j >= 0 {
        cur += m[code_index(q[i as usize])][code_index(s[j as usize])] + bias[i as usize] as i32;
        if cur > best_l {
            best_l = cur;
            start = i as usize;
        } else if cur < best_l - x {
            break;
        }
        i -= 1;
        j -= 1;
    }
    (best_r + best_l, start, end.max(qp))
}

/// First element of a compact, subject-sorted bucket with local subject >= s.
#[inline]
fn lower_bound_compact(v: &[u32], s: u32) -> usize {
    v.partition_point(|&p| (p & 0xFFFF) < s)
}

/// Write the indices of the non-zero entries of `counts` (ascending) into
/// `out`; returns how many were written. AVX2 when available.
fn nonzero_indices(counts: &[u32], out: &mut [u16]) -> usize {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("popcnt") && is_x86_feature_detected!("bmi1") {
            // SAFETY: features detected.
            return unsafe { nonzero_indices_avx2(counts, out) };
        }
    }
    let mut n = 0;
    for (i, &c) in counts.iter().enumerate() {
        out[n] = i as u16;
        n += (c != 0) as usize;
    }
    n
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt,bmi1")]
unsafe fn nonzero_indices_avx2(counts: &[u32], out: &mut [u16]) -> usize {
    use std::arch::x86_64::*;
    let mut n = 0usize;
    let chunks = counts.len() / 8;
    let zero = _mm256_setzero_si256();
    for k in 0..chunks {
        let v = _mm256_loadu_si256(counts.as_ptr().add(k * 8) as *const __m256i);
        let z = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(v, zero))) as u32;
        let mut m = !z & 0xFF;
        while m != 0 {
            let b = m.trailing_zeros() as usize;
            out[n] = (k * 8 + b) as u16;
            n += 1;
            m &= m - 1;
        }
    }
    for i in chunks * 8..counts.len() {
        out[n] = i as u16;
        n += (counts[i] != 0) as usize;
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{ScoreMatrix, ScoringSystem};

    #[test]
    fn profile_xdrop_matches_reference_xdrop() {
        let m = ScoreMatrix::new(ScoringSystem::Blosum62);
        let mut st = 12345u64;
        let mut rnd = |k: u64| {
            st = st.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (st >> 33) % k
        };
        for t in 0..3000 {
            let ql = 1 + rnd(300) as usize;
            let q: Vec<u8> = (0..ql).map(|_| rnd(24) as u8).collect();
            let mut s: Vec<u8> = (0..1 + rnd(300) as usize).map(|_| rnd(24) as u8).collect();
            if t % 2 == 0 {
                // related segment so extensions run long
                let off = rnd(s.len() as u64) as usize;
                for (k, c) in s[off..].iter_mut().enumerate() {
                    if k < ql && rnd(3) > 0 {
                        *c = q[k];
                    }
                }
            }
            let bias: Vec<i8> = (0..ql).map(|_| if t % 3 == 0 { rnd(7) as i8 - 3 } else { 0 }).collect();
            let mut prof = vec![0i16; ql * 32];
            for i in 0..ql {
                for c in 0..32 {
                    prof[i * 32 + c] = (m.matrix[code_index(q[i])][if c < 24 { c } else { crate::alphabet::AA_X as usize }] + bias[i] as i32) as i16;
                }
            }
            let x = [10, 20, 35][t % 3];
            for _ in 0..5 {
                let qp = rnd(ql as u64) as usize;
                let sp = rnd(s.len() as u64) as usize;
                assert_eq!(ungapped_xdrop_prof(&prof, ql, &s, qp, sp, x), ungapped_xdrop(&q, &s, qp, sp, &m.matrix, &bias, x), "t={t} qp={qp} sp={sp}");
            }
        }
    }
}
