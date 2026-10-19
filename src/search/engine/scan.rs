//! Batched database-scan schedule ("query tables").
//!
//! The per-query schedules look every query seed up in the database's seed
//! postings: the work of a query is proportional to the postings it touches,
//! and the residues behind every posting are fetched from random places of
//! the database. For many queries against a large database this schedule
//! inverts the roles. The seeds of a whole batch of queries go into one small
//! lookup table per seed shape, and the database is streamed once per batch,
//! subject by subject: every database position looks its seed codes up in
//! the tables, and every match is a hit of one query that is tested with the
//! window filter right away — against the query's residues (a few MB for the
//! whole batch, cache resident) while the subject is in L1. Survivors are
//! grouped by (subject, query) and go through exactly the per-pair logic of
//! the index schedules (diagonal coverage, X-drop extension, region
//! clustering, per-query candidate cap under the same total order), so for the
//! same seeds and parameters the candidate regions — and hence the reported
//! alignments — are identical; `tests::scan_matches_index_schedules` checks
//! this.
//!
//! The scan needs no seed postings, only the residues and the database seed
//! mask, so indexes built for it are small and quick to build, and the seed
//! shapes are chosen at search time.

use super::*;
use crate::alphabet::{SeedAlphabet, SeedShape};
use crate::search::hugebuf::HugeBuf;

/// Seeding scheme of the scan schedule.
#[derive(Clone, Debug)]
pub struct ScanSeeds {
    pub alphabet: SeedAlphabet,
    pub shapes: Vec<SeedShape>,
}

/// Padding (bytes) around every query and around the residue streams; at
/// least the window filter's half-width so window loads never leave a buffer.
const PAD: usize = 32;
/// Seed byte of a residue that is not in the seed alphabet.
const NONSEED: u8 = 0x7F;
/// Flag on the seed byte of a residue masked from seeding (the reduced code
/// still counts as an identity in the window filter).
const MASKED: u8 = 0x80;
/// Query sentinel (masked, ambiguous or padding): never equal to a database byte.
const QSENT: u8 = 0xFE;
/// Residues per scan work unit.
const CHUNK_RESIDUES: usize = 1 << 14;
/// Largest number of residues per scan round (see [`Engine::search_scan`]).
const SUPER_RESIDUES: usize = 1 << 20;
/// Survivors a round should hold at most: rounds shrink after one that
/// exceeded it (dense stretches of the database: self hits, close homologs)
/// and grow back after sparse ones.
const ROUND_SURVIVORS: usize = 1 << 21;

/// Default for `SearchParams::scan_pair_survivors`: a (query, subject) pair
/// projected to hold more window-filter survivors than this is scanned and
/// visited on its own, in query-position slices of about this many.
const PAIR_SURVIVORS: usize = 1 << 22;

/// A (query, subject) pair taken out of the round's work items because its
/// survivors would exceed the pair limit (`projected` = estimated total).
#[derive(Clone, Copy, Debug)]
struct PairJob {
    subject: usize,
    q: usize,
    projected: usize,
}
/// Code spaces up to this size use a dense offset table.
const DENSE_MAX_CODES: u64 = 1 << 22;

/// The database residues with PAD readable bytes before and after (the
/// vector loads of the window filter and X-drop run past subject ends; the
/// lanes outside a subject are masked, so the padding values do not matter).
enum RawResidues<'a> {
    /// Not needed (no exact-identity weight).
    None,
    /// A padded copy (the index storage has no slack around its residues).
    Own(HugeBuf<u8>),
    /// The index's own storage: PAD bytes before the residues .. PAD after.
    View(&'a [u8]),
}

impl std::ops::Deref for RawResidues<'_> {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        match self {
            RawResidues::None => &[],
            RawResidues::Own(h) => h,
            RawResidues::View(v) => v,
        }
    }
}

/// Database-side data of the scan schedule (built once per engine).
pub(crate) struct ScanDb<'a> {
    /// Seed byte of every residue code: seed-alphabet code, or NONSEED.
    seed_map: [u8; 256],
    /// Database seed mask: sorted, disjoint residue ranges `[a, b)`.
    mask_runs: Vec<(u64, u64)>,
    /// PAD + residues (alignment encoding) + PAD; only when the window
    /// filter counts exact identities.
    raw: RawResidues<'a>,
    /// Per shape: bitmap of over-frequent codes (empty = none).
    frequent: Vec<Vec<u64>>,
    /// Per shape: the occurrence cap.
    pub(crate) caps: Vec<u64>,
    /// Subject ranges `[lo, hi)` scanned as one work unit.
    chunks: Vec<(usize, usize)>,
    /// How the database seed mask was obtained.
    pub(crate) mask_source: &'static str,
    pub(crate) build_s: f64,
}

/// Seed bytes of one stretch of the database (one scan round): PAD +
/// seed-alphabet code | MASKED, or NONSEED, per residue + PAD. Built per
/// round instead of for the whole database (one byte per residue saved).
pub(crate) struct RoundSeeds {
    buf: Vec<u8>,
    /// First residue of the stretch: residue `r` is at `buf[PAD + r - base]`.
    base: usize,
}

/// Seed bytes of residues `[r0, r1)` (global residue offsets).
fn seed_bytes(seed_map: &[u8; 256], mask_runs: &[(u64, u64)], res: &[u8], r0: usize, r1: usize) -> RoundSeeds {
    let mut buf = vec![NONSEED; r1 - r0 + 2 * PAD];
    for (x, &c) in buf[PAD..PAD + r1 - r0].iter_mut().zip(&res[r0..r1]) {
        *x = seed_map[c as usize];
    }
    let k = mask_runs.partition_point(|&(_, b)| b as usize <= r0);
    for &(a, b) in &mask_runs[k..] {
        let (a, b) = (a as usize, b as usize);
        if a >= r1 {
            break;
        }
        for x in &mut buf[PAD + a.max(r0) - r0..PAD + b.min(r1) - r0] {
            if *x != NONSEED {
                *x |= MASKED;
            }
        }
    }
    RoundSeeds { buf, base: r0 }
}

impl ScanDb<'_> {
    /// Seed bytes of residues `[r0, r1)` (global residue offsets).
    fn round_seeds(&self, res: &[u8], r0: usize, r1: usize) -> RoundSeeds {
        seed_bytes(&self.seed_map, &self.mask_runs, res, r0, r1)
    }
}

/// One seed table: codes -> batch positions of the seed starts.
///
/// Two layouts. *Dense* (code spaces a batch fills well): an occupancy
/// bitmap and one entry offset per code. *Blocked* (large code spaces): code
/// `c` belongs to block `c mod nb` with tag `c / nb` (< 32); one 64-bit
/// record per block holds the block's occupancy bits and the offset of its
/// first entry, and the block's entries are tagged and sorted by tag. A lookup
/// touches one record and one short run of entries (instead of a bitmap, a
/// rank table, an offset table and the entries). Blocking by residue keeps
/// codes with a common prefix — similar frequencies — in different blocks, so
/// a frequent code does not share its block with its frequent siblings.
struct Table {
    /// Code of an invalid window (absent in both layouts).
    inv: u32,
    dense: bool,
    /// Blocked: number of blocks, and `c / nb == (c * magic) >> mshift`.
    nb: u32,
    magic: u64,
    mshift: u32,
    /// Dense: entry offset per code (+1).
    offs: HugeBuf<u32>,
    /// Dense: occupancy bitmap.
    bits: HugeBuf<u64>,
    /// Blocked: per block, occupancy << 32 | offset of the first entry; two
    /// sentinel records (no codes, offset = number of entries).
    recs: HugeBuf<u64>,
    /// Batch positions (blocked: position << 5 | code & 31), sorted by code,
    /// followed by `ENTRY_PAD` padding entries.
    entries: HugeBuf<u32>,
    /// Number of entries (without the padding).
    n: usize,
}

/// Tag bits of a blocked-layout entry (32 codes per block).
const BLK_SHIFT: u32 = 5;
const BLK_MASK: u32 = (1 << BLK_SHIFT) - 1;
/// Padding entries after a table's entries (vector reads of a block's run).
const ENTRY_PAD: usize = 8;

impl Table {
    fn len(&self) -> usize {
        self.n
    }

    /// Blocked layout: (block, tag) of a code.
    #[inline(always)]
    fn split(&self, code: u32) -> (usize, u32) {
        let t = ((code as u64 * self.magic) >> self.mshift) as u32;
        (code as usize - t as usize * self.nb as usize, t)
    }

    /// Dense layout from (code, batch position) lists in batch order.
    fn dense(ncodes: u64, per_q: &[(Vec<(u32, u32)>, u64)], nseeds: usize) -> Table {
        let mut offs = vec![0u32; ncodes as usize + 1];
        for (v, _) in per_q {
            for &(c, _) in v {
                offs[c as usize + 1] += 1;
            }
        }
        for i in 0..ncodes as usize {
            offs[i + 1] += offs[i];
        }
        let mut pos = offs.clone();
        let mut bits = vec![0u64; (ncodes / 64) as usize + 2];
        let mut entries = vec![0u32; nseeds + ENTRY_PAD];
        for (v, _) in per_q {
            for &(c, g) in v {
                let p = &mut pos[c as usize];
                entries[*p as usize] = g;
                *p += 1;
                bits[(c >> 6) as usize] |= 1u64 << (c & 63);
            }
        }
        Table {
            inv: ncodes as u32,
            dense: true,
            nb: 0,
            magic: 0,
            mshift: 0,
            offs: HugeBuf::from_vec(offs),
            bits: HugeBuf::from_vec(bits),
            recs: HugeBuf::default(),
            entries: HugeBuf::from_vec(entries),
            n: nseeds,
        }
    }

    /// Blocked layout from (code, batch position) lists in batch order.
    fn blocked(ncodes: u64, per_q: &[(Vec<(u32, u32)>, u64)], nseeds: usize) -> Table {
        // 32 * nb > ncodes, so every code (and the invalid code ncodes) has a
        // tag below 32. Exact division by multiply-shift: with
        // 2^mshift > 32 nb^2 the rounding error stays below 1 / nb.
        let nb = (ncodes + BLK_MASK as u64) / (BLK_MASK as u64 + 1) + 1;
        assert!(nb < 1 << 26, "code space too large for the blocked layout");
        let mshift = 64 - (32 * nb * nb).leading_zeros() + 1;
        let magic = (((1u128 << mshift) + nb as u128 - 1) / nb as u128) as u64;
        let mut t = Table {
            inv: ncodes as u32,
            dense: false,
            nb: nb as u32,
            magic,
            mshift,
            offs: HugeBuf::default(),
            bits: HugeBuf::default(),
            recs: HugeBuf::default(),
            entries: HugeBuf::default(),
            n: nseeds,
        };
        for c in [0u64, 1, nb - 1, nb, nb + 1, 31 * nb - 1, 31 * nb, ncodes.saturating_sub(1), ncodes] {
            if c <= ncodes {
                let (k, tag) = t.split(c as u32);
                debug_assert_eq!((k as u64, tag as u64), (c % nb, c / nb));
            }
        }
        let nbk = nb as usize;
        // Entries per block, a stable counting sort by block (batch order
        // within a block), then a stable sort by tag within each block.
        let mut start = vec![0u32; nbk + 2];
        let mut occ = vec![0u32; nbk + 2];
        for (v, _) in per_q {
            for &(c, _) in v {
                let (k, tag) = t.split(c);
                start[k + 1] += 1;
                occ[k] |= 1 << tag;
            }
        }
        for i in 0..nbk + 1 {
            start[i + 1] += start[i];
        }
        let mut pos = start.clone();
        let mut entries = vec![0u32; nseeds + ENTRY_PAD];
        for (v, _) in per_q {
            for &(c, g) in v {
                let (k, tag) = t.split(c);
                let p = &mut pos[k];
                entries[*p as usize] = g << BLK_SHIFT | tag;
                *p += 1;
            }
        }
        for k in 0..nbk {
            let run = &mut entries[start[k] as usize..start[k + 1] as usize];
            if run.len() > 1 {
                run.sort_by_key(|&e| e & BLK_MASK);
            }
        }
        let recs: Vec<u64> = (0..nbk + 2).map(|k| (occ[k] as u64) << 32 | start[k.min(nbk)] as u64).collect();
        t.recs = HugeBuf::from_vec(recs);
        t.entries = HugeBuf::from_vec(entries);
        t
    }
}

/// A batch of queries laid out for the scan.
struct Batch {
    /// Offset of residue 0 of every query in the concatenated arrays.
    base: Vec<usize>,
    lens: Vec<usize>,
    /// Window-filter bytes: the alignment code of every unmasked residue in
    /// the seed alphabet, QSENT for masked/ambiguous residues and padding (the
    /// seed-alphabet code is derived on the fly, so one window load per hit
    /// serves both identity counts).
    win: HugeBuf<u8>,
    /// Seed-alphabet map as byte-shuffle tables (codes 0..15 and 16..23).
    rtab_lo: [u8; 16],
    rtab_hi: [u8; 16],
    /// Alignment code | (composition bias as u8) << 8 per position (X-drop).
    qcb: HugeBuf<u16>,
    /// Substitution scores as i16 (vectorized X-drop).
    m16: [i16; 768],
    /// Query containing batch position `k << QCOARSE_SHIFT` (start of the
    /// search for the query of a position; see [`Batch::query_of`]).
    qcoarse: Vec<u32>,
    region_thr: Vec<i32>,
    tables: Vec<Table>,
}

/// Return freed heap memory to the system (glibc keeps the many small
/// buffers freed after building the seed tables).
fn release_free_memory() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    // SAFETY: malloc_trim only releases unused heap pages.
    unsafe {
        libc::malloc_trim(0);
    }
}

/// Debug: current and peak RSS (MB) at a stage of the scan (SABER_MEMTRACE).
fn memtrace(stage: &str) {
    if std::env::var_os("SABER_MEMTRACE").is_none() {
        return;
    }
    let s = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    let get = |k: &str| s.lines().find(|l| l.starts_with(k)).and_then(|l| l.split_whitespace().nth(1)).and_then(|v| v.parse::<u64>().ok()).unwrap_or(0) / 1024;
    eprintln!("[memtrace] {stage}: rss {} MB (anon {} MB, file {} MB) peak {} MB", get("VmRSS:"), get("RssAnon:"), get("RssFile:"), get("VmHWM:"));
}

/// Time-stamp counter (cycle-level stage accounting of the scan).
#[inline(always)]
fn tsc() -> u64 {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: rdtsc has no preconditions.
    unsafe {
        std::arch::x86_64::_rdtsc()
    }
    #[cfg(not(target_arch = "x86_64"))]
    0
}

/// Stable sort of survivors by batch position: insertion sort for short
/// runs, else two 11-bit (plus a third if needed) LSD radix passes.
fn sort_by_position(v: &mut [Surv], tmp: &mut Vec<Surv>) {
    let n = v.len();
    if n <= 32 {
        for i in 1..n {
            let x = v[i];
            let k = x.g;
            let mut j = i;
            while j > 0 && v[j - 1].g > k {
                v[j] = v[j - 1];
                j -= 1;
            }
            v[j] = x;
        }
        return;
    }
    let maxg = v.iter().map(|x| x.g).max().unwrap_or(0);
    let bits = 32 - maxg.leading_zeros();
    tmp.clear();
    tmp.resize(n, Surv::default());
    let mut shift = 0u32;
    let mut in_v = true;
    while shift < bits.max(1) {
        let (src, dst): (&[Surv], &mut [Surv]) = if in_v { (&*v, &mut tmp[..]) } else { (&tmp[..], &mut *v) };
        let mut cnt = [0u32; 2048];
        for x in src {
            cnt[((x.g >> shift) & 2047) as usize] += 1;
        }
        let mut acc = 0u32;
        for c in cnt.iter_mut() {
            let t = *c;
            *c = acc;
            acc += t;
        }
        for x in src {
            let b = ((x.g >> shift) & 2047) as usize;
            dst[cnt[b] as usize] = *x;
            cnt[b] += 1;
        }
        in_v = !in_v;
        shift += 11;
    }
    if !in_v {
        v.copy_from_slice(&tmp[..n]);
    }
}

/// Survivors between the prefetch of their query codes and their X-drop.
const XDROP_AHEAD: usize = 8;

/// Granularity of [`Batch::qcoarse`] (positions).
const QCOARSE_SHIFT: u32 = 8;

impl Batch {
    /// Query whose residues include batch position `g`.
    #[inline(always)]
    fn query_of(&self, g: usize) -> usize {
        // SAFETY: g is a batch position (< total), so g >> shift < qcoarse.len().
        let mut q = unsafe { *self.qcoarse.get_unchecked(g >> QCOARSE_SHIFT) } as usize;
        while q + 1 < self.base.len() && self.base[q + 1] <= g {
            q += 1;
        }
        q
    }
}

/// A candidate region as held during the scan (half the size of [`Region`]).
#[derive(Clone, Copy)]
struct RegionC {
    subject: u32,
    d_min: i32,
    d_max: i32,
    best: i32,
}

impl From<Region> for RegionC {
    fn from(r: Region) -> Self {
        RegionC { subject: r.subject as u32, d_min: r.d_min as i32, d_max: r.d_max as i32, best: r.best }
    }
}

impl From<RegionC> for Region {
    fn from(r: RegionC) -> Self {
        Region { subject: r.subject as u64, d_min: r.d_min as i64, d_max: r.d_max as i64, best: r.best }
    }
}

/// [`region_order`] for compact regions.
fn region_order_c(a: &RegionC, b: &RegionC) -> std::cmp::Ordering {
    b.best.cmp(&a.best).then(a.subject.cmp(&b.subject)).then(a.d_min.cmp(&b.d_min)).then(a.d_max.cmp(&b.d_max))
}

/// A window-filter survivor: batch position and chunk-local residue offset
/// (the subject is recovered from the offset when the chunk is visited).
#[derive(Clone, Copy, Default)]
struct Surv {
    g: u32,
    t: u32,
}

/// A unit of scan work: subjects `[lo, hi)` against batch positions
/// `[g_lo, g_hi)` (whole batch: `0..u32::MAX`).
#[derive(Clone, Copy, Debug)]
struct Piece {
    lo: usize,
    hi: usize,
    g_lo: u32,
    g_hi: u32,
}

impl Piece {
    fn whole(lo: usize, hi: usize) -> Self {
        Piece { lo, hi, g_lo: 0, g_hi: u32::MAX }
    }
    fn restricted(&self) -> bool {
        self.g_lo > 0 || self.g_hi < u32::MAX
    }
}

/// Pieces scanned together, shape by shape from `s0`, with their survivors
/// (each piece's survivors stay in memory until its last shape, then the
/// piece is visited).
struct Item {
    pieces: Vec<Piece>,
    surv: Vec<Vec<Surv>>,
    s0: usize,
}

/// Per-thread scan state.
struct ScanState {
    w: Worker,
    /// Candidate regions per query (compact: they are held for the whole scan).
    regions: Vec<Vec<RegionC>>,
    surv: Vec<Surv>,
    sorted: Vec<Surv>,
    counts: Vec<u32>,
    hits: Vec<(u32, u32)>,
    /// Positions of one subject with table matches: (position, entry range,
    /// code tag), then (position, range in `hl`, -).
    pend: Vec<(u32, u32, u32, u32)>,
    /// The same for the next subject (pipelined one subject ahead).
    pend2: Vec<(u32, u32, u32, u32)>,
    /// Batch positions of the hits of one subject, grouped by position.
    hl: Vec<u32>,
    /// (first, end, query) of every (subject, query) survivor group.
    groups: Vec<(u32, u32, u32)>,
    /// Positions of one subject whose seed code is present: (position, code).
    cand: Vec<(u32, u32)>,
    /// Seed code (or the invalid code) per position.
    codes: Vec<u32>,
    /// Radix sort scratch.
    tmp: Vec<Surv>,
    /// Visit scratch: chunk-local subject starts, the subject of every
    /// 32-residue block of the chunk, the subject of every survivor.
    starts: Vec<u32>,
    blk: Vec<u32>,
    sidx: Vec<u32>,
}

/// Seed codes of positions `0..np` of one subject (seed bytes `sr`):
/// `inv` for windows with a masked or non-seed residue.
fn codes_scalar(sr: &[u8], np: usize, offsets: &[u16], asize: u8, a: u32, inv: u32, codes: &mut Vec<u32>) {
    codes.clear();
    codes.extend((0..np).map(|p| {
        let (code, ok) = seed_code_flag(sr, offsets, p, asize, a);
        if ok {
            code
        } else {
            inv
        }
    }));
}

/// [`codes_scalar`] with AVX2, eight positions per vector. `base` is the
/// offset of the subject's first seed byte in `sred` (reads may run up to
/// 7 + span bytes past the subject, inside the buffer's padding or the next
/// subject; those lanes are not kept).
#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn codes_avx2(sred: &[u8], base: usize, np: usize, offsets: &[u16], asize: u8, inv: u32, codes: &mut Vec<u32>) {
    use std::arch::x86_64::*;
    codes.clear();
    codes.reserve(np + 8);
    let out = codes.as_mut_ptr();
    let av = _mm256_set1_epi32(asize as i32);
    let am1 = _mm256_set1_epi32(asize as i32 - 1);
    let invv = _mm256_set1_epi32(inv as i32);
    let ptr = sred.as_ptr().add(base);
    let mut p = 0usize;
    while p < np {
        let mut code = _mm256_setzero_si256();
        let mut bad = _mm256_setzero_si256();
        for &o in offsets {
            let v = _mm256_cvtepu8_epi32(_mm_loadl_epi64(ptr.add(p + o as usize) as *const __m128i));
            bad = _mm256_or_si256(bad, _mm256_cmpgt_epi32(v, am1));
            code = _mm256_add_epi32(_mm256_mullo_epi32(code, av), v);
        }
        _mm256_storeu_si256(out.add(p) as *mut __m256i, _mm256_blendv_epi8(code, invv, bad));
        p += 8;
    }
    codes.set_len(np);
}

/// Seed codes of positions `0..np` of seed bytes `sb` (seed-alphabet code;
/// any byte >= `asize` makes a window invalid, code `inv`). `sb` must be
/// readable 8 + span bytes past position `np` (the vector path reads whole
/// vectors). Used by the index builder's seed statistics.
pub(crate) fn seed_codes_into(sb: &[u8], np: usize, offsets: &[u16], asize: u8, inv: u32, out: &mut Vec<u32>) {
    #[cfg(target_arch = "x86_64")]
    if is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 detected; reads stay inside `sb` (caller's padding).
        unsafe { seed_codes_avx2(sb, np, offsets, asize, inv, out) };
        return;
    }
    codes_scalar(sb, np, offsets, asize, asize as u32, inv, out);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn seed_codes_avx2(sb: &[u8], np: usize, offsets: &[u16], asize: u8, inv: u32, out: &mut Vec<u32>) {
    codes_avx2(sb, 0, np, offsets, asize, inv, out);
}

/// Positions whose code is present in `tab`, as (position, code) in `cand`
/// (branch-free compaction); returns the count.
#[inline(always)]
fn present_positions(tab: &Table, codes: &[u32], cand: &mut Vec<(u32, u32)>) -> usize {
    cand.clear();
    cand.reserve(codes.len());
    let out = cand.as_mut_ptr();
    let mut n = 0usize;
    // SAFETY: n <= p < codes.len() <= capacity; codes are < inv or inv, and
    // both lookup arrays cover the invalid code (always absent).
    unsafe {
        if tab.dense {
            let bits = tab.bits.as_ptr();
            for (p, &c) in codes.iter().enumerate() {
                let w = *bits.add((c >> 6) as usize);
                *out.add(n) = (p as u32, c);
                n += ((w >> (c & 63)) & 1) as usize;
            }
        } else {
            let recs = tab.recs.as_ptr();
            for (p, &c) in codes.iter().enumerate() {
                let (k, tag) = tab.split(c);
                let r = *recs.add(k);
                *out.add(n) = (p as u32, c);
                n += ((r >> (32 + tag)) & 1) as usize;
            }
        }
        cand.set_len(n);
    }
    n
}

/// Batch positions of the entries `[ea, eb)` of a block that carry `tag`,
/// written to `out` (room for eb - ea values); returns their number.
///
/// # Safety
/// `out` must be valid for eb - ea writes and `ea..eb` inside `entries`.
#[inline(always)]
unsafe fn expand_run(entries: &[u32], ea: usize, eb: usize, tag: u32, out: *mut u32) -> usize {
    let mut n = 0usize;
    for k in ea..eb {
        let e = *entries.get_unchecked(k);
        *out.add(n) = e >> BLK_SHIFT;
        n += ((e & BLK_MASK) == tag) as usize;
    }
    n
}

/// [`expand_run`] with AVX2: a block's entries are sorted by tag, so the
/// entries of one tag are contiguous; eight entries per step, the matching
/// lanes moved to the front with one permute, stopping after the run.
///
/// # Safety
/// AVX2 available; `out` valid for eb - ea + 8 writes; `entries` readable up
/// to eb + 8 (the table's padding).
#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn expand_run_avx2(entries: *const u32, ea: usize, eb: usize, tag: u32, out: *mut u32) -> usize {
    use std::arch::x86_64::*;
    let tv = _mm256_set1_epi32(tag as i32);
    let tmask = _mm256_set1_epi32(BLK_MASK as i32);
    let lanes = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    let mut n = 0usize;
    let mut k = ea;
    while k < eb {
        let e = _mm256_loadu_si256(entries.add(k) as *const __m256i);
        let t = _mm256_and_si256(e, tmask);
        let valid = if eb - k >= 8 { 0xFFu32 } else { (1u32 << (eb - k)) - 1 };
        let eq = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(t, tv))) as u32 & valid;
        let gt = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(t, tv))) as u32 & valid;
        if eq != 0 {
            let first = eq.trailing_zeros() as i32;
            let moved = _mm256_permutevar8x32_epi32(_mm256_srli_epi32::<{ BLK_SHIFT as i32 }>(e), _mm256_add_epi32(lanes, _mm256_set1_epi32(first)));
            _mm256_storeu_si256(out.add(n) as *mut __m256i, moved);
            n += eq.count_ones() as usize;
        }
        if gt != 0 {
            break;
        }
        k += 8;
    }
    n
}

/// Branch-free seed code: (code, valid). An invalid window (a residue that
/// is masked or outside the seed alphabet) yields code 0.
#[inline(always)]
fn seed_code_flag(sr: &[u8], offsets: &[u16], p: usize, asize: u8, a: u32) -> (u32, bool) {
    let mut c = 0u32;
    let mut bad = false;
    for &o in offsets {
        // SAFETY: callers pass p <= len - span and every offset < span.
        let r = unsafe { *sr.get_unchecked(p + o as usize) };
        bad |= r >= asize;
        c = c.wrapping_mul(a).wrapping_add(r as u32);
    }
    (if bad { 0 } else { c }, !bad)
}

#[inline]
fn seed_code(sr: &[u8], offsets: &[u16], p: usize, asize: u8, a: u32) -> Option<u32> {
    let mut c = 0u32;
    for &o in offsets {
        // SAFETY: callers pass p <= len - span and every offset < span.
        let r = unsafe { *sr.get_unchecked(p + o as usize) };
        if r >= asize {
            return None;
        }
        c = c * a + r as u32;
    }
    Some(c)
}

/// Ungapped X-drop through (qp, sp) with the batch's query codes and bias;
/// identical results to [`ungapped_xdrop_prof`] for the same query. The
/// rightward and leftward extensions are independent and advance in lockstep
/// (two dependency chains in flight instead of one).
#[inline]
fn xdrop_cat(mflat: &[i32; 768], qcb: &[u16], s: &[u8], qp: usize, sp: usize, x: i32, steps: &mut u64) -> (i32, usize, usize) {
    let qlen = qcb.len();
    let nr = (qlen - qp).min(s.len() - sp);
    let nl = qp.min(sp);
    // SAFETY (all score reads): right reads qp + t < qlen and sp + t < s.len()
    // for t < nr; left reads qp - t and sp - t for 1 <= t <= nl; codes < 32
    // index the 24 x 32 table. qcb = code | (bias as u8) << 8.
    let sc = |i: usize, j: usize| -> i32 {
        unsafe {
            let v = *qcb.get_unchecked(i);
            *mflat.get_unchecked(((v & 0xFF) as usize) * 32 + (*s.get_unchecked(j) as usize & 31)) + (v >> 8) as u8 as i8 as i32
        }
    };
    let (mut best_r, mut cur_r, mut end_t, mut tr) = (0i32, 0i32, 0usize, 0usize);
    let (mut best_l, mut cur_l, mut st_t, mut tl) = (0i32, 0i32, 0usize, 1usize);
    let mut go_r = nr > 0;
    let mut go_l = nl > 0;
    while go_r && go_l {
        cur_r += sc(qp + tr, sp + tr);
        tr += 1;
        end_t = if cur_r > best_r { tr } else { end_t };
        best_r = best_r.max(cur_r);
        go_r = tr < nr && cur_r >= best_r - x;
        cur_l += sc(qp - tl, sp - tl);
        st_t = if cur_l > best_l { tl } else { st_t };
        best_l = best_l.max(cur_l);
        tl += 1;
        go_l = tl <= nl && cur_l >= best_l - x;
    }
    while go_r {
        cur_r += sc(qp + tr, sp + tr);
        tr += 1;
        end_t = if cur_r > best_r { tr } else { end_t };
        best_r = best_r.max(cur_r);
        go_r = tr < nr && cur_r >= best_r - x;
    }
    while go_l {
        cur_l += sc(qp - tl, sp - tl);
        st_t = if cur_l > best_l { tl } else { st_t };
        best_l = best_l.max(cur_l);
        tl += 1;
        go_l = tl <= nl && cur_l >= best_l - x;
    }
    *steps += (tr + tl - 1) as u64;
    (best_r + best_l, qp - st_t, qp + end_t)
}

/// 32-byte aligned 16 x u16 scratch.
#[repr(C, align(32))]
struct V16([u16; 16]);

/// Substitution scores as i16 for [`xdrop_avx2`] (same layout as `mflat`).
fn mflat16(mflat: &[i32; 768]) -> [i16; 768] {
    let mut m = [0i16; 768];
    for (d, &v) in m.iter_mut().zip(mflat.iter()) {
        *d = v.clamp(-1000, 1000) as i16;
    }
    m
}

/// One direction of [`xdrop_avx2`]: at most `n` steps from `q`/`s` (right:
/// offsets 0, 1, ...; left: offsets -1, -2, ...). Returns (best score, steps
/// to the best, steps taken).
///
/// Sixteen steps per vector: the scores of 16 residue pairs (scalar table
/// reads), their prefix sums and running maxima (relative to the score at the
/// block start, so 16-bit lanes cannot overflow), and the first lane where
/// the score falls more than `x` below the maximum. Lanes past `n` score
/// -1000, which ends the extension at the boundary like the scalar loop.
///
/// # Safety
/// AVX2 available; `q` and `s` readable 16 elements beyond the `n` steps in
/// the extension's direction (the padded batch and residue buffers).
#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn xdrop_dir_avx2<const RIGHT: bool>(m16: &[i16; 768], q: *const u16, s: *const u8, n: usize, x: i32) -> (i32, usize, usize) {
    use std::arch::x86_64::*;
    if n == 0 {
        return (0, 0, 0);
    }
    let lanes = _mm256_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
    let xv = _mm256_set1_epi16(x as i16);
    let bc7 = _mm256_setr_epi8(14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15);
    let mp = m16.as_ptr();
    let lo8 = _mm256_set1_epi16(0xFF);
    let m31 = _mm256_set1_epi16(31);
    // reverses the eight 16-bit lanes of each 128-bit half
    let rev = _mm256_setr_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);
    let (mut cur, mut best, mut end) = (0i32, 0i32, 0usize);
    let mut t0 = 0usize;
    let mut ibuf = V16([0u16; 16]);
    let mut buf = V16([0u16; 16]);
    loop {
        // Score table indices (query code * 32 + subject code) and the query
        // bias of 16 steps, in step order.
        let (qv, sv) = if RIGHT {
            (_mm256_loadu_si256(q.add(t0) as *const __m256i), _mm256_cvtepu8_epi16(_mm_loadu_si128(s.add(t0) as *const __m128i)))
        } else {
            let qv = _mm256_loadu_si256(q.sub(t0 + 16) as *const __m256i);
            let sv = _mm256_cvtepu8_epi16(_mm_loadu_si128(s.sub(t0 + 16) as *const __m128i));
            (
                _mm256_permute4x64_epi64::<0x4E>(_mm256_shuffle_epi8(qv, rev)),
                _mm256_permute4x64_epi64::<0x4E>(_mm256_shuffle_epi8(sv, rev)),
            )
        };
        let idx = _mm256_or_si256(_mm256_slli_epi16::<5>(_mm256_and_si256(qv, lo8)), _mm256_and_si256(sv, m31));
        _mm256_store_si256(ibuf.0.as_mut_ptr() as *mut __m256i, idx);
        // sixteen scalar table reads (a vector gather is slower here)
        for k in 0..16 {
            buf.0[k] = std::ptr::read_volatile(mp.add(*ibuf.0.get_unchecked(k) as usize)) as u16;
        }
        let mut sc = _mm256_add_epi16(_mm256_load_si256(buf.0.as_ptr() as *const __m256i), _mm256_srai_epi16::<8>(qv));
        let vk = n - t0;
        if vk < 16 {
            let inval = _mm256_cmpgt_epi16(lanes, _mm256_set1_epi16(vk as i16 - 1));
            sc = _mm256_blendv_epi8(sc, _mm256_set1_epi16(-1000), inval);
        }
        // prefix sums (in-lane, then the low half's total into the high half)
        let mut r = sc;
        r = _mm256_add_epi16(r, _mm256_slli_si256::<2>(r));
        r = _mm256_add_epi16(r, _mm256_slli_si256::<4>(r));
        r = _mm256_add_epi16(r, _mm256_slli_si256::<8>(r));
        r = _mm256_add_epi16(r, _mm256_shuffle_epi8(_mm256_permute2x128_si256::<0x08>(r, r), bc7));
        // running maxima, at least the best so far (relative); the zeros
        // shifted in are below it
        let base = best - cur;
        let mut m = r;
        m = _mm256_max_epi16(m, _mm256_slli_si256::<2>(m));
        m = _mm256_max_epi16(m, _mm256_slli_si256::<4>(m));
        m = _mm256_max_epi16(m, _mm256_slli_si256::<8>(m));
        m = _mm256_max_epi16(m, _mm256_shuffle_epi8(_mm256_permute2x128_si256::<0x08>(m, m), bc7));
        m = _mm256_max_epi16(m, _mm256_set1_epi16(base as i16));
        let stop = _mm256_movemask_epi8(_mm256_cmpgt_epi16(_mm256_sub_epi16(m, xv), r)) as u32;
        let last = if stop != 0 { (stop.trailing_zeros() / 2) as usize } else { 15 };
        let mut ma = [0i16; 16];
        _mm256_storeu_si256(ma.as_mut_ptr() as *mut __m256i, m);
        let bt = ma[last] as i32;
        if bt > base {
            // the best was first reached at the last strict improvement
            let eq = _mm256_movemask_epi8(_mm256_cmpeq_epi16(r, _mm256_set1_epi16(bt as i16))) as u32;
            end = t0 + (eq.trailing_zeros() / 2) as usize + 1;
            best = cur + bt;
        }
        if stop != 0 {
            return (best, end, t0 + last + 1);
        }
        let mut ra = [0i16; 16];
        _mm256_storeu_si256(ra.as_mut_ptr() as *mut __m256i, r);
        cur += ra[15] as i32;
        t0 += 16;
        if t0 >= n {
            return (best, end, t0);
        }
    }
}

/// [`xdrop_cat`] with AVX2 (identical results): `qcb` and `s` point at the
/// query's and the subject's first residue in buffers padded by at least 16
/// elements on both sides.
///
/// # Safety
/// AVX2 available; the padding as described; qp < qlen, sp < sl.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn xdrop_avx2(m16: &[i16; 768], qcb: *const u16, qlen: usize, s: *const u8, sl: usize, qp: usize, sp: usize, x: i32, steps: &mut u64) -> (i32, usize, usize) {
    let nr = (qlen - qp).min(sl - sp);
    let nl = qp.min(sp);
    let (best_r, end_t, tr) = xdrop_dir_avx2::<true>(m16, qcb.add(qp), s.add(sp), nr, x);
    let (best_l, st_t, tl) = xdrop_dir_avx2::<false>(m16, qcb.add(qp), s.add(sp), nl, x);
    *steps += (tr + tl) as u64;
    (best_r + best_l, qp - st_t, qp + end_t)
}

impl<'a> Engine<'a> {
    /// Use the scan schedule with this seeding scheme (protein databases).
    pub fn with_scan_seeds(mut self, seeds: ScanSeeds) -> Self {
        self.scan_seeds = Some(seeds);
        self
    }

    /// True when [`Self::search_scan`] can run with the current parameters.
    pub fn scan_available(&self) -> bool {
        self.scan_seeds.as_ref().map_or(false, |s| !s.shapes.is_empty())
            && self.is_protein
            && !self.params.exhaustive
            && self.params.two_hit_window == 0
            && self.params.seed_min_score == 0
            && self.nbr.iter().all(|g| g.is_none())
            && self.backend.is_none()
    }

    /// Over-frequent seed caps of the scan schedule (one per shape).
    pub fn scan_caps(&self) -> Vec<u64> {
        self.scan_db.get().map(|d| d.caps.clone()).unwrap_or_default()
    }

    /// (database seed-mask source, seconds spent building the scan data).
    pub fn scan_info(&self) -> (String, f64) {
        self.scan_db.get().map(|d| (d.mask_source.to_string(), d.build_s)).unwrap_or_default()
    }

    /// Database-side scan data (built on first use).
    pub(crate) fn scan_db(&self) -> &ScanDb<'a> {
        self.scan_db.get_or_init(|| self.build_scan_db())
    }

    fn build_scan_db(&self) -> ScanDb<'a> {
        let t0 = Instant::now();
        let seeds = self.scan_seeds.as_ref().expect("scan seeds configured");
        let idx = self.index;
        let res = idx.residues();
        let n = res.len();
        let map = seeds.alphabet.map;
        let mut seed_map = [NONSEED; 256];
        for (c, x) in seed_map.iter_mut().enumerate().take(24) {
            let r = map[c];
            *x = if r == 255 { NONSEED } else { r };
        }
        // Database seed mask: stored runs, or tantan recomputed exactly as the
        // index builder computes it.
        let masked_db = idx.kind() == SeqKind::Protein && idx.seed_config().masking == MaskingMode::Tantan;
        let mut mask_source = "none";
        let mut mask_runs: Vec<(u64, u64)> = Vec::new();
        if masked_db {
            match idx.mask_runs() {
                Some(runs) => {
                    mask_source = "index";
                    mask_runs = runs.chunks_exact(2).map(|r| (r[0], r[1])).collect();
                }
                None => {
                    mask_source = "recomputed";
                    let mm: crate::scoring::ScoringSystem = idx.meta().masking_matrix.parse().unwrap_or(crate::scoring::ScoringSystem::Blosum62);
                    let tt = Tantan::protein(&ScoreMatrix::new(mm));
                    let offs = idx.subject_offsets();
                    mask_runs = (0..idx.num_subjects())
                        .into_par_iter()
                        .map_init(
                            || (TantanScratch::default(), Vec::new()),
                            |(sc, m), i| {
                                let s = idx.subject_codes(i);
                                tt.mask(s, sc, m);
                                let base = offs[i];
                                let mut out = Vec::new();
                                let mut j = 0;
                                while j < m.len() {
                                    if m[j] {
                                        let st = j;
                                        while j < m.len() && m[j] {
                                            j += 1;
                                        }
                                        out.push((base + st as u64, base + j as u64));
                                    } else {
                                        j += 1;
                                    }
                                }
                                out
                            },
                        )
                        .flatten()
                        .collect();
                }
            }
            mask_runs.sort_unstable();
        }
        let raw = if self.params.hit_filter_exact_weight > 0 {
            match idx.residues_padded(PAD) {
                Some(v) => RawResidues::View(v),
                None => {
                    let mut v = HugeBuf::filled(n + 2 * PAD, 0xFFu8);
                    v[PAD..PAD + n].copy_from_slice(res);
                    RawResidues::Own(v)
                }
            }
        } else {
            RawResidues::None
        };
        // Over-frequent codes per shape (the index schedules' rule: skip query
        // seeds occurring more than the cap in the database).
        let offs = idx.subject_offsets();
        let ns = idx.num_subjects();
        let asize = seeds.alphabet.size;
        let a = asize as u32;
        let mut frequent = Vec::new();
        let mut caps = Vec::new();
        // Statistics stored in the index (built with the same seed alphabet
        // and shapes) save one database pass per shape.
        let icfg = idx.seed_config();
        // Statistics are stored per index shape: usable for every search
        // shape with the same pattern (same seed alphabet).
        let same_alphabet = icfg.alphabet.map == seeds.alphabet.map;
        for sh in seeds.shapes.iter() {
            let ncodes = (a as u64).pow(sh.weight() as u32);
            let stored_si = if same_alphabet { icfg.shapes.iter().position(|x| x.pattern == sh.pattern) } else { None };
            if let Some(isi) = stored_si {
                if let Some((total, list)) = idx.scan_freq(isi) {
                    let cap = if self.params.seed_max_occ > 0 {
                        self.params.seed_max_occ as u64
                    } else {
                        super::super::index::auto_occ_cap(total, ncodes)
                    };
                    let mut bits = Vec::new();
                    // The list holds the codes above the index's listing
                    // threshold: complete for caps at or above it.
                    let usable = cap >= idx.scan_freq_min();
                    for &x in list.iter().filter(|_| usable) {
                        let (code, n) = ((x >> 32) as usize, x & 0xFFFF_FFFF);
                        if n > cap && code < ncodes as usize {
                            if bits.is_empty() {
                                bits = vec![0u64; ((ncodes + 63) / 64) as usize];
                            }
                            bits[code >> 6] |= 1u64 << (code & 63);
                        }
                    }
                    if usable {
                        frequent.push(bits);
                        caps.push(cap);
                        continue;
                    }
                }
            }
            let (seed_map_ref, runs_ref) = (&seed_map, &mask_runs);
            let count_into = |lo: usize, hi: usize, cnt: &mut Vec<u32>| -> u64 {
                let mut tot = 0u64;
                let mut s0 = lo;
                while s0 < hi {
                    // blocks of about SUPER_RESIDUES residues
                    let mut s1 = s0 + 1;
                    while s1 < hi && (offs[s1 + 1] - offs[s0]) as usize <= SUPER_RESIDUES {
                        s1 += 1;
                    }
                    let rs = seed_bytes(seed_map_ref, runs_ref, res, offs[s0] as usize, offs[s1] as usize);
                    for s in s0..s1 {
                        let (b, e) = (offs[s] as usize - rs.base, offs[s + 1] as usize - rs.base);
                        if e - b < sh.span {
                            continue;
                        }
                        let sr = &rs.buf[PAD + b..PAD + e];
                        for p in 0..=sr.len() - sh.span {
                            if let Some(c) = seed_code(sr, &sh.offsets, p, asize, a) {
                                cnt[c as usize] += 1;
                                tot += 1;
                            }
                        }
                    }
                    s0 = s1;
                }
                tot
            };
            let nthreads = rayon::current_num_threads().max(1);
            let step = (ns + nthreads - 1) / nthreads.max(1);
            let parts: Vec<(Vec<u32>, u64)> = (0..nthreads)
                .into_par_iter()
                .map(|t| {
                    let lo = (t * step).min(ns);
                    let hi = ((t + 1) * step).min(ns);
                    let mut cnt = vec![0u32; ncodes as usize];
                    let tot = count_into(lo, hi, &mut cnt);
                    (cnt, tot)
                })
                .collect();
            let total: u64 = parts.iter().map(|p| p.1).sum();
            let cap = if self.params.seed_max_occ > 0 {
                self.params.seed_max_occ as u64
            } else {
                super::super::index::auto_occ_cap(total, ncodes)
            };
            let mut bits = Vec::new();
            if total > cap {
                let mut any = false;
                let mut bm = vec![0u64; ((ncodes + 63) / 64) as usize];
                for c in 0..ncodes as usize {
                    let tot: u64 = parts.iter().map(|p| p.0[c] as u64).sum();
                    if tot > cap {
                        bm[c >> 6] |= 1u64 << (c & 63);
                        any = true;
                    }
                }
                if any {
                    bits = bm;
                }
            }
            frequent.push(bits);
            caps.push(cap);
        }
        // Work units: consecutive subjects with about CHUNK_RESIDUES residues.
        let mut chunks = Vec::new();
        let mut lo = 0usize;
        while lo < ns {
            let mut hi = lo;
            let start = offs[lo];
            while hi < ns && (hi == lo || (offs[hi + 1] - start) as usize <= CHUNK_RESIDUES) {
                hi += 1;
            }
            chunks.push((lo, hi));
            lo = hi;
        }
        ScanDb { seed_map, mask_runs, raw, frequent, caps, chunks, mask_source, build_s: t0.elapsed().as_secs_f64() }
    }

    /// Search a batch of queries with the scan schedule. Results are
    /// identical to [`Self::search`] for the same seeds and parameters.
    pub fn search_scan(&self, queries: &[&[u8]]) -> Vec<(Vec<Hit>, Counters)> {
        let nq = queries.len();
        if nq == 0 {
            return Vec::new();
        }
        memtrace("scan start");
        let db = self.scan_db();
        memtrace("scan db built");
        let seeds = self.scan_seeds.as_ref().expect("scan seeds configured");
        let t_prep = Instant::now();
        // Per-query preparation (bias, mask, statistics).
        let preps: Vec<(QueryPrep, Counters)> = queries
            .par_iter()
            .map_init(Worker::new, |w, q| {
                let before = w.counters;
                w.counters.queries += 1;
                w.counters.query_residues += q.len() as u64;
                let p = if q.is_empty() {
                    QueryPrep { bias: Vec::new(), mask: Vec::new(), has_bias: false, qstats: self.stats.for_query(1), floor: i64::MAX, region_floor: i64::MAX }
                } else {
                    self.prepare(q, w)
                };
                let mut d = w.counters;
                d.sub(&before);
                (p, d)
            })
            .collect();
        let mut gctr = Counters::default();
        let b = self.build_batch(queries, &preps, seeds, db, &mut gctr);
        release_free_memory();
        memtrace("batch built");
        gctr.t_seed_ns += t_prep.elapsed().as_nanos() as u64;
        let cap = self.params.max_candidates.max(1);
        let max_q = b.lens.iter().copied().max().unwrap_or(0);
        let need = max_q + self.index.meta().max_subject_len as usize + 2;
        // Scan the database in rounds: a group of work units (about
        // SUPER_RESIDUES residues) is scanned shape by shape by all threads, so
        // that only one shape's seed table competes with the query windows for
        // the cache at any time; the survivors of the group (all shapes) are
        // then visited. One scratch state per thread.
        let nthreads = rayon::current_num_threads().max(1);
        let states: Vec<std::sync::Mutex<ScanState>> = (0..nthreads)
            .map(|_| {
                let mut w = Worker::new();
                w.diag_epoch.resize(need, 0);
                w.diag_cover.resize(need, 0);
                w.diag_last.resize(need, 0);
                std::sync::Mutex::new(ScanState {
                    w,
                    regions: vec![Vec::new(); nq],
                    surv: Vec::new(),
                    sorted: Vec::new(),
                    counts: Vec::new(),
                    hits: Vec::new(),
                    pend: Vec::new(),
                    pend2: Vec::new(),
                    hl: Vec::new(),
                    groups: Vec::new(),
                    cand: Vec::new(),
                    codes: Vec::new(),
                    tmp: Vec::new(),
                    starts: Vec::new(),
                    blk: Vec::new(),
                    sidx: Vec::new(),
                })
            })
            .collect();
        if self.index.num_subjects() > 0 && b.tables.iter().any(|t| t.len() > 0) {
            let offs = self.index.subject_offsets();
            let nch = db.chunks.len();
            let state = |stv: &[std::sync::Mutex<ScanState>]| -> usize { rayon::current_thread_index().unwrap_or(0) % stv.len() };
            let mut c0 = 0usize;
            // Start small (the survivor density is unknown) and grow.
            let mut round_res = 2 * CHUNK_RESIDUES;
            let budget = ROUND_SURVIVORS;
            let pair_limit = if self.params.scan_pair_survivors > 0 { self.params.scan_pair_survivors } else { PAIR_SURVIVORS };
            // survivors a piece may hold before its scan of a shape stops at
            // the next subject boundary
            let piece_limit = (2 * budget).min(pair_limit);
            let nshapes = seeds.shapes.len();
            let mut splits = 0u64;
            while c0 < nch {
                let start = offs[db.chunks[c0].0];
                let mut c1 = c0 + 1;
                while c1 < nch && (offs[db.chunks[c1].1] - start) as usize <= round_res {
                    c1 += 1;
                }
                // seed bytes of the round's residues
                let rs = db.round_seeds(self.index.residues(), start as usize, offs[db.chunks[c1 - 1].1] as usize);
                // Survivors are held until their piece has seen every shape.
                // An item whose projected survivors exceed the budget is cut
                // (by pieces, then by subjects, then by queries): the first
                // part is finished and visited, the others resume later with
                // their survivors. Every (subject, query) pair stays in one
                // piece, so the results do not depend on the cuts.
                let mut work = vec![Item {
                    pieces: (c0..c1).map(|w| Piece::whole(db.chunks[w].0, db.chunks[w].1)).collect(),
                    surv: (c0..c1).map(|_| Vec::new()).collect(),
                    s0: 0,
                }];
                let mut round_n = 0usize;
                let mut jobs: Vec<PairJob> = Vec::new();
                while let Some(mut it) = work.pop() {
                    debug_assert!(it.s0 < nshapes);
                    for si in it.s0..nshapes {
                        let pieces = &it.pieces;
                        let stops: Vec<Option<usize>> = it
                            .surv
                            .par_iter_mut()
                            .zip(pieces.par_iter())
                            .map(|(sv, pc)| {
                                let mut st = states[state(&states)].lock().expect("scan state");
                                let t0 = Instant::now();
                                let stop = self.scan_chunk_shapes(db, &rs, &b, seeds, pc, si..si + 1, &mut st, sv, piece_limit);
                                st.w.counters.t_extend_ns += t0.elapsed().as_nanos() as u64;
                                stop
                            })
                            .collect();
                        // A piece that stopped early keeps the subjects it
                        // scanned; the rest becomes a new item that resumes
                        // at this shape, with its survivors of the earlier
                        // shapes (offsets rebased to its first subject; no
                        // pair is split, subjects are whole).
                        for (k, stop) in stops.into_iter().enumerate() {
                            if let Some(sx) = stop {
                                let pc = it.pieces[k];
                                if sx < pc.hi {
                                    let shift = (offs[sx] - offs[pc.lo]) as u32;
                                    let sv = std::mem::take(&mut it.surv[k]);
                                    let (mut keep, mut rest) = (Vec::with_capacity(sv.len()), Vec::new());
                                    for x in sv {
                                        if x.t < shift {
                                            keep.push(x);
                                        } else {
                                            rest.push(Surv { g: x.g, t: x.t - shift });
                                        }
                                    }
                                    it.surv[k] = keep;
                                    work.push(Item { pieces: vec![Piece { lo: sx, ..pc }], surv: vec![rest], s0: si });
                                    it.pieces[k].hi = sx;
                                    splits += 1;
                                }
                            }
                        }
                        if si + 1 < nshapes {
                            splits += self.split_item(&mut it, si + 1, nshapes, budget, pair_limit, &b, &mut work, &mut jobs);
                        }
                    }
                    round_n += it.surv.iter().map(|v| v.len()).sum::<usize>();
                    let pieces = &it.pieces;
                    it.surv.par_iter_mut().zip(pieces.par_iter()).for_each(|(sv, pc)| {
                        let mut st = states[state(&states)].lock().expect("scan state");
                        let t0 = Instant::now();
                        self.visit_chunk(db, &b, pc.lo, pc.hi, cap, sv, &mut st);
                        st.w.counters.t_extend_ns += t0.elapsed().as_nanos() as u64;
                    });
                }
                // pairs too large for one piece: each scanned and visited on
                // its own, in query-position slices
                jobs.par_iter().for_each(|job| {
                    let mut st = states[state(&states)].lock().expect("scan state");
                    let t0 = Instant::now();
                    self.stream_pair(db, &rs, &b, seeds, *job, cap, pair_limit, &mut st);
                    st.w.counters.t_extend_ns += t0.elapsed().as_nanos() as u64;
                });
                // size the next round by this one's density
                if round_n > budget {
                    round_res = (round_res / 2).max(CHUNK_RESIDUES);
                } else if round_n < budget / 4 {
                    round_res = (round_res * 2).min(SUPER_RESIDUES);
                }
                c0 = c1;
            }
            gctr.scan_splits += splits;
        }
        // The seed tables and query arrays are not needed for scoring.
        drop(b);
        release_free_memory();
        memtrace("scan done");
        let states: Vec<ScanState> = states.into_iter().map(|m| m.into_inner().expect("scan state")).collect();
        // Merge per-query candidates (top `cap` under the global total order).
        let mut per_query: Vec<Vec<Region>> = (0..nq).map(|_| Vec::new()).collect();
        for st in states {
            gctr.add(&st.w.counters);
            for (q, r) in st.regions.into_iter().enumerate() {
                if !r.is_empty() {
                    per_query[q].extend(r.into_iter().map(Region::from));
                }
            }
        }
        let mut out: Vec<(Vec<Hit>, Counters)> = per_query
            .into_par_iter()
            .zip(queries.par_iter())
            .zip(preps.into_par_iter())
            .map_init(Worker::new, |w, ((mut regions, q), (prep, pctr))| {
                let before = w.counters;
                if q.is_empty() {
                    return (Vec::new(), pctr);
                }
                regions.sort_unstable_by(region_order);
                if regions.len() > cap {
                    w.counters.regions_dropped_by_cap += (regions.len() - cap) as u64;
                    regions.truncate(cap);
                }
                w.counters.regions += regions.len() as u64;
                let t2 = Instant::now();
                let windows = self.plan_windows(&regions, q.len());
                let scores = self.score_windows(q, &prep, &regions, &windows, w);
                w.counters.t_score_ns += t2.elapsed().as_nanos() as u64;
                let hits = self.finish(q, &prep, &regions, &windows, &scores, w);
                let mut d = w.counters;
                d.sub(&before);
                d.add(&pctr);
                (hits, d)
            })
            .collect();
        out[0].1.add(&gctr);
        memtrace("scored");
        out
    }

    /// Lay out a batch: concatenated query arrays and one seed table per shape.
    fn build_batch(&self, queries: &[&[u8]], preps: &[(QueryPrep, Counters)], seeds: &ScanSeeds, db: &ScanDb, ctr: &mut Counters) -> Batch {
        let nq = queries.len();
        let mut base = Vec::with_capacity(nq);
        let mut lens = Vec::with_capacity(nq);
        let mut cur = PAD;
        for q in queries {
            base.push(cur);
            lens.push(q.len());
            cur += q.len() + PAD;
        }
        let total = cur;
        let mut red = vec![QSENT; total];
        let mut win = HugeBuf::filled(total, QSENT);
        let mut qcb = HugeBuf::filled(total, crate::alphabet::AA_X as u16);
        let map = &seeds.alphabet.map;
        for (qi, q) in queries.iter().enumerate() {
            let o = base[qi];
            let prep = &preps[qi].0;

            for (i, &c) in q.iter().enumerate() {
                qcb[o + i] = c as u16 | ((prep.bias.get(i).copied().unwrap_or(0) as u8 as u16) << 8);
                if !prep.mask.get(i).copied().unwrap_or(false) {
                    let r = map[(c as usize).min(23)];
                    if r != 255 {
                        red[o + i] = r;
                        win[o + i] = c;
                    }
                }
            }
        }
        let region_thr: Vec<i32> = preps.iter().map(|(p, _)| self.hsp_threshold_for(p)).collect();
        // Seed tables.
        let asize = seeds.alphabet.size;
        let a = asize as u32;
        let mut tables = Vec::with_capacity(seeds.shapes.len());
        for (si, sh) in seeds.shapes.iter().enumerate() {
            let ncodes = (a as u64).pow(sh.weight() as u32);
            let freq = &db.frequent[si];
            // (code, batch position) of every kept seed, query by query
            // (batch positions increase, so entries end up sorted per code).
            let per_q: Vec<(Vec<(u32, u32)>, u64)> = queries
                .par_iter()
                .enumerate()
                .map(|(qi, q)| {
                    let mut v = Vec::new();
                    let mut skipped = 0u64;
                    if q.len() >= sh.span {
                        let o = base[qi];
                        let r = &red[o..o + q.len()];
                        for p in 0..=q.len() - sh.span {
                            if let Some(c) = seed_code(r, &sh.offsets, p, asize, a) {
                                if !freq.is_empty() && (freq[(c >> 6) as usize] >> (c & 63)) & 1 == 1 {
                                    skipped += 1;
                                    continue;
                                }
                                v.push((c, (o + p) as u32));
                            }
                        }
                    }
                    (v, skipped)
                })
                .collect();
            let nseeds: usize = per_q.iter().map(|x| x.0.len()).sum();
            ctr.seed_positions += nseeds as u64;
            ctr.seeds_frequent_skipped += per_q.iter().map(|x| x.1).sum::<u64>();
            // Dense offsets when the batch fills a good part of the code space;
            // otherwise the blocked layout keeps the table cache sized.
            let table = if ncodes <= DENSE_MAX_CODES && ncodes <= 2 * nseeds as u64 {
                Table::dense(ncodes, &per_q, nseeds)
            } else {
                Table::blocked(ncodes, &per_q, nseeds)
            };
            tables.push(table);
        }
        let mut rtab_lo = [0xFFu8; 16];
        let mut rtab_hi = [0xFFu8; 16];
        for cc in 0..24usize {
            let r = map[cc];
            if cc < 16 {
                rtab_lo[cc] = r;
            } else {
                rtab_hi[cc - 16] = r;
            }
        }
        let qcoarse: Vec<u32> = (0..(total >> QCOARSE_SHIFT) + 1)
            .map(|k| (base.partition_point(|&x| x <= k << QCOARSE_SHIFT).max(1) - 1) as u32)
            .collect();
        Batch { base, lens, win, rtab_lo, rtab_hi, qcb, m16: mflat16(&self.mflat), qcoarse, region_thr, tables }
    }

    /// Seed lookups and window filter of subjects `[lo, hi)` for the given
    /// shapes; survivors are appended to `surv`. Returns (hits, hits filtered).
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::too_many_arguments)]
    fn scan_hits<const AVX: bool>(
        &self,
        db: &ScanDb,
        rs: &RoundSeeds,
        b: &Batch,
        seeds: &ScanSeeds,
        pc: &Piece,
        shapes: std::ops::Range<usize>,
        st: &mut ScanState,
        surv: &mut Vec<Surv>,
        limit: usize,
    ) -> (u64, u64, Option<usize>) {
        let (lo, hi) = (pc.lo, pc.hi);
        let offs = self.index.subject_offsets();
        let asize = seeds.alphabet.size;
        let a = asize as u32;
        let min_ids = self.params.hit_filter_min_ids;
        let kx = self.params.hit_filter_exact_weight;
        let c = seeds.shapes.iter().map(|s| s.span).min().unwrap_or(1) / 2;
        let (mut hits, mut filtered) = (0u64, 0u64);
        let mut stop: Option<usize> = None;
        for (si, sh) in seeds.shapes.iter().enumerate().skip(shapes.start).take(shapes.len()) {
            let tab = &b.tables[si];
            if tab.len() == 0 {
                continue;
            }
            // Two-stage software pipeline over the subjects. Stage A of a
            // subject — seed codes, presence (one bitmap word or block record
            // per position), entry ranges with a prefetch of the entries —
            // runs one subject ahead of its stage B — the hits' batch
            // positions from the entries, window filter — so that the entries
            // (a random access each) arrive while the previous subject is
            // filtered.
            let mut prev: Option<(usize, usize, usize)> = None;
            let mut s = lo;
            loop {
                let t_a = tsc();
                let mut next: Option<(usize, usize, usize)> = None;
                while s < hi {
                    let (sb, se) = (offs[s] as usize, offs[s + 1] as usize);
                    let sl = se - sb;
                    s += 1;
                    if sl < sh.span {
                        continue;
                    }
                    let np = sl - sh.span + 1;
                    #[cfg(target_arch = "x86_64")]
                    if AVX {
                        // SAFETY: AVX = true only from the AVX2 wrapper; the
                        // loads stay inside the padded residue buffer.
                        unsafe { codes_avx2(&rs.buf, PAD + sb - rs.base, np, &sh.offsets, asize, tab.inv, &mut st.codes) };
                    } else {
                        codes_scalar(&rs.buf[PAD + sb - rs.base..PAD + se - rs.base], np, &sh.offsets, asize, a, tab.inv, &mut st.codes);
                    }
                    #[cfg(not(target_arch = "x86_64"))]
                    codes_scalar(&rs.buf[PAD + sb - rs.base..PAD + se - rs.base], np, &sh.offsets, asize, a, tab.inv, &mut st.codes);
                    let ncand = present_positions(tab, &st.codes, &mut st.cand);
                    st.w.counters.scan_candidates += ncand as u64;
                    let t_r = tsc();
                    st.w.counters.c_scan_cand += t_r - t_a;
                    // Entry ranges.
                    st.pend2.clear();
                    if tab.dense {
                        #[cfg(target_arch = "x86_64")]
                        for &(_, code) in &st.cand {
                            // SAFETY: prefetch hint inside the offset table.
                            unsafe { std::arch::x86_64::_mm_prefetch(tab.offs.as_ptr().add(code as usize) as *const i8, std::arch::x86_64::_MM_HINT_T0) };
                        }
                        for &(p, code) in &st.cand {
                            // SAFETY: present codes are < ncodes; offs has ncodes + 1 entries.
                            let (ea, eb) = unsafe { (*tab.offs.get_unchecked(code as usize), *tab.offs.get_unchecked(code as usize + 1)) };
                            st.pend2.push((p, ea, eb, 0));
                            #[cfg(target_arch = "x86_64")]
                            // SAFETY: prefetch hint inside the entries.
                            unsafe {
                                std::arch::x86_64::_mm_prefetch(tab.entries.as_ptr().add(ea as usize) as *const i8, std::arch::x86_64::_MM_HINT_T0)
                            };
                        }
                    } else {
                        for &(p, code) in &st.cand {
                            let (k, tag) = tab.split(code);
                            // SAFETY: blocks are < nb and recs has nb + 2 records.
                            let (ea, eb) = unsafe { (*tab.recs.get_unchecked(k) as u32, *tab.recs.get_unchecked(k + 1) as u32) };
                            st.pend2.push((p, ea, eb, tag));
                            #[cfg(target_arch = "x86_64")]
                            {
                                // The whole run of the block (all its lines).
                                let ep = tab.entries.as_ptr();
                                let mut k = ea as usize & !15;
                                while k < eb as usize {
                                    // SAFETY: prefetch hint inside the entries.
                                    unsafe { std::arch::x86_64::_mm_prefetch(ep.add(k) as *const i8, std::arch::x86_64::_MM_HINT_T0) };
                                    k += 16;
                                }
                            }
                        }
                    }
                    st.w.counters.c_scan_range += tsc() - t_r;
                    next = Some((sb, sl, s - 1));
                    break;
                }
                if let Some((sb, sl, ps)) = prev {
                    // Batch positions of the hits, grouped by subject position.
                    let t_x = tsc();
                    st.hl.clear();
                    if pc.restricted() {
                        // a piece of the batch: a code's entries are in batch
                        // order (within its tag's run in a block), so the
                        // positions in [g_lo, g_hi) are found by binary search
                        for it in st.pend.iter_mut() {
                            let (ea, eb) = (it.1 as usize, it.2 as usize);
                            let h0 = st.hl.len();
                            let run = &tab.entries[ea..eb];
                            if tab.dense {
                                let a = run.partition_point(|&g| g < pc.g_lo);
                                let z = run.partition_point(|&g| g < pc.g_hi);
                                st.hl.extend_from_slice(&run[a..z.max(a)]);
                            } else {
                                let tag = it.3;
                                let ta = run.partition_point(|&e| (e & BLK_MASK) < tag);
                                let tb = run.partition_point(|&e| (e & BLK_MASK) <= tag);
                                let sub = &run[ta..tb.max(ta)];
                                let a = sub.partition_point(|&e| (e >> BLK_SHIFT) < pc.g_lo);
                                let z = sub.partition_point(|&e| (e >> BLK_SHIFT) < pc.g_hi);
                                st.hl.extend(sub[a..z.max(a)].iter().map(|&e| e >> BLK_SHIFT));
                            }
                            *it = (it.0, h0 as u32, st.hl.len() as u32, 0);
                        }
                    } else {
                        for it in st.pend.iter_mut() {
                            let (ea, eb) = (it.1 as usize, it.2 as usize);
                            let h0 = st.hl.len();
                            if tab.dense {
                                st.hl.extend_from_slice(&tab.entries[ea..eb]);
                            } else {
                                let tag = it.3;
                                st.hl.reserve(eb - ea + 8);
                                // SAFETY: capacity reserved for eb - ea + 8 more (the
                                // vector path stores whole vectors); ea..eb + 8 is
                                // inside the padded entries.
                                unsafe {
                                    let out = st.hl.as_mut_ptr().add(h0);
                                    #[cfg(target_arch = "x86_64")]
                                    let n = if AVX { expand_run_avx2(tab.entries.as_ptr(), ea, eb, tag, out) } else { expand_run(&tab.entries, ea, eb, tag, out) };
                                    #[cfg(not(target_arch = "x86_64"))]
                                    let n = expand_run(&tab.entries, ea, eb, tag, out);
                                    st.hl.set_len(h0 + n);
                                }
                            }
                            *it = (it.0, h0 as u32, st.hl.len() as u32, 0);
                        }
                    }
                    if shapes.len() == 1 && ps > lo && surv.len() + st.hl.len() > limit {
                        // the piece is full: stop before this subject (not
                        // its first, so that the piece keeps at least one;
                        // the caller scans the rest, from it, separately)
                        stop = Some(ps);
                        st.w.counters.c_scan_expand += tsc() - t_x;
                        break;
                    }
                    hits += st.hl.len() as u64;
                    let t_c = tsc();
                    // chunk-local offset of the subject's first residue
                    let t_base = (sb - offs[lo] as usize) as u32;
                    if min_ids == 0 {
                        for &(p, h0, h1, _) in &st.pend {
                            for &g in &st.hl[h0 as usize..h1 as usize] {
                                surv.push(Surv { g, t: t_base + p });
                            }
                        }
                    } else {
                        #[cfg(target_arch = "x86_64")]
                        if AVX {
                            // SAFETY: AVX = true only from the AVX2 wrapper; window
                            // loads stay inside the padded buffers.
                            unsafe { self.filter_avx2(db, rs, b, &st.hl, &st.pend, sb, sl, c, kx, min_ids, t_base, surv, &mut filtered) };
                        } else {
                            self.filter_scalar(db, rs, b, seeds, &st.hl, &st.pend, sb, sl, c, kx, min_ids, t_base, surv, &mut filtered);
                        }
                        #[cfg(not(target_arch = "x86_64"))]
                        self.filter_scalar(db, rs, b, seeds, &st.hl, &st.pend, sb, sl, c, kx, min_ids, t_base, surv, &mut filtered);
                    }
                    let t_d = tsc();
                    st.w.counters.c_scan_expand += t_c - t_x;
                    st.w.counters.c_scan_filter += t_d - t_c;
                }
                let Some(nx) = next else { break };
                if let Some((_, _, ps)) = prev {
                    // Too many survivors held: stop after this subject; the
                    // caller rescans the rest of the piece (from the subject
                    // whose stage A just ran) as a separate item.
                    if surv.len() > limit && shapes.len() == 1 {
                        stop = Some(ps + 1);
                        break;
                    }
                }
                std::mem::swap(&mut st.pend, &mut st.pend2);
                prev = Some(nx);
            }
        }
        (hits, filtered, stop)
    }

    /// Window filter (portable): seed-alphabet identities (+ weighted exact
    /// identities) between every hit's query window and its subject window.
    #[allow(clippy::too_many_arguments)]
    fn filter_scalar(
        &self,
        db: &ScanDb,
        rs: &RoundSeeds,
        b: &Batch,
        seeds: &ScanSeeds,
        hl: &[u32],
        pend: &[(u32, u32, u32, u32)],
        sb: usize,
        sl: usize,
        c: usize,
        kx: u32,
        min_ids: u32,
        t_base: u32,
        surv: &mut Vec<Surv>,
        filtered: &mut u64,
    ) {
        for &(p, h0, h1, _) in pend {
            let ss = p as i64 + c as i64 - 16;
            let lo_l = (-ss).max(0) as u32;
            let hi_l = (sl as i64 - ss).min(32) as u32;
            let lane = (if hi_l >= 32 { u32::MAX } else { (1u32 << hi_l) - 1 }) & !((1u32 << lo_l) - 1);
            let g0 = (PAD as i64 + sb as i64 + ss) as usize;
            for &g in &hl[h0 as usize..h1 as usize] {
                let qo = g as usize + c - 16;
                let mut r = 0u32;
                let mut e = 0u32;
                for t in 0..32 {
                    if (lane >> t) & 1 == 0 {
                        continue;
                    }
                    let qc = b.win[qo + t];
                    let qr = if qc < 24 { seeds.alphabet.map[qc as usize] } else { QSENT };
                    if qr != 255 && qr == rs.buf[g0 - rs.base + t] & 0x7F {
                        r += 1;
                        if kx > 0 && qc == db.raw[g0 + t] {
                            e += 1;
                        }
                    }
                }
                if r + kx * e >= min_ids {
                    surv.push(Surv { g, t: t_base + p });
                } else {
                    *filtered += 1;
                }
            }
        }
    }

    /// [`Self::scan_hits`] compiled with AVX2 (the window filter inlined).
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2,popcnt")]
    #[allow(clippy::too_many_arguments)]
    unsafe fn scan_hits_avx2(
        &self,
        db: &ScanDb,
        rs: &RoundSeeds,
        b: &Batch,
        seeds: &ScanSeeds,
        pc: &Piece,
        shapes: std::ops::Range<usize>,
        st: &mut ScanState,
        surv: &mut Vec<Surv>,
        limit: usize,
    ) -> (u64, u64, Option<usize>) {
        self.scan_hits::<true>(db, rs, b, seeds, pc, shapes, st, surv, limit)
    }

    /// After `done` of `nshapes` shapes, cut an item whose projected
    /// survivors (so far x nshapes / done) exceed `budget`: into groups of
    /// pieces, else (one piece) into subject ranges, else (one subject) into
    /// query ranges. The parts after the first become work items, with their
    /// survivors, that resume at shape `done`. A query whose survivors alone
    /// are projected above `pair_limit` becomes a [`PairJob`] (its survivors
    /// are dropped and its positions excluded from the item; the pair is
    /// scanned again by `stream_pair`). Returns the number of cuts and jobs.
    #[allow(clippy::too_many_arguments)]
    fn split_item(&self, it: &mut Item, done: usize, nshapes: usize, budget: usize, pair_limit: usize, b: &Batch, work: &mut Vec<Item>, jobs: &mut Vec<PairJob>) -> u64 {
        let proj = |n: usize| n.saturating_mul(nshapes) / done.max(1);
        let total: usize = it.surv.iter().map(|v| v.len()).sum();
        if proj(total) <= budget {
            return 0;
        }
        // cut points of a sequence of counts: each part's projection stays
        // within the budget unless a single element exceeds it
        let cuts_of = |counts: &mut dyn Iterator<Item = usize>| -> Vec<usize> {
            let mut cuts = Vec::new();
            let mut acc = 0usize;
            for (k, x) in counts.enumerate() {
                let x = proj(x);
                if acc > 0 && acc + x > budget {
                    cuts.push(k);
                    acc = 0;
                }
                acc += x;
            }
            cuts
        };
        if it.pieces.len() > 1 {
            let cuts = cuts_of(&mut it.surv.iter().map(|v| v.len()));
            for &a in cuts.iter().rev() {
                let pieces = it.pieces.split_off(a);
                let surv = it.surv.split_off(a);
                work.push(Item { pieces, surv, s0: done });
            }
            return cuts.len() as u64;
        }
        let pc = it.pieces[0];
        let sv = std::mem::take(&mut it.surv[0]);
        if pc.hi - pc.lo > 1 {
            // by subjects (survivor offsets are relative to the piece's first residue)
            let offs = self.index.subject_offsets();
            let starts: Vec<u32> = offs[pc.lo..=pc.hi].iter().map(|&o| (o - offs[pc.lo]) as u32).collect();
            let subj = |t: u32| starts.partition_point(|&x| x <= t) - 1;
            let mut cnt = vec![0usize; pc.hi - pc.lo];
            for x in &sv {
                cnt[subj(x.t)] += 1;
            }
            let cuts = cuts_of(&mut cnt.into_iter());
            if cuts.is_empty() {
                it.surv[0] = sv;
                return 0;
            }
            let mut parts: Vec<Vec<Surv>> = (0..=cuts.len()).map(|_| Vec::new()).collect();
            for x in sv {
                let k = cuts.partition_point(|&c| c <= subj(x.t));
                let shift = if k == 0 { 0 } else { starts[cuts[k - 1]] };
                parts[k].push(Surv { g: x.g, t: x.t - shift });
            }
            for k in (1..=cuts.len()).rev() {
                let lo = pc.lo + cuts[k - 1];
                let hi = if k < cuts.len() { pc.lo + cuts[k] } else { pc.hi };
                work.push(Item { pieces: vec![Piece { lo, hi, ..pc }], surv: vec![parts.pop().expect("part")], s0: done });
            }
            it.pieces[0].hi = pc.lo + cuts[0];
            it.surv[0] = parts.pop().expect("first part");
            return cuts.len() as u64;
        }
        // one subject: by queries (cut at query starts)
        let mut qs: Vec<u32> = sv.iter().map(|x| b.query_of(x.g as usize) as u32).collect();
        qs.sort_unstable();
        let mut runs: Vec<(u32, usize)> = Vec::new();
        for &q in &qs {
            match runs.last_mut() {
                Some((lq, n)) if *lq == q => *n += 1,
                _ => runs.push((q, 1)),
            }
        }
        if runs.iter().any(|r| proj(r.1) > pair_limit) {
            // Parts over [g_lo, g_hi): every oversized query is a part of
            // its own (streamed), the queries between are grouped within the
            // budget as below.
            let mut parts_r: Vec<(u32, u32, Option<(usize, usize)>)> = Vec::new();
            let mut from = pc.g_lo;
            let mut acc = 0usize;
            for &(q, n) in &runs {
                let q = q as usize;
                let gs = (b.base[q] as u32).max(from);
                if proj(n) > pair_limit {
                    if gs > from {
                        parts_r.push((from, gs, None));
                    }
                    let ge = ((b.base[q] + b.lens[q]) as u32).min(pc.g_hi).max(gs);
                    parts_r.push((gs, ge, Some((q, proj(n)))));
                    from = ge;
                    acc = 0;
                } else {
                    let x = proj(n);
                    if acc > 0 && acc + x > budget && gs > from {
                        parts_r.push((from, gs, None));
                        from = gs;
                        acc = 0;
                    }
                    acc += x;
                }
            }
            if from < pc.g_hi {
                parts_r.push((from, pc.g_hi, None));
            }
            let mut bins: Vec<Vec<Surv>> = parts_r.iter().map(|_| Vec::new()).collect();
            for x in sv {
                let k = parts_r.partition_point(|r| r.1 <= x.g);
                if k < parts_r.len() && parts_r[k].2.is_none() {
                    bins[k].push(x);
                }
            }
            let mut keep: Option<(u32, u32, Vec<Surv>)> = None;
            let mut n_cut = 0u64;
            for ((from, to, job), part) in parts_r.into_iter().zip(bins) {
                match job {
                    Some((q, projected)) => {
                        jobs.push(PairJob { subject: pc.lo, q, projected });
                        n_cut += 1;
                    }
                    None if keep.is_none() => keep = Some((from, to, part)),
                    None => {
                        work.push(Item { pieces: vec![Piece { g_lo: from, g_hi: to, ..pc }], surv: vec![part], s0: done });
                        n_cut += 1;
                    }
                }
            }
            let (from, to, part) = keep.unwrap_or((pc.g_lo, pc.g_lo, Vec::new()));
            it.pieces[0].g_lo = from;
            it.pieces[0].g_hi = to;
            it.surv[0] = part;
            return n_cut;
        }
        let cuts = cuts_of(&mut runs.iter().map(|r| r.1));
        if cuts.is_empty() {
            it.surv[0] = sv;
            return 0;
        }
        let gcut: Vec<u32> = cuts.iter().map(|&k| b.base[runs[k].0 as usize] as u32).collect();
        let mut parts: Vec<Vec<Surv>> = (0..=cuts.len()).map(|_| Vec::new()).collect();
        for x in sv {
            parts[gcut.partition_point(|&g| g <= x.g)].push(x);
        }
        for k in (1..=cuts.len()).rev() {
            let g_lo = gcut[k - 1];
            let g_hi = if k < cuts.len() { gcut[k] } else { pc.g_hi };
            work.push(Item { pieces: vec![Piece { g_lo, g_hi, ..pc }], surv: vec![parts.pop().expect("part")], s0: done });
        }
        it.pieces[0].g_hi = gcut[0];
        it.surv[0] = parts.pop().expect("first part");
        cuts.len() as u64
    }


    /// Scan a piece: seed lookups and window filter for the given shapes.
    /// With one shape, the scan stops after the first subject at which the
    /// piece holds more than `limit` survivors and returns the next subject
    /// (the rest of the piece is to be scanned separately).
    #[allow(clippy::too_many_arguments)]
    #[inline(never)]
    fn scan_chunk_shapes(
        &self,
        db: &ScanDb,
        rs: &RoundSeeds,
        b: &Batch,
        seeds: &ScanSeeds,
        pc: &Piece,
        shapes: std::ops::Range<usize>,
        st: &mut ScanState,
        surv: &mut Vec<Surv>,
        limit: usize,
    ) -> Option<usize> {
        #[cfg(target_arch = "x86_64")]
        let (hits, filtered, stop) = if self.avx2 {
            // SAFETY: AVX2 (and popcnt) detected at engine construction.
            unsafe { self.scan_hits_avx2(db, rs, b, seeds, pc, shapes, st, surv, limit) }
        } else {
            self.scan_hits::<false>(db, rs, b, seeds, pc, shapes, st, surv, limit)
        };
        #[cfg(not(target_arch = "x86_64"))]
        let (hits, filtered, stop) = self.scan_hits::<false>(db, rs, b, seeds, pc, shapes, st, surv, limit);
        st.w.counters.postings_visited += hits;
        st.w.counters.hits += hits;
        st.w.counters.hits_filtered += filtered;
        stop
    }

    /// Survivors of subjects `[lo, hi)` (all shapes) grouped by (subject,
    /// query) and visited: coverage check, X-drop, regions.
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    fn visit_chunk(&self, db: &ScanDb, b: &Batch, lo: usize, hi: usize, cap: usize, surv: &mut Vec<Surv>, st: &mut ScanState) {
        if surv.is_empty() {
            return;
        }
        let t0 = tsc();
        self.visit_chunk_inner(db, b, lo, hi, cap, surv, st);
        st.w.counters.c_scan_visit += tsc() - t0;
    }

    #[allow(clippy::too_many_arguments)]
    fn visit_chunk_inner(&self, db: &ScanDb, b: &Batch, lo: usize, hi: usize, cap: usize, surv: &mut Vec<Surv>, st: &mut ScanState) {
        let offs = self.index.subject_offsets();
        // The vectorized X-drop reads the padded residue copy.
        let vec_xdrop = self.avx2 && !db.raw.is_empty();
        std::mem::swap(&mut st.surv, surv);
        let mut t_sort = 0u64;
        // Group survivors by subject, then order each subject's survivors by
        // (query position, shape, subject position): the order in which the
        // index schedules visit the hits of one (query, subject) pair.
        let nsub = hi - lo;
        // Chunk-local subject starts and the subject holding the first
        // residue of every 32-residue block (subjects are rarely shorter).
        let base0 = offs[lo];
        st.starts.clear();
        st.starts.extend(offs[lo..=hi].iter().map(|&o| (o - base0) as u32));
        let clen = st.starts[nsub] as usize;
        st.blk.clear();
        {
            let mut s = 0usize;
            for k in 0..=(clen >> 5) {
                let r = (k << 5) as u32;
                while s + 1 < nsub && st.starts[s + 1] <= r {
                    s += 1;
                }
                st.blk.push(s as u32);
            }
        }
        st.counts.clear();
        st.counts.resize(nsub + 1, 0);
        st.sidx.clear();
        st.sidx.reserve(st.surv.len());
        for x in &st.surv {
            // SAFETY: t < clen, so t >> 5 indexes blk and the walk stops at a
            // subject < nsub (starts[nsub] = clen > t).
            let mut s = unsafe { *st.blk.get_unchecked((x.t >> 5) as usize) } as usize;
            while unsafe { *st.starts.get_unchecked(s + 1) } <= x.t {
                s += 1;
            }
            st.sidx.push(s as u32);
            st.counts[s + 1] += 1;
        }
        for i in 0..nsub {
            st.counts[i + 1] += st.counts[i];
        }
        st.sorted.clear();
        st.sorted.resize(st.surv.len(), Surv::default());
        {
            let mut pos = st.counts.clone();
            for (x, &si) in st.surv.iter().zip(&st.sidx) {
                let p = &mut pos[si as usize];
                st.sorted[*p as usize] = *x;
                *p += 1;
            }
        }
        let sorted = std::mem::take(&mut st.sorted);
        let mut sorted = sorted;
        let mut groups = std::mem::take(&mut st.groups);
        for si in 0..nsub {
            let (a0, a1) = (st.counts[si] as usize, st.counts[si + 1] as usize);
            if a0 == a1 {
                continue;
            }
            let grp = &mut sorted[a0..a1];
            let t0s = st.starts[si];
            let ts0 = tsc();
            // Survivors arrive in (shape, subject position) order, so a stable
            // sort by batch position gives (position, shape, subject position).
            sort_by_position(grp, &mut st.tmp);
            let subject = (lo + si) as u64;
            let subj = self.index.subject_codes(subject as usize);
            // SAFETY: the residue copy holds every subject after PAD bytes.
            let sraw = if vec_xdrop { Some(unsafe { db.raw.as_ptr().add(PAD + offs[subject as usize] as usize) }) } else { None };
            // (subject, query) groups; queries come in increasing order.
            groups.clear();
            let mut i = 0;
            while i < grp.len() {
                let q = b.query_of(grp[i].g as usize);
                let qend = (b.base[q] + b.lens[q]) as u32;
                let mut j = i + 1;
                while j < grp.len() && grp[j].g < qend {
                    j += 1;
                }
                groups.push((i as u32, j as u32, q as u32));
                i = j;
            }
            t_sort += tsc() - ts0;
            #[cfg(target_arch = "x86_64")]
            {
                // The X-drop reads the subject's residues and the query codes
                // around every survivor: request them before the visits.
                let mut o = 0;
                while o < subj.len() {
                    // SAFETY: prefetch is a hint inside the residue buffer.
                    unsafe { std::arch::x86_64::_mm_prefetch(subj.as_ptr().add(o) as *const i8, std::arch::x86_64::_MM_HINT_T0) };
                    o += 64;
                }
            }
            // The query codes around the survivors (the first block of each
            // X-drop direction), requested XDROP_AHEAD survivors ahead.
            let mut pf = 0usize;
            for gi in 0..groups.len() {
                let (ga, gb, q) = groups[gi];
                #[cfg(target_arch = "x86_64")]
                {
                    let target = (gb as usize + XDROP_AHEAD).min(grp.len());
                    while pf < target {
                        let g = grp[pf].g as usize;
                        // SAFETY: prefetch hints inside the padded batch (PAD >= 16 positions each side).
                        unsafe {
                            let p0 = b.qcb.as_ptr().add(g - 16) as *const i8;
                            std::arch::x86_64::_mm_prefetch(p0, std::arch::x86_64::_MM_HINT_T0);
                            std::arch::x86_64::_mm_prefetch(p0.add(62), std::arch::x86_64::_MM_HINT_T0);
                        }
                        pf += 1;
                    }
                }
                let _ = (ga, &mut pf);
                let q = q as usize;
                st.hits.clear();
                for x in &grp[ga as usize..gb as usize] {
                    st.hits.push((x.g - b.base[q] as u32, x.t - t0s));
                }
                self.begin_subject(&mut st.w);
                let hits = std::mem::take(&mut st.hits);
                self.visit_scan(b, q, subj, sraw, &hits, &mut st.w);
                st.hits = hits;
                let list = &mut st.regions[q];
                self.finish_subject_into(subject, &mut st.w, b.region_thr[q], list);
                if list.len() > 2 * cap {
                    list.sort_unstable_by(region_order_c);
                    st.w.counters.regions_dropped_by_cap += (list.len() - cap) as u64;
                    list.truncate(cap);
                }
            }
        }
        st.groups = groups;
        st.sorted = sorted;
        st.w.counters.c_scan_sort += t_sort;
        std::mem::swap(&mut st.surv, surv);
        // release the buffer: one per work unit would otherwise keep its peak size
        *surv = Vec::new();
    }

    /// Window filter (AVX2): one 32-byte load per hit, the query's alignment
    /// codes, reduced to the seed alphabet with two byte shuffles; the windows
    /// of the hits 16 ahead are requested into L1. Not inlined: a small
    /// function keeps its constants in registers.
    #[cfg(target_arch = "x86_64")]
    #[inline(never)]
    #[target_feature(enable = "avx2,popcnt")]
    #[allow(clippy::too_many_arguments)]
    unsafe fn filter_avx2(
        &self,
        db: &ScanDb,
        rs: &RoundSeeds,
        b: &Batch,
        hl: &[u32],
        pend: &[(u32, u32, u32, u32)],
        sb: usize,
        sl: usize,
        c: usize,
        kx: u32,
        min_ids: u32,
        t_base: u32,
        surv: &mut Vec<Surv>,
        filtered: &mut u64,
    ) {
        use std::arch::x86_64::*;
        let tlo = _mm256_broadcastsi128_si256(_mm_loadu_si128(b.rtab_lo.as_ptr() as *const __m128i));
        let thi = _mm256_broadcastsi128_si256(_mm_loadu_si128(b.rtab_hi.as_ptr() as *const __m128i));
        let nib = _mm256_set1_epi8(0x0f);
        let fifteen = _mm256_set1_epi8(15);
        let m7f = _mm256_set1_epi8(0x7F);
        let wp = b.win.as_ptr();
        let hp = hl.as_ptr();
        let nh = hl.len();
        for &(p, h0, h1, _) in pend {
            if h0 == h1 {
                continue;
            }
            let ss = p as i64 + c as i64 - 16;
            let lo_l = (-ss).max(0) as u32;
            let hi_l = (sl as i64 - ss).min(32) as u32;
            let lane = (if hi_l >= 32 { u32::MAX } else { (1u32 << hi_l) - 1 }) & !((1u32 << lo_l) - 1);
            let g0 = (PAD as i64 + sb as i64 + ss) as usize;
            let sv = _mm256_and_si256(_mm256_loadu_si256(rs.buf.as_ptr().add(g0 - rs.base) as *const __m256i), m7f);
            let svr = if kx > 0 { _mm256_loadu_si256(db.raw.as_ptr().add(g0) as *const __m256i) } else { _mm256_setzero_si256() };
            for k in h0 as usize..h1 as usize {
                if k + 16 < nh {
                    // both lines a window may span
                    let wq = wp.add(*hp.add(k + 16) as usize + c - 16);
                    _mm_prefetch(wq as *const i8, _MM_HINT_T0);
                    _mm_prefetch(wq.add(31) as *const i8, _MM_HINT_T0);
                }
                let g = *hp.add(k);
                let w = _mm256_loadu_si256(wp.add(g as usize + c - 16) as *const __m256i);
                // codes 0..23 -> seed-alphabet code; sentinels (bit 7 set) -> 0xFF
                let idx = _mm256_and_si256(w, nib);
                let red = _mm256_blendv_epi8(_mm256_shuffle_epi8(tlo, idx), _mm256_shuffle_epi8(thi, idx), _mm256_cmpgt_epi8(w, fifteen));
                let red = _mm256_or_si256(red, _mm256_cmpgt_epi8(_mm256_setzero_si256(), w));
                let eq = _mm256_movemask_epi8(_mm256_cmpeq_epi8(red, sv)) as u32 & lane;
                let mut v = eq.count_ones();
                if kx > 0 {
                    let ex = _mm256_movemask_epi8(_mm256_cmpeq_epi8(w, svr)) as u32 & eq;
                    v += kx * ex.count_ones();
                }
                if v >= min_ids {
                    surv.push(Surv { g, t: t_base + p });
                } else {
                    *filtered += 1;
                }
            }
        }
    }

    /// One (query, subject) pair scanned and visited on its own: the query's
    /// positions in slices of about `slice` projected survivors, each slice
    /// scanned with every shape (a piece restricted to the slice), its
    /// survivors ordered by (query position, shape, subject position) and
    /// visited, with one visit state for the whole pair — the order and the
    /// decisions of an unsliced visit, so the same regions.
    #[allow(clippy::too_many_arguments)]
    fn stream_pair(&self, db: &ScanDb, rs: &RoundSeeds, b: &Batch, seeds: &ScanSeeds, job: PairJob, cap: usize, slice: usize, st: &mut ScanState) {
        let (s, q) = (job.subject, job.q);
        let g0 = b.base[q] as u32;
        let qlen = b.lens[q];
        if qlen == 0 {
            return;
        }
        let per_pos = (job.projected / qlen).max(1);
        let step = (slice / per_pos).clamp(1, qlen) as u32;
        let offs = self.index.subject_offsets();
        let subj = self.index.subject_codes(s);
        let vec_xdrop = self.avx2 && !db.raw.is_empty();
        // SAFETY: the residue copy holds every subject after PAD bytes.
        let sraw = if vec_xdrop { Some(unsafe { db.raw.as_ptr().add(PAD + offs[s] as usize) }) } else { None };
        let mut t_visit = 0u64;
        self.begin_subject(&mut st.w);
        let mut surv: Vec<Surv> = Vec::new();
        let mut hits: Vec<(u32, u32)> = Vec::new();
        let g1 = g0 + qlen as u32;
        let mut ga = g0;
        while ga < g1 {
            let gb = ga.saturating_add(step).min(g1);
            let pc = Piece { lo: s, hi: s + 1, g_lo: ga, g_hi: gb };
            surv.clear();
            for si in 0..seeds.shapes.len() {
                self.scan_chunk_shapes(db, rs, b, seeds, &pc, si..si + 1, st, &mut surv, usize::MAX);
            }
            // (shape, subject position) order per shape; a stable sort by
            // batch position gives the visit order
            let tv = tsc();
            sort_by_position(&mut surv, &mut st.tmp);
            hits.clear();
            hits.extend(surv.iter().map(|x| (x.g - g0, x.t)));
            self.visit_scan(b, q, subj, sraw, &hits, &mut st.w);
            t_visit += tsc() - tv;
            ga = gb;
        }
        let list = &mut st.regions[q];
        self.finish_subject_into(s as u64, &mut st.w, b.region_thr[q], list);
        if list.len() > 2 * cap {
            list.sort_unstable_by(region_order_c);
            st.w.counters.regions_dropped_by_cap += (list.len() - cap) as u64;
            list.truncate(cap);
        }
        st.w.counters.scan_pairs_streamed += 1;
        st.w.counters.c_scan_visit += t_visit;
    }

    /// Diagonal coverage check and X-drop extension of the (window-filtered)
    /// hits of one (query, subject) pair; same decisions as
    /// [`Self::visit_hits_fast`] without the filter.
    #[inline(never)]
    fn visit_scan(&self, b: &Batch, q: usize, subj: &[u8], sraw: Option<*const u8>, hits: &[(u32, u32)], w: &mut Worker) {
        let qlen = b.lens[q];
        let o = b.base[q];
        let qcb = &b.qcb[o..o + qlen];
        let sl = subj.len();
        let epoch = w.epoch;
        let thr = self.params.hsp_threshold;
        let x = self.params.x_drop_ungapped;
        let (mut avoided, mut ext) = (0u64, 0u64);
        let mut steps = 0u64;
        let de = &mut w.diag_epoch[..];
        let dc = &mut w.diag_cover[..];
        let dl = &mut w.diag_last[..];
        let hsps = &mut w.hsps;
        for &(qp, sp) in hits {
            let d = sp as i64 - qp as i64;
            let di = (d + qlen as i64) as usize;
            if di >= de.len() || sp as usize >= sl || qp as usize >= qlen {
                continue;
            }
            if de[di] != epoch {
                de[di] = epoch;
                dc[di] = 0;
                dl[di] = u32::MAX;
            }
            if qp < dc[di] {
                avoided += 1;
                continue;
            }
            ext += 1;
            let (score, lo, hi) = match sraw {
                #[cfg(target_arch = "x86_64")]
                // SAFETY: AVX2 detected (sraw is only set then); the batch
                // and residue buffers are padded by PAD >= 16 on both sides.
                Some(sptr) => unsafe { xdrop_avx2(&b.m16, b.qcb.as_ptr().add(o), qlen, sptr, sl, qp as usize, sp as usize, x, &mut steps) },
                _ => xdrop_cat(&self.mflat, qcb, subj, qp as usize, sp as usize, x, &mut steps),
            };
            dc[di] = hi as u32;
            dl[di] = qp;
            if score >= thr {
                hsps.push(Hsp { diag: d, q_lo: lo as u32, q_hi: hi as u32, score });
            }
        }
        w.counters.extensions_avoided += avoided;
        w.counters.ungapped_extensions += ext;
        w.counters.xdrop_steps += steps;
        w.counters.scan_groups += 1;
    }

    /// [`Self::finish_subject`] for one query of a batch: the subject's
    /// regions are appended to `out` (the per-query cap is applied later
    /// under the same total order).
    #[inline(never)]
    fn finish_subject_into(&self, subject: u64, w: &mut Worker, region_thr: i32, out: &mut Vec<RegionC>) {
        if w.hsps.is_empty() {
            return;
        }
        w.counters.hsps += w.hsps.len() as u64;
        w.hsps.sort_unstable_by_key(|h| h.diag);
        let join = 2 * self.params.window_slack;
        let mut regs: Vec<Region> = Vec::new();
        let mut evidence: Vec<(u32, i32)> = Vec::new();
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
        if region_thr > self.params.hsp_threshold {
            let before = regs.len();
            let mut k = 0;
            regs.retain(|r| {
                let (n, sum) = evidence[k];
                k += 1;
                r.best >= region_thr || (n >= 2 && sum >= region_thr)
            });
            w.counters.regions_weak_dropped += (before - regs.len()) as u64;
        }
        regs.sort_unstable_by(|a, b| b.best.cmp(&a.best).then(a.d_min.cmp(&b.d_min)));
        regs.truncate(self.params.max_regions_per_subject.max(1));
        out.extend(regs.into_iter().map(RegionC::from));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{ScoreMatrix, ScoringSystem};

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self, m: u64) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 33) % m
        }
    }

    /// Random (code, batch position) seeds in batch order, as build_batch makes them.
    fn random_seeds(r: &mut Lcg, ncodes: u64, nq: usize) -> Vec<(Vec<(u32, u32)>, u64)> {
        let mut g = 32u32;
        (0..nq)
            .map(|_| {
                let len = 1 + r.next(60) as u32;
                let v: Vec<(u32, u32)> = (0..len)
                    .map(|p| {
                        // skewed codes: some frequent, most rare
                        let c = if r.next(4) == 0 { r.next(50) } else { r.next(ncodes) };
                        (c as u32, g + p)
                    })
                    .collect();
                g += len + 32;
                (v, 0)
            })
            .collect()
    }

    #[test]
    fn blocked_and_dense_tables_return_exactly_the_seeds_of_a_code() {
        let mut r = Lcg(5);
        for ncodes in [1000u64, 100_003, 10_000_000] {
            let per_q = random_seeds(&mut r, ncodes, 300);
            let n: usize = per_q.iter().map(|x| x.0.len()).sum();
            let mut want: std::collections::BTreeMap<u32, Vec<u32>> = Default::default();
            for (v, _) in &per_q {
                for &(c, g) in v {
                    want.entry(c).or_default().push(g);
                }
            }
            let tabs = [Table::dense(ncodes, &per_q, n), Table::blocked(ncodes, &per_q, n)];
            for tab in &tabs {
                let probe: Vec<u32> = want.keys().copied().chain((0..200).map(|_| r.next(ncodes) as u32)).chain([tab.inv]).collect();
                let mut cand = Vec::new();
                present_positions(tab, &probe, &mut cand);
                let present: std::collections::BTreeSet<u32> = cand.iter().map(|&(_, c)| c).collect();
                for &c in &probe {
                    assert_eq!(present.contains(&c), want.contains_key(&c), "ncodes {ncodes} dense {} code {c}", tab.dense);
                }
                for (&c, gs) in &want {
                    let mut got: Vec<u32> = if tab.dense {
                        tab.entries[tab.offs[c as usize] as usize..tab.offs[c as usize + 1] as usize].to_vec()
                    } else {
                        let (k, tag) = tab.split(c);
                        let (ea, eb) = (tab.recs[k] as u32 as usize, tab.recs[k + 1] as u32 as usize);
                        let mut out = vec![0u32; eb - ea + 8];
                        // SAFETY: room for eb - ea + 8 values; entries padded.
                        let k1 = unsafe { expand_run(&tab.entries, ea, eb, tag, out.as_mut_ptr()) };
                        #[cfg(target_arch = "x86_64")]
                        if is_x86_feature_detected!("avx2") {
                            let mut o2 = vec![0u32; eb - ea + 8];
                            // SAFETY: AVX2 checked; as above.
                            let k2 = unsafe { expand_run_avx2(tab.entries.as_ptr(), ea, eb, tag, o2.as_mut_ptr()) };
                            assert_eq!(&o2[..k2], &out[..k1]);
                        }
                        out.truncate(k1);
                        out
                    };
                    got.sort_unstable();
                    assert_eq!(&got, gs, "ncodes {ncodes} dense {} code {c}", tab.dense);
                }
            }
        }
    }

    #[test]
    fn survivor_sort_is_stable_by_position() {
        let mut r = Lcg(9);
        for n in [0usize, 1, 5, 32, 33, 500, 5000] {
            let mut v: Vec<Surv> = (0..n)
                .map(|i| Surv { g: r.next(if n > 100 { 1 << 22 } else { 40 }) as u32, t: i as u32 })
                .collect();
            let mut want = v.clone();
            want.sort_by_key(|x| x.g); // stable
            let mut tmp = Vec::new();
            sort_by_position(&mut v, &mut tmp);
            assert!(v.iter().zip(&want).all(|(a, b)| a.g == b.g && a.t == b.t), "n {n}");
        }
    }

    #[test]
    fn round_seed_bytes_apply_the_mask_and_pad() {
        let mut map = [NONSEED; 256];
        for c in 0..20 {
            map[c] = (c % 10) as u8;
        }
        let res: Vec<u8> = (0..100u8).map(|i| i % 24).collect();
        let runs = vec![(5u64, 9u64), (40, 41), (95, 120)];
        let rs = seed_bytes(&map, &runs, &res, 30, 97);
        assert_eq!(rs.base, 30);
        assert_eq!(rs.buf.len(), 67 + 2 * PAD);
        assert!(rs.buf[..PAD].iter().chain(&rs.buf[PAD + 67..]).all(|&x| x == NONSEED));
        for r0 in 30..97usize {
            let x = rs.buf[PAD + r0 - 30];
            let base = map[res[r0] as usize];
            let masked = r0 == 40 || r0 >= 95;
            let want = if base == NONSEED { NONSEED } else if masked { base | MASKED } else { base };
            assert_eq!(x, want, "residue {r0}");
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn vector_xdrop_matches_scalar() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let mx = ScoreMatrix::new(ScoringSystem::Blosum62);
        let mut mflat = [0i32; 768];
        for a in 0..24 {
            for c in 0..32 {
                mflat[a * 32 + c] = mx.matrix[a][if c < 24 { c } else { crate::alphabet::AA_X as usize }];
            }
        }
        let m16 = mflat16(&mflat);
        let mut r = Lcg(11);
        for trial in 0..400 {
            let ql = 1 + r.next(300) as usize;
            let sl = 1 + r.next(400) as usize;
            let pad = 32usize;
            let mut q = vec![crate::alphabet::AA_X as u16; ql + 2 * pad];
            for i in 0..ql {
                let bias = if trial % 3 == 0 { 0 } else { r.next(7) as i8 - 3 };
                q[pad + i] = r.next(24) as u16 | ((bias as u8 as u16) << 8);
            }
            let mut sv = vec![0xFFu8; sl + 2 * pad];
            for j in 0..sl {
                sv[pad + j] = r.next(24) as u8;
            }
            // related stretches: long extensions across several blocks
            if trial % 2 == 0 {
                let off = r.next(sl as u64) as usize;
                for t in 0..ql.min(sl - off) {
                    if r.next(5) > 0 {
                        sv[pad + off + t] = (q[pad + t] & 0xFF) as u8;
                    }
                }
            }
            for _ in 0..50 {
                let qp = r.next(ql as u64) as usize;
                let sp = r.next(sl as u64) as usize;
                for x in [10, 20, 35] {
                    let (mut s1, mut s2) = (0u64, 0u64);
                    let want = xdrop_cat(&mflat, &q[pad..pad + ql], &sv[pad..pad + sl], qp, sp, x, &mut s1);
                    // SAFETY: AVX2 checked; both buffers are padded by 32.
                    let got = unsafe { xdrop_avx2(&m16, q.as_ptr().add(pad), ql, sv.as_ptr().add(pad), sl, qp, sp, x, &mut s2) };
                    assert_eq!(got, want, "trial {trial} qp {qp} sp {sp} x {x} ql {ql} sl {sl}");
                }
            }
        }
    }
}

