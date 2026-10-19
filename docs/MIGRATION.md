# Migrating to SABER 1.3

SABER 1.3 keeps 1.2's engine, output formats and index format (v4) and
changes how protein searches are scheduled and tuned. Command lines keep
working; results change because the presets and the E-value model changed.
Part I lists the changes from 1.2; part II (unchanged) the changes from
1.0/1.1 to 1.2.

## Part I — from 1.2 to 1.3

### I.1 Results change: presets, E-values, scoring band

| | 1.2 | 1.3 |
|---|---|---|
| presets | one or two short spaced seeds (`110111`, `1111011`, …) over Murphy-10, window filter 11–13 seed-alphabet identities | DIAMOND's 16 weight-7 `--very-sensitive` seed shapes over Murphy-10 (`very-fast`: the first 4, `fast`: the first 8, `default`: all 16), `sensitive`: eight weight-6 shapes, `more-sensitive`: all 24; window filter value = seed-alphabet identities + 2 × exact identities, threshold 28 (`very-fast`) or 27. Region threshold 0.7 (`very-fast`, `fast`) or 0.65 × the reporting floor (all five presets). `legacy` and `exhaustive` are unchanged |
| E-values | one search space for the database (BLAST) | `--evalue-model per-subject` (default): search space of the query and the hit's subject (BLAST length adjustment for the pair) × the number of database sequences — a score is less significant in a long subject than in a short one, as in FASTA/SSEARCH and HMMER per-sequence E-values. `--evalue-model database` restores 1.2's E-values. Bit scores, alignments and the ranking within a subject are the same under both |
| scoring window | the subject window of `--band` (64) around the region's diagonals; long pairs in a band | each candidate region is scored inside its diagonals ± 48 (`--score-band`; chains of strong regions within 1,024 diagonals are joined, as 1.2 did for long pairs), first in a ±16 band; only regions reaching 0.75 × a lower bound of their reporting floor are rescored in the full band. `--score-band 0` restores 1.2's windows |
| frequency cap (`--seed-max-occ 0`) | 64 × mean occurrences, at least 4,096 | 64 × mean occurrences, at least **256** |
| stringent `-E` / `--min-bitscore` | the region threshold (a fraction of the reporting floor) grew with the stringency, so a stringent search dropped candidate regions an `-E 10` search reported with E below the cut-off | region threshold and two-stage scoring threshold relative to the floor at E = 10 or `-E`, whichever is more permissive: a stringent search reports exactly the `-E 10` alignments that pass its cut-off |

Every 1.3 setting can be set back individually (`--seed-shape`,
`--hit-filter-ids`, `--hit-filter-exact-weight 0`, `--score-band 0`,
`--seed-max-occ 4096`, `--evalue-model database`), but 1.2's presets
themselves are not kept as named presets.

### I.2 Scheduling: the scan schedule

Protein searches (blastp, blastx, tblastn) with single-hit seeding now use
the **scan schedule** by default (`--batch-order auto` = `scan`): seed tables
are built for a whole batch of queries and the database is streamed past
them once per batch. It reports exactly the alignments of the query-major and
partition-major schedules (`--batch-order query|partition`, tested), which
remain available. Two-hit seeding (`--index-min-hits 2`, `--neighborhood-t`),
database minimizers (`--minimizer-window` > 1), blastn, GPU scoring and
`--exhaustive` use the query-major schedule as before.

* `--makedb` no longer stores seed postings, which the scan schedule does not
  need: the index holds the residues, the seed mask and, per seed shape, the
  list of codes occurring more than 128 times (for the frequency cap). The
  scale benchmark's index shrinks from 183 MB (1.2) to 44 MB. To search an
  index with `--batch-order query|partition`, build it with
  `--index-postings` — with the 16-shape `default` preset that stores 16
  posting lists, about 95 bytes per database residue (3.4 GB instead of
  43 MB for the 36 M-residue ×4 MMseqs database of the index benchmark).
  SABER 1.2 reads 1.3 indexes built with
  `--index-postings` (with their seeds and 1.2's parameters); without them it
  stops with `index is missing section 'off.0.0'`.
* 1.2 indexes are read by 1.3. Their frequent-code lists start at 4,096
  occurrences, below 1.3's automatic cap (256), so 1.3 recounts the seeds of
  the database once per search (a few seconds on 38 M residues); rebuild the
  index to avoid it. The same recount happens when a preset's seed shapes are
  not those the index was built with (e.g. an index built with the `default`
  preset searched with `--preset sensitive`).
* With the scan schedule, a `--db` search seeds with the search's own seed
  alphabet and shapes (the preset's, or `--seed-alphabet` / `--seed-shape`);
  1.2 always used the seeds stored in the index and ignored those flags with
  a note. The query-major and partition-major schedules still use the
  index's seeds (and print the note). The database's seed mask is always the
  one stored in the index.
* Memory follows the query batch: about 64 bytes per query residue for the
  16-shape presets (seed tables) plus about 40 MB of fixed table records.
  Batches hold at most max(512 K, database residues / 4) query residues (up
  to 8 M; `--memory-budget` lowers this), so small databases are searched in
  several passes with small tables.
* Every query batch streams the whole database, so **single-query searches
  of a large index are slower than in 1.2** (a long protein against a
  36 M-residue index: about 2 s instead of 0.1 s; DIAMOND `--very-sensitive`
  11 s). Search many queries in one run, and do not lower `--query-batch`
  without need. If a workload consists of many separate single-query
  searches, build the index with `--index-postings` and search it with
  `--batch-order query` (same alignments, about 0.5 s per query, but about
  95 bytes per residue of index).

### I.3 New and changed options

| Option | |
|---|---|
| `--evalue-model per-subject\|database` | new (default `per-subject`) |
| `--batch-order scan` | new value; `auto` now chooses `scan` for eligible searches |
| `--index-postings` | new: store seed postings (with `--makedb`, or in an in-memory build) |
| `--hit-filter-exact-weight N` | new: weight of exact identities in the window filter value (presets: 2; `legacy`: 0) |
| `--score-band N` | new: scoring band half-width around the region's diagonals (presets: 48; `legacy`: 0) |
| `--band N` | now applies only when `--score-band` is 0 |
| `--hit-filter-ids N` | same meaning, in units of the new filter value (presets: 27/28; 1.2: 11–13 seed-alphabet identities) |
| `--seed-max-occ 0` | automatic cap floor 4,096 → 256 |
| `--stats` | new counters `scan_candidates`, `scan_groups`, `scan_splits`, `scan_pairs_streamed`, `xdrop_steps`, `band_rescored` and per-stage cycle counts `c_scan_*` |

Advanced, hidden options: `--score-prefilter-band N` (first-pass band, presets
16; 0 = one pass), `--score-prefilter-frac F` (0.75), `--scan-batch-residues N`
(query residues per scan batch), `--scan-pair-survivors N` (survivors above
which a (query, subject) pair is scanned on its own, in slices; 4 M).

### I.4 Library API

`SearchParams` gained `hit_filter_exact_weight`, `score_band`,
`score_prefilter_band`, `score_prefilter_frac`; `SearchStatistics::with_model`
and `EvalueModel` select the E-value model (`SearchStatistics::new` keeps the
database model); `Engine::search_scan` runs the scan schedule for a batch of
queries; `Counters` gained the fields listed above. `presets::Preset` gained
`exact_weight`, `score_band` and `prefilter_band`.

---

## Part II — from 1.0 / 1.1 to 1.2

SABER 1.2 replaces the search pipeline of 1.0/1.1 with a single engine shared
by every mode, the CLI and the library. Most command lines keep working
unchanged; this part lists everything that behaved differently in 1.2 and why.

### II.1 Persistent indexes must be rebuilt

| | 1.0.x / 1.1.x | 1.2 |
|---|---|---|
| File format | v2/v3 (`KmerIndex::save`) | **v4** (`search::index`) |
| Ambiguous residues | stored as sentinel 255, later scored as `*` | stored as X; never seeded |
| Seeding scheme | k only | alphabet, seed shapes, minimizer window, masking — all stored |
| Layout | one global posting list | partitions of ≤ 65,535 subjects, compact 32-bit postings per partition (64-bit only where needed); bucket offsets dense, or sparse (bitmap + ranks) when the seed code space is much larger than the partition |
| Integrity | none | header + metadata + per-section CRC32, file-size check, structural checks on open (offsets, ranks, bounds), atomic write (temp file + rename) |

Opening an old index fails with a message that names its format version and
the rebuild command. Rebuild once:

```bash
saber --makedb db.sdx -d db.fasta                     # default preset seeds
saber --makedb db.sdx -d db.fasta --preset fast       # seeds of another preset
saber --mode pn --makedb db_tr.sdx -d genomes.fna     # translated (tblastn) index
saber -q q.faa --db db.sdx --verify-index             # optional full checksum pass
```

In 1.2 the seeding configuration is fixed when the index is built. Search-time
flags that only affect seeding (`--seed-alphabet`, `--seed-shape`,
`--index-k`, `--minimizer-window`) are ignored with `--db` (a note is
printed); thresholds such as `--hit-filter-ids`, `--hsp-threshold`, `-E`,
`-k` still apply. (In 1.3 this holds for the query-major and partition-major
schedules only; the scan schedule seeds with the search's own seeds, see
I.2.)

### II.2 CLI behaviour changes

| Option / output | 1.0/1.1 | 1.2 |
|---|---|---|
| default seeding | contiguous full-alphabet 4-mers | `--preset default`: Murphy-10 reduced alphabet, spaced seed `110111`, window identity filter. `--preset legacy` restores 4-mer seeding (now bounded, masked and composition-corrected) |
| `--sensitive` | neighborhood words T=13 + two-hit | alias for `--preset sensitive`; the old behaviour is `--neighborhood-t 13` |
| `--no-index` | full Smith–Waterman per subject | same meaning (also `--exhaustive`, `--preset exhaustive`), now through the SIMD engine |
| `--prefilter*` | k-mer prefilter for the no-index path | accepted, ignored with a warning |
| `--index-diag-bucket`, `--keep-ascii` | tuning/debug | accepted, no effect |
| `--band N` | slice the subject around the best diagonal | subject window slack around the seed diagonals (0 = whole subject); traceback uses genuine banded DP that widens until it reproduces the window score |
| `--index-max-candidates` | cap applied after materialising every hit | bounded heap of candidate regions; drops are counted and reported |
| `-g/-e` | SABER convention (first gap residue costs `-g`) | unchanged by default; `--gap-convention blast` accepts BLAST/DIAMOND values (BLAST 11/1 = SABER 12/1) |
| nn gap default | 11/1 (no published blastn statistics) | 7/2 (= BLAST 5/2, the blastn default for +2/−3) when `-g/-e` are not given |
| E-values / bit scores | **ungapped** λ, K for gapped alignments (E-values too small by orders of magnitude) | gapped λ, K, H, α, β from NCBI `blast_stat.c` with BLAST length adjustment; translated searches use the nucleotide database length / 3 |
| unsupported matrix/gap pairs | silently used unrelated parameters | rejected; `--allow-uncalibrated-evalues` runs them with E-values flagged UNCALIBRATED (m0 banner, `evalue_calibrated=false` in JSON/stats) |
| composition | none | Hauser et al. 2016 local composition-bias correction (as in MMseqs2/DIAMOND); `--comp-bias off` gives plain substitution scores |
| masking | none | tantan for **seeding only** (query and database); alignments always use the original residues; `--masking none` disables |
| `--min-bitscore` | documented as overriding `-E`, implemented as both | both filters apply (documented) |
| `-k` | alignments per query | subjects per query; `--max-hsps` (default 1) alignments per subject — identical to before with the default |
| m8/m9 column 6 | gap **columns** | gap **opens** (BLAST semantics); JSON has both `gaps` and `gap_opens` |
| CIGAR | `I`/`D` swapped | SAM semantics: `I` = query residue against a gap, `D` = subject residue against a gap |
| coverage | alignment columns / length (could exceed 100 % with gaps) | residues consumed on each sequence / its length |
| traceback window | capped at 2,048 residues (longer alignments silently truncated) | no cap; adaptive widening, with an explicit memory cap that is counted (`band_capped`) and reported as `exact=false` |
| nx / pn coordinates | protein (frame) coordinates | nucleotide coordinates, BLAST style (minus strand: start > end); IDs keep the `|frame±N` suffix; JSON adds `query_frame`/`subject_frame` |
| nx / pn engine | brute force per frame (no index) | the indexed engine; pn uses a translated index (`--mode pn --makedb`) |
| output order | per query in input order | unchanged; ties broken deterministically (score, subject index, coordinates), identical for any thread count or `--batch-order` |
| `--mismatch -3` | rejected by the parser (only `--mismatch=-3` worked) | accepted (negative numbers as separate arguments) |
| `--genetic-code archaeal` | NCBI table 15 (Blepharisma nuclear: TAG = Q) | NCBI table 11, which archaea use; table 15 is `blepharisma` or `15` |
| `--outfmt bam` | wrote SAM text without saying so | still writes SAM text, now with a note on stderr (convert with `samtools view -b`) |
| SAM `NM` tag | mismatches only | edit distance: mismatches + inserted and deleted residues |
| `--quiet`, diagnostics | `--quiet` still printed the timer lines (on stderr) | `--quiet` prints errors only; otherwise a one-line run summary goes to stderr (results only on stdout, as before); `--stats FILE` writes a JSON run summary (`hits_filtered` counts every seed hit rejected by the window filter; `extensions_avoided` counts hits that passed it but lay on an already extended stretch of their diagonal) |

New options: `--preset`, `--gap-convention`, `--allow-uncalibrated-evalues`,
`--masking`, `--comp-bias`, `--seed-alphabet`, `--seed-shape`,
`--seed-max-occ`, `--hit-filter-ids`, `--max-hsps`, `--query-cover`,
`--subject-cover`, `--memory-budget`, `--tmpdir`, `--query-batch`,
`--batch-order`, `--stats`, `--exhaustive`, `--verify-index`.

### II.3 Library API changes

* `algorithm::index::{KmerIndex, SeedHit, Candidate}` were removed: they
  materialised every seed hit before capping candidates. Use
  `search::api::Searcher` (same engine as the CLI) or
  `search::index::SeedIndex` + `search::engine::Engine` directly.
  `algorithm::{DEFAULT_K_PROTEIN, DEFAULT_K_NUCLEOTIDE}` remain.
* `db::DatabaseSearch` keeps its signature but now runs the engine in
  exhaustive mode with gapped statistics; the `EValue` argument is retained for
  source compatibility and no longer used.
* `gpu::GpuExecutor::score_batch` takes **ASCII** sequences for every
  implementation (1.1's CUDA protein path expected pre-encoded bytes while the
  CPU fallback expected ASCII). The new `gpu::ScoreBackend` trait takes
  canonical codes and must return exact `i32` scores.
* `AlignmentResult` gained `gap_opens`, `positives`, `query_length`,
  `subject_length`, `query_frame`, `subject_frame`, `evalue_calibrated` and
  `exact` (all `#[serde(default)]`).
* `saber_rs::VERSION` now comes from `Cargo.toml`.
* `TranslationTable::Archaeal` now translates with NCBI table 11 (as
  `Bacterial`); the new `TranslationTable::Blepharisma` is NCBI table 15, which
  `Archaeal` used before.

## Toolchain

Rust 1.75 or newer (`rust-version = "1.75"`); all dependency pins of 1.1 are
kept. 1.3 was built and tested with rustc 1.75.0 (`cargo test --release`:
128 unit, 23 integration and 1 doc test); `cargo clippy` of rustc 1.97 reports
style lints only.
