# Known limitations and untested paths (SABER 1.3.0)

## Not validated in this release

* **GPU (`--features gpu-cuda`, `--gpu`)** — no CUDA device was available.
  What *was* verified (1.2.0): the CUDA kernel source, compiled as C++ by
  `tools/gpu_emulate/emulate.cpp`, reproduces the independent reference scores
  for 437 windows (protein, nucleotide, PAM30, bias, scores above 16 bits,
  out-of-range codes); the Rust host code type-checks against cudarc 0.9.15
  when the CUDA headers (`cuda.h`) are installed — without a CUDA toolkit,
  `cargo check --features gpu-cuda` fails.
  Not verified: NVRTC compilation of the kernel, launches, transfers, speed.
  The backend scores one query's candidate windows per launch with fresh device
  allocations; it does not batch across queries, keep database blocks resident
  or overlap transfers with compute, so it is not expected to beat the CPU path
  as shipped. With `--gpu` the query-major schedule is used (the scan schedule
  does not drive the GPU backend). `saber --doctor --gpu` compares GPU and CPU
  scores on a device.
* **More than two threads** — the development machine has 2 vCPUs. The scan
  schedule parallelises over database chunks within each round (one barrier
  per seed shape and round) and the scoring stage over queries; scaling beyond
  2 threads has not been measured.
* **Databases larger than RAM** — the external-memory `--makedb`, the
  memory-budget spill of in-memory builds and the partition-major schedule are
  implemented and tested for identical results, but no database larger than
  the 8 GB of RAM of the test machine was searched. The scan schedule streams
  the database residues from the memory-mapped index and holds the seed tables
  of one query batch.
* **Non-x86 / no-AVX2 hosts** — scalar fallbacks produce identical results
  (tested on x86 by forcing the scalar paths, including the scan schedule's
  table lookups, window filter and X-drop: `engine_consistency::
  scan_schedule_scalar_and_avx2_paths_match_the_index_schedule`), but their
  speed was not benchmarked. AVX-512 is not used (its intrinsics are not
  stable in Rust 1.75).
* **Huge pages** — the seed tables and query arrays request transparent huge
  pages (`MADV_HUGEPAGE`, whole 2 MiB blocks only). The benchmarks ran with
  THP in `madvise` mode; with THP disabled (`never`, or `SABER_NO_HUGEPAGES=1`)
  the table lookups take more TLB misses (not measured).
* **nn (blastn) performance** — correctness is tested; speed and sensitivity
  were not benchmarked against blastn. Nucleotide seeds are contiguous 11-mers
  (no spaced seeds, no DUST masking), searched with the query-major schedule.
* **nn searches the query's plus strand only** (as in 1.0–1.2; blastn searches
  both strands by default). To find minus-strand matches, add the reverse
  complement of each query as a second query.

## Heuristics and their limits

* Every preset except `exhaustive` is a seed-and-extend heuristic: an
  alignment is examined only if a seed hit passes the window filter and its
  ungapped extension reaches `--hsp-threshold`. The 1.3 window filter is
  stricter than 1.2's: a 32-residue window around the seed must reach 27 (28
  for `very-fast`) = seed-alphabet identities + 2 × exact identities
  (`--hit-filter-ids`, `--hit-filter-exact-weight`). Candidate regions are
  capped per query (`--index-max-candidates`) and seeds that occur more often
  than `--seed-max-occ` in the database are skipped — automatically above
  64 × the mean occurrences or 256, whichever is larger (1.2: 4,096), so
  members of large, highly similar families can lose seeds that 1.2 kept.
  When a cap or skip applies, the run summary (`--stats`, stderr) says so with
  counts.
* **Scoring band.** Each candidate region is scored inside its diagonals ± 48
  (`--score-band`); a subject's strong regions whose diagonal ranges lie
  within 1,024 of each other are scored as one band, so an indel of up to
  ~1,000 residues between two seeded parts stays inside it
  (`cli::long_alignment_is_not_truncated`: 6,000 residues across a 260-residue
  net indel). An optimal alignment that leaves the band is scored inside it
  (the score is exact for the band, not for the whole subject); a part that
  deviates further from every seeded diagonal yields a separate or shorter
  alignment. `--score-band 0` scores the whole subject window as 1.2 did.
* **Two-stage scoring.** A region is first scored inside ± 16 diagonals; only
  if that score reaches 0.75 × a lower bound of the query's reporting floor at
  E = 10 (or at `-E`, if more permissive) is it rescored in the ± 48 band. An
  alignment whose score comes mostly from outside the narrow band can be
  missed (`--score-prefilter-band 0` disables the first stage).
* **Region threshold.** All five heuristic presets drop candidate regions
  none of whose ungapped segments — nor the sum of two or more — reaches
  0.65 (`very-fast`, `fast`: 0.7) × the query's reporting score floor at
  E = 10 (or at `-E`, if more permissive; 1.2: at `-E`, for `default`,
  `fast` and `very-fast`, at 0.45–0.5). Both thresholds were tuned at
  `-E 10`; being relative to that floor, a stringent `-E` or `--min-bitscore`
  only filters the output of an `-E 10` search
  (`cli::stringent_evalue_only_filters_an_evalue_10_search`); it saves only
  the tracebacks of the alignments it drops.
* tantan masks periodic sequence for seeding, which can hide homologues whose
  shared residues are periodic: in a synthetic test, homologues with an A at
  every third position (69 % identity over up to 344 residues; `tests/cli.rs`
  `protein_case`) are masked as period-3 repeats and found by no heuristic
  preset, but by `--masking none` or `--preset exhaustive`. Alignments always
  use the unmasked sequence; only seeding is affected.
* Banded traceback starts on the candidate region's seed diagonals and widens
  until it reproduces the exact band score. If that band needs more than
  32 MB of direction matrix (`--traceback-max-mb`), the traceback restarts from
  a band around the optimum's end cell (found by a score-only pass); only if
  that band also exceeds the cap is the best alignment found so far reported,
  with `exact=false` (JSON) and counted as `band_capped`. m8 output does not
  show the `exact` flag.
* **Scan-schedule memory.** The seed tables of a query batch take about
  4 bytes per query residue and seed shape plus 2.5 MB per weight-7 shape over
  Murphy-10 (≈ 64 bytes per residue + 40 MB for the 16-shape presets).
  Window-filter survivors are held until their database chunk has seen every
  shape. Work is cut by chunks, subjects and queries so that about 2 M
  survivors (16 MB) are held at a time; the scan of a shape stops at a
  subject boundary once a piece holds 4 M; and a (query, subject) pair
  projected above 4 M survivors — repeat-rich sequences with masking and the
  frequency cap both disabled — is scanned and visited on its own, in slices
  of its query positions (`scan_pairs_streamed`; same results). What is not
  bounded is the survivors of one subject for one seed shape (its length ×
  the batch positions sharing a seed with it). Adversarial test (long
  low-complexity queries against 2,000 low-complexity sequences, no masking,
  no frequency cap, 16 seed shapes): the first 4 of the benchmark's 10
  queries take 61.9 s at a peak of 564 MB; all 10 ran for 30 minutes without
  finishing, at a peak of 685 MB (`bench/results/checks/`). Before the pair
  bounds, the 10 were killed by the out-of-memory killer at 6.0 GB after
  448 s. On such input memory is bounded, time is not: with the `legacy`
  preset's seeds (contiguous 4-mers, full alphabet) all 10 take 27.4 s at
  868 MB.
* **Few queries against a large database.** The scan schedule streams the
  whole database once per query batch, so a batch of one query costs about
  as much as a full pass over the database: a long protein searched alone
  against a prebuilt index of 36 M residues takes 1.9–2.4 s (DIAMOND
  `--very-sensitive` 11–12 s; SABER 1.2.0, whose index held seed postings,
  0.1 s), and 20 such proteins searched one per batch (`--query-batch 1`)
  take 50.6 s instead of 7.4 s in one batch (1.2.0: 9.8 s). Keep queries in
  batches (the default). For many separate single-query searches of one
  database, an index built with `--index-postings` and searched with
  `--batch-order query` gives the same alignments in 0.5–0.6 s per query,
  but that index takes about 95 bytes per residue (3.4 GB for these 36 M
  residues; 90 s to build) and the search maps most of it (peak RSS 3.3 GB)
  (`bench/results/checks/few_queries_large_database.jsonl`,
  `single_query_postings_index.jsonl`).

## Statistics

* **Per-subject E-values** (`--evalue-model per-subject`, the 1.3 default):
  the search space is the query–subject pair's (BLAST length adjustment for
  the pair) × the number of database sequences. Hits in long subjects get
  larger E-values than under BLAST's database search space, hits in short
  subjects smaller ones; the ranking of hits across subjects changes
  accordingly (on the benchmarks this ranks true homologues above false
  positives more often: `IMPLEMENTATION_CHECKLIST.md` A.3). E-values are
  therefore not comparable with BLAST's, DIAMOND's or SABER 1.2's;
  `--evalue-model database` restores 1.2's. Calibration against shuffled
  databases: `bench/results/summary.md`, "Null behaviour". Much of 1.3's
  ranking gain at scale comes from this model: with database E-values the
  same alignments reach a pooled recall of 0.077 at 0.01 false positives per
  query instead of 0.083 (DIAMOND `--very-sensitive`: 0.080). At a fixed
  cut-off it changes which hits pass: on the two-domain benchmark
  1.3 reports 1,655 of the stronger domains' true positives at E ≤ 1e-3
  (1.2: 1,671) and, with database E-values, 1,720.
* Gapped Karlin–Altschul parameters exist only for the published
  matrix/gap combinations (NCBI `blast_stat.c`: BLOSUM45/50/62/80/90,
  PAM30/70/250 and blastn reward/penalty tables). Other combinations are
  rejected unless `--allow-uncalibrated-evalues`, in which case E-values are
  approximate and flagged.
* Composition correction is the per-position bias correction of Hauser et al.
  (2016), as used by MMseqs2 and DIAMOND by default — not BLAST's
  composition-based matrix adjustment.
* For translated searches the database length is the nucleotide length / 3
  (BLAST convention), queries of nx use the frame length.

## Benchmark scope

* One 2-vCPU virtual machine (Xeon @ 2.80 GHz) in one session. It is a shared
  host: the CPU time of identical runs varied by up to ~20 % over the day, so
  the headline pair was run alternately (5 runs each) and every other main
  configuration 3 times; translated, two-domain, long-protein, scaling,
  decoy and sweep runs are single runs.
* Without a preparation step (SABER searching the FASTA, DIAMOND running
  `makedb` + `blastp`) the margin at scale is small: 32.7 vs 33.8 s at the
  median, SABER faster in each of the 5 alternating pairs (0.91–0.97×) —
  less than the machine's drift over a day, so it may not hold on other
  hardware. With prebuilt indexes the margin is 16 %.
* The at-scale benchmark's negatives are reversed proteins: they keep length
  and composition but not the local structure of real non-homologues, and they
  come in groups of four similar sequences (the ×4 database holds each
  original three more times at 15 % substitution). SCOP40's labels contain
  some cross-fold pairs with genuine sequence similarity, which every tool
  ranks as top "false positives".
* The presets were chosen on the tune splits of the two SCOP40-based
  benchmarks (and the `sensitive` seeds on tune-split homologous pairs); the
  test-split results are held out from that choice, but both splits come from
  the same SCOP40 collection, so they say nothing about other kinds of
  proteins. The `default` preset uses DIAMOND's very-sensitive seed shapes.
* The MMseqs2 example set has no labels; recall is measured against SABER's
  own exhaustive Smith–Waterman search, a reference for the alignment problem,
  not for biological homology. It shares SABER's scoring, composition
  correction and E-value statistics (per-subject in 1.3), so near the
  E ≤ 1e-3 cut-off it favours SABER over DIAMOND; the recall of the reference
  pairs with E ≤ 1e-10 is reported alongside. 74 of the 500 queries are also
  database members, and self-hits count for every tool.
* **Memory on small databases.** The scan schedule's seed tables have a fixed
  part (2.5 MB per seed shape) and grow with the query batch, so on small
  databases 1.3 uses more memory than 1.2 and than DIAMOND: SCOP40
  all-against-all `very-fast` 58 MB vs DIAMOND `--fast` 29 MB, `sensitive` 83 MB
  vs `--very-sensitive` 59 MB, `more-sensitive` 161 MB vs `--ultra-sensitive`
  117 MB (`default` 137 MB; SABER 1.2.0 `default` 29 MB); on the MMseqs2
  example's 25 % and 50 % subsets `default` takes 84 and 89 MB,
  `--very-sensitive` 44 and 73 MB.
* DIAMOND 2.2.8 only; BLAST, MMseqs2 and SWORD were not run.

## Not implemented

* A sort-merge ("double-indexed") seed join à la DIAMOND. The scan schedule
  looks every database position up in hashed seed tables of a query batch
  instead (one random access per position and shape).
* GPU items of checklist B.9 beyond correctness (allocation reuse, cross-query
  batching, resident database blocks, transfer overlap, minimum batch size).
* Little-endian index files only (a byte-order mismatch is detected and asks for
  a rebuild).
* BAM output: `--outfmt bam` writes SAM text (as in 1.0–1.2, with a note on
  stderr). The SAM header carries no `@SQ` lines; samtools needs the reference
  lengths (`samtools view -t ref.fai`, see `examples/example_sam_output.sh`).
