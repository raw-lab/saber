# Known limitations and untested paths (SABER 1.4.0)

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
  ungapped extension reaches `--hsp-threshold`. The window filter (since
  1.3) is stricter than 1.2's: a 32-residue window around the seed must
  reach 27 (28 for `very-fast`) = seed-alphabet identities + 2 × exact
  identities (`--hit-filter-ids`, `--hit-filter-exact-weight`). Candidate
  regions are capped per query (`--index-max-candidates`) and seeds that
  occur more often than `--seed-max-occ` in the database are skipped —
  automatically above 64 × the mean occurrences or 256 × max(1, database
  residues / 40 M), whichever is larger (1.2: 4,096), so members of large,
  highly similar families can lose seeds that 1.2 kept.
  When a cap or skip applies, the run summary (`--stats`, stderr) says so with
  counts.
* **The frequency cap below 40 M residues.** Query seeds whose code occurs
  more often in the database than the cap are skipped; the automatic cap is
  64 × the mean occurrences per code, at least 256 × max(1, database
  residues / 40 M). Above 40 M residues the floor grows with the database
  (1.4), so the seeds used no longer depend on its size; below, the fixed
  floor of 256 applies, so the larger the database, the more seeds are
  skipped (the 11,206 SCOP40 queries with `default`: 0.13 M of 27 M seeds
  skipped against 11 M residues, 0.48 M against 20 M, 1.37 M against 38 M).
  Together with the region thresholds, which also still rise with the
  database below 40 M residues, this makes SABER lose true homologues as the
  database grows from 11 to 38 M residues, and its sensitivity falls faster
  than DIAMOND's there (`IMPLEMENTATION_CHECKLIST.md` A.2–A.3).
  `--seed-max-occ` sets a fixed cap. No cap at all (`--seed-max-occ
  1000000000`) took 19 % more time at scale (23 % more CPU) for +0.0035 /
  +0.0040 sensitivity on the test / tune split (third machine, 2 alternating
  runs each; `bench/results/checks_extra/`).
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
  the tracebacks of the alignments it drops. Since 1.4 the floor is that of a
  database of at most 40 M residues (`--ref-db-residues`): on larger
  databases the thresholds stay at that level instead of rising with the
  database, so no candidate region is dropped only because the database is
  large (`cli::region_heuristics_stop_tightening_above_the_reference_size`).
  The price is time: the number of candidate regions band-scored grows in
  proportion to the database beyond that size. The two-stage scoring
  threshold follows the same rule. Below 40 M residues the thresholds still
  rise with the database, as in 1.3.
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
  queries take 60.7 s at a peak of 558 MB (1.4.0 run, machine B;
  `bench/results/lowc/`). In the 1.3.0 run all 10 ran for 30 minutes without
  finishing, at a peak of 685 MB, and before the pair bounds (1.3
  development) the 10 were killed by the out-of-memory killer at 5.8 GB
  after 448 s (`bench/results/v1.3.0/checks/`). On such input memory is
  bounded, time is not: with the `legacy` preset's seeds (contiguous 4-mers,
  full alphabet) all 10 take 27.3 s at 867 MB (1.4.0 run).
* **Few queries against a large database.** The scan schedule streams the
  whole database once per query batch, so a batch of one query costs a full
  pass over the database. Since 1.4, batches of up to 262,144 seeds per shape
  carry a presence sketch that rules out most database positions cheaply:
  one SCOP40 domain against the prebuilt 38.2 M-residue at-scale index takes
  1.2 s (1.3.0: 2.1 s; DIAMOND `--very-sensitive`: 9.5 s), 100 domains 1.5 s
  (machine A, `qscale` section). Twenty long proteins searched one per batch
  (`--query-batch 1`) against the 36 M-residue ×4 database take 24.3 s with
  2 threads on machine B (`threads` section). In the 1.3.0 run, on a machine
  of type A, the same workload took 50.6 s, and 7.4 s in one batch. Keep
  queries in batches (the default). For many separate single-query searches
  of one database, an index built with `--index-postings` and searched with
  `--batch-order query` gives the same alignments in 0.5–0.6 s per query
  (measured with 1.3), but that index takes about 95 bytes per residue
  (3.2 GB for these 36 M residues; 82 s to build on machine B) and the search
  maps most of it (peak RSS 3.2 GB) (`bench/results/summary.md`, "Index
  construction"; `bench/results/v1.3.0/checks/few_queries_large_database.jsonl`,
  `single_query_postings_index.jsonl`).

## Statistics

* **Per-subject E-values** (`--evalue-model per-subject`, the default since
  1.3): the search space is the query–subject pair's (BLAST length adjustment
  for the pair) × the number of database sequences. Hits in long subjects get
  larger E-values than under BLAST's database search space, hits in short
  subjects smaller ones; the ranking of hits across subjects changes
  accordingly (on the benchmarks this ranks true homologues above false
  positives more often: `IMPLEMENTATION_CHECKLIST.md` B.3). E-values are
  therefore not comparable with BLAST's, DIAMOND's or SABER 1.2's;
  `--evalue-model database` restores 1.2's. Calibration against shuffled
  databases: `bench/results/summary.md`, "Null behaviour". Much of the
  ranking gain at scale comes from this model: with database E-values the
  same alignments reach a pooled recall of 0.077 at 0.01 false positives per
  query instead of 0.083 (DIAMOND `--very-sensitive`: 0.080). The gain is
  not only a length effect: against decoys cut to SCOP40 domain lengths
  (`lenmatch` section) `default` still reaches 0.082 (database E-values:
  0.077; DIAMOND: 0.080). At a fixed cut-off the model changes which hits
  pass: on the two-domain benchmark 1.4 reports 1,655 of the stronger
  domains' true positives at E ≤ 1e-3 (1.2: 1,671) and, with database
  E-values, 1,720 (measured with 1.3, whose output is the same).
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

* 2-vCPU virtual machines in one session (`bench/README.md`). The virtual
  machine was replaced twice: machine A (Xeon @ 2.80 GHz, 33 MiB L3) ran the
  at-scale, scaling, length-matched, query-count and SCOP40 sections, machine
  B (Xeon @ 2.10 GHz, 260 MiB L3) the others and a repeat of the headline
  pairs, and a third machine of A's CPU type the one-off check. These are
  shared virtual machines whose speed varies: the third machine ran the same
  one-off searches 30–45 % slower than machine A had. The headline pairs were
  therefore run alternately (5 runs each at scale, 3 at every scaling size)
  and the main configurations of each section 3 times; `summary.md` marks
  the single runs (among them the other tiers of the scaling section, DIAMOND
  `--ultra-sensitive` and SABER 1.3.0 at scale, `exhaustive`, and the
  translated, two-domain, long-protein, decoy, low-complexity and sweep
  runs). Comparisons between tools are made only within one machine.
* **The speed ratio depends on the machine.** SABER `default` against DIAMOND
  `--very-sensitive` with prebuilt indexes took 0.80× of DIAMOND's time at
  scale on machine A and 0.95× on machine B; at 146.8 M residues 0.95× and
  0.99× (identical outputs on both). On other hardware SABER's lead may be
  smaller or absent, and on databases much larger than 147 M residues DIAMOND
  may be faster (an extrapolation puts the crossover near 300 M residues on
  machine A).
* **Without a preparation step** (SABER searching the FASTA, DIAMOND running
  `makedb` + `blastp`) the result at scale changes with the machine: SABER
  was slower on machine B (28.5 vs 27.8 s, all 5 alternating pairs), faster on
  the third machine (39.8 vs 41.9 s, all 5 alternating runs; 1.3.0 39.7 s)
  and in machine A's single run (27.5 vs 32.1 s). At 146.8 M residues it was
  slower (100.0 vs 91.4 s, single runs): SABER's index construction takes
  3–4.4× as long as DIAMOND's `makedb` on machine A at 11–147 M residues
  (7.4× at 147 M residues on machine B).
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
  correction and E-value statistics (per-subject since 1.3), so near the
  E ≤ 1e-3 cut-off it favours SABER over DIAMOND; the recall of the reference
  pairs with E ≤ 1e-10 is reported alongside. 74 of the 500 queries are also
  database members, and self-hits count for every tool.
* **Memory on small databases.** The scan schedule's seed tables have a fixed
  part (2.5 MB per seed shape) and grow with the query batch, so on small
  databases SABER (since 1.3) uses more memory than 1.2 and than DIAMOND:
  - SCOP40 all-against-all: `very-fast` 59 MB vs DIAMOND `--fast` 27 MB,
    `sensitive` 88 MB vs `--very-sensitive` 58 MB, `more-sensitive` 163 MB
    vs `--ultra-sensitive` 114 MB (`default` 137 MB; SABER 1.2.0 `default`
    29 MB in the 1.3.0 run);
  - the MMseqs2 example's 25 % and 50 % subsets: `default` 92 and 97 MB,
    `--very-sensitive` 50 and 64 MB;
  - the 11.0 and 20.1 M-residue scaling databases: `default` 243 and 264 MB,
    `--very-sensitive` 132 and 217 MB. From the 38.2 M-residue database on,
    SABER uses less (280 vs 385 MB; at 146.8 M residues 503 vs 1,400 MB).

  1.4's presence sketches add up to 512 KB per seed shape for batches of at
  most 262,144 seeds (8 MB for the 16-shape presets): the MMseqs2 example
  takes 104 MB, against 96 MB for 1.3.0.
* DIAMOND 2.2.8 only; BLAST, MMseqs2 and SWORD were not run.

## Not implemented

* A sort-merge ("double-indexed") seed join à la DIAMOND. The scan schedule
  looks every database position up in hashed seed tables of a query batch
  instead (one random access per position and shape). Per additional
  database residue the two approaches cost about the same on the benchmark
  (1.05–1.09 CPU seconds per million residues for SABER, 1.06–1.07 for
  DIAMOND between 38 and 147 M residues). SABER's lead is therefore a fixed
  offset, and its time ratio approaches 1 as the database grows. Keeping a
  lead on much larger databases needs a cheaper per-position path; the stage
  times of both tools at 147 M residues are in `IMPLEMENTATION_CHECKLIST.md`
  A.2.
* GPU items of checklist C.9 beyond correctness (allocation reuse, cross-query
  batching, resident database blocks, transfer overlap, minimum batch size).
* Little-endian index files only (a byte-order mismatch is detected and asks for
  a rebuild).
* BAM output: `--outfmt bam` writes SAM text (as in 1.0–1.2, with a note on
  stderr). The SAM header carries no `@SQ` lines; samtools needs the reference
  lengths (`samtools view -t ref.fai`, see `examples/example_sam_output.sh`).
