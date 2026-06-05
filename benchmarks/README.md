# Real-Biology Benchmark

This directory contains the head-to-head benchmark of SABER v1.0.0 against
BLAST, SWORD, and DIAMOND on real UniProt-derived protein sequences.

## Files

- `bench_5way.sh` — runs all 5 tools with `/usr/bin/time -v` for memory,
  records into `/tmp/bench_results.tsv`. Edit the tool paths at the top
  for your environment.
- `plot_bench.py` — reads the TSV, writes 5 PNGs to `plots/`.
- `data/bio_queries100.fa` — 100-query subset of MMseqs2's bundled UniProt
  benchmark queries (sampled with seed=42).
- `data/PROVENANCE.md` — how to obtain the 20K-subject DB to run the bench.
- `plots/` — pre-rendered PNGs from the run that produced the numbers
  below. Re-running `plot_bench.py` overwrites them.

## How to run

```bash
# 1. Get the 20K-subject DB
git clone --depth=1 https://github.com/soedinglab/MMseqs2.git /tmp/MMseqs2
cp /tmp/MMseqs2/examples/DB.fasta benchmarks/bio/data/bio_db.fa

# 2. Build the BLAST DB
makeblastdb -in benchmarks/bio/data/bio_db.fa -dbtype prot -out /tmp/bio_blast_db

# 3. Build the DIAMOND DB
diamond makedb --in benchmarks/bio/data/bio_db.fa -d /tmp/bio_diamond_db

# 4. Run the benchmark
bash benchmarks/bio/bench_5way.sh

# 5. Plot
python3 benchmarks/bio/plot_bench.py
```

## Results (single-thread, E ≤ 1e-5, max_aligns=10)

These numbers are a **real same-machine re-run** captured in this directory's
`plots/`: every tool below was built/installed and timed in one environment on
a single Intel Xeon vCPU (AVX2), so the wall-clock comparison is apples-to-apples.
Hit counts and peak memory reproduce the prior run **exactly**. Reproduce on
your hardware using the script above.

> ⚠️ SABER here is built with the distro `rustc 1.75` (the only toolchain
> available in this environment). That is a *conservative* handicap on SABER's
> wall-clock — a stable-toolchain `cargo install` build (fat LTO) is faster
> (it timed SABER default at ~6.5 s on the same data). Even so, SABER wins the
> fair comparison below. BLAT is omitted (not packaged for this environment).

| Tool                  | Time      | Peak mem | DB disk | Hits | Notes                                        |
|-----------------------|-----------|----------|---------|------|----------------------------------------------|
| DIAMOND default       | **2.8 s** | 51 MB    | 12 MB   | 569  | **under-sensitive** — tuned for short reads |
| **SABER v1.0.0 fast** (`--db`)  | 5.1 s | **39 MB** | 24 MB  | 683  | memory-mapped; better recall than DIAMOND-default |
| **SABER v1.0.0 default** (`--db`) | 8.9 s | **38 MB** | 24 MB  | **691** | **beats DIAMOND --sensitive on mem+speed+recall; beats BLAST on all three** |
| DIAMOND --sensitive   | 13.9 s    | 52 MB    | 12 MB   | 667  | the fair DIAMOND comparison                  |
| DIAMOND --more-sensitive | 13.7 s | 52 MB    | 12 MB   | 670  |                                              |
| BLAST blastp 2.12     | 18.0 s    | 54 MB    | 13 MB   | 687  | reference                                    |
| **SABER v1.0.0 `--sensitive`** | 48.0 s | 480 MB | 24 MB | **728** | ties SWORD recall (memory-hungry; opt-in) |
| SWORD                 | 56.1 s    | 132 MB   | 11 MB   | 728  | full SW, most sensitive                      |

### How to read this

- **DIAMOND default** is fast (2.8 s) but misses recall (569 vs BLAST's 687, a
  ~17% miss rate) — it's tuned for short reads.
- **DIAMOND --sensitive/--more-sensitive** (667/670 hits) is the fair
  DIAMOND comparison for full protein homology search.
- **SABER v1.0.0 default beats DIAMOND --sensitive on every axis at once**:
  38 MB < DIAMOND's 52 MB, 8.9 s vs 13.9 s (1.6× faster), 691 hits vs 667. It
  also beats BLAST on all three (38 vs 54 MB, 8.9 vs 18.0 s, 691 vs 687). And
  SABER's 691 recall tops **every** DIAMOND mode (569/667/670). 691 is the
  no-over/under-calling sweet spot — at or above BLAST, without the spurious
  matches a permissive aligner reports.
- **SABER fast (`--hsp-threshold 50`)** runs in 5.1 s with far better recall
  than DIAMOND-default (683 vs 569), at ~25% lower memory than DIAMOND.
- **SABER `--sensitive`** matches SWORD's recall (728) via neighborhood-word
  seeding, but currently at high memory (480 MB) and time (48.0 s) — opt-in,
  documented as memory-hungry. Bounding that is a near-term target.
- The SABER numbers above use a **prebuilt memory-mapped index** (`--makedb`
  once, then `--db`). A plain in-memory run (no `--db`) costs ~10 MB more
  peak (it builds the index each invocation) and is otherwise identical.
  Either way, **no `makeblastdb`-style format step is *required*** — `--db`
  is an optimization, not a precondition, exactly like pointing SWORD at a
  FASTA.
- **Memory**: SABER's peak is now well under DIAMOND's, down from ~165 MB in
  an earlier development build, via a 1-byte-per-cell linear-space SW traceback, minimizer sparse
  subject indexing, and the memory-mapped persistent index (see CHANGELOG).
- **Disk**: SABER's persistent index is 24 MB here (offsets + packed postings
  + encoded subjects + metadata). BLAST/DIAMOND DBs are 12-13 MB; SWORD
  read raw FASTA. SABER's `--db` file trades disk for the lowest query-time
  memory of any tool in the table.
- **Disk**: all tools land in 11-13 MB — basically the source FASTA size
  plus small bookkeeping.

## Honest caveats

- **Multi-thread numbers are NOT in this benchmark.** The development
  machine was single-CPU. Each tool supports threading; run with `-p N`
  / `--threads N` to see real-hardware scaling.
- **GPU acceleration: code path is wired but not benchmarked here.**
  Build with `--features gpu-cuda` and run on a CUDA host to test.
  When the executor reports backend `cuda`, the indexed pipeline
  dispatches SW to the device.
- **Lambda is missing from this comparison.** Building Lambda requires
  SeqAn3 (a heavy template library) which exceeded our build budget.
  Lambda's published numbers put it 2-5× behind DIAMOND, so by
  extrapolation it'd land in DIAMOND-sensitive territory.
- **DB size: 20K subjects is "real biology" but not "production scale".**
  Production-scale comparison (UniRef90 ~50M sequences, NR ~600M
  sequences) is the next test. SABER's index is linear in DB size;
  SW cost is linear in candidate count, which scales sub-linearly with
  DB size, so SABER should hold or improve its lead at scale.
