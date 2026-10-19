> **Historical (SABER ≤ 1.1.0).** These are the earlier releases' own benchmark
> notes and numbers, kept for reference; they were not re-measured for 1.2.0
> or 1.3.0. The current benchmark suite, protocol and results are in `bench/`
> and `bench/results/summary.md` (1.2.0's in `bench/results/v1.2.0/`).

# SABER performance benchmark record

Hardware: 1 core, Intel Xeon @ 2.80 GHz, AVX2, 3 GB RAM. rustc 1.75, `--release`.
Data: synthetic protein sets, UniProt-shaped log-normal lengths (mean 314 aa),
Swiss-Prot background composition, planted homolog families at 25-90% identity.
Generator: `benchmarks/gen_bench.py`. Timings are best-of-N via wait4 rusage.

## SABER before vs after

| Dataset | Metric | Baseline v1.0.2 | Optimized | Change |
|---|---|---|---|---|
| 20K seqs / 6.3M res | wall | 5.72 s | 2.10 s | **2.7x faster** |
| | peak RSS | 32.9 MB | 35.8 MB | +2.9 MB (fixed cost) |
| | alignments | 1101 | 1101 | identical |
| 60K seqs / 18.8M res | wall | 8.07 s | 3.49 s | **2.3x faster** |
| | peak RSS | 86.0 MB | 87.4 MB | +1.4 MB |
| | alignments | 1034 | 1036 | +2 |

The RSS delta is fixed rayon/thread-pool overhead, not data-proportional:
+9% at 33 MB, +1.6% at 87 MB.

Stage-level, windowed index build: 421 ms -> 283 ms (1.5x).

## vs DIAMOND 2.1.9 (same core, same data)

### 20 000 sequences
| Tool | wall | peak RSS | hits | planted homologs found |
|---|---|---|---|---|
| **SABER optimized** | **2.10 s** | **35.9 MB** | **1101** | **827** |
| DIAMOND default | 1.77 s | 43.3 MB | 832 | 824 |
| DIAMOND --sensitive | 9.55 s | 43.5 MB | 842 | 825 |
| DIAMOND makedb (prep) | +2.93 s | 29.8 MB | - | - |

### 60 000 sequences
| Tool | wall | peak RSS | hits |
|---|---|---|---|
| **SABER optimized** | **3.49 s** | **87.4 MB** | **1036** |
| DIAMOND default | 5.01 s | 95.5 MB | 832 |
| DIAMOND --sensitive | 26.54 s | 96.1 MB | 859 |
| DIAMOND makedb (prep) | +1.19 s | 44.3 MB | - |

vs DIAMOND --sensitive: 4.5x faster at 20K, 7.6x at 60K, lower memory, more hits.
vs DIAMOND default: ahead on memory and hit count at both sizes, ahead on time
at 60K. DIAMOND search alone is faster at 20K, but SABER needs no makedb step.
DIAMOND search time grew 2.8x from 20K->60K; SABER's grew 1.7x.

## Kernel vs two Rust alignment libraries

1500 real pipeline pairs, 2.26e8 DP cells. Gap penalties translated so all
three price a length-k gap identically at -(10+k): rust-bio and block-aligner
charge `open + k*extend`, SABER and BLAST charge `open + (k-1)*extend`.

| Engine | sec | GCUPS | exact vs rust-bio full DP |
|---|---|---|---|
| **SABER (SIMD AVX2)** | 0.074 | 3.06 | **1500/1500 (100.0%)** |
| block-aligner 0.5.0 (SIMD) | 0.071 | 3.19 | 1191/1500 (79.4%) |
| rust-bio 1.5.0 (scalar SW) | 5.333 | 0.04 | ground truth |

SABER is 72x faster than rust-bio's scalar Smith-Waterman and within 4% of
block-aligner, but block-aligner under-reports 309 of 1500 scores because its
adaptive band drops alignments that wander outside it. SABER is the only SIMD
engine here that reproduces exact Smith-Waterman.

## Equivalence checks

Baseline vs optimized, medium dataset, alignments lost / gained:

| Configuration | baseline | optimized | lost | gained |
|---|---|---|---|---|
| default | 1101 | 1101 | 0 | 0 |
| `--minimizer-window 4` | 1101 | 1101 | 0 | 0 |
| `--minimizer-window 8` | 1011 | 1011 | 0 | 0 |
| `--index-k 5` | 897 | 897 | 0 | 0 |
| `--max-aligns 50 -E 0.001` | 817 | 817 | 0 | 0 |
| `--mode nn` (nucleotide) | 3 | 3 | 0 | 0 |
| `--neighborhood-t 11 --index-min-hits 2` | 1941 | 1966 | 167 | 192 |

The neighborhood configuration is the one behavioural difference. It comes
from the unified alignment window producing more accurate (higher) scores,
which reshuffles marginal hits at the bottom of the per-query top-10. The
gained hits outscore the lost ones (30.1 vs 28.7 mean bits) and ground-truth
recall improves (872 vs 870 planted homologs), so it is a marginal-ranking
improvement rather than a regression.

## Parallel-path correctness

Speedup from the parallel query loop could not be measured on a 1-core
sandbox, but correctness was validated under thread oversubscription. Every
thread count produces byte-identical, run-to-run stable output:

| threads | vs threads=1 | run-to-run stable | alignments |
|---|---|---|---|
| 2 | identical | yes | 1101 |
| 4 | identical | yes | 1101 |
| 8 | identical | yes | 1101 |
| 16 | identical | yes | 1101 |

This confirms the ordered chunked write preserves input order, the
thread-local traceback scratch does not leak state between queries, and
parallelism introduces no nondeterminism. It does NOT confirm speedup.
Run `benchmarks/validate_scaling.sh <query.fa> <db.fa>` on real multi-core
hardware for that number.

## Caveats

- Multi-core scaling of the parallel query loop is UNVALIDATED (1-core sandbox).
- Index-build memory fixes matter at k>=5 and high thread counts; neither could
  be exercised here.
- Data is synthetic. Planted families give clean ground truth but no
  low-complexity regions or compositional bias, so candidate-set sizes are
  optimistic versus real UniProt and the DIAMOND hit-count gap may narrow.
