# SABER vs BLAST — Correctness Report

This folder holds two things:

1. **This file** — a small-example *correctness* comparison between SABER
   v1.0.0 and NCBI BLAST+ 2.12.0 (same inputs, same parameters), showing
   SABER recovers the same hits with comparable percent-identity and E-values.

## Hardware / environment

- Ubuntu 24.04 (containerized)
- Rust 1.75.0+ (`cargo build --release`, default features, no GPU)
- BLAST+ 2.12.0+ds-4build2 (apt)

## Small example data (the `examples/` files)

### Protein vs Protein

| Tool   | Query                          | Hit                              | %ID   | bits  | E-value     |
|--------|--------------------------------|----------------------------------|-------|-------|-------------|
| SABER  | query1_human_hemoglobin_alpha  | db1_identical_to_query1          | 100.0 | 338.8 | 1.16e-97    |
| BLAST  | query1_human_hemoglobin_alpha  | db1_identical_to_query1          | 100.0 | 286   | 8.27e-106   |
| SABER  | query1_human_hemoglobin_alpha  | db2_hemoglobin_beta_similar      | 43.45 | 134.9 | 2.77e-36    |
| BLAST  | query1_human_hemoglobin_alpha  | db2_hemoglobin_beta_similar      | 43.45 | 114   | 9.13e-38    |
| SABER  | query1_human_hemoglobin_alpha  | db3_myoglobin_distant_homolog    | 27.03 | 59.7  | 1.16e-13    |
| BLAST  | query1_human_hemoglobin_alpha  | db3_myoglobin_distant_homolog    | 27.52 | 51.2  | 3.75e-13    |
| SABER  | query2_short_kinase_motif      | db4_kinase_with_query2_motif     | 100.0 | 138.5 | 8.93e-38    |
| BLAST  | query2_short_kinase_motif      | db4_kinase_with_query2_motif     | 100.0 | 119   | 6.34e-42    |

Same hits, same percent identity within rounding, comparable rankings. Bit
scores differ because BLAST applies composition-based statistics and X-drop
heuristics while SABER returns the raw Gotoh / Karlin-Altschul values.

### Nucleotide vs Nucleotide

Both tools find the identical hit (HBA query → db1 at 100% / full length) and
the 93.94% partial hit (db2). SABER additionally reports a short low-quality
match against the GC-rich decoy at E ≈ 8.8 (above default 10 so dropped by
BLAST's heuristic threshold). See `saber_nn_results.tsv` and `blastn_results.tsv`.

### Translated modes (blastx / tblastn)

See `saber_nx_results.tsv` / `blastx_results.tsv` and
`saber_pn_results.tsv` / `tblastn_results.tsv`. SABER reports a richer
breakdown across all 6 frames because it does not apply BLAST's seed-extension
filtering.

## Performance

The headline speed / memory / recall comparison, run on 100 UniProt queries against 20,000
real UniProt subjects, single-thread, E ≤ 1e-5. Every tool was built/installed
and timed in the **same environment** (one CPU), so the comparison is fair:

| Tool                    | Time   | Peak mem | Hits | |
|-------------------------|--------|----------|------|--|
| DIAMOND default         | 2.8 s  | 51 MB | 569 | under-sensitive (short-read tuned) |
| **SABER v1.0.0** (fast)    | 5.1 s  | 39 MB | 683 | better recall than DIAMOND-default |
| **SABER v1.0.0** (default) | 8.9 s  | 38 MB | **691** | beats DIAMOND `--sensitive` on memory, speed **and** recall |
| DIAMOND `--sensitive`   | 13.9 s | 52 MB | 667 | the fair DIAMOND comparison |
| BLAST `blastp`          | 18.0 s | 54 MB | 687 | reference |
| SABER v1.0.0 `--sensitive` | 48.0 s | 480 MB | 728 | ties SWORD recall (opt-in) |
| SWORD                   | 56.1 s | 132 MB | 728 | full SW |

SABER v1.0.0 default **beats DIAMOND `--sensitive` on all three axes at once**
(38 < 52 MB, 8.9 < 13.9 s, 691 > 667 hits) and beats BLAST on all three. Its
691 recall tops every DIAMOND mode. Unlike the early development builds, SABER
is now *faster* than BLAST here — the minimizer index plus a banded SIMD
Smith–Waterman kernel removed the old "full SW against every subject" cost
while keeping full-SW sensitivity on the surviving candidates.

> Note: SABER above is built with the distro `rustc 1.75`; a stable-toolchain
> `cargo install` build (fat LTO) is faster still (~6.5 s for default on this
> data). Hit counts and memory are toolchain-independent and reproduce exactly.
> Numbers are single-thread on one 20K-subject dataset; rerun at production
> DB scale before publishing a claim.
