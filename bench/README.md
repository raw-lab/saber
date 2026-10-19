# SABER benchmarks

Everything needed to reproduce the SABER 1.3.0 / 1.2.0 / DIAMOND comparison
in `bench/results/` (summary: `results/summary.md`, machine-readable:
`results/summary.json`). The 1.2.0 release's suite and results — SABER
1.0.2 / 1.1.0 / 1.2.0 / DIAMOND, three machines — are kept unchanged in
`results/v1.2.0/` (with its README). Nothing here is synthetic except the
derived decoy and stress sets, which are generated deterministically from the
public data.

## Data

| File | Source | SHA-256 |
|---|---|---|
| `scop40.fa`, `scop40.labels.tsv` | SCOPe 2.08 40 % domains, converted from `test_data/scop40.cal.gz` of https://github.com/rcedgar/reseek (commit `9cb33e1`, `scop40.cal.gz` sha256 `69808731…`) with `scripts/cal_to_fasta.py` | `93268cef…` / `bbbcf71d…` |
| `mmseqs_QUERY.fasta`, `mmseqs_DB.fasta` | `examples/QUERY.fasta`, `examples/DB.fasta` of https://github.com/soedinglab/MMseqs2 (commit `76de5d8`) — 500 queries, 20,000 UniProt sequences, 9.06 M residues | `c99bc94a…` / `55d48bb7…` |

Full SHA-256 checksums of the inputs, the derived sets and the binaries are
recorded in `results/env/environment.txt` (the at-scale database's in
`results/scale/dataset.sha256`).

```bash
python3 -I bench/scripts/cal_to_fasta.py scop40.cal.gz DATA/scop40.fa DATA/scop40.labels.tsv
python3 -I bench/scripts/make_datasets.py DATA DERIVED     # decoys, back-translations, stress sets
V13=... V12=... DMD=... DATA=... DER=... OUT=... bench/scripts/run_all.sh   # every section, ≈ 1.5 h on 2 cores
bench/scripts/run_all.sh scale scale_extra scale_algo       # only the at-scale comparison and its checks
python3 -I bench/scripts/summarize.py OUT                   # tables + plots
bench/scripts/collect_results.sh OUT bench/results          # machine-readable files only
```

Derived sets (`make_datasets.py`, fixed seed): composition-preserving decoys
(global shuffle; shuffle within 10-residue windows), back-translated SCOP40
domains with random synonymous codons, random flanks and random strand (blastx /
tblastn), 40-residue fragments, 400 two-domain queries (domains from different
SCOP classes joined by `GSGSGS`), low-complexity / tandem-repeat stress queries
and database, MMseqs database subsets (25 %, 50 %) and a ×4 database (the
original and three 15 %-mutated copies), long-protein queries (> 1,500 aa).

**At-scale database** (`scripts/make_scale_db.py`, run by the `scale`
section): the 11,206 SCOP40 domains followed by 80,000 reversed proteins
`REV_0 … REV_79999` — REV_i is the i-th sequence of the ×4 database written
backwards (91,206 sequences, 38.2 M residues). Reversed proteins keep the
length and composition of real proteins but have no homologs among the SCOP40
domains, so every hit to a `REV_` sequence is a false positive, like a hit to a
domain of another fold. The ×4 database holds each original sequence four
times (three mutated copies), so the decoys come in groups of four similar
sequences.

## Tools

* SABER 1.3.0 (this release) and SABER 1.2.0 (previous release), both built
  with rustc 1.75.0, `cargo build --release` (opt-level 3, thin LTO, one
  codegen unit, no `target-cpu`).
* DIAMOND 2.2.8 (official Linux binary), `blastp`/`blastx`, `--threads 2`.

## Machine and protocol

One 2-vCPU virtual machine (Intel Xeon @ 2.80 GHz, Cascade Lake, 33 MiB L3,
8 GB RAM, transparent huge pages in `madvise` mode), one session; every timing
record carries the command, the CPU model and a timestamp.

* Every timed run goes through `scripts/runbench.py`: wall time, user+system
  CPU and the child's peak RSS from `wait4()` (reported in MB of 2^20 bytes).
  The main configurations are repeated 3× (median reported, min–max given);
  single runs are marked.
* The machine is a shared virtual machine: the CPU time of identical runs
  varied by up to ~20 % over a day. The headline comparison (SABER 1.3
  `default` vs DIAMOND `--very-sensitive` at scale) is therefore run
  alternately, 5 times each, so that drift affects both tools alike.
  Two DIAMOND options that trade memory or strategy for speed are checked
  at scale: one seed-index chunk (`-c1`, section `scale_extra`, which also
  runs SABER with the 1.2 E-value model) and the query-indexed seed search
  (`--algo 1`, section `scale_algo`, alternating with SABER). DIAMOND's
  precomputed seed index (`makeidx`, searched with `--target-indexed -c1`)
  was tried once by hand: 58 s to build, and the search took 69.9 s instead
  of 34 s, so it is not part of the suite.
* A one-off search at scale with no preparation step (section `oneoff`):
  SABER searches the FASTA database directly (index built in memory and
  included in its time), DIAMOND runs `makedb` and then `blastp`
  (`--very-sensitive`) as one timed command; 5 alternating runs each.
* SCOP40 is searched with `-E 10` and no per-query limit (`-k 100000` /
  DIAMOND `-k 0`), so the evaluation sees complete ranked lists; "practical"
  timings use `-E 1e-3 -k 25`.
* Index construction is timed separately from search (SABER `--makedb`,
  DIAMOND `makedb`). At scale and in the index section both tools search
  prebuilt indexes; in the SCOP40 all-vs-all and MMseqs2 sections SABER builds
  its index in memory on every run (included in its time) and DIAMOND uses a
  prebuilt `.dmnd`.
* Cold-cache runs evict the index file with `posix_fadvise(DONTNEED)` (no root
  needed; the query file stays cached).
* **Parameter choices** (the 1.3 presets) were made on the **tune** splits only
  — superfamilies with even MD5 parity — of SCOP40 all-against-all and of the
  at-scale benchmark; the **test** split is held out. The weight-6 seeds of
  `sensitive` were chosen greedily for coverage of 18,102 homologous pairs of
  tune-split queries. The `sweep` section re-runs the candidate
  configurations with the release build (tune splits only).

## Metrics (`scripts/eval_scop.py`)

Pairs are classified with SCOP labels (Foldseek/MMseqs2 convention):
same family or superfamily = true positive, different fold = false positive,
same fold but different superfamily = ignored, self hits ignored. At scale
(`--decoy-prefix REV_ --decoys 80000`) every hit to a reversed protein is a
false positive too, and the 80,000 decoys count among the negatives.

* **sensitivity to the first FP** (per query, family / superfamily level,
  averaged over queries with ≥ 1 possible TP): fraction of possible TPs ranked
  above the query's first FP (by E-value).
* **pooled ROC**: hits of all queries sorted by E-value. Recall = TP / all
  same-superfamily pairs; precision = TP / (TP + FP); false-discovery
  proportion FDP = FP / (TP + FP); FPR = FP / (FP + TN) using all negative
  pairs; EPQ = FP per query.
* reported: recall at EPQ 0.01 / 0.1 / 1, recall at FDP 1 / 5 / 10 %, FDP at
  fixed recall, TP/FP/FDP/EPQ at E-value cut-offs (calibration).
* decoys: hits per query at E-value cut-offs when real queries are searched
  against shuffled databases (every hit is a false positive; a calibrated
  E-value predicts ≈ E hits per query).
* MMseqs workload (no labels): recall of SABER's exhaustive Smith–Waterman
  result (`--preset exhaustive`, E ≤ 1e-3 query–subject pairs), and recall of
  its E ≤ 1e-10 pairs (less dependent on each tool's E-value statistics).
* two-domain queries (`scripts/eval_multi.py`): fraction of each constituent
  domain's single-query true positives that are still found when the domain is
  part of a two-domain query.

E-values are not comparable between tools (and SABER 1.3's per-subject
E-values differ from 1.2's), so the threshold-free measures — sensitivity to
the first FP, recall at a given number of false positives per query, FDP at
matched recall — are the primary accuracy measures; counts at fixed E-value
cut-offs are reported alongside.

## Checks run by hand (`results/checks/`)

* `diamond_makeidx_target_indexed.jsonl`: DIAMOND's precomputed seed index
  at scale (`makeidx --very-sensitive`, then `blastp --very-sensitive -c1
  --target-indexed`; the first attempt without `-c1` is rejected by DIAMOND).
* `diamond_very-sensitive_verbose.log`: DIAMOND's own stage times at scale
  (`--verbose`): reference masking 1.1 s, reference histograms 3.3 s, then
  16 seed shapes.
* `masker_before_after.txt`, `index_build_*.jsonl`: the AVX2 tantan masker
  against the previous build (speed, identical masks, byte-identical
  indexes, identical search outputs) and the cost of index construction by
  component.
* `oneoff_trial_new_build.jsonl`: a first one-off trial (2 alternating runs)
  before the `oneoff` section was added.
* `lowc10_default_bounded_build_30min.jsonl`: all 10 adversarial
  low-complexity queries with the `default` seeds, no masking and no
  frequency cap (the `lowc` section runs the first 4): stopped by `timeout`
  after 30 minutes at a peak of 685 MB. `lowc10_default_before_pair_bounds_oom.jsonl`:
  the same run in an earlier run of the suite, with the build before the
  pair bounds (sha256 `b0066185…`, then installed as `saber-1.3.0`): killed
  by the out-of-memory killer at 6.0 GB after 448 s.

Checks after the final run, `scripts/extra_checks.sh` (the files here were
written by the same commands run by hand):

* `two_domain_database_evalues.metrics.json` (+ `_times.jsonl`): the
  two-domain benchmark with `--evalue-model database`, single-domain
  reference searched the same way.
* `few_queries_large_database.jsonl`, `single_query_stats.json`: 20 long
  proteins against the ×4 database (36 M residues) in one batch and one per
  batch, SABER 1.3 and 1.2.0; three of them as single-query searches of
  prebuilt indexes (SABER 1.3, SABER 1.2.0 with its postings index, DIAMOND
  `--very-sensitive`), and the stage counters of one.
* `single_query_postings_index.jsonl`: the same three queries against an
  index with seed postings (`--index-postings`, 3.4 GB): the query-major
  schedule (`--batch-order query`) and the scan schedule give the same
  alignments.

## Masking validation

`scripts/validate_masking.sh OUT FASTA...` runs the reference `tantan -p -f 1`
(tantan 49) and SABER's implementation (`cargo run --example mask_probs`) on the
same proteins and compares the per-residue repeat probabilities
(`scripts/compare_masking.py`). The results of the 1.2.0 validation
(`results/v1.2.0/masking_vs_tantan.txt`, 764,596 residues) apply unchanged:
1.3's AVX2 masker computes bit-identical probabilities
(`masking::tests::avx2_masker_is_bit_identical_to_the_scalar_body`; identical
masks for all 38.2 M residues of the at-scale database,
`results/checks/masker_before_after.txt`).

## 1.2.0 scripts

`scripts/defects.py` (defect reproductions of 1.0.2 / 1.1.0, inputs in
`bench/defects/`), `scripts/recheck.sh`, `scripts/translated_matched.sh` and
`scripts/devcheck.sh` belong to the 1.2.0 suite; their results are in
`results/v1.2.0/`.
