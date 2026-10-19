# SABER benchmarks (1.2.0)

> **SABER 1.2.0 benchmark (kept for reference).** This README and the
> results in this directory belong to the 1.2.0 release. Paths `results/…`
> below refer to this directory (`bench/results/v1.2.0/`); `scripts/…` to
> `bench/scripts/`, except `run_all.sh`, `summarize.py` and
> `collect_results.sh`, whose 1.2.0 versions are in `scripts/` here (the
> 1.3.0 suite replaced them). The 1.3.0 benchmark is described in
> `bench/README.md`.

Everything needed to reproduce the 1.0.2 / 1.1.0 / 1.2.0 / DIAMOND comparison
in `bench/results/`. Nothing here is synthetic except the derived decoy and
stress sets, which are generated deterministically from the public data.

## Data

| File | Source | SHA-256 |
|---|---|---|
| `scop40.fa`, `scop40.labels.tsv` | SCOPe 2.08 40 % domains, converted from `test_data/scop40.cal.gz` of https://github.com/rcedgar/reseek (commit `9cb33e1`, `scop40.cal.gz` sha256 `69808731…`) with `scripts/cal_to_fasta.py` | `93268cef…` / `bbbcf71d…` |
| `mmseqs_QUERY.fasta`, `mmseqs_DB.fasta` | `examples/QUERY.fasta`, `examples/DB.fasta` of https://github.com/soedinglab/MMseqs2 (commit `76de5d8`) — 500 queries, 20,000 UniProt sequences, 9.06 M residues | `c99bc94a…` / `55d48bb7…` |

Full SHA-256 checksums of the inputs are recorded in `results/env/environment_host*.txt`, those of the five binaries (SABER 1.0.2, 1.1.0, the benchmarked and the delivered 1.2.0 build, DIAMOND) in `results/env/final_binary.txt`.

```bash
python3 -I bench/scripts/cal_to_fasta.py scop40.cal.gz DATA/scop40.fa DATA/scop40.labels.tsv
python3 -I bench/scripts/make_datasets.py DATA DERIVED     # decoys, back-translations, stress sets
bench/scripts/run_all.sh                                    # all sections (≈ 3–4 h on 2 cores)
python3 -I bench/scripts/summarize.py RESULTS               # tables + plots
```

Derived sets (`make_datasets.py`, fixed seed): composition-preserving decoys
(global shuffle; shuffle within 10-residue windows), back-translated SCOP40
domains with random synonymous codons, random flanks and random strand (blastx /
tblastn), 40-residue fragments, 400 two-domain queries (domains from different
SCOP classes joined by `GSGSGS`), low-complexity / tandem-repeat stress queries
and database, MMseqs database subsets (25 %, 50 %) and a ×4 database (three
15 %-mutated copies), long-protein queries (> 1,500 aa).

## Tools

* SABER 1.0.2 (`saber-1.0.2-u64-offsets.tar.gz`) and 1.1.0
  (`saber-optimized-v1.1.0.tar.gz`), unmodified, and SABER 1.2.0 — all built
  with rustc 1.75.0, `cargo build --release` (opt-level 3, thin LTO, one
  codegen unit, no `target-cpu`). The 1.1.0 package did not bump its version
  string, so its binary prints `saber 1.0.2`; `results/env/` identifies the
  binaries by SHA-256.
* DIAMOND 2.2.8 (official Linux binary), `blastp`/`blastx`, `--threads 2`.

## Machines

The session's container moved twice during the final run, so the published
results come from three 2-vCPU virtual machines with 8 GB RAM: **A1** (Intel
Xeon @ 2.80 GHz, 33 MiB L3) for the SCOP40, decoy, MMseqs2 and index
sections; **B** (Intel Xeon @ 2.10 GHz, 260 MiB L3; 1.7–2.1× faster on SCOP40
all-vs-all: SABER 1.2 1.7×, DIAMOND 1.8×, SABER 1.1.0 2.1×) for the
two-domain, thread-scaling and preset-sweep sections; **A2** (same CPU model
as A1, a different VM instance) for the translated and low-complexity
sections. Every timing record carries `cpu` and `host` fields, `summary.md`
names the host under each section heading, and timings are compared only
within a section.

The machine labels are partly inferred. `runbench.py` began recording the CPU
model during the final run; earlier records (all of A1's, and B's before
21:34) were labelled from their timestamps (container moves at 21:15 and
21:45). In the A1 sections the DIAMOND and SABER 1.0.2 / 1.1.0 rows were
measured between 14:34 and 17:17 and reused with `REUSE=1`, the SABER 1.2 rows
between 20:10 and 20:37. A move to another instance of the same CPU model in
between would not have been detected; `results/recheck/` re-measures two A1
configurations later on a machine of that model.

## Protocol

* Every timed run goes through `scripts/runbench.py`: wall time, user+system
  CPU and the child's peak RSS from `wait4()`, plus the CPU model. The main
  SCOP40, MMseqs2, index and thread configurations are repeated 3× (median
  reported, min–max given). Single runs: SABER 1.0.2 / 1.1.0 where slow,
  exhaustive searches, the decoy, translated, low-complexity, two-domain,
  long-protein, scaling and preset-sweep runs (`summary.md` shows no range for
  them).
* `REUSE=1 run_all.sh …` keeps successful baseline records (DIAMOND, SABER
  1.0.2 / 1.1.0) from a previous invocation instead of re-measuring them; SABER
  1.2 runs are always re-measured.
* SCOP40 is searched all-against-all with `-E 10` and no per-query limit
  (`-k 100000` / DIAMOND `-k 0`), so the evaluation sees complete ranked lists;
  a second "practical" timing uses `-E 1e-3 -k 25`.
* Index construction is timed separately from search (SABER `--makedb`, DIAMOND
  `makedb`); DIAMOND searches always use a prebuilt `.dmnd`; SABER is shown both
  building its index in memory on every run (its normal mode) and with `--db`.
  Cold-cache runs evict the index file with `posix_fadvise(DONTNEED)` (no root
  needed; the FASTA query file stays cached).
* Parameter choices (presets) were made on the SCOP40 **tune** split only
  (superfamilies with even MD5 parity); the **test** split is held out.

## Follow-up runs (after the main run, one session each)

* `scripts/recheck.sh` — the SABER-vs-DIAMOND timing comparisons of the A1
  sections (SCOP40 general tier, MMseqs2 in-memory and `--db`, ×4 database,
  long proteins) re-timed with both tools in one session, with the delivered
  build; output identity is checked against the main run.
* `scripts/translated_matched.sh` — the 34-query blastx / tblastn subsets
  searched by SABER 1.1.0 (full six-frame dynamic programming), SABER 1.2
  `default`, `more-sensitive` and `exhaustive`, and DIAMOND on one machine, so
  translated speed can be compared at matched sensitivity.
* `scripts/devcheck.sh` — a development build from before the performance
  changes against the benchmarked build: time, peak RSS, stage times and
  output differences (`scripts/compare_m8.py`).
* `results/masking_f32_vs_f64.txt` — single- vs double-precision tantan on the
  masking validation set.

`summary.md` reports them in their own sections ("Re-check", "Translated
subsets at matched sensitivity", "Before / after").

## Metrics (`scripts/eval_scop.py`)

Pairs are classified with SCOP labels (Foldseek/MMseqs2 convention):
same family or superfamily = true positive, different fold = false positive,
same fold but different superfamily = ignored, self hits ignored.

* **sensitivity to the first FP** (per query, family / superfamily level,
  averaged over queries with ≥ 1 possible TP): fraction of possible TPs ranked
  above the query's first FP (by E-value).
* **pooled ROC**: hits of all queries sorted by E-value. Recall = TP / all
  same-superfamily pairs; precision = TP / (TP + FP); false-discovery
  proportion FDP = FP / (TP + FP); FPR = FP / (FP + TN) using all different-fold
  pairs as negatives; EPQ = FP per query.
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

## Masking validation

`scripts/validate_masking.sh OUT FASTA...` runs the reference `tantan -p -f 1`
(tantan 49) and SABER's implementation (`cargo run --example mask_probs`) on the
same proteins and compares the per-residue repeat probabilities
(`scripts/compare_masking.py`): agreement of the ≥ 0.5 mask calls used for
seeding, Pearson r, mean absolute difference. `results/masking_vs_tantan.txt`
was produced on the first 3,000 domains of `scop40.fa` plus the 500 MMseqs2
example queries (764,596 residues; commands in the script header).

## Defect reproductions

`scripts/defects.py V102 V110 V12 defects > results/defects.json` runs the
three SABER versions with identical arguments on the small inputs in
`bench/defects/` (plus the extra runs listed in the script's `EXTRA`) and
records score, E-value, ranges, CIGAR consistency, gap fields, coverage and
stray stdout lines:

| Case | Input | Shows |
|---|---|---|
| `amb` | subject with 23 X/B/Z | reported score unaffected by the 255 sentinel in 1.x (963); 1.2 gives 963 with `--comp-bias off`, 950 with its default composition correction |
| `amb2` | 120-residue pair identical except for an X at every fifth position | 1.0.2 / 1.1.0 report no hit (their ungapped extension stops at the sentinel; 1.1.0 `--no-index` finds the pair, E = 5e-73); 1.2 reports it |
| `gap` | alignment with a 5- and a 3-residue gap | 1.x CIGAR swaps `I`/`D` (consumes 152 query residues for a 150-residue range), column 6 counts gap columns, E-value from ungapped λ, K (score 809: E = 7.9e-109 in 1.x, 1.8e-83 in 1.2 with `--comp-bias off`) |
| `long` | 3,000-residue identical pair | 1.x truncates the alignment to 2,048 columns |
| `sat` | 4,000-residue poly-W pair (score 44,000) | 1.x reports 22,528 (2,048-column cap); 1.2 reports 44,000 with masking and composition correction disabled (its default search does not report a pure homopolymer pair) |
