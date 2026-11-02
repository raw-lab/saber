# SABER benchmarks

Everything needed to reproduce the SABER 1.4.0 / 1.3.0 / DIAMOND comparison
in `bench/results/` (summary: `results/summary.md`, machine-readable:
`results/summary.json`, plots: `results/*.png`). The suites and results of
the earlier releases are kept unchanged with their READMEs: 1.3.0 (SABER
1.3.0 / 1.2.0 / DIAMOND) in `results/v1.3.0/`, 1.2.0 (SABER 1.0.2 / 1.1.0 /
1.2.0 / DIAMOND, three machines) in `results/v1.2.0/`. Nothing here is
synthetic except the derived decoy and stress sets, which are generated
deterministically from the public data.

## Data

| File | Source | SHA-256 |
|---|---|---|
| `scop40.fa`, `scop40.labels.tsv` | SCOPe 2.08 40 % domains, converted from `test_data/scop40.cal.gz` of https://github.com/rcedgar/reseek (commit `9cb33e1`, `scop40.cal.gz` sha256 `69808731…`) with `scripts/cal_to_fasta.py` | `93268cef…` / `bbbcf71d…` |
| `mmseqs_QUERY.fasta`, `mmseqs_DB.fasta` | `examples/QUERY.fasta`, `examples/DB.fasta` of https://github.com/soedinglab/MMseqs2 (commit `76de5d8`) — 500 queries, 20,000 UniProt sequences, 9.06 M residues | `c99bc94a…` / `55d48bb7…` |

Full SHA-256 checksums of the inputs, the derived sets (including the scaling
databases) and the binaries are recorded in `results/env/environment.txt`
(also per section: `results/scale/dataset.sha256`,
`results/scaling/dataset.sha256`, `results/lenmatch/dataset.sha256`).

```bash
python3 -I bench/scripts/cal_to_fasta.py scop40.cal.gz DATA/scop40.fa DATA/scop40.labels.tsv
python3 -I bench/scripts/make_datasets.py DATA DERIVED     # decoys, back-translations, stress sets
VN=... VO=... DMD=... DATA=... DER=... OUT=... bench/scripts/run_all.sh   # every section, ≈ 2.5 h on 2 vCPUs
bench/scripts/run_all.sh scale scale_extra scale_algo       # only the at-scale comparison and its checks
SCALING_K="1 4 16" bench/scripts/run_all.sh scaling        # database-size scaling at three sizes
bench/scripts/cross_machine_check.sh                        # headline pairs again (another machine), ≈ 15 min
bench/scripts/oneoff_check.sh                               # one-off search: SABER 1.4, 1.3.0, DIAMOND, ≈ 10 min
bench/scripts/extra_checks_14.sh                            # 1.3.0 counters, the frequency cap, 1.3 vs 1.4 outputs
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

**Scaling databases** (`scripts/make_scaling_dbs.py`, run by the `scaling`
section): the 11,206 SCOP40 domains followed by K × 20,000 reversed proteins
for K = 1, 2, 4, 8, 16 — 31,206 to 331,206 sequences, 11.0 to 146.8 M
residues. The decoy sources are, in order, the ×4 database and then further
copies of the MMseqs2 database at 15 % substitution (fixed seed), so
`scale_k4.fa` is the at-scale database (identical file) and every database
contains the decoys of the smaller ones; each original protein gives K
similar decoys.

**Length-matched decoys** (`scripts/make_lenmatched_db.py`, run by the
`lenmatch` section): the at-scale database's 80,000 reversed proteins cut
into consecutive pieces with the lengths of SCOP40 domains (lengths drawn
from SCOP40 with a fixed seed; a protein's last piece is dropped if shorter
than 30 residues) — 244,143 decoys, 255,349 sequences, 37.2 M residues. True
and false subjects then have the same length distribution, which removes the
advantage that per-subject E-values get from decoys longer than the domains.

## Tools

* SABER 1.4.0 (this release) and SABER 1.3.0 (previous release), both built
  with rustc 1.75.0, `cargo build --release` (opt-level 3, thin LTO, one
  codegen unit, no `target-cpu`). In the `scaling` and `qscale` sections
  SABER 1.3.0 searches the indexes built by 1.4.0 (the same index format;
  1.3.0's own index of the at-scale database gives the same output).
* DIAMOND 2.2.8 (official Linux binary), `blastp`/`blastx`, `--threads 2`.

## Machine and protocol

One session on 2-vCPU virtual machines with 8 GB RAM and transparent huge
pages in `madvise` mode. Every timing record carries the command, the CPU
model and a timestamp, and `summary.md` names the machine of every section
(machines are named A, B, … by CPU type, in the order they first appear):

* **Machine A** (Intel Xeon @ 2.80 GHz, 33 MiB L3): sections
  `env`, `scale`, `scaling`, `lenmatch`, `qscale` and `scop` (the 1.4.0 run
  took 2 h 23 min in all, from 20:54 to 23:17).
* **Machine B** (Intel Xeon @ 2.10 GHz, 260 MiB L3): the virtual
  machine was replaced during the `decoy` section (22:33). That section's
  partial output was removed, and it and every later section ran on the new
  machine (`run.log` records both events). `cross_machine_check.sh` then
  repeated the headline pairs on machine B: SABER `default` and DIAMOND
  `--very-sensitive` alternately at scale (5 runs each, prebuilt indexes of
  the `scale` section) and at K = 16 (indexes built again, 3 runs each). It
  also compared their outputs with machine A's (all identical).
* A further replacement gave a new virtual machine of machine A's CPU type,
  on which `oneoff_check.sh` ran the one-off search at scale: SABER 1.4,
  SABER 1.3.0 and DIAMOND `makedb` + `blastp --very-sensitive` alternately,
  5 runs each. A machine of a type seen earlier is not necessarily the same
  host.
* Comparisons between tools are made only within one machine; the
  SABER/DIAMOND time ratio differs between A and B (see `summary.md`,
  "Cross-machine check").

* Every timed run goes through `scripts/runbench.py`: wall time, user+system
  CPU and the child's peak RSS from `wait4()` (reported in MB of 2^20 bytes).
  The main configurations are repeated 3× (median reported, min–max given);
  single runs are marked.
* These are shared virtual machines whose speed varies (in earlier sessions
  the CPU time of identical runs varied by up to ~20 % over a day; here the
  third machine ran the one-off searches 30–45 % slower than machine A). The
  headline comparison (SABER 1.4 `default` vs DIAMOND `--very-sensitive` at
  scale) is therefore run alternately, 5 times each, so that drift affects
  both tools alike.
  Two DIAMOND options that trade memory or strategy for speed are checked
  at scale: one seed-index chunk (`-c1`, section `scale_extra`, which also
  runs SABER with the 1.2 E-value model) and the query-indexed seed search
  (`--algo 1`, section `scale_algo`, alternating with SABER). DIAMOND's
  precomputed seed index (`makeidx`, searched with `--target-indexed -c1`)
  was tried once by hand in the 1.3 session: 58 s to build, and the search
  took 69.9 s instead of 33 s, so it is not part of the suite.
* A one-off search at scale with no preparation step (section `oneoff`):
  SABER searches the FASTA database directly (index built in memory and
  included in its time), DIAMOND runs `makedb` and then `blastp`
  (`--very-sensitive`) as one timed command; 5 alternating runs each.
* **Scaling with database size** (section `scaling`, K = 1, 2, 4, 8, 16):
  per size, index construction by both tools (single runs); SABER 1.4
  `default` and DIAMOND `--very-sensitive` alternately, 3 runs each (median
  and range reported); SABER `very-fast`, `more-sensitive`, SABER 1.3.0
  `default`, DIAMOND default and `--ultra-sensitive` (limited to 1 hour) one
  run each; the one-off pair (SABER on the FASTA, DIAMOND `makedb` +
  `blastp`) one run each. Outputs are evaluated with `--decoys` = K × 20,000.
  The section checks that the in-memory and prebuilt-index outputs are
  identical (`scaling/outputs_check.txt`) and records the database sizes
  (`scaling/dbsizes.tsv`). The sizes can be run in separate invocations
  (`SCALING_K`); the published results ran all five in one invocation.
* **Length-matched decoys** (section `lenmatch`): SABER `default` and DIAMOND
  `--very-sensitive` alternately, 3 runs each; SABER `default` with database
  E-values, `more-sensitive` and DIAMOND `--ultra-sensitive` one run each.
* **Scaling with the number of queries** (section `qscale`): 1, 10, 100,
  1,000 and all 11,206 SCOP40 domains (every N-th) against the at-scale
  indexes; SABER 1.4 `default`, SABER 1.3.0 `default` and DIAMOND
  `--very-sensitive` in turn, 3 runs each.
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
* **Parameter choices** (the 1.3 presets, unchanged in 1.4) were made on the **tune** splits only
  — superfamilies with even MD5 parity — of SCOP40 all-against-all and of the
  at-scale benchmark; the **test** split is held out. The weight-6 seeds of
  `sensitive` were chosen greedily for coverage of 18,102 homologous pairs of
  tune-split queries. The `sweep` section re-runs the candidate
  configurations with the release build (tune splits only). 1.4's one new
  parameter, the reference database size of the region heuristics (40 M
  residues), was not tuned: it is about the size of the at-scale database
  the presets were tuned on (38.2 M residues).

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

* `diamond_very-sensitive_k16_verbose.log`: DIAMOND's own stage times on the
  K = 16 scaling database (`--verbose`, same binary and settings as the
  suite; run during development): reference masking 3.5 s, reference
  histograms 10.7 s, then per seed shape about 1.0 s for the reference seed
  array, 0.8 s for the hash join and 2.1 s for the ungapped search, 8.0 s of
  gapped alignment; 87.5 s in all.
* `v13_v14_outputs.txt` (`scripts/compare_scaling_outputs.py OUT
  scop40.labels.tsv`, run on the suite's outputs): SABER 1.4 and 1.3.0
  `default` outputs compared per scaling database, and the true positives
  each tool loses between K = 4 and K = 16 (and K = 4 → 8, K = 1 → 16)
  although their E-value, scaled to the larger database, stays at most 10.
  `v13_v14_outputs_build_ff519871.txt` is the same comparison for a first
  run of the `scale`, `scaling`, `lenmatch` and `qscale` sections with the
  build before the scaled frequency-cap floor (sha256 `ff519871…`); that run
  was superseded by the release build's, and its file is kept because it
  shows the effect of the region change alone.
* `scop40_v13_v14_ab.jsonl`: SABER 1.4 and 1.3.0 `default` on SCOP40
  all-against-all, alternately, 5 runs each, on machine B (the `scop`
  section runs the two versions one after the other). Timings only; the
  outputs went to `/dev/null`.
* `first_run_ff519871/`: timings and `--stats` of the superseded first run
  (see `v13_v14_outputs_build_ff519871.txt`); its `README.txt` explains why
  its counters show the query seeds that 1.3's cap rule skips.
* `dev_ab/`: the development A/B measurements behind the 1.4 speed items —
  the adopted window-filter change and the variants that were measured and
  not adopted (entry prefetch distances, larger scan rounds, gathers in the
  X-drop, range-partitioned seed counting at index construction); its
  `README.txt` lists the binaries and files.

Scripted checks run after the suite, each with its own `environment.txt`:

* `results/checks_machine/` (`scripts/cross_machine_check.sh`, machine B):
  the headline pairs at scale and at K = 16, and their outputs compared with
  machine A's.
* `results/checks_oneoff/` (`scripts/oneoff_check.sh`, the third machine):
  the one-off search at scale with SABER 1.4, SABER 1.3.0 and DIAMOND.
* `results/checks_extra/` (`scripts/extra_checks_14.sh`, the third
  machine): SABER 1.3.0's own `--stats` on the 74.4 and 146.8 M-residue
  scaling databases (candidate regions, skipped query seeds; outputs
  compared with the scaling section's); SABER 1.4 `default` at scale with
  and without the automatic frequency cap, alternately, with accuracy on
  both splits; and byte comparisons of the suite's 1.4 and 1.3.0 outputs.
  The length-matched database's size (`results/lenmatch/datasets.txt`) was
  counted afterwards with the `awk` command in that script.

The 1.3.0 checks (DIAMOND's `makeidx`, the adversarial low-complexity runs,
few-query searches, two-domain E-value models, …) are in
`results/v1.3.0/checks/`.

## Masking validation

`scripts/validate_masking.sh OUT FASTA...` runs the reference `tantan -p -f 1`
(tantan 49) and SABER's implementation (`cargo run --example mask_probs`) on the
same proteins and compares the per-residue repeat probabilities
(`scripts/compare_masking.py`). The results of the 1.2.0 validation
(`results/v1.2.0/masking_vs_tantan.txt`, 764,596 residues) apply unchanged:
1.3's AVX2 masker computes bit-identical probabilities
(`masking::tests::avx2_masker_is_bit_identical_to_the_scalar_body`; identical
masks for all 38.2 M residues of the at-scale database,
`results/v1.3.0/checks/masker_before_after.txt`); 1.4 did not change the
masker.

## Scripts of earlier suites

`scripts/defects.py` (defect reproductions of 1.0.2 / 1.1.0, inputs in
`bench/defects/`), `scripts/recheck.sh`, `scripts/translated_matched.sh` and
`scripts/devcheck.sh` belong to the 1.2.0 suite; their results are in
`results/v1.2.0/`. `scripts/extra_checks.sh` belongs to the 1.3.0 suite
(results in `results/v1.3.0/checks/`).
