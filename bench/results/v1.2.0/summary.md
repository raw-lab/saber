## SCOP40 all-vs-all (11,206 domains, 2 threads, E ≤ 10, no hit limit)

Measured on: host A1 (Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB)).

Sensitivity columns are on the **held-out test split**; tune-split sensitivity in brackets. Time = wall seconds, median of runs (range); RSS = peak resident memory.

| Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP (superfamily) | Sens. (family) | Recall @ EPQ 0.01 | Recall @ FDP 5 % | TP at E ≤ 1e-3 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|---|---|---|
| saber_v1.0.2 | 107.0 (105.7–108.6) | 149.9 | 27 | 0.158 [0.150] | 0.307 | 0.0645 | 0.0723 | 12397 | 63 |
| saber_v1.1.0 | 43.5 (42.4–46.1) | 77.0 | 32 | 0.158 [0.150] | 0.307 | 0.0645 | 0.0723 | 12397 | 63 |
| saber_v1.2_very-fast | 4.2 (4.2–4.5) | 7.4 | 28 | 0.110 [0.111] | 0.220 | 0.0456 | 0.0487 | 7766 | 5 |
| saber_v1.2_fast | 5.6 (5.5–5.8) | 9.9 | 28 | 0.153 [0.150] | 0.299 | 0.0625 | 0.0728 | 10946 | 10 |
| saber_v1.2_default | 10.0 (9.9–10.4) | 18.5 | 30 | 0.197 [0.187] | 0.380 | 0.0805 | 0.0977 | 13911 | 10 |
| saber_v1.2_sensitive | 13.1 (12.7–13.2) | 24.0 | 30 | 0.203 [0.192] | 0.389 | 0.0825 | 0.1003 | 14217 | 10 |
| saber_v1.2_more-sensitive | 20.5 (19.7–20.5) | 38.7 | 37 | 0.225 [0.206] | 0.426 | 0.0888 | 0.1103 | 15331 | 11 |
| saber_v1.2_legacy | 8.5 (8.2–8.6) | 15.8 | 26 | 0.150 [0.144] | 0.300 | 0.0624 | 0.0698 | 10836 | 7 |
| saber_v1.2_exhaustive | 391.9 | 768.2 | 36 | 0.265 [0.245] | 0.496 | 0.0943 | 0.1254 | 16754 | 15 |
| diamond_fast | 1.7 (1.7–1.8) | 3.3 | 27 | 0.013 [0.012] | 0.026 | 0.0036 | 0.0036 | 673 | 0 |
| diamond_default | 3.0 (2.7–3.1) | 5.5 | 32 | 0.045 [0.045] | 0.092 | 0.0155 | 0.0157 | 2775 | 1 |
| diamond_mid-sensitive | 7.5 (7.3–8.4) | 14.2 | 33 | 0.099 [0.098] | 0.202 | 0.0402 | 0.0419 | 7148 | 6 |
| diamond_sensitive | 14.8 (14.2–14.8) | 27.5 | 41 | 0.143 [0.141] | 0.291 | 0.0656 | 0.0690 | 11857 | 14 |
| diamond_more-sensitive | 16.5 (15.4–16.7) | 30.7 | 48 | 0.143 [0.142] | 0.292 | 0.0661 | 0.0694 | 11940 | 14 |
| diamond_very-sensitive | 13.0 (12.2–13.3) | 24.5 | 58 | 0.188 [0.176] | 0.370 | 0.0851 | 0.0973 | 15433 | 17 |
| diamond_ultra-sensitive | 45.3 (42.3–46.1) | 86.1 | 115 | 0.218 [0.200] | 0.423 | 0.0959 | 0.1150 | 17395 | 22 |

TP / FP at E ≤ 1e-3 = same-superfamily / different-fold pairs reported by the test-split queries at that cut-off. E-values are not interchangeable between tools: the same cut-off yields more true and more false positives from DIAMOND than from the SABER preset of similar sensitivity to the first FP.

### False-discovery proportion at matched recall (test split)

Pooled over all test-split queries: hits sorted by E-value, FDP = FP / (TP + FP) at the most stringent cut-off that reaches the recall (TP / all same-superfamily pairs). 'not reached' = the tool never reaches that recall at E ≤ 10.

| Tool / preset | FDP @ recall 0.06 | FDP @ recall 0.08 | FDP @ recall 0.09 | FDP @ recall 0.1 | FDP @ recall 0.11 |
|---|---|---|---|---|---|
| saber_v1.0.2 | 0.0008 | 0.3678 | not reached | not reached | not reached |
| saber_v1.1.0 | 0.0008 | 0.3691 | not reached | not reached | not reached |
| saber_v1.2_very-fast | not reached | not reached | not reached | not reached | not reached |
| saber_v1.2_fast | 0.0021 | not reached | not reached | not reached | not reached |
| saber_v1.2_default | 0.0001 | 0.0035 | 0.0138 | 0.0710 | 0.2844 |
| saber_v1.2_sensitive | 0.0001 | 0.0021 | 0.0099 | 0.0458 | 0.1871 |
| saber_v1.2_more-sensitive | 0.0001 | 0.0007 | 0.0045 | 0.0130 | 0.0458 |
| saber_v1.2_legacy | 0.0020 | not reached | not reached | not reached | not reached |
| saber_v1.2_exhaustive | 0.0000 | 0.0005 | 0.0016 | 0.0058 | 0.0136 |
| diamond_fast | not reached | not reached | not reached | not reached | not reached |
| diamond_default | not reached | not reached | not reached | not reached | not reached |
| diamond_mid-sensitive | not reached | not reached | not reached | not reached | not reached |
| diamond_sensitive | 0.0005 | not reached | not reached | not reached | not reached |
| diamond_more-sensitive | 0.0004 | not reached | not reached | not reached | not reached |
| diamond_very-sensitive | 0.0000 | 0.0007 | 0.0073 | not reached | not reached |
| diamond_ultra-sensitive | 0.0000 | 0.0003 | 0.0009 | 0.0048 | 0.0187 |

### Practical settings (E ≤ 1e-3, 25 targets per query), SCOP40 all-vs-all

| Run | Time s | RSS MB |
|---|---|---|
| saber_v1.2_fast | 4.3 (4.2–4.3) | 26 |
| saber_v1.2_default | 9.0 (9.0–9.4) | 26 |
| saber_v1.1.0 | 31.1 | 31 |
| diamond_default | 2.8 (2.7–2.8) | 29 |
| diamond_very-sensitive | 13.0 (12.7–13.2) | 57 |
| diamond_ultra-sensitive | 43.8 (43.4–44.2) | 115 |

### SCOP40 index build and persistent-index search

| Run | Time s | RSS MB |
|---|---|---|
| diamond_makedb | 0.1 (0.1–0.1) | 19 |
| saber_v1.2_makedb_default | 0.4 (0.4–0.4) | 17 |
| saber_v1.2_default_db_warm | 9.1 (8.8–9.3) | 27 |

## Null behaviour: real SCOP40 queries against shuffled databases

Measured on: host A1 (Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB)).

Real queries against composition-preserving decoys (`shuf`: whole-sequence shuffle; `wshuf10`: shuffle inside 10-residue windows). Hits to decoys derived from a sequence of the query's own fold are excluded (a windowed shuffle of a homolog keeps part of its similarity); every remaining hit is a false positive. A calibrated E-value gives at most about E false hits per query (EPQ ≤ E).

| Run | EPQ at E ≤ 1e-3 | E ≤ 0.01 | E ≤ 0.1 | E ≤ 1 | E ≤ 10 |
|---|---|---|---|---|---|
| diamond_default_shuf | 0.0000 | 0.0001 | 0.0008 | 0.0029 | 0.0103 |
| diamond_default_wshuf10 | 0.0000 | 0.0001 | 0.0010 | 0.0034 | 0.0106 |
| diamond_ultra-sensitive_shuf | 0.0012 | 0.0064 | 0.0242 | 0.0920 | 0.2509 |
| diamond_ultra-sensitive_wshuf10 | 0.0007 | 0.0073 | 0.0327 | 0.1116 | 0.2996 |
| diamond_very-sensitive_shuf | 0.0006 | 0.0031 | 0.0126 | 0.0431 | 0.1070 |
| diamond_very-sensitive_wshuf10 | 0.0004 | 0.0046 | 0.0180 | 0.0535 | 0.1258 |
| saber_v1.1.0_shuf | 0.0026 | 0.0099 | 0.0620 | 0.3585 | 1.8990 |
| saber_v1.1.0_wshuf10 | 0.0034 | 0.0115 | 0.0663 | 0.3843 | 2.0949 |
| saber_v1.2_default_shuf | 0.0003 | 0.0014 | 0.0162 | 0.1872 | 1.7155 |
| saber_v1.2_default_wshuf10 | 0.0002 | 0.0025 | 0.0207 | 0.2204 | 1.9460 |
| saber_v1.2_fast_shuf | 0.0003 | 0.0012 | 0.0081 | 0.0894 | 0.7970 |
| saber_v1.2_fast_wshuf10 | 0.0000 | 0.0009 | 0.0095 | 0.1045 | 0.9105 |
| saber_v1.2_more-sensitive_shuf | 0.0006 | 0.0029 | 0.0270 | 0.2963 | 2.7860 |
| saber_v1.2_more-sensitive_wshuf10 | 0.0004 | 0.0038 | 0.0311 | 0.3475 | 3.1271 |

## MMseqs2 example workload (500 queries × 20,000 UniProt sequences, E ≤ 1e-3)

Measured on: host A1 (Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB)).

Recall is measured against SABER's own exhaustive Smith–Waterman search (`--preset exhaustive`: 16542 query–subject pairs at E ≤ 1e-3). That reference shares SABER's scoring, composition correction and E-value statistics, so near the cut-off it favours SABER; 'Recall, E ≤ 1e-10 pairs' uses only the 12797 reference pairs with E ≤ 1e-10, which a tool that finds them reports at E ≤ 1e-3 whatever its statistics. 'Extra' = pairs a tool reports that the reference does not (mostly borderline E-values under different statistics). 74 of the 500 queries are also database members; self-hits count for every tool.

| Tool / preset | Time s | RSS MB | Pairs | Recall of exhaustive SW | Recall, E ≤ 1e-10 pairs | Extra pairs |
|---|---|---|---|---|---|---|
| saber_v1.0.2 | 68.0 | 54 | 14906 | 0.7143 | 0.8211 | 3090 |
| saber_v1.1.0 | 30.1 (29.5–31.2) | 68 | 14905 | 0.7143 | 0.8211 | 3089 |
| saber_v1.2_very-fast | 5.3 (5.0–5.3) | 104 | 12878 | 0.7784 | 0.8842 | 2 |
| saber_v1.2_fast | 5.0 (4.8–5.5) | 81 | 14399 | 0.8703 | 0.9434 | 2 |
| saber_v1.2_default | 7.4 (7.3–8.1) | 106 | 15607 | 0.9434 | 0.9794 | 2 |
| saber_v1.2_default_db | 6.0 (5.9–6.1) | 104 | 15607 | 0.9434 | 0.9794 | 2 |
| saber_v1.2_default_nofloor | 13.5 | 107 | 15793 | 0.9546 | 0.9834 | 2 |
| saber_v1.2_sensitive | 14.4 (14.3–15.0) | 108 | 15878 | 0.9597 | 0.9853 | 2 |
| saber_v1.2_exhaustive | 196.0 | 140 | 16542 | 1.0000 | 1.0000 | 0 |
| diamond_fast | 1.5 (1.5–1.5) | 47 | 4323 | 0.2610 | 0.3314 | 5 |
| diamond_default | 2.1 (2.1–2.3) | 47 | 6801 | 0.4091 | 0.5064 | 34 |
| diamond_mid-sensitive | 5.8 (5.7–6.7) | 56 | 10561 | 0.6323 | 0.7355 | 101 |
| diamond_sensitive | 10.9 (10.7–11.9) | 61 | 12742 | 0.7592 | 0.8615 | 184 |
| diamond_very-sensitive | 7.1 (7.0–7.7) | 113 | 15425 | 0.9143 | 0.9623 | 300 |

### Scaling (single run each)

| Run | Time s | RSS MB |
|---|---|---|
| scale_db25_saber_v1.1.0 | 9.2 | 25 |
| scale_db25_saber_v1.2_default | 2.6 | 47 |
| scale_db25_diamond_very-sensitive | 2.5 | 38 |
| scale_db25_diamond_default | 0.7 | 35 |
| scale_db50_saber_v1.1.0 | 17.4 | 39 |
| scale_db50_saber_v1.2_default | 4.2 | 72 |
| scale_db50_diamond_very-sensitive | 4.2 | 70 |
| scale_db50_diamond_default | 1.3 | 37 |
| scale_db400_saber_v1.1.0 | 75.7 | 256 |
| scale_db400_saber_v1.2_default | 26.4 | 251 |
| scale_db400_diamond_very-sensitive | 27.5 | 362 |
| scale_db400_diamond_default | 8.1 | 133 |
| scaleq_50_saber_v1.1.0 | 3.0 | 53 |
| scaleq_50_saber_v1.2_default | 2.2 | 63 |
| scaleq_50_diamond_very-sensitive | 4.5 | 96 |
| scaleq_100_saber_v1.1.0 | 6.1 | 54 |
| scaleq_100_saber_v1.2_default | 2.7 | 66 |
| scaleq_100_diamond_very-sensitive | 5.0 | 109 |
| scaleq_250_saber_v1.1.0 | 14.7 | 67 |
| scaleq_250_saber_v1.2_default | 4.2 | 75 |
| scaleq_250_diamond_very-sensitive | 5.6 | 112 |

### Long proteins (99 queries > 1,500 aa against the 20,000-sequence database)

The queries are every fifth database sequence longer than 1,500 residues, so each also finds itself (self-hits count for every tool). Single runs. The SABER 1.2 row is the benchmarked build, whose traceback memory cap lowered the score of 9 of these alignments; the delivered build restores them (sections 'Re-check' and 'Before / after').

| Run | Time s | RSS MB | Pairs at E ≤ 1e-3 |
|---|---|---|---|
| saber_v1.1.0 | 31.8 | 76 | 8697 |
| saber_v1.2_default | 8.4 | 115 | 4712 |
| diamond_very-sensitive | 8.1 | 105 | 4672 |

## Index construction, size and persistent-index search (MMseqs database and ×4 database)

Measured on: host A1 (Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB)).

| Run | Time s | RSS MB |
|---|---|---|
| makedb_saber_v1.2_mmseqs_DB | 2.1 (1.9–2.1) | 60 |
| makedb_diamond_mmseqs_DB | 0.5 (0.4–0.5) | 28 |
| makedb_saber_v1.1.0_mmseqs_DB | 0.5 | 44 |
| db_warm_saber_v1.2_mmseqs_DB | 6.1 (5.9–6.2) | 104 |
| db_warm_diamond_very-sensitive_mmseqs_DB | 7.2 (7.0–7.7) | 109 |
| db_cold_saber_v1.2_mmseqs_DB | 6.2 (5.9–6.5) | 104 |
| db_cold_diamond_very-sensitive_mmseqs_DB | 13.8 (7.6–13.8) | 108 |
| makedb_saber_v1.2_mmseqs_DB_p16M | 1.9 | 37 |
| db_order_query_mmseqs_DB | 5.8 | 105 |
| db_order_partition_mmseqs_DB | 5.7 | 102 |
| makedb_saber_v1.2_db400 | 7.9 (7.5–7.9) | 151 |
| makedb_diamond_db400 | 1.6 (1.5–2.0) | 58 |
| makedb_saber_v1.1.0_db400 | 2.3 | 210 |
| db_warm_saber_v1.2_db400 | 20.1 (20.1–20.2) | 248 |
| db_warm_diamond_very-sensitive_db400 | 24.3 (23.8–24.4) | 364 |
| db_cold_saber_v1.2_db400 | 20.9 (20.6–21.0) | 248 |
| db_cold_diamond_very-sensitive_db400 | 23.6 (23.2–23.9) | 362 |
| makedb_saber_v1.2_db400_p16M | 7.6 | 55 |
| db_order_query_db400 | 19.7 | 249 |
| db_order_partition_db400 | 18.6 | 249 |

File sizes / checks:
```
-rw-r--r-- 1 root root 11755080 Oct 18 17:09 /home/claude/bench/final/index/mmseqs_DB.dmnd
-rw-r--r-- 1 root root 43091108 Oct 18 20:33 /home/claude/bench/final/index/mmseqs_DB.sdx
-rw-r--r-- 1 root root 25671747 Oct 18 20:33 /home/claude/bench/final/index/mmseqs_DB.v110.sdx
mmseqs_DB: query-major == partition-major output
-rw-r--r-- 1 root root  39883840 Oct 18 17:12 /home/claude/bench/final/index/db400.dmnd
-rw-r--r-- 1 root root 175243442 Oct 18 20:34 /home/claude/bench/final/index/db400.sdx
-rw-r--r-- 1 root root 156758132 Oct 18 20:34 /home/claude/bench/final/index/db400.v110.sdx
db400: query-major == partition-major output
```

## Translated searches

Measured on: host A2 (Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB)).

blastx (`nx`): back-translated SCOP40 domains (random synonymous codons, random flanks and strand) against SCOP40 proteins; tblastn (`pn`): a fixed 1-in-37 subset of SCOP40 proteins (303 queries) against the back-translated set. `nxsub_*` / `pnsub_*`: a fixed 1-in-333 query subset (34 queries), used to compare with SABER 1.1.0, which searched translated data by full dynamic programming (its sensitivity numbers on 34 queries are indicative only). Metrics are over the queries of each run (both splits): superfamily-level sensitivity to the first FP, pooled recall at 0.01 false positives per query, and false positives at E ≤ 1e-3.

| Run | Queries | Time s | RSS MB | Sens. to 1st FP (SF) | Recall @ EPQ 0.01 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|
| nx_diamond_default | 11206 | 3.8 | 56 | 0.045 | 0.0123 | 4 |
| nx_diamond_ultra-sensitive | 11206 | 54.5 | 148 | 0.213 | 0.0718 | 41 |
| nx_diamond_very-sensitive | 11206 | 14.5 | 107 | 0.185 | 0.0641 | 33 |
| nx_saber_v1.2_default | 11206 | 25.1 | 33 | 0.189 | 0.0609 | 25 |
| nx_saber_v1.2_fast | 11206 | 12.9 | 30 | 0.150 | 0.0485 | 19 |
| nx_saber_v1.2_more-sensitive | 11206 | 52.6 | 41 | 0.210 | 0.0659 | 33 |
| nxsub_diamond_very-sensitive | 34 | 1.0 | 29 | 0.163 | 0.0713 | 0 |
| nxsub_saber_v1.1.0 | 34 | 332.3 | 12 | 0.228 | 0.0771 | 2 |
| nxsub_saber_v1.2_default | 34 | 0.4 | 18 | 0.137 | 0.0684 | 0 |
| pn_saber_v1.2_default | 303 | 2.5 | 62 | 0.194 | 0.0579 | 0 |
| pn_saber_v1.2_more-sensitive | 303 | 3.5 | 106 | 0.215 | 0.0558 | 2 |
| pnsub_saber_v1.1.0 | 34 | 330.4 | 32 | 0.235 | 0.0765 | 0 |
| pnsub_saber_v1.2_default | 34 | 2.0 | 60 | 0.150 | 0.0661 | 0 |

### Translated subsets at matched sensitivity (follow-up run)

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

The 34-query blastx / tblastn subsets above, searched again in one session on one machine after the main run (`bench/scripts/translated_matched.sh`): SABER 1.1.0 (full dynamic programming over six frames), SABER 1.2 `default`, `more-sensitive` and `exhaustive`, and DIAMOND (blastx only), so speed can be compared at similar sensitivity. Metrics as above (34 queries, indicative only); SABER 1.1.0 single runs, the others median of 3.

| Run | Time s | RSS MB | Sens. to 1st FP (SF) | Recall @ EPQ 0.01 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|
| nxsub_saber_v1.1.0 | 339.2 | 12 | 0.228 | 0.0771 | 2 |
| nxsub_saber_v1.2_default | 0.3 (0.3–0.3) | 18 | 0.137 | 0.0684 | 0 |
| nxsub_saber_v1.2_more-sensitive | 0.4 (0.4–0.5) | 27 | 0.172 | 0.0759 | 0 |
| nxsub_saber_v1.2_exhaustive | 8.7 (8.7–9.5) | 23 | 0.224 | 0.0783 | 0 |
| nxsub_diamond_very-sensitive | 0.9 (0.9–1.0) | 30 | 0.163 | 0.0713 | 0 |
| nxsub_diamond_ultra-sensitive | 3.0 (2.9–3.3) | 47 | 0.201 | 0.0788 | 0 |
| pnsub_saber_v1.1.0 | 329.7 | 33 | 0.235 | 0.0765 | 0 |
| pnsub_saber_v1.2_default | 1.6 (1.5–1.6) | 60 | 0.150 | 0.0661 | 0 |
| pnsub_saber_v1.2_more-sensitive | 2.0 (1.9–2.0) | 102 | 0.181 | 0.0730 | 0 |
| pnsub_saber_v1.2_exhaustive | 9.9 (9.8–10.1) | 95 | 0.234 | 0.0788 | 0 |

## Re-check of the A1 SABER-vs-DIAMOND timings (follow-up run)

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

In the A1 sections the DIAMOND rows were measured between 14:34 and 17:17 and the SABER 1.2 rows between 20:10 and 20:37. These comparisons were re-measured with the same commands, both tools in one session on one machine of A1's CPU model (`bench/scripts/recheck.sh`, SABER = the delivered build), next to the medians of the original sections. Absolute times differ between the sessions (another virtual machine instance or host load); compare the two tools within a column. The delivered build's output equals the benchmarked build's except on the long-protein set, where its traceback fallback changes 9 of 4,712 alignments (all to higher scores).

| Section: run | Time s (re-check) | RSS MB | Time s (original section) | RSS MB |
|---|---|---|---|---|
| scop: diamond_very-sensitive | 9.6 (9.3–10.2) | 61 | 13.0 (12.2–13.3) | 58 |
| scop: saber_v1.2_default | 8.0 (7.8–8.3) | 30 | 10.0 (9.9–10.4) | 30 |
| mmseqs: diamond_very-sensitive | 5.9 (5.8–6.1) | 104 | 7.1 (7.0–7.7) | 113 |
| mmseqs: saber_v1.2_default | 6.1 (6.1–6.3) | 106 | 7.4 (7.3–8.1) | 106 |
| mmseqs: saber_v1.2_default_db | 5.0 (5.0–5.1) | 104 | 6.0 (5.9–6.1) | 104 |
| index: db_warm_diamond_very-sensitive_db400 | 17.7 (17.6–18.8) | 363 | 24.3 (23.8–24.4) | 364 |
| index: db_warm_saber_v1.2_db400 | 17.8 (17.1–18.2) | 248 | 20.1 (20.1–20.2) | 248 |
| mmseqs: long_diamond_very-sensitive | 6.8 (6.7–7.4) | 114 | 8.1 | 105 |
| mmseqs: long_saber_v1.2_default | 7.2 (7.0–8.3) | 110 | 8.4 | 115 |

Outputs: scop/diamond_very-sensitive: output identical; scop/saber_v1.2_default: output identical; mmseqs/diamond_very-sensitive: output identical; mmseqs/saber_v1.2_default: output identical; mmseqs/saber_v1.2_default_db: output identical; mmseqs/long_diamond_very-sensitive: output identical; mmseqs/long_saber_v1.2_default: output DIFFERENT.

## Before / after the performance changes (follow-up run)

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

A development build from before the index-construction, sparse-offset, 16-bit traceback, filter-first and X-drop changes (`_old`, sha256 63442e7d…, full hash in `env/final_binary.txt`) against the benchmarked build (`_new`) on one machine in one session; `--stats` stage times are summed over threads.

| Run | Time s | RSS MB | traceback stage s | seed stage s |
|---|---|---|---|---|
| scop_very-fast_old | 3.7 | 102 | 1.20 | 0.55 |
| scop_default_old | 10.8 | 32 | 2.52 | 2.50 |
| scop_more-sensitive_old | 20.3 | 48 | 3.18 | 5.49 |
| mm_default_old | 7.5 | 157 | 2.84 | 1.64 |
| mm_long_old | 9.5 | 221 | 5.29 | 1.62 |
| nx_default_old | 37.6 | 35 | 3.69 | 10.35 |
| scop_very-fast_new | 2.7 | 28 | 0.58 | 0.24 |
| scop_default_new | 7.5 | 29 | 1.17 | 0.28 |
| scop_more-sensitive_new | 16.7 | 36 | 1.50 | 2.63 |
| mm_default_new | 6.3 | 107 | 1.41 | 1.42 |
| mm_long_new | 7.2 | 117 | 2.61 | 1.25 |
| nx_default_new | 24.9 | 33 | 1.66 | 0.87 |

Outputs of the development build compared with the released outputs (first = development build):

* SCOP40 very-fast: 36311 vs 36311 lines; pairs only in first 0, only in second 0; 6 common pairs differ (6 with the same E-value and bit score, 0 higher / 0 lower bit score in the second); columns that change: pident 6, length 6, mismatch 6, qstart 3, sstart 3, gapopen 3, qend 1, send 1
* SCOP40 default: 80652 vs 80652 lines; pairs only in first 0, only in second 0; 6 common pairs differ (6 with the same E-value and bit score, 0 higher / 0 lower bit score in the second); columns that change: pident 6, length 6, mismatch 6, qstart 4, sstart 4, gapopen 3, qend 2, send 2
* SCOP40 more-sensitive: 106215 vs 106215 lines; pairs only in first 0, only in second 0; 6 common pairs differ (6 with the same E-value and bit score, 0 higher / 0 lower bit score in the second); columns that change: pident 6, length 6, mismatch 6, qstart 4, sstart 4, gapopen 3, qend 2, send 2
* MMseqs2 example default: identical (15607 lines)
* blastx default: 115378 vs 115378 lines; pairs only in first 0, only in second 0; 7 common pairs differ (7 with the same E-value and bit score, 0 higher / 0 lower bit score in the second); columns that change: pident 7, length 7, mismatch 7, gapopen 4, qstart 4, sstart 4, qend 3, send 2
* long proteins, benchmarked build: 4712 vs 4712 lines; pairs only in first 0, only in second 0; 11 common pairs differ (2 with the same E-value and bit score, 0 higher / 9 lower bit score in the second); columns that change: pident 11, length 11, mismatch 11, qstart 11, qend 11, sstart 11, gapopen 9, send 9, evalue 9, bitscore 9
* long proteins, delivered build: 4712 vs 4712 lines; pairs only in first 0, only in second 0; 2 common pairs differ (2 with the same E-value and bit score, 0 higher / 0 lower bit score in the second); columns that change: pident 2, length 2, mismatch 2, qstart 2, qend 2, sstart 2, send 2, gapopen 1
* SCOP40 default, NEW re-run: identical (80652 lines)
* blastx default, NEW re-run: identical (115378 lines)

## Adversarial low-complexity / repeat inputs

Measured on: host A2 (Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB)).

100 homopolymer / tandem-repeat queries (some with an embedded real domain) against the MMseqs database plus 2,000 low-complexity sequences. `mmq_*`: the 500 real MMseqs queries against the same polluted database — hits to the artificial `lowcdb*` sequences are false positives. SABER 1.0.2/1.1.0 ran under a 6 GB address-space limit; the adversarial `lowc10` runs and the 1.x runs under a 10-minute time limit (exit 124 = stopped by the time limit; exit −6 = aborted after a failed memory allocation, e.g. `memory allocation of 4294967296 bytes failed`).

| Run | Time s | RSS MB | Exit | Hits (E ≤ 1e-3) | Hits to artificial low-complexity sequences |
|---|---|---|---|---|---|
| lowc_saber_v1.2_default | 1.5 | 62 | 0 | 72 | 0 |
| lowc_saber_v1.2_default_nomask | 27.1 | 123 | 0 | 10448 | 5647 |
| lowc10_saber_v1.2_legacy_nomask_nocap | 11.2 | 103 | 0 | 2751 | 1786 |
| lowc10_saber_v1.1.0 | 116.0 | 3791 | -6 | – (aborted) | – |
| lowc_diamond_default | 1.4 | 44 | 0 | 21 | 0 |
| lowc_diamond_very-sensitive | 3.3 | 109 | 0 | 71 | 0 |
| lowc_saber_v1.1.0 | 21.2 | 4205 | -6 | – (aborted) | – |
| lowc_saber_v1.0.2 | 395.1 | 3674 | -6 | – (aborted) | – |
| mmq_saber_v1.2_default | 6.6 | 109 | 0 | 15526 | 0 |
| mmq_diamond_very-sensitive | 6.8 | 105 | 0 | 15368 | 0 |
| mmq_saber_v1.1.0 | 31.1 | 224 | 0 | 18600 | 3883 |

## Two-domain queries and short fragments

Measured on: host B (Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB)).

Retention = fraction of a domain's true positives from its single-domain search (each tool's own, E ≤ 1e-3) that are still found when the domain is part of a two-domain query; the TP columns give the absolute numbers found in the two-domain searches. Single runs.

| Tool / preset | Retention of the weaker domain's TPs | Retention (stronger) | Weaker-domain TPs (two-domain / alone) | Stronger-domain TPs (two-domain / alone) | FP in two-domain searches (E ≤ 1e-3) | Fragments: sens. to 1st FP (SF, test) | Fragment search time s |
|---|---|---|---|---|---|---|---|
| saber_v1.1.0 | 0.948 | 0.962 | 184 / 194 | 1499 / 1558 | 7 | 0.051 | 2.6 |
| saber_v1.2_default | 0.963 | 0.980 | 183 / 190 | 1671 / 1705 | 0 | 0.061 | 0.8 |
| saber_v1.2_more-sensitive | 0.959 | 0.973 | 210 / 219 | 1783 / 1832 | 0 | 0.075 | 1.7 |
| diamond_very-sensitive | 0.938 | 0.947 | 212 / 226 | 1749 / 1846 | 0 | 0.043 | 1.5 |
| diamond_ultra-sensitive | 0.925 | 0.941 | 234 / 253 | 1937 / 2058 | 0 | 0.047 | 5.7 |

## Thread scaling (2-core machine)

Measured on: host B (Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB)).

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| t1_saber_v1.2_default | 11.4 (11.2–11.6) | 11.4 | 24 |
| t1_diamond_very-sensitive | 14.2 (14.2–14.3) | 14.2 | 47 |
| t1_saber_v1.1.0 | 38.8 | 38.4 | 28 |
| t1_saber_v1.2_intraquery | 11.2 (11.0–13.4) | 11.0 | 203 |
| t2_saber_v1.2_default | 6.0 (5.9–6.0) | 11.2 | 30 |
| t2_diamond_very-sensitive | 7.2 (7.2–7.2) | 14.0 | 55 |
| t2_saber_v1.1.0 | 21.1 | 38.7 | 31 |
| t2_saber_v1.2_intraquery | 7.1 (7.1–7.3) | 11.4 | 208 |

## Preset selection sweep (SCOP40 **tune** split only)

Measured on: host B (Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB)).

Single runs, 2 threads. These numbers chose the presets; the test-split results above are the held-out evaluation.

| Configuration | Time s | Sens. to 1st FP (SF, tune) | Sens. (family) | Recall @ EPQ 0.01 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|
| full_1111_twohit40 | 1.8 | 0.044 | 0.096 | 0.0110 | 0 |
| full_11111_hf0 | 1.8 | 0.089 | 0.188 | 0.0238 | 10 |
| m10_11101111_hf12 | 1.9 | 0.111 | 0.231 | 0.0290 | 5 |
| full_1111_hf10 | 1.9 | 0.118 | 0.249 | 0.0307 | 15 |
| full_110111_hf0 | 1.9 | 0.100 | 0.210 | 0.0262 | 12 |
| m10_1111011_hf13 | 2.3 | 0.144 | 0.300 | 0.0370 | 8 |
| full_1111_hf8 | 2.6 | 0.137 | 0.288 | 0.0362 | 15 |
| m10_1111011_hf12 | 2.8 | 0.150 | 0.312 | 0.0384 | 8 |
| m10_111111_hf11 | 3.0 | 0.140 | 0.289 | 0.0362 | 12 |
| m10_1111011_hf11 | 3.5 | 0.153 | 0.317 | 0.0392 | 10 |
| m10_1111011+1101111_hf12 | 4.1 | 0.168 | 0.345 | 0.0427 | 10 |
| m10_110111_hf13 | 4.6 | 0.178 | 0.360 | 0.0461 | 6 |
| m10_11111_hf12 | 5.4 | 0.174 | 0.358 | 0.0453 | 8 |
| full_1111_hf0 | 5.5 | 0.144 | 0.302 | 0.0379 | 15 |
| m10_110111_hf12 | 5.9 | 0.187 | 0.377 | 0.0477 | 10 |
| m10_110111_hf11 | 8.2 | 0.192 | 0.386 | 0.0484 | 14 |
| m10_1101011+111011_hf12 | 9.9 | 0.200 | 0.404 | 0.0505 | 13 |
| m10_1101011+111011_hf11 | 13.2 | 0.206 | 0.412 | 0.0510 | 17 |
| full_1111_nbT13_twohit | 53.7 | 0.198 | 0.404 | 0.0515 | 11 |

