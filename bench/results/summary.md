## At scale: SCOP40 against SCOP40 + 80,000 reversed proteins (91,206 sequences, 38.2 M residues)

**Measured on more than one machine — compare rows only within one machine** (marked after the name): A = Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB); B = Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

11,206 SCOP40 queries, 2 threads, E ≤ 10, no hit limit; every hit to a reversed protein (`REV_*`) or to a domain of another fold is a false positive. Prebuilt indexes for both tools (`saber --makedb`, `diamond makedb`; built once, timed below). SABER 1.4 `default` and DIAMOND `--very-sensitive` were run alternately (5 runs each) so that machine drift affects both alike; the other rows are median of 3 (single runs: `diamond_ultra-sensitive`, `saber_v1.4_default_dbmodel`, `saber_v1.3_default`). Accuracy columns are on the **held-out test split** (tune split in brackets): sensitivity = fraction of a query's same-superfamily domains ranked above its first false positive, averaged over queries; recall @ EPQ x = pooled fraction of same-superfamily pairs found before x false positives per query; TP / FP = same-superfamily / different-fold (or decoy) pairs reported at E ≤ 1e-3. Time = wall seconds, median (range) of the runs; RSS = peak resident memory (median).

| Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP (superfamily) | Recall @ EPQ 0.01 | Recall @ EPQ 0.1 | Recall @ FDP 5 % | TP at E ≤ 1e-3 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|---|---|---|
| saber_v1.4_very-fast [A] | 6.2 (6.1–6.4) | 11.5 | 111 | 0.1460 [0.1419] | 0.0639 [0.0373] | 0.0693 | 0.0697 | 10008 | 1 |
| saber_v1.4_fast [A] | 12.1 (11.9–12.6) | 22.6 | 161 | 0.1713 [0.1597] | 0.0741 [0.0422] | 0.0811 | 0.0826 | 11367 | 1 |
| saber_v1.4_default [A] | 24.2 (23.9–24.9) | 45.0 | 280 | 0.1929 [0.1770] | 0.0826 [0.0467] | 0.0917 | 0.0941 | 12453 | 1 |
| saber_v1.4_sensitive [A] | 32.0 (31.9–32.3) | 60.1 | 246 | 0.2052 [0.1930] | 0.0874 [0.0482] | 0.0978 | 0.1009 | 13023 | 2 |
| saber_v1.4_more-sensitive [A] | 55.7 (55.6–56.3) | 103.4 | 420 | 0.2143 [0.1983] | 0.0898 [0.0495] | 0.1011 | 0.1049 | 13306 | 2 |
| diamond_fast [A] | 4.5 (4.3–4.6) | 8.7 | 143 | 0.0125 [0.0122] | 0.0036 [0.0026] | 0.0036 | 0.0036 | 664 | 0 |
| diamond_default [A] | 6.7 (6.6–6.8) | 12.8 | 162 | 0.0448 [0.0445] | 0.0152 [0.0099] | 0.0154 | 0.0154 | 2693 | 0 |
| diamond_mid-sensitive [A] | 19.7 (19.5–19.8) | 37.6 | 154 | 0.0978 [0.0965] | 0.0387 [0.0234] | 0.0406 | 0.0404 | 6799 | 1 |
| diamond_sensitive [A] | 37.4 (36.9–40.8) | 72.1 | 163 | 0.1421 [0.1398] | 0.0626 [0.0352] | 0.0674 | 0.0677 | 10939 | 5 |
| diamond_more-sensitive [A] | 38.4 (37.8–38.8) | 73.1 | 155 | 0.1424 [0.1407] | 0.0630 [0.0363] | 0.0679 | 0.0682 | 11015 | 5 |
| diamond_very-sensitive [A] | 30.3 (29.6–30.8) | 58.3 | 385 | 0.1832 [0.1729] | 0.0801 [0.0459] | 0.0884 | 0.0905 | 13720 | 5 |
| diamond_ultra-sensitive [A] | 116.8 | 227.0 | 903 | 0.2075 [0.1916] | 0.0882 [0.0501] | 0.0997 | 0.1032 | 15013 | 5 |
| diamond_very-sensitive_c1 [B] | 27.9 (27.3–30.3) | 53.1 | 385 | 0.1832 [0.1729] | 0.0801 [0.0459] | 0.0884 | 0.0905 | 13720 | 5 |
| saber_v1.4_default_dbmodel [B] | 25.8 | 47.9 | 281 | 0.1867 [0.1723] | 0.0769 [0.0438] | 0.0859 | 0.0879 | 12219 | 1 |
| saber_v1.3_default [A] | 24.4 | 45.7 | 304 | 0.1929 [0.1770] | 0.0826 [0.0467] | 0.0917 | 0.0941 | 12453 | 1 |

`diamond_very-sensitive_c1` = DIAMOND `--very-sensitive -c1` (one seed-index chunk, DIAMOND's speed/memory option); `saber_v1.4_default_dbmodel` = SABER 1.4 `default` with `--evalue-model database` (database-wide E-values, as in 1.2: same alignments and bit scores, different E-values and therefore a different ranking across subjects). Both are run by `run_all.sh scale_extra`.

### False-discovery proportion at matched recall (test split, at scale)

| Tool / preset | FDP @ recall 0.06 | FDP @ recall 0.07 | FDP @ recall 0.08 | FDP @ recall 0.09 | FDP @ recall 0.1 |
|---|---|---|---|---|---|
| saber_v1.4_very-fast | 0.0006 | 0.0563 | not reached | not reached | not reached |
| saber_v1.4_fast | 0.0001 | 0.0007 | 0.0245 | not reached | not reached |
| saber_v1.4_default | 0.0000 | 0.0002 | 0.0018 | 0.0216 | 0.1467 |
| saber_v1.4_sensitive | 0.0000 | 0.0001 | 0.0005 | 0.0062 | 0.0422 |
| saber_v1.4_more-sensitive | 0.0000 | 0.0001 | 0.0005 | 0.0034 | 0.0233 |
| diamond_fast | not reached | not reached | not reached | not reached | not reached |
| diamond_default | not reached | not reached | not reached | not reached | not reached |
| diamond_mid-sensitive | not reached | not reached | not reached | not reached | not reached |
| diamond_sensitive | 0.0011 | not reached | not reached | not reached | not reached |
| diamond_more-sensitive | 0.0009 | not reached | not reached | not reached | not reached |
| diamond_very-sensitive | 0.0000 | 0.0003 | 0.0034 | 0.0447 | not reached |
| diamond_ultra-sensitive | 0.0000 | 0.0000 | 0.0005 | 0.0049 | 0.0302 |
| diamond_very-sensitive_c1 | 0.0000 | 0.0003 | 0.0034 | 0.0447 | not reached |
| saber_v1.4_default_dbmodel | 0.0001 | 0.0007 | 0.0082 | 0.0780 | 0.3813 |
| saber_v1.3_default | 0.0000 | 0.0002 | 0.0018 | 0.0216 | 0.1467 |

Alternating runs, SABER 1.4 `default` / DIAMOND `--very-sensitive` wall time: 24.9 / 30.4 s, 24.2 / 30.3 s, 23.9 / 29.6 s, 24.1 / 30.8 s, 24.3 / 30.2 s (ratio 0.78–0.82); CPU time 45 / 58 s, 45 / 59 s, 45 / 57 s, 45 / 60 s, 46 / 58 s.

### Index construction (scale database)

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| makedb_saber_v1.4 | 3.6 (3.5–3.9) | 6.6 | 73 |
| makedb_diamond | 1.2 (1.1–1.2) | 2.1 | 59 |
| makedb_saber_v1.3 | 3.8 | 6.8 | 73 |

```
-rw-r--r-- 1 root root 40690880 Oct 19 20:54 /home/claude/bench/final14b/scale/scale.dmnd
-rw-r--r-- 1 root root 43983827 Oct 19 20:54 /home/claude/bench/final14b/scale/scale.sdx
-rw-r--r-- 1 root root 43983827 Oct 19 20:54 /home/claude/bench/final14b/scale/scale.vo.sdx
```

### Practical settings at scale (E ≤ 1e-3, 25 targets per query)

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| practical_saber_v1.4_default | 24.5 (24.0–24.5) | 45.5 | 280 |
| practical_diamond_very-sensitive | 30.6 (29.5–30.6) | 59.1 | 385 |

### DIAMOND's seed-search algorithm at scale

Measured on machine B: Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

DIAMOND chooses between its double-indexed and query-indexed seed search automatically; here `--very-sensitive --algo 1` (query-indexed) is forced, run alternately with SABER 1.4 `default` (whose accuracy is in the table above; same output). Accuracy columns are on the **held-out test split** (tune split in brackets): sensitivity = fraction of a query's same-superfamily domains ranked above its first false positive, averaged over queries; recall @ EPQ x = pooled fraction of same-superfamily pairs found before x false positives per query; TP / FP = same-superfamily / different-fold (or decoy) pairs reported at E ≤ 1e-3. Time = wall seconds, median (range) of the runs; RSS = peak resident memory (median).

| Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP (superfamily) | Recall @ EPQ 0.01 | Recall @ EPQ 0.1 | Recall @ FDP 5 % | TP at E ≤ 1e-3 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|---|---|---|
| saber_v1.4_default | 25.4 (25.1–25.9) | 47.1 | 280 | – [–] | – [–] | – | – | – | – |
| diamond_very-sensitive_algo1 | 59.9 (56.8–60.6) | 114.9 | 430 | 0.2076 [0.1913] | 0.0879 [0.0499] | 0.0992 | 0.1028 | 15007 | 5 |

Alternating runs, SABER / DIAMOND `--algo 1` wall time: 25.9 / 59.9 s, 25.1 / 60.6 s, 25.4 / 56.8 s.

### One-off search at scale, index construction included

Measured on machine B: Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

No preparation step: SABER 1.4 `default` searches the FASTA database directly (`-d`, index built in memory, included in its time); DIAMOND runs `makedb` and then `blastp --very-sensitive` (one timed command). Alternating runs.

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| saber_v1.4_default_fasta | 28.5 (27.8–30.8) | 51.8 | 310 |
| diamond_very-sensitive_makedb_blastp | 27.8 (27.3–27.9) | 53.5 | 385 |

Alternating runs, SABER / DIAMOND wall time: 30.8 / 27.3 s, 27.8 / 27.7 s, 28.0 / 27.8 s, 28.5 / 27.9 s, 29.5 / 27.9 s.

Output check: in-memory index output == prebuilt index output.

## SCOP40 all-against-all (11,206 domains, 2 threads, E ≤ 10, no hit limit)

Measured on machine A: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

SABER builds its index in memory on every run (`-d`; `default_db` = prebuilt `--db`), DIAMOND uses a prebuilt `.dmnd`. Accuracy columns are on the **held-out test split** (tune split in brackets): sensitivity = fraction of a query's same-superfamily domains ranked above its first false positive, averaged over queries; recall @ EPQ x = pooled fraction of same-superfamily pairs found before x false positives per query; TP / FP = same-superfamily / different-fold (or decoy) pairs reported at E ≤ 1e-3. Time = wall seconds, median (range) of the runs; RSS = peak resident memory (median).

| Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP (superfamily) | Recall @ EPQ 0.01 | Recall @ EPQ 0.1 | Recall @ FDP 5 % | TP at E ≤ 1e-3 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|---|---|---|
| saber_v1.4_very-fast | 2.4 (2.4–2.6) | 4.1 | 59 | 0.1611 [0.1563] | 0.0700 [0.0418] | 0.0780 | 0.0792 | 11473 | 3 |
| saber_v1.4_fast | 3.8 (3.6–5.2) | 6.2 | 94 | 0.1895 [0.1755] | 0.0802 [0.0470] | 0.0921 | 0.0944 | 13098 | 3 |
| saber_v1.4_default | 6.1 (6.1–6.4) | 10.7 | 137 | 0.2105 [0.1925] | 0.0878 [0.0511] | 0.1023 | 0.1058 | 14242 | 4 |
| saber_v1.4_sensitive | 5.4 (5.3–5.4) | 9.8 | 88 | 0.2229 [0.2076] | 0.0929 [0.0536] | 0.1089 | 0.1135 | 14909 | 5 |
| saber_v1.4_more-sensitive | 10.3 (10.2–10.3) | 18.2 | 163 | 0.2340 [0.2150] | 0.0958 [0.0551] | 0.1137 | 0.1189 | 15308 | 5 |
| saber_v1.4_legacy | 6.1 (6.0–6.1) | 11.6 | 37 | 0.1504 [0.1433] | 0.0640 [0.0383] | 0.0701 | 0.0707 | 10587 | 4 |
| saber_v1.4_exhaustive | 153.7 | 304.2 | 35 | 0.2724 [0.2514] | 0.1023 [0.0581] | 0.1248 | 0.1345 | 16145 | 6 |
| diamond_fast | 1.2 (1.2–1.2) | 2.3 | 27 | 0.0125 [0.0123] | 0.0036 [0.0026] | 0.0036 | 0.0036 | 673 | 0 |
| diamond_default | 1.9 (1.9–2.0) | 3.7 | 29 | 0.0451 [0.0451] | 0.0155 [0.0101] | 0.0157 | 0.0157 | 2775 | 1 |
| diamond_mid-sensitive | 5.2 (5.1–5.2) | 10.0 | 34 | 0.0989 [0.0978] | 0.0402 [0.0243] | 0.0421 | 0.0419 | 7148 | 6 |
| diamond_sensitive | 9.5 (9.4–9.5) | 18.2 | 38 | 0.1429 [0.1407] | 0.0656 [0.0374] | 0.0690 | 0.0690 | 11857 | 14 |
| diamond_more-sensitive | 9.7 (9.7–10.0) | 18.9 | 39 | 0.1432 [0.1415] | 0.0661 [0.0386] | 0.0694 | 0.0694 | 11940 | 14 |
| diamond_very-sensitive | 8.4 (8.4–8.5) | 16.5 | 58 | 0.1880 [0.1763] | 0.0851 [0.0498] | 0.0960 | 0.0973 | 15433 | 17 |
| diamond_ultra-sensitive | 29.5 (28.7–29.9) | 57.7 | 114 | 0.2182 [0.1996] | 0.0959 [0.0552] | 0.1117 | 0.1150 | 17395 | 22 |
| saber_v1.3_default | 5.7 (5.7–6.0) | 10.2 | 135 | 0.2105 [0.1925] | 0.0878 [0.0511] | 0.1023 | 0.1058 | 14242 | 4 |
| saber_v1.3_more-sensitive | 10.2 (10.2–10.3) | 18.1 | 160 | 0.2340 [0.2150] | 0.0958 [0.0551] | 0.1137 | 0.1189 | 15308 | 5 |
| saber_v1.4_default_db | 5.6 (5.6–5.8) | 10.0 | 133 | 0.2105 [0.1925] | 0.0878 [0.0511] | 0.1023 | 0.1058 | 14242 | 4 |

### False-discovery proportion at matched recall (test split)

Pooled over all test-split queries: hits sorted by E-value, FDP = FP / (TP + FP) at the most stringent cut-off that reaches the recall. 'not reached' = the tool never reaches that recall at E ≤ 10.

| Tool / preset | FDP @ recall 0.06 | FDP @ recall 0.08 | FDP @ recall 0.09 | FDP @ recall 0.1 | FDP @ recall 0.11 |
|---|---|---|---|---|---|
| saber_v1.4_very-fast | 0.0003 | 0.0632 | not reached | not reached | not reached |
| saber_v1.4_fast | 0.0001 | 0.0034 | 0.0200 | 0.1487 | not reached |
| saber_v1.4_default | 0.0000 | 0.0007 | 0.0049 | 0.0187 | 0.0915 |
| saber_v1.4_sensitive | 0.0000 | 0.0004 | 0.0015 | 0.0086 | 0.0304 |
| saber_v1.4_more-sensitive | 0.0000 | 0.0003 | 0.0009 | 0.0056 | 0.0165 |
| saber_v1.4_legacy | 0.0009 | not reached | not reached | not reached | not reached |
| saber_v1.4_exhaustive | 0.0000 | 0.0003 | 0.0009 | 0.0019 | 0.0069 |
| diamond_fast | not reached | not reached | not reached | not reached | not reached |
| diamond_default | not reached | not reached | not reached | not reached | not reached |
| diamond_mid-sensitive | not reached | not reached | not reached | not reached | not reached |
| diamond_sensitive | 0.0005 | not reached | not reached | not reached | not reached |
| diamond_more-sensitive | 0.0004 | not reached | not reached | not reached | not reached |
| diamond_very-sensitive | 0.0000 | 0.0007 | 0.0073 | not reached | not reached |
| diamond_ultra-sensitive | 0.0000 | 0.0003 | 0.0009 | 0.0048 | 0.0187 |
| saber_v1.3_default | 0.0000 | 0.0007 | 0.0049 | 0.0187 | 0.0915 |
| saber_v1.3_more-sensitive | 0.0000 | 0.0003 | 0.0009 | 0.0056 | 0.0165 |
| saber_v1.4_default_db | 0.0000 | 0.0007 | 0.0049 | 0.0187 | 0.0915 |

### Practical settings (E ≤ 1e-3, 25 targets per query), SCOP40 all-vs-all

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| practical_saber_v1.4_fast | 3.3 (3.2–3.4) | 5.9 | 92 |
| practical_saber_v1.4_default | 5.6 (5.5–5.6) | 9.9 | 136 |
| practical_diamond_default | 1.9 (1.9–2.0) | 3.7 | 30 |
| practical_diamond_very-sensitive | 8.2 (8.2–8.4) | 16.1 | 61 |
| practical_diamond_ultra-sensitive | 28.9 (28.1–29.4) | 56.6 | 114 |

### SCOP40 index build and persistent-index search

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| diamond_makedb | 0.1 (0.1–0.1) | 0.1 | 19 |
| saber_v1.4_makedb_default | 0.4 (0.4–0.5) | 0.6 | 28 |
| saber_v1.4_default_db_warm | 5.6 (5.6–5.8) | 10.0 | 133 |

## Null behaviour: real SCOP40 queries against shuffled databases

Measured on machine B: Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

Real queries against composition-preserving decoys (`shuf`: whole-sequence shuffle; `wshuf10`: shuffle inside 10-residue windows). Hits to decoys derived from a sequence of the query's own fold are excluded; every remaining hit is a false positive. A calibrated E-value gives at most about E false hits per query (EPQ ≤ E). `dbmodel` = SABER 1.4 with `--evalue-model database` (database-wide E-values, as in 1.2).

| Run | EPQ at E ≤ 1e-3 | E ≤ 0.01 | E ≤ 0.1 | E ≤ 1 | E ≤ 10 |
|---|---|---|---|---|---|
| diamond_default_shuf | 0.0000 | 0.0001 | 0.0008 | 0.0029 | 0.0103 |
| diamond_default_wshuf10 | 0.0000 | 0.0001 | 0.0010 | 0.0034 | 0.0106 |
| diamond_ultra-sensitive_shuf | 0.0012 | 0.0064 | 0.0242 | 0.0920 | 0.2509 |
| diamond_ultra-sensitive_wshuf10 | 0.0007 | 0.0073 | 0.0327 | 0.1116 | 0.2996 |
| diamond_very-sensitive_shuf | 0.0006 | 0.0031 | 0.0126 | 0.0431 | 0.1070 |
| diamond_very-sensitive_wshuf10 | 0.0004 | 0.0046 | 0.0180 | 0.0535 | 0.1258 |
| saber_v1.3_default_shuf | 0.0001 | 0.0006 | 0.0043 | 0.0637 | 0.7152 |
| saber_v1.3_default_wshuf10 | 0.0002 | 0.0006 | 0.0073 | 0.0763 | 0.8099 |
| saber_v1.4_default_dbmodel_shuf | 0.0005 | 0.0020 | 0.0202 | 0.2311 | 2.0698 |
| saber_v1.4_default_dbmodel_wshuf10 | 0.0004 | 0.0030 | 0.0229 | 0.2680 | 2.3158 |
| saber_v1.4_default_shuf | 0.0001 | 0.0006 | 0.0043 | 0.0637 | 0.7152 |
| saber_v1.4_default_wshuf10 | 0.0002 | 0.0006 | 0.0073 | 0.0763 | 0.8099 |
| saber_v1.4_fast_shuf | 0.0000 | 0.0004 | 0.0032 | 0.0451 | 0.5033 |
| saber_v1.4_fast_wshuf10 | 0.0001 | 0.0004 | 0.0053 | 0.0533 | 0.5671 |
| saber_v1.4_more-sensitive_shuf | 0.0001 | 0.0009 | 0.0077 | 0.0949 | 1.0573 |
| saber_v1.4_more-sensitive_wshuf10 | 0.0002 | 0.0006 | 0.0099 | 0.1110 | 1.1934 |

## MMseqs2 example workload (500 queries × 20,000 UniProt sequences, E ≤ 1e-3)

Measured on machine B: Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

Recall is measured against SABER 1.4's exhaustive Smith–Waterman search (`--preset exhaustive`: 16344 query–subject pairs at E ≤ 1e-3). That reference shares SABER's scoring, composition correction and E-value statistics, so near the cut-off it favours SABER; 'Recall, E ≤ 1e-10 pairs' uses only the 12600 reference pairs with E ≤ 1e-10, which a tool that finds them reports at E ≤ 1e-3 whatever its statistics. 'Extra' = pairs a tool reports that the reference does not. SABER builds its index in memory on every run; DIAMOND uses a prebuilt `.dmnd`. 74 of the 500 queries are also database members.

| Tool / preset | Time s | RSS MB | Pairs | Recall of exhaustive SW | Recall, E ≤ 1e-10 pairs | Extra pairs |
|---|---|---|---|---|---|---|
| saber_v1.4_very-fast | 2.0 (1.8–2.0) | 63 | 14565 | 0.8912 | 0.9577 | 0 |
| saber_v1.4_fast | 2.6 (2.5–2.6) | 73 | 15075 | 0.9224 | 0.9696 | 0 |
| saber_v1.4_default | 4.0 (4.0–4.2) | 104 | 15510 | 0.9490 | 0.9806 | 0 |
| saber_v1.4_sensitive | 3.7 (3.7–3.8) | 69 | 15632 | 0.9564 | 0.9827 | 0 |
| saber_v1.4_more-sensitive | 6.1 (6.0–6.2) | 125 | 15779 | 0.9654 | 0.9857 | 0 |
| saber_v1.4_exhaustive | 76.7 | 100 | 16344 | 1.0000 | 1.0000 | 0 |
| diamond_fast | 1.1 (1.1–1.1) | 46 | 4323 | 0.2644 | 0.3367 | 2 |
| diamond_default | 1.5 (1.5–1.5) | 54 | 6801 | 0.4141 | 0.5137 | 33 |
| diamond_mid-sensitive | 4.0 (3.9–4.2) | 57 | 10561 | 0.6387 | 0.7431 | 122 |
| diamond_sensitive | 7.3 (7.1–7.4) | 63 | 12742 | 0.7658 | 0.8690 | 226 |
| diamond_very-sensitive | 5.1 (5.0–5.2) | 108 | 15425 | 0.9216 | 0.9647 | 362 |
| diamond_ultra-sensitive | 17.1 (16.8–17.2) | 118 | 15800 | 0.9419 | 0.9684 | 406 |
| saber_v1.3_default | 4.1 (4.0–4.2) | 96 | 15510 | 0.9490 | 0.9806 | 0 |

### Database and query-count scaling (single runs)

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| scale_db25_saber_v1.4_default | 1.5 | 2.8 | 92 |
| scale_db25_saber_v1.3_default | 1.3 | 2.4 | 84 |
| scale_db25_diamond_very-sensitive | 1.7 | 3.2 | 50 |
| scale_db25_diamond_default | 0.5 | 0.9 | 31 |
| scale_db50_saber_v1.4_default | 2.5 | 4.5 | 97 |
| scale_db50_saber_v1.3_default | 2.2 | 4.1 | 89 |
| scale_db50_diamond_very-sensitive | 3.0 | 5.7 | 64 |
| scale_db50_diamond_default | 0.8 | 1.6 | 39 |
| scale_db400_saber_v1.4_default | 11.9 | 22.3 | 141 |
| scale_db400_saber_v1.3_default | 12.2 | 23.0 | 134 |
| scale_db400_diamond_very-sensitive | 16.3 | 31.4 | 372 |
| scale_db400_diamond_default | 4.8 | 9.2 | 134 |
| scaleq_50_saber_v1.4_default | 1.6 | 2.5 | 86 |
| scaleq_50_diamond_very-sensitive | 3.0 | 5.7 | 97 |
| scaleq_100_saber_v1.4_default | 1.6 | 3.0 | 92 |
| scaleq_100_diamond_very-sensitive | 3.2 | 6.0 | 106 |
| scaleq_250_saber_v1.4_default | 2.8 | 4.5 | 96 |
| scaleq_250_diamond_very-sensitive | 3.6 | 6.9 | 113 |

### Long proteins (99 queries > 1,500 aa against the 20,000-sequence database)

Every fifth database sequence longer than 1,500 residues, so each also finds itself. Single runs.

| Run | Time s | RSS MB | Pairs at E ≤ 1e-3 |
|---|---|---|---|
| saber_v1.4_default | 5.1 | 107 | 4744 |
| saber_v1.3_default | 5.1 | 98 | 4744 |
| diamond_very-sensitive | 5.9 | 113 | 4672 |

## Index construction, size and persistent-index search (MMseqs database and ×4 database)

Measured on machine B: Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| makedb_saber_v1.4_mmseqs_DB | 1.0 (0.8–1.2) | 1.4 | 37 |
| makedb_diamond_mmseqs_DB | 0.2 (0.2–0.2) | 0.4 | 28 |
| makedb_saber_v1.3_mmseqs_DB | 0.8 | 1.5 | 37 |
| db_warm_saber_v1.4_mmseqs_DB | 3.4 (3.2–3.5) | 6.4 | 97 |
| db_warm_diamond_very-sensitive_mmseqs_DB | 5.1 (4.9–5.2) | 9.7 | 109 |
| db_cold_saber_v1.4_mmseqs_DB | 3.3 (3.3–3.4) | 6.2 | 94 |
| db_cold_diamond_very-sensitive_mmseqs_DB | 5.1 (4.9–5.1) | 9.6 | 108 |
| makedb_saber_v1.4_mmseqs_DB_postings_p4M | 20.1 | 19.6 | 355 |
| db_order_scan_mmseqs_DB | 3.5 | 6.5 | 406 |
| db_order_query_mmseqs_DB | 6.4 | 12.2 | 838 |
| db_order_partition_mmseqs_DB | 6.4 | 11.9 | 859 |
| makedb_saber_v1.4_db400 | 3.0 (2.9–3.2) | 5.3 | 76 |
| makedb_diamond_db400 | 1.0 (0.9–1.1) | 1.7 | 58 |
| makedb_saber_v1.3_db400 | 3.4 | 5.9 | 76 |
| db_warm_saber_v1.4_db400 | 8.9 (8.9–9.6) | 16.8 | 117 |
| db_warm_diamond_very-sensitive_db400 | 16.0 (15.9–16.8) | 30.9 | 363 |
| db_cold_saber_v1.4_db400 | 9.3 (9.1–9.4) | 17.4 | 116 |
| db_cold_diamond_very-sensitive_db400 | 16.2 (15.8–16.2) | 30.9 | 362 |
| makedb_saber_v1.4_db400_postings_p4M | 82.3 | 81.1 | 1304 |
| db_order_scan_db400 | 9.9 | 18.2 | 1370 |
| db_order_query_db400 | 20.9 | 40.4 | 3321 |
| db_order_partition_db400 | 18.8 | 35.5 | 3341 |

File sizes / checks:
```
-rw-r--r-- 1 root root 11755080 Oct 19 22:42 /home/claude/bench/final14b/index/mmseqs_DB.dmnd
-rw-r--r-- 1 root root 10324348 Oct 19 22:42 /home/claude/bench/final14b/index/mmseqs_DB.sdx
-rw-r--r-- 1 root root 10324348 Oct 19 22:42 /home/claude/bench/final14b/index/mmseqs_DB.vo.sdx
-rw-r--r-- 1 root root 833367820 Oct 19 22:43 /home/claude/bench/final14b/index/mmseqs_DB.p4M.sdx
mmseqs_DB: scan == query-major output
mmseqs_DB: scan == partition-major output
-rw-r--r-- 1 root root 39883840 Oct 19 22:44 /home/claude/bench/final14b/index/db400.dmnd
-rw-r--r-- 1 root root 42965141 Oct 19 22:44 /home/claude/bench/final14b/index/db400.sdx
-rw-r--r-- 1 root root 42965141 Oct 19 22:44 /home/claude/bench/final14b/index/db400.vo.sdx
-rw-r--r-- 1 root root 3423358448 Oct 19 22:48 /home/claude/bench/final14b/index/db400.p4M.sdx
db400: scan == query-major output
db400: scan == partition-major output
```

## Translated searches

Measured on machine B: Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

blastx (`nx`): back-translated SCOP40 domains (random synonymous codons, random flanks and strand) against SCOP40 proteins; tblastn (`pn`): a fixed 1-in-37 subset of SCOP40 proteins (303 queries) against the back-translated set (DIAMOND has no tblastn mode). Metrics over the queries of each run (both splits). Single runs.

| Run | Queries | Time s | RSS MB | Sens. to 1st FP (SF) | Recall @ EPQ 0.01 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|
| nx_diamond_default | 11206 | 3.0 | 56 | 0.045 | 0.0123 | 4 |
| nx_diamond_ultra-sensitive | 11206 | 45.6 | 154 | 0.213 | 0.0718 | 41 |
| nx_diamond_very-sensitive | 11206 | 12.7 | 101 | 0.185 | 0.0641 | 33 |
| nx_saber_v1.3_default | 11206 | 18.5 | 140 | 0.199 | 0.0655 | 12 |
| nx_saber_v1.4_default | 11206 | 18.8 | 147 | 0.199 | 0.0655 | 12 |
| nx_saber_v1.4_fast | 11206 | 11.5 | 93 | 0.181 | 0.0602 | 11 |
| nx_saber_v1.4_more-sensitive | 11206 | 29.6 | 214 | 0.219 | 0.0710 | 20 |
| pn_saber_v1.3_default | 303 | 2.4 | 86 | 0.206 | 0.0610 | 1 |
| pn_saber_v1.4_default | 303 | 2.1 | 94 | 0.206 | 0.0610 | 1 |
| pn_saber_v1.4_more-sensitive | 303 | 3.0 | 101 | 0.228 | 0.0658 | 2 |

## Adversarial low-complexity / repeat inputs

Measured on machine B: Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

100 homopolymer / tandem-repeat queries (some with an embedded real domain) against the MMseqs database plus 2,000 low-complexity sequences. `lowc10_*` / `lowc4_*`: every 10th query (10 queries) / the first 4 of them, with masking, composition correction and the frequency cap all off (every repeat hit enumerated; 20-minute time limit, exit 124 = stopped by it). `mmq_*`: the 500 real MMseqs queries against the same polluted database — hits to the artificial `lowcdb*` sequences are false positives.

| Run | Time s | RSS MB | Exit | Hits (E ≤ 1e-3) | Hits to artificial low-complexity sequences |
|---|---|---|---|---|---|
| lowc_saber_v1.4_default | 1.3 | 85 | 0 | 73 | 0 |
| lowc_saber_v1.4_default_nomask | 19.9 | 176 | 0 | 8155 | 4218 |
| lowc4_saber_v1.4_default_nomask_nocap | 60.7 | 558 | 0 | 2332 | 1618 |
| lowc10_saber_v1.4_legacy_nomask_nocap | 27.3 | 867 | 0 | 2723 | 1786 |
| lowc_diamond_default | 1.2 | 45 | 0 | 21 | 0 |
| lowc_diamond_very-sensitive | 2.9 | 110 | 0 | 71 | 0 |
| mmq_saber_v1.4_default | 4.0 | 106 | 0 | 15469 | 0 |
| mmq_diamond_very-sensitive | 5.4 | 114 | 0 | 15368 | 0 |

## Two-domain queries and short fragments

Measured on machine B: Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

Retention = fraction of a domain's true positives from its single-domain search (each tool's own, E ≤ 1e-3) that are still found when the domain is part of a two-domain query. Fragments: 40-residue pieces of SCOP40 domains. Single runs.

| Tool / preset | Retention of the weaker domain's TPs | Retention (stronger) | Weaker-domain TPs (two-domain / alone) | Stronger-domain TPs (two-domain / alone) | FP in two-domain searches (E ≤ 1e-3) | Fragments: sens. to 1st FP (SF, test) | Fragment search time s |
|---|---|---|---|---|---|---|---|
| saber_v1.4_default | 0.960 | 0.961 | 190 / 198 | 1655 / 1723 | 0 | 0.071 | 1.4 |
| saber_v1.4_more-sensitive | 0.963 | 0.958 | 211 / 219 | 1742 / 1818 | 0 | 0.082 | 2.0 |
| diamond_very-sensitive | 0.938 | 0.947 | 212 / 226 | 1749 / 1846 | 0 | 0.043 | 1.7 |
| diamond_ultra-sensitive | 0.925 | 0.941 | 234 / 253 | 1937 / 2058 | 0 | 0.047 | 6.4 |

## Thread scaling (2-core machine)

Measured on machine B: Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

`onequery`: 20 long proteins searched one query per batch (`--query-batch 1`) against the ×4 database (index built in memory); every batch streams the whole database, in chunks shared by the threads (few-query searches: `docs/LIMITATIONS.md`).

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| t1_saber_v1.4_default | 10.1 (10.1–10.2) | 10.0 | 117 |
| t1_diamond_very-sensitive | 15.9 (15.5–15.9) | 15.8 | 47 |
| t1_saber_v1.4_onequery | 44.6 (44.5–44.9) | 43.4 | 101 |
| t2_saber_v1.4_default | 6.2 (6.1–6.2) | 10.7 | 137 |
| t2_diamond_very-sensitive | 8.2 (8.2–8.5) | 15.9 | 57 |
| t2_saber_v1.4_onequery | 24.3 (24.2–24.3) | 43.2 | 102 |

## Preset selection (tune splits only)

Measured on machine B: Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

Configurations considered for the 1.3 presets (unchanged in 1.4), on the tune splits of the at-scale benchmark (`sc`) and of SCOP40 all-against-all (`av`). Single runs, 2 threads; the presets' held-out results are in the sections above. Names: seed set (DS4/DS8/DS16 = the first 4/8/16 of DIAMOND's very-sensitive weight-7 seeds over Murphy-10; W8 = eight weight-6 seeds; DS16W8 = both), window filter threshold T, region threshold fraction f, two-stage scoring prefilter band pf.

| Configuration | Scale: time s | Scale: sens. (tune) | Scale: recall @ EPQ 0.01 (tune) | All-vs-all: time s | All-vs-all: sens. (tune) | All-vs-all: recall @ EPQ 0.01 (tune) |
|---|---|---|---|---|---|---|
| ds4_T28_f70 | 8.4 | 0.1419 | 0.0373 | 2.4 | 0.1563 | 0.0418 |
| ds8_T28_f70 | 12.0 | 0.1577 | 0.0418 | 3.6 | 0.1730 | 0.0467 |
| ds8_T27_f70 | 12.9 | 0.1597 | 0.0422 | 3.6 | 0.1755 | 0.0470 |
| ds16_T28_f70 | 22.6 | 0.1721 | 0.0460 | 5.7 | 0.1897 | 0.0508 |
| ds16_T27_f70 | 24.2 | 0.1742 | 0.0461 | 6.2 | 0.1925 | 0.0511 |
| ds16_T27_f65_pf16 | 26.0 | 0.1770 | 0.0467 | 6.3 | 0.1925 | 0.0511 |
| ds16_T27_f60 | 29.8 | 0.1781 | 0.0468 | 5.9 | 0.1925 | 0.0511 |
| ds16_T26_f60_pf16 | 31.9 | 0.1800 | 0.0469 | 6.5 | 0.1952 | 0.0515 |
| w8_T27_f65_pf16 | 35.6 | 0.1930 | 0.0482 | 5.4 | 0.2076 | 0.0536 |
| ds16_T25_f55_pf16 | 35.8 | 0.1809 | 0.0473 | 6.9 | 0.1961 | 0.0519 |
| ds16w8_T27_f65_pf16 | 60.8 | 0.1983 | 0.0495 | 10.5 | 0.2150 | 0.0551 |

## Scaling with database size

Measured on machine A: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

SCOP40 (11,206 domains) against SCOP40 + K × 20,000 reversed proteins (`make_scaling_dbs.py`; K = 4 is the at-scale database above). The decoy sources are the MMseqs2 example database and copies of it at 15 % substitution, so each original protein gives K similar decoys. Prebuilt indexes; SABER 1.4 `default` and DIAMOND `--very-sensitive` run alternately (median of the runs, range in brackets), the other rows single runs. Accuracy on the held-out test split (tune split in brackets); every hit to a decoy or another fold is a false positive.

| K | Database | Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP | Recall @ EPQ 0.01 | TP / FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|---|---|
| 1 | 31,206 seq., 11.0 M res. | saber_v1.4_very-fast | 3.7 | 6.8 | 102 | 0.1557 [0.1522] | 0.0677 [0.0396] | 10891 / 2 |
| 1 |  | saber_v1.4_default | 12.3 (12.2–12.4) | 22.3 | 243 | 0.2044 [0.1866] | 0.0859 [0.0487] | 13432 / 2 |
| 1 |  | saber_v1.4_more-sensitive | 23.0 | 41.6 | 362 | 0.2253 [0.2073] | 0.0932 [0.0523] | 14348 / 2 |
| 1 |  | saber_v1.3_default | 12.8 | 23.4 | 243 | 0.2044 [0.1866] | 0.0859 [0.0487] | 13432 / 2 |
| 1 |  | diamond_default | 3.3 | 6.2 | 59 | 0.0449 [0.0447] | 0.0153 [0.0100] | 2730 / 0 |
| 1 |  | diamond_very-sensitive | 14.9 (14.8–15.2) | 29.0 | 132 | 0.1858 [0.1749] | 0.0818 [0.0473] | 14418 / 9 |
| 1 |  | diamond_ultra-sensitive | 57.4 | 111.2 | 329 | 0.2135 [0.1966] | 0.0917 [0.0525] | 15974 / 11 |
| 2 | 51,206 seq., 20.1 M res. | saber_v1.4_very-fast | 5.3 | 9.8 | 113 | 0.1521 [0.1490] | 0.0664 [0.0389] | 10519 / 1 |
| 2 |  | saber_v1.4_default | 17.4 (16.8–17.9) | 32.3 | 264 | 0.2002 [0.1827] | 0.0844 [0.0477] | 12996 / 1 |
| 2 |  | saber_v1.4_more-sensitive | 35.2 | 64.5 | 385 | 0.2213 [0.2039] | 0.0916 [0.0510] | 13870 / 2 |
| 2 |  | saber_v1.3_default | 18.1 | 33.6 | 274 | 0.2002 [0.1827] | 0.0844 [0.0477] | 12996 / 1 |
| 2 |  | diamond_default | 4.6 | 8.6 | 97 | 0.0449 [0.0446] | 0.0153 [0.0099] | 2713 / 0 |
| 2 |  | diamond_very-sensitive | 20.5 (20.5–20.6) | 39.8 | 217 | 0.1846 [0.1741] | 0.0809 [0.0467] | 14081 / 5 |
| 2 |  | diamond_ultra-sensitive | 79.7 | 154.6 | 519 | 0.2114 [0.1949] | 0.0903 [0.0516] | 15532 / 8 |
| 4 | 91,206 seq., 38.2 M res. | saber_v1.4_very-fast | 6.0 | 11.2 | 111 | 0.1460 [0.1419] | 0.0639 [0.0373] | 10008 / 1 |
| 4 |  | saber_v1.4_default | 23.5 (23.0–24.9) | 44.3 | 280 | 0.1929 [0.1770] | 0.0826 [0.0467] | 12453 / 1 |
| 4 |  | saber_v1.4_more-sensitive | 53.7 | 99.9 | 420 | 0.2143 [0.1983] | 0.0898 [0.0495] | 13306 / 2 |
| 4 |  | saber_v1.3_default | 25.4 | 47.5 | 305 | 0.1929 [0.1770] | 0.0826 [0.0467] | 12453 / 1 |
| 4 |  | diamond_default | 6.6 | 12.7 | 146 | 0.0448 [0.0445] | 0.0152 [0.0099] | 2693 / 0 |
| 4 |  | diamond_very-sensitive | 29.6 (29.6–30.2) | 57.5 | 385 | 0.1832 [0.1729] | 0.0801 [0.0459] | 13720 / 5 |
| 4 |  | diamond_ultra-sensitive | 117.2 | 228.1 | 903 | 0.2075 [0.1916] | 0.0882 [0.0501] | 15013 / 5 |
| 8 | 171,206 seq., 74.4 M res. | saber_v1.4_very-fast | 10.9 | 20.5 | 152 | 0.1432 [0.1402] | 0.0630 [0.0369] | 9756 / 1 |
| 8 |  | saber_v1.4_default | 43.4 (43.3–43.7) | 82.3 | 349 | 0.1888 [0.1745] | 0.0804 [0.0456] | 12042 / 1 |
| 8 |  | saber_v1.4_more-sensitive | 101.5 | 190.2 | 536 | 0.2096 [0.1953] | 0.0873 [0.0487] | 12834 / 2 |
| 8 |  | saber_v1.3_default | 40.8 | 77.1 | 357 | 0.1865 [0.1725] | 0.0796 [0.0453] | 11966 / 1 |
| 8 |  | diamond_default | 10.9 | 21.0 | 261 | 0.0447 [0.0445] | 0.0151 [0.0098] | 2671 / 0 |
| 8 |  | diamond_very-sensitive | 49.4 (49.2–50.5) | 96.2 | 722 | 0.1819 [0.1717] | 0.0786 [0.0452] | 13355 / 7 |
| 8 |  | diamond_ultra-sensitive | 193.4 | 375.9 | 1694 | 0.2040 [0.1883] | 0.0850 [0.0489] | 14486 / 8 |
| 16 | 331,206 seq., 146.8 M res. | saber_v1.4_very-fast | 19.6 | 37.5 | 261 | 0.1419 [0.1387] | 0.0617 [0.0365] | 9487 / 1 |
| 16 |  | saber_v1.4_default | 84.9 (84.8–85.1) | 161.6 | 503 | 0.1852 [0.1714] | 0.0787 [0.0448] | 11607 / 1 |
| 16 |  | saber_v1.4_more-sensitive | 198.8 | 372.8 | 790 | 0.2047 [0.1910] | 0.0840 [0.0468] | 12321 / 2 |
| 16 |  | saber_v1.3_default | 76.4 | 145.0 | 465 | 0.1814 [0.1681] | 0.0776 [0.0442] | 11493 / 1 |
| 16 |  | diamond_default | 19.6 | 37.6 | 491 | 0.0447 [0.0444] | 0.0150 [0.0098] | 2658 / 0 |
| 16 |  | diamond_very-sensitive | 88.9 (87.1–90.2) | 172.7 | 1400 | 0.1799 [0.1695] | 0.0763 [0.0441] | 12952 / 2 |
| 16 |  | diamond_ultra-sensitive | 342.1 | 665.6 | 3189 | 0.1988 [0.1844] | 0.0828 [0.0477] | 13961 / 3 |

### Index construction and the one-off search (no preparation step) by database size

One-off: SABER 1.4 `default` searching the FASTA (index built in memory, included) and DIAMOND `makedb` + `blastp --very-sensitive` as one timed command; single runs. The prebuilt-index times are the medians above.

| K | SABER makedb s / MB | DIAMOND makedb s / MB | SABER one-off s / MB | DIAMOND makedb + blastp s / MB | SABER prebuilt s | DIAMOND prebuilt s |
|---|---|---|---|---|---|---|
| 1 | 1.2 / 38 | 0.3 / 29 | 13.8 / 245 | 15.6 / 133 | 12.3 | 14.9 |
| 2 | 2.2 / 60 | 0.6 / 38 | 18.5 / 277 | 20.6 / 217 | 17.4 | 20.5 |
| 4 | 4.2 / 73 | 1.1 / 59 | 27.5 / 308 | 32.1 / 385 | 23.5 | 29.6 |
| 8 | 8.2 / 120 | 2.1 / 101 | 51.3 / 415 | 52.9 / 723 | 43.4 | 49.4 |
| 16 | 17.7 / 214 | 4.0 / 183 | 100.0 / 637 | 91.4 / 1401 | 84.9 | 88.9 |

Output checks: K=1: in-memory index output == prebuilt index output; K=2: in-memory index output == prebuilt index output; K=4: in-memory index output == prebuilt index output; K=8: in-memory index output == prebuilt index output; K=16: in-memory index output == prebuilt index output.

## Length-matched decoys

Measured on machine A: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

The at-scale database's decoys are whole proteins, longer than the SCOP40 domains, and per-subject E-values make a hit in a long subject less significant — which favours SABER there. Here the same 80,000 reversed proteins are cut into pieces with the lengths of SCOP40 domains (`make_lenmatched_db.py`), so true and false subjects have the same length distribution. `saber_v1.4_default_dbmodel` = database-wide E-values. Accuracy columns are on the **held-out test split** (tune split in brackets): sensitivity = fraction of a query's same-superfamily domains ranked above its first false positive, averaged over queries; recall @ EPQ x = pooled fraction of same-superfamily pairs found before x false positives per query; TP / FP = same-superfamily / different-fold (or decoy) pairs reported at E ≤ 1e-3. Time = wall seconds, median (range) of the runs; RSS = peak resident memory (median).

```
/home/claude/bench/derived/scaling/scale_k4L.fa: 255349 sequences (244143 decoys), 37213650 residues
```

| Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP (superfamily) | Recall @ EPQ 0.01 | Recall @ EPQ 0.1 | Recall @ FDP 5 % | TP at E ≤ 1e-3 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|---|---|---|
| saber_v1.4_default | 25.4 (25.0–25.4) | 46.2 | 293 | 0.1894 [0.1733] | 0.0822 [0.0463] | 0.0904 | 0.0924 | 11828 | 1 |
| saber_v1.4_default_dbmodel | 25.4 | 47.0 | 291 | 0.1897 [0.1740] | 0.0769 [0.0446] | 0.0867 | 0.0889 | 12401 | 1 |
| saber_v1.4_more-sensitive | 55.8 | 103.1 | 437 | 0.2082 [0.1929] | 0.0884 [0.0497] | 0.0994 | 0.1024 | 12544 | 2 |
| diamond_very-sensitive | 32.1 (29.3–32.8) | 56.7 | 374 | 0.1835 [0.1733] | 0.0799 [0.0460] | 0.0886 | 0.0908 | 13733 | 6 |
| diamond_ultra-sensitive | 116.8 | 221.3 | 831 | 0.2086 [0.1924] | 0.0883 [0.0502] | 0.0998 | 0.1036 | 15033 | 8 |

| Tool / preset | FDP @ recall 0.06 | FDP @ recall 0.07 | FDP @ recall 0.08 | FDP @ recall 0.09 | FDP @ recall 0.1 |
|---|---|---|---|---|---|
| saber_v1.4_default | 0.0000 | 0.0003 | 0.0021 | 0.0276 | not reached |
| saber_v1.4_default_dbmodel | 0.0001 | 0.0006 | 0.0069 | 0.0622 | 0.3209 |
| saber_v1.4_more-sensitive | 0.0000 | 0.0001 | 0.0007 | 0.0042 | 0.0317 |
| diamond_very-sensitive | 0.0000 | 0.0004 | 0.0037 | 0.0415 | not reached |
| diamond_ultra-sensitive | 0.0000 | 0.0001 | 0.0007 | 0.0051 | 0.0293 |

## Scaling with the number of queries

Measured on machine A: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

1 to 11,206 SCOP40 domains (every N-th) against the at-scale database (38.2 M residues), prebuilt indexes; median of the runs. SABER 1.4 `default`, SABER 1.3.0 `default`, DIAMOND `--very-sensitive`.

| Queries | SABER 1.4 s | SABER 1.3.0 s | DIAMOND s | SABER 1.4 MB | DIAMOND MB |
|---|---|---|---|---|---|
| 1 | 1.2 (1.1–1.2) | 2.1 (1.9–2.2) | 9.5 (9.4–9.6) | 76 | 365 |
| 10 | 1.2 (1.1–1.2) | 1.9 (1.9–2.1) | 11.5 (11.1–11.6) | 78 | 367 |
| 100 | 1.5 (1.5–1.6) | 2.1 (2.0–2.4) | 11.9 (11.9–11.9) | 81 | 382 |
| 1,000 | 4.5 (4.3–4.7) | 4.3 (4.2–4.5) | 13.8 (13.7–14.3) | 105 | 380 |
| 11,206 | 23.1 (23.0–23.3) | 23.6 (23.3–24.0) | 28.8 (28.3–29.4) | 282 | 385 |

## Cross-machine check: the headline pairs on a second machine

Measured on machine B: Intel(R) Xeon(R) Processor @ 2.10GHz (cache 266240 KB).

`cross_machine_check.sh`: SABER 1.4 `default` and DIAMOND `--very-sensitive` alternately on the machine that ran the later sections, at scale (prebuilt indexes of the scale section; 5 runs each) and on the K = 16 scaling database (indexes built again; 3 runs each). Medians; the first machine's numbers are those of the scale and scaling sections.

| Database | Machine | SABER s | DIAMOND s | ratio | SABER MB | DIAMOND MB |
|---|---|---|---|---|---|---|
| at scale (38.2 M residues) | A | 24.2 (23.9–24.9) | 30.3 (29.6–30.8) | 0.80 | 280 | 385 |
| at scale (38.2 M residues) | B | 26.0 (25.7–26.6) | 27.4 (26.9–28.3) | 0.95 | 280 | 385 |
| K = 16 (146.8 M residues) | A | 84.9 (84.8–85.1) | 88.9 (87.1–90.2) | 0.95 | 503 | 1400 |
| K = 16 (146.8 M residues) | B | 84.7 (84.6–86.9) | 85.7 (84.4–87.7) | 0.99 | 502 | 1400 |

Outputs: scale_saber_v1.4_default.m8: identical to the main run's output; scale_diamond_very-sensitive.m8: identical to the main run's output; k16_saber_v1.4_default.m8: identical to the main run's output; k16_diamond_very-sensitive.m8: identical to the main run's output.

## One-off search at scale: SABER 1.4, SABER 1.3.0 and DIAMOND, alternating

Measured on a new virtual machine of type A: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

`oneoff_check.sh`: no preparation step, at scale (38.2 M residues). SABER 1.4 and SABER 1.3.0 `default` search the FASTA database (index built in memory, included in the time); DIAMOND runs `makedb` and then `blastp --very-sensitive` as one timed command. The three were run alternately. Machines are named by CPU type: this check ran on a new virtual machine (the second replacement of the session), so a machine of a type seen earlier is not necessarily the same host.

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| saber_v1.4_default_fasta | 39.8 (38.3–40.7) | 72.9 | 310 |
| saber_v1.3_default_fasta | 39.7 (37.9–41.4) | 73.5 | 307 |
| diamond_very-sensitive_makedb_blastp | 41.9 (41.6–44.2) | 80.4 | 388 |

Alternating runs, SABER 1.4 / SABER 1.3.0 / DIAMOND wall time: 40.3 / 38.8 / 41.7 s, 38.3 / 37.9 / 42.6 s, 38.7 / 41.4 / 44.2 s, 40.7 / 39.7 / 41.9 s, 39.8 / 40.2 / 41.6 s.

Outputs: saber_v1.4_default_fasta.m8: identical to the scale section's prebuilt-index output; saber_v1.3_default_fasta.m8: identical to SABER 1.4's output.

## Additional checks: SABER 1.3.0's counters, the frequency cap, 1.3.0 vs 1.4 outputs

Measured on a new virtual machine of type A: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

`extra_checks_14.sh`. Counters from `--stats` do not depend on the machine: candidate regions (`regions`) and query seeds skipped by the frequency cap (`seeds_frequent_skipped`) of SABER 1.3.0 (run here on indexes built by 1.4) and of SABER 1.4 (the scaling section's runs).

| Database | SABER 1.3.0 regions | SABER 1.4 regions | ratio | SABER 1.3.0 seeds skipped | SABER 1.4 seeds skipped |
|---|---|---|---|---|---|
| K = 8 (74.4 M residues) | 3,230,710 | 5,172,373 | 1.60 | 1,686,817 | 1,443,280 |
| K = 16 (146.8 M residues) | 3,903,270 | 10,245,641 | 2.62 | 1,667,699 | 1,432,680 |

The frequency cap at scale (38.2 M residues): SABER 1.4 `default` with the automatic cap and without it (`--seed-max-occ 1000000000`), alternately; accuracy on both splits (the capped run's accuracy is the scale section's, its output is identical).

| Run | Time s | CPU s | RSS MB | Sens. test [tune] | Recall @ EPQ 0.01 test [tune] |
|---|---|---|---|---|---|
| scale_saber_v1.4_default | 34.4 (32.3–36.4) | 61.3 | 280 | 0.1929 [0.1770] | 0.0826 [0.0467] |
| scale_saber_v1.4_default_nocap | 40.8 (40.4–41.2) | 75.5 | 283 | 0.1964 [0.1810] | 0.0833 [0.0471] |

Outputs: k8 1.3.0 output vs the scaling section's 1.3.0 output: identical (52386 lines); k16 1.3.0 output vs the scaling section's 1.3.0 output: identical (50118 lines); at scale, 1.4 default here vs the scale section's output: identical (55121 lines); scop/saber_v1.4_default.m8 vs scop/saber_v1.3_default.m8: identical (64092 lines); scop/saber_v1.4_more-sensitive.m8 vs scop/saber_v1.3_more-sensitive.m8: identical (76585 lines); mmseqs/saber_v1.4_default.m8 vs mmseqs/saber_v1.3_default.m8: identical (15510 lines); mmseqs/long_saber_v1.4_default.m8 vs mmseqs/long_saber_v1.3_default.m8: identical (4744 lines); scale/saber_v1.4_default.m8 vs scale/saber_v1.3_default.m8: identical (55121 lines); translated/nx_saber_v1.4_default.m8 vs translated/nx_saber_v1.3_default.m8: identical (76624 lines); translated/pn_saber_v1.4_default.m8 vs translated/pn_saber_v1.3_default.m8: identical (1779 lines).

