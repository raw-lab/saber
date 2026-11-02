## At scale: SCOP40 against SCOP40 + 80,000 reversed proteins (91,206 sequences, 38.2 M residues)

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

11,206 SCOP40 queries, 2 threads, E ≤ 10, no hit limit; every hit to a reversed protein (`REV_*`) or to a domain of another fold is a false positive. Prebuilt indexes for both tools (`saber --makedb`, `diamond makedb`; built once, timed below). SABER 1.3 `default` and DIAMOND `--very-sensitive` were run alternately (5 runs each) so that machine drift affects both alike; the other rows are median of 3 (single runs: `diamond_ultra-sensitive`, `saber_v1.3_default_dbmodel`, `saber_v1.2_default`). Accuracy columns are on the **held-out test split** (tune split in brackets): sensitivity = fraction of a query's same-superfamily domains ranked above its first false positive, averaged over queries; recall @ EPQ x = pooled fraction of same-superfamily pairs found before x false positives per query; TP / FP = same-superfamily / different-fold (or decoy) pairs reported at E ≤ 1e-3. Time = wall seconds, median (range) of the runs; RSS = peak resident memory (median).

| Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP (superfamily) | Recall @ EPQ 0.01 | Recall @ EPQ 0.1 | Recall @ FDP 5 % | TP at E ≤ 1e-3 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|---|---|---|
| saber_v1.3_very-fast | 7.0 (6.7–7.3) | 13.1 | 129 | 0.1460 [0.1419] | 0.0639 [0.0373] | 0.0693 | 0.0697 | 10008 | 1 |
| saber_v1.3_fast | 13.5 (13.4–13.5) | 25.0 | 180 | 0.1713 [0.1597] | 0.0741 [0.0422] | 0.0811 | 0.0826 | 11367 | 1 |
| saber_v1.3_default | 27.9 (27.5–28.2) | 51.3 | 303 | 0.1929 [0.1770] | 0.0826 [0.0467] | 0.0917 | 0.0941 | 12453 | 1 |
| saber_v1.3_sensitive | 39.0 (37.7–39.4) | 72.8 | 269 | 0.2052 [0.1930] | 0.0874 [0.0482] | 0.0978 | 0.1009 | 13023 | 2 |
| saber_v1.3_more-sensitive | 64.8 (63.1–65.6) | 117.8 | 449 | 0.2143 [0.1983] | 0.0898 [0.0495] | 0.1011 | 0.1049 | 13306 | 2 |
| diamond_fast | 4.7 (4.7–4.9) | 8.9 | 144 | 0.0125 [0.0122] | 0.0036 [0.0026] | 0.0036 | 0.0036 | 664 | 0 |
| diamond_default | 7.4 (7.1–7.6) | 13.9 | 147 | 0.0448 [0.0445] | 0.0152 [0.0099] | 0.0154 | 0.0154 | 2693 | 0 |
| diamond_mid-sensitive | 21.9 (20.9–23.0) | 41.3 | 155 | 0.0978 [0.0965] | 0.0387 [0.0234] | 0.0406 | 0.0404 | 6799 | 1 |
| diamond_sensitive | 41.7 (39.5–41.7) | 78.5 | 161 | 0.1421 [0.1398] | 0.0626 [0.0352] | 0.0674 | 0.0677 | 10939 | 5 |
| diamond_more-sensitive | 40.2 (39.7–41.7) | 76.2 | 156 | 0.1424 [0.1407] | 0.0630 [0.0363] | 0.0679 | 0.0682 | 11015 | 5 |
| diamond_very-sensitive | 33.3 (32.5–34.3) | 63.9 | 386 | 0.1832 [0.1729] | 0.0801 [0.0459] | 0.0884 | 0.0905 | 13720 | 5 |
| diamond_ultra-sensitive | 120.7 | 233.3 | 907 | 0.2075 [0.1916] | 0.0882 [0.0501] | 0.0997 | 0.1032 | 15013 | 5 |
| diamond_very-sensitive_c1 | 32.9 (32.5–33.3) | 62.4 | 386 | 0.1832 [0.1729] | 0.0801 [0.0459] | 0.0884 | 0.0905 | 13720 | 5 |
| saber_v1.3_default_dbmodel | 28.4 | 52.7 | 304 | 0.1867 [0.1723] | 0.0769 [0.0438] | 0.0859 | 0.0879 | 12219 | 1 |
| saber_v1.2_default | 146.6 | 282.6 | 198 | 0.1802 [0.1716] | 0.0735 [0.0428] | 0.0830 | 0.0849 | 11870 | 2 |

`diamond_very-sensitive_c1` = DIAMOND `--very-sensitive -c1` (one seed-index chunk, DIAMOND's speed/memory option); `saber_v1.3_default_dbmodel` = SABER 1.3 `default` with `--evalue-model database` (1.2's E-values: same alignments and bit scores, different E-values and therefore a different ranking across subjects). Both are run by `run_all.sh scale_extra`.

### False-discovery proportion at matched recall (test split, at scale)

| Tool / preset | FDP @ recall 0.06 | FDP @ recall 0.07 | FDP @ recall 0.08 | FDP @ recall 0.09 | FDP @ recall 0.1 |
|---|---|---|---|---|---|
| saber_v1.3_very-fast | 0.0006 | 0.0563 | not reached | not reached | not reached |
| saber_v1.3_fast | 0.0001 | 0.0007 | 0.0245 | not reached | not reached |
| saber_v1.3_default | 0.0000 | 0.0002 | 0.0018 | 0.0216 | 0.1467 |
| saber_v1.3_sensitive | 0.0000 | 0.0001 | 0.0005 | 0.0062 | 0.0422 |
| saber_v1.3_more-sensitive | 0.0000 | 0.0001 | 0.0005 | 0.0034 | 0.0233 |
| diamond_fast | not reached | not reached | not reached | not reached | not reached |
| diamond_default | not reached | not reached | not reached | not reached | not reached |
| diamond_mid-sensitive | not reached | not reached | not reached | not reached | not reached |
| diamond_sensitive | 0.0011 | not reached | not reached | not reached | not reached |
| diamond_more-sensitive | 0.0009 | not reached | not reached | not reached | not reached |
| diamond_very-sensitive | 0.0000 | 0.0003 | 0.0034 | 0.0447 | not reached |
| diamond_ultra-sensitive | 0.0000 | 0.0000 | 0.0005 | 0.0049 | 0.0302 |
| diamond_very-sensitive_c1 | 0.0000 | 0.0003 | 0.0034 | 0.0447 | not reached |
| saber_v1.3_default_dbmodel | 0.0001 | 0.0007 | 0.0082 | 0.0780 | 0.3813 |
| saber_v1.2_default | 0.0002 | 0.0017 | 0.0170 | 0.1389 | not reached |

Alternating runs, SABER 1.3 `default` / DIAMOND `--very-sensitive` wall time: 27.9 / 32.6 s, 28.1 / 34.3 s, 27.5 / 33.3 s, 27.6 / 33.9 s, 28.2 / 32.5 s (ratio 0.81–0.87); CPU time 51 / 62 s, 52 / 66 s, 51 / 64 s, 51 / 65 s, 52 / 63 s.

### Index construction (scale database)

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| makedb_saber_v1.3 | 4.0 (4.0–4.1) | 7.1 | 73 |
| makedb_diamond | 1.2 (1.2–1.3) | 2.1 | 59 |
| makedb_saber_v1.2 | 5.8 | 8.5 | 141 |

```
-rw-r--r-- 1 root root  40690880 Oct 19 14:26 /home/claude/bench/final13/scale/scale.dmnd
-rw-r--r-- 1 root root  43983827 Oct 19 14:26 /home/claude/bench/final13/scale/scale.sdx
-rw-r--r-- 1 root root 183298099 Oct 19 14:26 /home/claude/bench/final13/scale/scale.v12.sdx
```

### Practical settings at scale (E ≤ 1e-3, 25 targets per query)

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| practical_saber_v1.3_default | 26.3 (25.3–26.4) | 48.6 | 304 |
| practical_diamond_very-sensitive | 32.9 (32.3–32.9) | 62.9 | 386 |

### DIAMOND's seed-search algorithm at scale

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

DIAMOND chooses between its double-indexed and query-indexed seed search automatically; here `--very-sensitive --algo 1` (query-indexed) is forced, run alternately with SABER 1.3 `default` (whose accuracy is in the table above; same output). Accuracy columns are on the **held-out test split** (tune split in brackets): sensitivity = fraction of a query's same-superfamily domains ranked above its first false positive, averaged over queries; recall @ EPQ x = pooled fraction of same-superfamily pairs found before x false positives per query; TP / FP = same-superfamily / different-fold (or decoy) pairs reported at E ≤ 1e-3. Time = wall seconds, median (range) of the runs; RSS = peak resident memory (median).

| Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP (superfamily) | Recall @ EPQ 0.01 | Recall @ EPQ 0.1 | Recall @ FDP 5 % | TP at E ≤ 1e-3 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|---|---|---|
| saber_v1.3_default | 28.9 (28.9–30.8) | 53.9 | 304 | – [–] | – [–] | – | – | – | – |
| diamond_very-sensitive_algo1 | 72.4 (71.9–73.0) | 137.5 | 434 | 0.2076 [0.1913] | 0.0879 [0.0499] | 0.0992 | 0.1028 | 15007 | 5 |

Alternating runs, SABER / DIAMOND `--algo 1` wall time: 30.8 / 71.9 s, 28.9 / 72.4 s, 28.9 / 73.0 s.

### One-off search at scale, index construction included

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

No preparation step: SABER 1.3 `default` searches the FASTA database directly (`-d`, index built in memory, included in its time); DIAMOND runs `makedb` and then `blastp --very-sensitive` (one timed command). Alternating runs.

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| saber_v1.3_default_fasta | 32.7 (31.0–34.2) | 59.6 | 308 |
| diamond_very-sensitive_makedb_blastp | 33.8 (33.4–35.4) | 65.0 | 386 |

Alternating runs, SABER / DIAMOND wall time: 32.7 / 33.7 s, 33.2 / 34.2 s, 34.2 / 35.4 s, 31.0 / 33.8 s, 32.3 / 33.4 s.

Output check: in-memory index output == prebuilt index output.

## SCOP40 all-against-all (11,206 domains, 2 threads, E ≤ 10, no hit limit)

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

SABER builds its index in memory on every run (`-d`; `default_db` = prebuilt `--db`), DIAMOND uses a prebuilt `.dmnd`. Accuracy columns are on the **held-out test split** (tune split in brackets): sensitivity = fraction of a query's same-superfamily domains ranked above its first false positive, averaged over queries; recall @ EPQ x = pooled fraction of same-superfamily pairs found before x false positives per query; TP / FP = same-superfamily / different-fold (or decoy) pairs reported at E ≤ 1e-3. Time = wall seconds, median (range) of the runs; RSS = peak resident memory (median).

| Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP (superfamily) | Recall @ EPQ 0.01 | Recall @ EPQ 0.1 | Recall @ FDP 5 % | TP at E ≤ 1e-3 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|---|---|---|
| saber_v1.3_very-fast | 2.7 (2.5–2.7) | 4.7 | 58 | 0.1611 [0.1563] | 0.0700 [0.0418] | 0.0780 | 0.0792 | 11473 | 3 |
| saber_v1.3_fast | 4.3 (4.1–4.6) | 7.5 | 92 | 0.1895 [0.1755] | 0.0802 [0.0470] | 0.0921 | 0.0944 | 13098 | 3 |
| saber_v1.3_default | 7.1 (6.9–7.3) | 12.3 | 137 | 0.2105 [0.1925] | 0.0878 [0.0511] | 0.1023 | 0.1058 | 14242 | 4 |
| saber_v1.3_sensitive | 7.0 (6.9–7.2) | 12.4 | 83 | 0.2229 [0.2076] | 0.0929 [0.0536] | 0.1089 | 0.1135 | 14909 | 5 |
| saber_v1.3_more-sensitive | 12.3 (11.8–12.8) | 21.5 | 161 | 0.2340 [0.2150] | 0.0958 [0.0551] | 0.1137 | 0.1189 | 15308 | 5 |
| saber_v1.3_legacy | 7.0 (7.0–7.5) | 12.9 | 37 | 0.1504 [0.1433] | 0.0640 [0.0383] | 0.0701 | 0.0707 | 10587 | 4 |
| saber_v1.3_exhaustive | 176.9 | 345.4 | 37 | 0.2724 [0.2514] | 0.1023 [0.0581] | 0.1248 | 0.1345 | 16145 | 6 |
| diamond_fast | 1.6 (1.5–1.6) | 2.9 | 29 | 0.0125 [0.0123] | 0.0036 [0.0026] | 0.0036 | 0.0036 | 673 | 0 |
| diamond_default | 2.3 (2.3–2.4) | 4.4 | 32 | 0.0451 [0.0451] | 0.0155 [0.0101] | 0.0157 | 0.0157 | 2775 | 1 |
| diamond_mid-sensitive | 6.6 (6.1–6.8) | 12.6 | 31 | 0.0989 [0.0978] | 0.0402 [0.0243] | 0.0421 | 0.0419 | 7148 | 6 |
| diamond_sensitive | 11.1 (11.0–11.6) | 21.1 | 46 | 0.1429 [0.1407] | 0.0656 [0.0374] | 0.0690 | 0.0690 | 11857 | 14 |
| diamond_more-sensitive | 11.5 (11.4–11.8) | 21.5 | 49 | 0.1432 [0.1415] | 0.0661 [0.0386] | 0.0694 | 0.0694 | 11940 | 14 |
| diamond_very-sensitive | 9.3 (9.3–9.8) | 18.0 | 59 | 0.1880 [0.1763] | 0.0851 [0.0498] | 0.0960 | 0.0973 | 15433 | 17 |
| diamond_ultra-sensitive | 33.5 (33.4–35.7) | 64.7 | 117 | 0.2182 [0.1996] | 0.0959 [0.0552] | 0.1117 | 0.1150 | 17395 | 22 |
| saber_v1.2_default | 8.0 (7.6–8.3) | 15.0 | 29 | 0.1969 [0.1869] | 0.0805 [0.0477] | 0.0943 | 0.0977 | 13911 | 10 |
| saber_v1.2_more-sensitive | 17.5 (17.1–17.6) | 32.6 | 37 | 0.2245 [0.2057] | 0.0888 [0.0510] | 0.1053 | 0.1103 | 15331 | 11 |
| saber_v1.3_default_db | 6.9 (6.8–8.0) | 11.8 | 132 | 0.2105 [0.1925] | 0.0878 [0.0511] | 0.1023 | 0.1058 | 14242 | 4 |

### False-discovery proportion at matched recall (test split)

Pooled over all test-split queries: hits sorted by E-value, FDP = FP / (TP + FP) at the most stringent cut-off that reaches the recall. 'not reached' = the tool never reaches that recall at E ≤ 10.

| Tool / preset | FDP @ recall 0.06 | FDP @ recall 0.08 | FDP @ recall 0.09 | FDP @ recall 0.1 | FDP @ recall 0.11 |
|---|---|---|---|---|---|
| saber_v1.3_very-fast | 0.0003 | 0.0632 | not reached | not reached | not reached |
| saber_v1.3_fast | 0.0001 | 0.0034 | 0.0200 | 0.1487 | not reached |
| saber_v1.3_default | 0.0000 | 0.0007 | 0.0049 | 0.0187 | 0.0915 |
| saber_v1.3_sensitive | 0.0000 | 0.0004 | 0.0015 | 0.0086 | 0.0304 |
| saber_v1.3_more-sensitive | 0.0000 | 0.0003 | 0.0009 | 0.0056 | 0.0165 |
| saber_v1.3_legacy | 0.0009 | not reached | not reached | not reached | not reached |
| saber_v1.3_exhaustive | 0.0000 | 0.0003 | 0.0009 | 0.0019 | 0.0069 |
| diamond_fast | not reached | not reached | not reached | not reached | not reached |
| diamond_default | not reached | not reached | not reached | not reached | not reached |
| diamond_mid-sensitive | not reached | not reached | not reached | not reached | not reached |
| diamond_sensitive | 0.0005 | not reached | not reached | not reached | not reached |
| diamond_more-sensitive | 0.0004 | not reached | not reached | not reached | not reached |
| diamond_very-sensitive | 0.0000 | 0.0007 | 0.0073 | not reached | not reached |
| diamond_ultra-sensitive | 0.0000 | 0.0003 | 0.0009 | 0.0048 | 0.0187 |
| saber_v1.2_default | 0.0001 | 0.0035 | 0.0138 | 0.0710 | 0.2844 |
| saber_v1.2_more-sensitive | 0.0001 | 0.0007 | 0.0045 | 0.0130 | 0.0458 |
| saber_v1.3_default_db | 0.0000 | 0.0007 | 0.0049 | 0.0187 | 0.0915 |

### Practical settings (E ≤ 1e-3, 25 targets per query), SCOP40 all-vs-all

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| practical_saber_v1.3_fast | 4.3 (4.0–4.3) | 7.3 | 91 |
| practical_saber_v1.3_default | 6.6 (6.4–7.3) | 11.4 | 134 |
| practical_diamond_default | 2.2 (2.2–2.5) | 4.2 | 32 |
| practical_diamond_very-sensitive | 9.9 (9.8–10.1) | 19.0 | 60 |
| practical_diamond_ultra-sensitive | 34.1 (34.1–34.2) | 65.8 | 118 |

### SCOP40 index build and persistent-index search

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| diamond_makedb | 0.1 (0.1–0.1) | 0.1 | 20 |
| saber_v1.3_makedb_default | 0.4 (0.4–0.4) | 0.7 | 28 |
| saber_v1.3_default_db_warm | 6.9 (6.8–8.0) | 11.8 | 132 |

## Null behaviour: real SCOP40 queries against shuffled databases

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

Real queries against composition-preserving decoys (`shuf`: whole-sequence shuffle; `wshuf10`: shuffle inside 10-residue windows). Hits to decoys derived from a sequence of the query's own fold are excluded; every remaining hit is a false positive. A calibrated E-value gives at most about E false hits per query (EPQ ≤ E). `dbmodel` = SABER 1.3 with `--evalue-model database` (the 1.2 E-values).

| Run | EPQ at E ≤ 1e-3 | E ≤ 0.01 | E ≤ 0.1 | E ≤ 1 | E ≤ 10 |
|---|---|---|---|---|---|
| diamond_default_shuf | 0.0000 | 0.0001 | 0.0008 | 0.0029 | 0.0103 |
| diamond_default_wshuf10 | 0.0000 | 0.0001 | 0.0010 | 0.0034 | 0.0106 |
| diamond_ultra-sensitive_shuf | 0.0012 | 0.0064 | 0.0242 | 0.0920 | 0.2509 |
| diamond_ultra-sensitive_wshuf10 | 0.0007 | 0.0073 | 0.0327 | 0.1116 | 0.2996 |
| diamond_very-sensitive_shuf | 0.0006 | 0.0031 | 0.0126 | 0.0431 | 0.1070 |
| diamond_very-sensitive_wshuf10 | 0.0004 | 0.0046 | 0.0180 | 0.0535 | 0.1258 |
| saber_v1.2_default_shuf | 0.0003 | 0.0014 | 0.0162 | 0.1872 | 1.7155 |
| saber_v1.2_default_wshuf10 | 0.0002 | 0.0025 | 0.0207 | 0.2204 | 1.9460 |
| saber_v1.3_default_dbmodel_shuf | 0.0005 | 0.0020 | 0.0202 | 0.2311 | 2.0698 |
| saber_v1.3_default_dbmodel_wshuf10 | 0.0004 | 0.0030 | 0.0229 | 0.2680 | 2.3158 |
| saber_v1.3_default_shuf | 0.0001 | 0.0006 | 0.0043 | 0.0637 | 0.7152 |
| saber_v1.3_default_wshuf10 | 0.0002 | 0.0006 | 0.0073 | 0.0763 | 0.8099 |
| saber_v1.3_fast_shuf | 0.0000 | 0.0004 | 0.0032 | 0.0451 | 0.5033 |
| saber_v1.3_fast_wshuf10 | 0.0001 | 0.0004 | 0.0053 | 0.0533 | 0.5671 |
| saber_v1.3_more-sensitive_shuf | 0.0001 | 0.0009 | 0.0077 | 0.0949 | 1.0573 |
| saber_v1.3_more-sensitive_wshuf10 | 0.0002 | 0.0006 | 0.0099 | 0.1110 | 1.1934 |

## MMseqs2 example workload (500 queries × 20,000 UniProt sequences, E ≤ 1e-3)

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

Recall is measured against SABER 1.3's exhaustive Smith–Waterman search (`--preset exhaustive`: 16344 query–subject pairs at E ≤ 1e-3). That reference shares SABER's scoring, composition correction and E-value statistics, so near the cut-off it favours SABER; 'Recall, E ≤ 1e-10 pairs' uses only the 12600 reference pairs with E ≤ 1e-10, which a tool that finds them reports at E ≤ 1e-3 whatever its statistics. 'Extra' = pairs a tool reports that the reference does not. SABER builds its index in memory on every run; DIAMOND uses a prebuilt `.dmnd`. 74 of the 500 queries are also database members.

| Tool / preset | Time s | RSS MB | Pairs | Recall of exhaustive SW | Recall, E ≤ 1e-10 pairs | Extra pairs |
|---|---|---|---|---|---|---|
| saber_v1.3_very-fast | 2.2 (2.1–2.2) | 62 | 14565 | 0.8912 | 0.9577 | 0 |
| saber_v1.3_fast | 2.9 (2.8–2.9) | 69 | 15075 | 0.9224 | 0.9696 | 0 |
| saber_v1.3_default | 4.6 (4.4–4.6) | 96 | 15510 | 0.9490 | 0.9806 | 0 |
| saber_v1.3_sensitive | 4.3 (4.0–4.3) | 73 | 15632 | 0.9564 | 0.9827 | 0 |
| saber_v1.3_more-sensitive | 6.9 (6.6–6.9) | 114 | 15779 | 0.9654 | 0.9857 | 0 |
| saber_v1.3_exhaustive | 80.9 | 141 | 16344 | 1.0000 | 1.0000 | 0 |
| diamond_fast | 1.2 (1.1–1.2) | 48 | 4323 | 0.2644 | 0.3367 | 2 |
| diamond_default | 1.7 (1.7–1.7) | 47 | 6801 | 0.4141 | 0.5137 | 33 |
| diamond_mid-sensitive | 4.6 (4.4–4.6) | 55 | 10561 | 0.6387 | 0.7431 | 122 |
| diamond_sensitive | 8.5 (8.1–8.8) | 63 | 12742 | 0.7658 | 0.8690 | 226 |
| diamond_very-sensitive | 6.1 (5.7–6.4) | 108 | 15425 | 0.9216 | 0.9647 | 362 |
| diamond_ultra-sensitive | 19.5 (18.9–19.8) | 119 | 15800 | 0.9419 | 0.9684 | 406 |
| saber_v1.2_default | 6.3 (6.3–6.6) | 106 | 15607 | 0.9470 | 0.9802 | 129 |

### Database and query-count scaling (single runs)

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| scale_db25_saber_v1.3_default | 1.7 | 3.0 | 84 |
| scale_db25_saber_v1.2_default | 2.2 | 4.3 | 44 |
| scale_db25_diamond_very-sensitive | 2.0 | 3.7 | 44 |
| scale_db25_diamond_default | 0.6 | 1.0 | 44 |
| scale_db50_saber_v1.3_default | 2.7 | 4.9 | 89 |
| scale_db50_saber_v1.2_default | 3.8 | 7.3 | 72 |
| scale_db50_diamond_very-sensitive | 3.5 | 6.4 | 73 |
| scale_db50_diamond_default | 1.0 | 1.9 | 33 |
| scale_db400_saber_v1.3_default | 13.6 | 25.6 | 135 |
| scale_db400_saber_v1.2_default | 21.3 | 40.3 | 251 |
| scale_db400_diamond_very-sensitive | 18.7 | 35.6 | 365 |
| scale_db400_diamond_default | 5.5 | 10.5 | 136 |
| scaleq_50_saber_v1.3_default | 1.8 | 3.2 | 82 |
| scaleq_50_diamond_very-sensitive | 3.3 | 6.2 | 98 |
| scaleq_100_saber_v1.3_default | 2.1 | 3.8 | 84 |
| scaleq_100_diamond_very-sensitive | 3.9 | 7.2 | 109 |
| scaleq_250_saber_v1.3_default | 3.0 | 5.6 | 89 |
| scaleq_250_diamond_very-sensitive | 4.3 | 8.1 | 114 |

### Long proteins (99 queries > 1,500 aa against the 20,000-sequence database)

Every fifth database sequence longer than 1,500 residues, so each also finds itself. Single runs.

| Run | Time s | RSS MB | Pairs at E ≤ 1e-3 |
|---|---|---|---|
| saber_v1.3_default | 5.6 | 100 | 4744 |
| saber_v1.2_default | 7.5 | 110 | 4712 |
| diamond_very-sensitive | 6.1 | 104 | 4672 |

## Index construction, size and persistent-index search (MMseqs database and ×4 database)

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| makedb_saber_v1.3_mmseqs_DB | 0.9 (0.9–1.1) | 1.7 | 37 |
| makedb_diamond_mmseqs_DB | 0.3 (0.3–0.3) | 0.5 | 29 |
| makedb_saber_v1.2_mmseqs_DB | 1.4 | 2.1 | 60 |
| db_warm_saber_v1.3_mmseqs_DB | 3.6 (3.6–4.0) | 6.8 | 89 |
| db_warm_diamond_very-sensitive_mmseqs_DB | 5.6 (5.6–5.7) | 10.7 | 115 |
| db_cold_saber_v1.3_mmseqs_DB | 3.7 (3.6–3.8) | 6.8 | 88 |
| db_cold_diamond_very-sensitive_mmseqs_DB | 5.8 (5.5–6.0) | 10.8 | 115 |
| makedb_saber_v1.3_mmseqs_DB_postings_p4M | 16.8 | 16.7 | 352 |
| db_order_scan_mmseqs_DB | 3.6 | 6.5 | 399 |
| db_order_query_mmseqs_DB | 6.0 | 11.5 | 838 |
| db_order_partition_mmseqs_DB | 6.1 | 11.5 | 858 |
| makedb_saber_v1.3_db400 | 3.6 (3.4–3.8) | 6.4 | 76 |
| makedb_diamond_db400 | 1.1 (1.0–1.3) | 2.0 | 59 |
| makedb_saber_v1.2_db400 | 5.1 | 7.8 | 151 |
| db_warm_saber_v1.3_db400 | 10.2 (9.8–10.7) | 19.2 | 125 |
| db_warm_diamond_very-sensitive_db400 | 18.5 (18.0–18.6) | 35.2 | 368 |
| db_cold_saber_v1.3_db400 | 10.2 (9.8–10.3) | 19.2 | 123 |
| db_cold_diamond_very-sensitive_db400 | 18.3 (18.0–19.0) | 35.0 | 365 |
| makedb_saber_v1.3_db400_postings_p4M | 80.7 | 82.0 | 1303 |
| db_order_scan_db400 | 11.2 | 20.7 | 1382 |
| db_order_query_db400 | 19.6 | 37.5 | 3321 |
| db_order_partition_db400 | 19.0 | 35.9 | 3340 |

File sizes / checks:
```
-rw-r--r-- 1 root root 11755080 Oct 19 15:13 /home/claude/bench/final13/index/mmseqs_DB.dmnd
-rw-r--r-- 1 root root 10324348 Oct 19 15:13 /home/claude/bench/final13/index/mmseqs_DB.sdx
-rw-r--r-- 1 root root 43091108 Oct 19 15:13 /home/claude/bench/final13/index/mmseqs_DB.v12.sdx
-rw-r--r-- 1 root root 833367820 Oct 19 15:14 /home/claude/bench/final13/index/mmseqs_DB.p4M.sdx
mmseqs_DB: scan == query-major output
mmseqs_DB: scan == partition-major output
-rw-r--r-- 1 root root  39883840 Oct 19 15:14 /home/claude/bench/final13/index/db400.dmnd
-rw-r--r-- 1 root root  42965141 Oct 19 15:14 /home/claude/bench/final13/index/db400.sdx
-rw-r--r-- 1 root root 175243442 Oct 19 15:14 /home/claude/bench/final13/index/db400.v12.sdx
-rw-r--r-- 1 root root 3423358448 Oct 19 15:19 /home/claude/bench/final13/index/db400.p4M.sdx
db400: scan == query-major output
db400: scan == partition-major output
```

## Translated searches

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

blastx (`nx`): back-translated SCOP40 domains (random synonymous codons, random flanks and strand) against SCOP40 proteins; tblastn (`pn`): a fixed 1-in-37 subset of SCOP40 proteins (303 queries) against the back-translated set (DIAMOND has no tblastn mode). Metrics over the queries of each run (both splits). Single runs.

| Run | Queries | Time s | RSS MB | Sens. to 1st FP (SF) | Recall @ EPQ 0.01 | FP at E ≤ 1e-3 |
|---|---|---|---|---|---|---|
| nx_diamond_default | 11206 | 3.6 | 58 | 0.045 | 0.0123 | 4 |
| nx_diamond_ultra-sensitive | 11206 | 53.1 | 154 | 0.213 | 0.0718 | 41 |
| nx_diamond_very-sensitive | 11206 | 14.8 | 107 | 0.185 | 0.0641 | 33 |
| nx_saber_v1.2_default | 11206 | 25.7 | 34 | 0.189 | 0.0609 | 25 |
| nx_saber_v1.3_default | 11206 | 20.3 | 141 | 0.199 | 0.0655 | 12 |
| nx_saber_v1.3_fast | 11206 | 12.9 | 89 | 0.181 | 0.0602 | 11 |
| nx_saber_v1.3_more-sensitive | 11206 | 32.5 | 213 | 0.219 | 0.0710 | 20 |
| pn_saber_v1.2_default | 303 | 2.3 | 62 | 0.194 | 0.0579 | 0 |
| pn_saber_v1.3_default | 303 | 2.8 | 86 | 0.206 | 0.0610 | 1 |
| pn_saber_v1.3_more-sensitive | 303 | 3.5 | 89 | 0.228 | 0.0658 | 2 |

## Adversarial low-complexity / repeat inputs

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

100 homopolymer / tandem-repeat queries (some with an embedded real domain) against the MMseqs database plus 2,000 low-complexity sequences. `lowc10_*` / `lowc4_*`: every 10th query (10 queries) / the first 4 of them, with masking, composition correction and the frequency cap all off (every repeat hit enumerated; 20-minute time limit, exit 124 = stopped by it). `mmq_*`: the 500 real MMseqs queries against the same polluted database — hits to the artificial `lowcdb*` sequences are false positives.

| Run | Time s | RSS MB | Exit | Hits (E ≤ 1e-3) | Hits to artificial low-complexity sequences |
|---|---|---|---|---|---|
| lowc_saber_v1.3_default | 2.0 | 85 | 0 | 73 | 0 |
| lowc_saber_v1.3_default_nomask | 23.8 | 153 | 0 | 8155 | 4218 |
| lowc4_saber_v1.3_default_nomask_nocap | 61.9 | 564 | 0 | 2332 | 1618 |
| lowc10_saber_v1.3_legacy_nomask_nocap | 27.4 | 868 | 0 | 2723 | 1786 |
| lowc_diamond_default | 1.7 | 47 | 0 | 21 | 0 |
| lowc_diamond_very-sensitive | 3.6 | 111 | 0 | 71 | 0 |
| mmq_saber_v1.3_default | 5.1 | 98 | 0 | 15469 | 0 |
| mmq_diamond_very-sensitive | 6.4 | 116 | 0 | 15368 | 0 |

## Two-domain queries and short fragments

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

Retention = fraction of a domain's true positives from its single-domain search (each tool's own, E ≤ 1e-3) that are still found when the domain is part of a two-domain query. Fragments: 40-residue pieces of SCOP40 domains. Single runs.

| Tool / preset | Retention of the weaker domain's TPs | Retention (stronger) | Weaker-domain TPs (two-domain / alone) | Stronger-domain TPs (two-domain / alone) | FP in two-domain searches (E ≤ 1e-3) | Fragments: sens. to 1st FP (SF, test) | Fragment search time s |
|---|---|---|---|---|---|---|---|
| saber_v1.3_default | 0.960 | 0.961 | 190 / 198 | 1655 / 1723 | 0 | 0.071 | 1.5 |
| saber_v1.3_more-sensitive | 0.963 | 0.958 | 211 / 219 | 1742 / 1818 | 0 | 0.082 | 2.5 |
| diamond_very-sensitive | 0.938 | 0.947 | 212 / 226 | 1749 / 1846 | 0 | 0.043 | 1.9 |
| diamond_ultra-sensitive | 0.925 | 0.941 | 234 / 253 | 1937 / 2058 | 0 | 0.047 | 7.7 |

## Thread scaling (2-core machine)

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

`onequery`: 20 long proteins searched one query per batch (`--query-batch 1`) against the ×4 database (index built in memory); every batch streams the whole database, in chunks shared by the threads (few-query searches: `docs/LIMITATIONS.md`).

| Run | Time s | CPU s | RSS MB |
|---|---|---|---|
| t1_saber_v1.3_default | 11.9 (11.7–12.6) | 11.5 | 114 |
| t1_diamond_very-sensitive | 19.2 (18.6–20.2) | 18.9 | 49 |
| t1_saber_v1.3_onequery | 90.0 (87.1–92.5) | 87.7 | 101 |
| t2_saber_v1.3_default | 7.3 (7.2–7.4) | 12.3 | 134 |
| t2_diamond_very-sensitive | 10.6 (10.0–11.9) | 20.3 | 55 |
| t2_saber_v1.3_onequery | 50.6 (50.3–52.3) | 92.5 | 102 |

## Preset selection (tune splits only)

Measured on: Intel(R) Xeon(R) Processor @ 2.80GHz (cache 33792 KB).

Configurations considered for the 1.3 presets, on the tune splits of the at-scale benchmark (`sc`) and of SCOP40 all-against-all (`av`). Single runs, 2 threads; the presets' held-out results are in the sections above. Names: seed set (DS4/DS8/DS16 = the first 4/8/16 of DIAMOND's very-sensitive weight-7 seeds over Murphy-10; W8 = eight weight-6 seeds; DS16W8 = both), window filter threshold T, region threshold fraction f, two-stage scoring prefilter band pf.

| Configuration | Scale: time s | Scale: sens. (tune) | Scale: recall @ EPQ 0.01 (tune) | All-vs-all: time s | All-vs-all: sens. (tune) | All-vs-all: recall @ EPQ 0.01 (tune) |
|---|---|---|---|---|---|---|
| ds4_T28_f70 | 7.7 | 0.1419 | 0.0373 | 2.8 | 0.1563 | 0.0418 |
| ds8_T28_f70 | 13.0 | 0.1577 | 0.0418 | 4.0 | 0.1730 | 0.0467 |
| ds8_T27_f70 | 14.4 | 0.1597 | 0.0422 | 4.2 | 0.1755 | 0.0470 |
| ds16_T28_f70 | 24.4 | 0.1721 | 0.0460 | 6.6 | 0.1897 | 0.0508 |
| ds16_T27_f65_pf16 | 28.0 | 0.1770 | 0.0467 | 7.5 | 0.1925 | 0.0511 |
| ds16_T27_f70 | 29.1 | 0.1742 | 0.0461 | 6.9 | 0.1925 | 0.0511 |
| ds16_T27_f60 | 33.1 | 0.1781 | 0.0468 | 7.7 | 0.1925 | 0.0511 |
| ds16_T26_f60_pf16 | 35.7 | 0.1800 | 0.0469 | 7.4 | 0.1952 | 0.0515 |
| w8_T27_f65_pf16 | 38.6 | 0.1930 | 0.0482 | 6.6 | 0.2076 | 0.0536 |
| ds16_T25_f55_pf16 | 39.6 | 0.1809 | 0.0473 | 7.7 | 0.1961 | 0.0519 |
| ds16w8_T27_f65_pf16 | 63.6 | 0.1983 | 0.0495 | 12.2 | 0.2150 | 0.0551 |

