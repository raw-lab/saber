Development A/B measurements behind the 1.4.0 speed items (IMPLEMENTATION_CHECKLIST.md A.4).
Same machine as the suite, the same day, before the suite ran; at-scale benchmark (SCOP40 queries
against SCOP40 + 80,000 reversed proteins, prebuilt index), 2 threads, runs alternating between the
variants. Timing records: runbench.py JSON lines (wall_s, cpu_s, maxrss_kb, the command). Every
variant's output was compared with the reference output (byte-identical in all cases).

Binaries (sha256 in binaries.sha256; dev builds of the 1.4.0 source with one change each):
  ref       1.4.0 before the window-filter repacking (sha256 8d0abc50..., then installed as saber-1.4.0)
  new1/f1   + packed window bytes, survivors written without branches (room for all hits of a subject)
  f3        + room reserved per database position instead (sha256 ff519871...; the released
            1.4.0, 72c2d9ec..., is f3 plus the scaled seed-frequency cap floor below)
  pf0/1/2   entry prefetch variants: whole runs into L1 one subject ahead (pf0 = adopted) /
            into L2 one subject ahead plus L1 8 candidates ahead (pf1) / only 8 candidates ahead (pf2)
  big, mid  larger scan rounds: survivor budget 8 M / 4 M and rounds of up to 4 M / 2 M residues
            (adopted: 2 M / 1 M)
  xg0/xg1   X-drop table reads: 16 scalar reads (xg0 = adopted) / two AVX2 gathers (xg1)
  f2        index construction: seed codes counted in cache-sized code ranges (radix partitioned)
            instead of directly with prefetching (adopted)
  The frequency-cap check (the run that led to the scaled cap floor of the released build):
  cap_k16   f3 (ff519871) on the K = 16 scaling database with --seed-max-occ 939
            (= 256 x 146.8 M / 40 M, what the scaled floor gives); f3's own automatic cap is
            the suite's first run (64 x mean occurrences: 840-873 per shape)
  cap_k8    the released build (72c2d9ec, scaled floor: cap 477) on the K = 8 database

Files:
  ab1.jsonl   ref vs new1 (2 alternating runs each)
  f3.jsonl    ref vs f3
  pf.jsonl    pf0 vs pf1 vs pf2
  big.jsonl   f1 vs big;  mid.jsonl  f1 vs mid vs big
  xg.jsonl    xg0 vs xg1
  mk.jsonl    --makedb of the at-scale database (38.2 M residues): ref vs f2 (index files identical)
  mk16.jsonl  --makedb of the K = 16 scaling database (146.8 M residues): ref vs f2 (index files identical)
  cap_k16.jsonl, cap_k16_939.*   timing, --stats and test/tune metrics of the cap_k16 run
  cap_k8.jsonl, cap_k8_new.*     the same for the cap_k8 run
  *.stats.json  --stats of the last run of each variant (TSC cycles per scan stage: c_scan_*)
