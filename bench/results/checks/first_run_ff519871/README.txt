First run of the scale, scaling, lenmatch and qscale sections of the 1.4.0
suite, with the build before the scaled seed-frequency cap floor (sha256
ff519871..., recorded in env/environment.txt). That build already had the
region-threshold change but kept 1.3's cap rule (64 x the mean occurrences per
code, at least 256), so its --stats count the query seeds that 1.3's rule
skips. The run was superseded by the release build's (72c2d9ec...); its files
are kept as evidence for IMPLEMENTATION_CHECKLIST.md A.1 and A.2:

  scale/times.jsonl                     timings at scale (machine A; SABER 1.4
                                        default / DIAMOND --very-sensitive
                                        alternately, 5 runs each)
  scaling/times.jsonl                   timings of the scaling section
  scaling/k*.saber_v1.4_default.stats.json
                                        --stats per database size
                                        (counters.seeds_frequent_skipped,
                                        counters.regions)
  ../v13_v14_outputs_build_ff519871.txt homologues lost with this build
