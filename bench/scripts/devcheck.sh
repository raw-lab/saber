#!/usr/bin/env bash
# Before / after the 1.2 performance changes: run a development build from
# before them (OLD) and the benchmarked build (NEW) on the same workloads in
# one session, then compare outputs with the outputs evaluated in the main run
# (scripts/compare_m8.py) and record time, peak RSS and --stats stage times.
#
#   OLD=/path/to/dev-build NEW=/path/to/saber bench/scripts/devcheck.sh
#   (DATA, DER, OUT as in run_all.sh; the delivered build's long-protein output
#    is taken from OUT/recheck, see recheck.sh)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OLD=${OLD:-/home/claude/bench/binaries/saber-1.2.0-dev-1638}
NEW=${NEW:-/home/claude/work/target175/release/saber}
DATA=${DATA:-/home/claude/bench/data}
DER=${DER:-/home/claude/bench/derived}
OUT=${OUT:-/home/claude/bench/final}
THREADS=${THREADS:-2}
SEC=$OUT/devcheck
mkdir -p "$SEC"
rm -f "$SEC/times.jsonl"
RB="python3 -I $HERE/runbench.py"
SC="-q $DATA/scop40.fa -d $DATA/scop40.fa -k 100000 -E 10 -t $THREADS --quiet"
MM="-q $DATA/mmseqs_QUERY.fasta -d $DATA/mmseqs_DB.fasta -k 100000 -E 1e-3 -t $THREADS --quiet"
LG="-q $OUT/mmseqs/long_queries.fa -d $DATA/mmseqs_DB.fasta -k 100000 -E 1e-3 -t $THREADS --quiet"
NX="--mode nx -q $DER/scop40.bt.fna -d $DATA/scop40.fa -k 100000 -E 10 -t $THREADS --quiet"
for v in old new; do
  B=$OLD; [ $v = new ] && B=$NEW
  for p in very-fast default more-sensitive; do
    $RB --label scop_${p}_$v --out "$SEC/times.jsonl" --stdout "$SEC/scop_${p}_$v.m8" -- $B $SC --preset $p --stats "$SEC/scop_${p}_$v.stats.json"
  done
  $RB --label mm_default_$v --out "$SEC/times.jsonl" --stdout "$SEC/mm_default_$v.m8" -- $B $MM --preset default --stats "$SEC/mm_default_$v.stats.json"
  $RB --label mm_long_$v --out "$SEC/times.jsonl" --stdout "$SEC/mm_long_$v.m8" -- $B $LG --stats "$SEC/mm_long_$v.stats.json"
  $RB --label nx_default_$v --out "$SEC/times.jsonl" --stdout "$SEC/nx_default_$v.m8" -- $B $NX --preset default --stats "$SEC/nx_default_$v.stats.json"
done
C="python3 -I $HERE/compare_m8.py"
{
  # development build vs the outputs evaluated in the main run
  $C "$SEC/scop_very-fast_old.m8" "$OUT/scop/saber_v1.2_very-fast.m8" "SCOP40 very-fast"
  $C "$SEC/scop_default_old.m8" "$OUT/scop/saber_v1.2_default.m8" "SCOP40 default"
  $C "$SEC/scop_more-sensitive_old.m8" "$OUT/scop/saber_v1.2_more-sensitive.m8" "SCOP40 more-sensitive"
  $C "$SEC/mm_default_old.m8" "$OUT/mmseqs/saber_v1.2_default.m8" "MMseqs2 example default"
  $C "$SEC/nx_default_old.m8" "$OUT/translated/nx_saber_v1.2_default.m8" "blastx default"
  $C "$SEC/mm_long_old.m8" "$OUT/mmseqs/long_saber_v1.2_default.m8" "long proteins, benchmarked build"
  if [ -f "$OUT/recheck/mmseqs_long_saber_v1.2_default.m8" ]; then
    $C "$SEC/mm_long_old.m8" "$OUT/recheck/mmseqs_long_saber_v1.2_default.m8" "long proteins, delivered build"
  fi
  # the NEW build re-run here must reproduce the main run's outputs
  $C "$SEC/scop_default_new.m8" "$OUT/scop/saber_v1.2_default.m8" "SCOP40 default, NEW re-run"
  $C "$SEC/nx_default_new.m8" "$OUT/translated/nx_saber_v1.2_default.m8" "blastx default, NEW re-run"
} > "$SEC/outputs_check.txt"
cat "$SEC/outputs_check.txt"
