#!/usr/bin/env bash
# Checks behind three statements of the 1.4.0 documentation, run after the
# suite (same variables as run_all.sh; needs the scale section's indexes and
# the scaling section's outputs):
#
#   1. SABER 1.3.0's own counters on the 74.4 and 146.8 M-residue scaling
#      databases (`--stats`): candidate regions scored and query seeds skipped
#      by the frequency cap, to set against 1.4's (scaling/k*.stats.json). The
#      outputs are compared with the scaling section's 1.3.0 outputs.
#   2. The frequency cap at scale: SABER 1.4 `default` with and without the
#      automatic cap (`--seed-max-occ 1000000000`), alternately, REPS2 runs
#      each, with accuracy on both splits.
#   3. SABER 1.4 against 1.3.0 outputs where both ran in the suite (SCOP40,
#      MMseqs2 example, long proteins, at scale, translated): byte comparison.
#   4. The length-matched database's size, appended to OUT/lenmatch/datasets.txt
#      (for the 1.4.0 run this step was run by hand after steps 1-3).
#
#   bench/scripts/extra_checks_14.sh
#
# Writes OUT/checks_extra/: environment.txt, times.jsonl, *.stats.json,
# *.metrics.json and outputs_check.txt.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VN=${VN:-/home/claude/bench/binaries/saber-1.4.0}
VO=${VO:-/home/claude/bench/binaries/saber-1.3.0}
DATA=${DATA:-/home/claude/bench/data}
DER=${DER:-/home/claude/bench/derived}
OUT=${OUT:-/home/claude/bench/final14}
THREADS=${THREADS:-2}
REPS2=${REPS2:-2}
PY="python3 -I"
RB="$PY $HERE/runbench.py"
EVAL="$PY $HERE/eval_scop.py"
LAB=$DATA/scop40.labels.tsv
SEC=$OUT/checks_extra
mkdir -p "$SEC"
log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$OUT/run.log"; }
timed() {
  local label=$1 reps=$2 out=$3; shift 4
  local so=/dev/null
  [ "$out" != "-" ] && so=$out
  $RB --label "$label" --reps "$reps" --out "$SEC/times.jsonl" --stdout "$so" -- "$@" 2>>"$SEC/runbench.err" || log "FAILED: $label (exit $?)"
}
{
  echo "date: $(date -Is)"
  echo "host: $(uname -a)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
  echo "cache: $(grep -m1 'cache size' /proc/cpuinfo | cut -d: -f2)"
  echo "cores: $(nproc)"
  echo "$VN: $($VN --version 2>&1 | tail -1) sha256 $(sha256sum $VN | cut -d' ' -f1)"
  echo "$VO: $($VO --version 2>&1 | tail -1) sha256 $(sha256sum $VO | cut -d' ' -f1)"
} > "$SEC/environment.txt"
Q=$DATA/scop40.fa
: > $SEC/outputs_check.txt
same() {  # same LABEL FILE1 FILE2
  if cmp -s "$2" "$3"; then echo "$1: identical ($(wc -l < "$2") lines)" >> $SEC/outputs_check.txt
  else echo "$1: DIFFERENT" >> $SEC/outputs_check.txt; fi
}

log "extra checks: SABER 1.3.0 counters on the K=8 and K=16 scaling databases"
for k in 8 16; do
  DB=$DER/scaling/scale_k$k.fa
  [ -s $DB ] || { log "extra checks: $DB missing"; exit 1; }
  timed k${k}_makedb_saber_v1.4 1 - -- $VN --makedb $SEC/k$k.sdx -d $DB --quiet
  timed k${k}_saber_v1.3_default 1 $SEC/k$k.saber_v1.3_default.m8 -- $VO -q $Q --db $SEC/k$k.sdx -k 100000 -E 10 -t $THREADS --quiet --stats $SEC/k$k.saber_v1.3_default.stats.json
  same "k$k 1.3.0 output vs the scaling section's 1.3.0 output" $SEC/k$k.saber_v1.3_default.m8 $OUT/scaling/k$k.saber_v1.3_default.m8
  rm -f $SEC/k$k.sdx $SEC/k$k.saber_v1.3_default.m8
done

log "extra checks: the frequency cap at scale, SABER 1.4 default with and without it, interleaved"
IDX=$OUT/scale/scale.sdx
[ -s $IDX ] || { log "extra checks: run the scale section first"; exit 1; }
for r in $(seq $REPS2); do
  timed scale_saber_v1.4_default 1 $SEC/scale_saber_v1.4_default.m8 -- $VN -q $Q --db $IDX -k 100000 -E 10 -t $THREADS --quiet --stats $SEC/scale_saber_v1.4_default.stats.json
  timed scale_saber_v1.4_default_nocap 1 $SEC/scale_saber_v1.4_default_nocap.m8 -- $VN -q $Q --db $IDX -k 100000 -E 10 -t $THREADS --seed-max-occ 1000000000 --quiet --stats $SEC/scale_saber_v1.4_default_nocap.stats.json
done
same "at scale, 1.4 default here vs the scale section's output" $SEC/scale_saber_v1.4_default.m8 $OUT/scale/saber_v1.4_default.m8
for n in scale_saber_v1.4_default_nocap; do
  for sub in tune test; do
    $EVAL $SEC/$n.m8 $LAB --subset $sub --name $n --decoy-prefix REV_ --decoys 80000 > $SEC/$n.$sub.metrics.json
  done
done
rm -f $SEC/scale_saber_v1.4_default.m8 $SEC/scale_saber_v1.4_default_nocap.m8

log "extra checks: 1.4 vs 1.3.0 outputs of the suite"
for p in "scop/saber_v1.4_default.m8 scop/saber_v1.3_default.m8" "scop/saber_v1.4_more-sensitive.m8 scop/saber_v1.3_more-sensitive.m8" \
         "mmseqs/saber_v1.4_default.m8 mmseqs/saber_v1.3_default.m8" "mmseqs/long_saber_v1.4_default.m8 mmseqs/long_saber_v1.3_default.m8" \
         "scale/saber_v1.4_default.m8 scale/saber_v1.3_default.m8" "translated/nx_saber_v1.4_default.m8 translated/nx_saber_v1.3_default.m8" \
         "translated/pn_saber_v1.4_default.m8 translated/pn_saber_v1.3_default.m8"; do
  set -- $p
  [ -s $OUT/$1 ] && [ -s $OUT/$2 ] && same "$1 vs $2" $OUT/$1 $OUT/$2
done
cat $SEC/outputs_check.txt

# 4. The size of the length-matched database (the lenmatch section records
#    only its checksum when the database already exists).
DBL=$DER/scaling/scale_k4L.fa
if [ -s $DBL ] && [ -d $OUT/lenmatch ] && ! grep -qs ' residues$' $OUT/lenmatch/datasets.txt; then
  awk '/^>/{n++; if ($0 ~ /^>REV_/) d++; next} {r+=length($0)} END{printf "%s: %d sequences (%d decoys), %d residues\n", FILENAME, n, d, r}' $DBL >> $OUT/lenmatch/datasets.txt
fi
log "extra checks: done"
