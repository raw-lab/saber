#!/usr/bin/env bash
# Headline comparisons repeated on the current machine, for a suite whose
# sections did not all run on one machine (the 1.4.0 run: the virtual machine
# was replaced after the scop section).
#
#   bench/scripts/cross_machine_check.sh            (same variables as run_all.sh)
#
# Writes OUT/checks_machine/: environment.txt (this machine), times.jsonl and
# outputs_check.txt. At scale (the scale section's prebuilt indexes): SABER
# 1.4 default and DIAMOND --very-sensitive alternately, KEYREPS runs each. At
# K = 16 (the scaling section's largest database): both indexes built again,
# REPS alternating runs each. The outputs are compared with the main run's.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VN=${VN:-/home/claude/bench/binaries/saber-1.4.0}
DMD=${DMD:-/home/claude/tools/d228/diamond}
DATA=${DATA:-/home/claude/bench/data}
DER=${DER:-/home/claude/bench/derived}
OUT=${OUT:-/home/claude/bench/final14}
THREADS=${THREADS:-2}
REPS=${REPS:-3}
KEYREPS=${KEYREPS:-5}
PY="python3 -I"
RB="$PY $HERE/runbench.py"
SEC=$OUT/checks_machine
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
  echo "mem: $(grep MemTotal /proc/meminfo)"
  echo "transparent huge pages: $(cat /sys/kernel/mm/transparent_hugepage/enabled 2>/dev/null)"
  echo "diamond: $($DMD version 2>&1 | head -1) sha256 $(sha256sum $DMD | cut -d' ' -f1)"
  echo "$VN: $($VN --version 2>&1 | tail -1) sha256 $(sha256sum $VN | cut -d' ' -f1)"
} > "$SEC/environment.txt"
Q=$DATA/scop40.fa
IDX=$OUT/scale
[ -s $IDX/scale.sdx ] && [ -s $IDX/scale.dmnd ] || { log "cross-machine check: run the scale section first"; exit 1; }
log "cross-machine check: at scale, SABER 1.4 default vs DIAMOND --very-sensitive, interleaved"
for r in $(seq $KEYREPS); do
  timed scale_saber_v1.4_default 1 $SEC/scale_saber_v1.4_default.m8 -- $VN -q $Q --db $IDX/scale.sdx -k 100000 -E 10 -t $THREADS --quiet
  timed scale_diamond_very-sensitive 1 $SEC/scale_diamond_very-sensitive.m8 -- $DMD blastp -q $Q -d $IDX/scale.dmnd --threads $THREADS -k 0 -e 10 --very-sensitive --quiet -o /dev/stdout
done
DB=$DER/scaling/scale_k16.fa
log "cross-machine check: K=16"
timed k16_makedb_saber_v1.4 1 - -- $VN --makedb $SEC/k16.sdx -d $DB --quiet
timed k16_makedb_diamond 1 - -- $DMD makedb --in $DB -d $SEC/k16 --threads $THREADS --quiet
for r in $(seq $REPS); do
  timed k16_saber_v1.4_default 1 $SEC/k16_saber_v1.4_default.m8 -- $VN -q $Q --db $SEC/k16.sdx -k 100000 -E 10 -t $THREADS --quiet
  timed k16_diamond_very-sensitive 1 $SEC/k16_diamond_very-sensitive.m8 -- $DMD blastp -q $Q -d $SEC/k16.dmnd --threads $THREADS -k 0 -e 10 --very-sensitive --quiet -o /dev/stdout
done
rm -f $SEC/k16.sdx $SEC/k16.dmnd
: > $SEC/outputs_check.txt
for p in "scale_saber_v1.4_default.m8 $OUT/scale/saber_v1.4_default.m8" "scale_diamond_very-sensitive.m8 $OUT/scale/diamond_very-sensitive.m8" \
         "k16_saber_v1.4_default.m8 $OUT/scaling/k16.saber_v1.4_default.m8" "k16_diamond_very-sensitive.m8 $OUT/scaling/k16.diamond_very-sensitive.m8"; do
  set -- $p
  if cmp -s $SEC/$1 $2; then echo "$1: identical to the main run's output" >> $SEC/outputs_check.txt; else echo "$1: DIFFERS from the main run's output" >> $SEC/outputs_check.txt; fi
done
cat $SEC/outputs_check.txt
log "cross-machine check: done"
