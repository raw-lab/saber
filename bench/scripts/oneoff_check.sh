#!/usr/bin/env bash
# One-off searches at scale (no preparation step), repeated on the current
# machine: SABER 1.4 and SABER 1.3.0 searching the FASTA database directly
# (index built in memory, included in the time) and DIAMOND `makedb` +
# `blastp --very-sensitive` as one timed command, run alternately KEYREPS
# times each. Added for the 1.4.0 run, whose oneoff section ran on machine B
# only (see bench/README.md).
#
#   bench/scripts/oneoff_check.sh            (same variables as run_all.sh)
#
# Writes OUT/checks_oneoff/: environment.txt (this machine), times.jsonl and
# outputs_check.txt. Needs the scale section's FASTA database
# (DER/scop40.scale.fa) and its SABER 1.4 output (OUT/scale/).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VN=${VN:-/home/claude/bench/binaries/saber-1.4.0}
VO=${VO:-/home/claude/bench/binaries/saber-1.3.0}
DMD=${DMD:-/home/claude/tools/d228/diamond}
DATA=${DATA:-/home/claude/bench/data}
DER=${DER:-/home/claude/bench/derived}
OUT=${OUT:-/home/claude/bench/final14}
THREADS=${THREADS:-2}
KEYREPS=${KEYREPS:-5}
PY="python3 -I"
RB="$PY $HERE/runbench.py"
SEC=$OUT/checks_oneoff
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
  echo "$VO: $($VO --version 2>&1 | tail -1) sha256 $(sha256sum $VO | cut -d' ' -f1)"
  echo "database: $DER/scop40.scale.fa sha256 $(sha256sum $DER/scop40.scale.fa | cut -d' ' -f1)"
} > "$SEC/environment.txt"
Q=$DATA/scop40.fa
DB=$DER/scop40.scale.fa
[ -s $DB ] && [ -s $OUT/scale/saber_v1.4_default.m8 ] || { log "one-off check: run the scale section first"; exit 1; }
log "one-off check: SABER 1.4 and 1.3.0 on the FASTA vs DIAMOND makedb + blastp --very-sensitive, interleaved"
for r in $(seq $KEYREPS); do
  timed saber_v1.4_default_fasta 1 $SEC/saber_v1.4_default_fasta.m8 -- $VN -q $Q -d $DB -k 100000 -E 10 -t $THREADS --quiet
  timed saber_v1.3_default_fasta 1 $SEC/saber_v1.3_default_fasta.m8 -- $VO -q $Q -d $DB -k 100000 -E 10 -t $THREADS --quiet
  timed diamond_very-sensitive_makedb_blastp 1 - -- bash -c "$DMD makedb --in $DB -d $SEC/oneoff --threads $THREADS --quiet && $DMD blastp -q $Q -d $SEC/oneoff.dmnd --threads $THREADS -k 0 -e 10 --very-sensitive --quiet -o /dev/null"
done
rm -f $SEC/oneoff.dmnd
: > $SEC/outputs_check.txt
if cmp -s $SEC/saber_v1.4_default_fasta.m8 $OUT/scale/saber_v1.4_default.m8; then
  echo "saber_v1.4_default_fasta.m8: identical to the scale section's prebuilt-index output" >> $SEC/outputs_check.txt
else
  echo "saber_v1.4_default_fasta.m8: DIFFERS from the scale section's prebuilt-index output" >> $SEC/outputs_check.txt
fi
if cmp -s $SEC/saber_v1.3_default_fasta.m8 $SEC/saber_v1.4_default_fasta.m8; then
  echo "saber_v1.3_default_fasta.m8: identical to SABER 1.4's output" >> $SEC/outputs_check.txt
else
  echo "saber_v1.3_default_fasta.m8: DIFFERS from SABER 1.4's output" >> $SEC/outputs_check.txt
fi
cat $SEC/outputs_check.txt
log "one-off check: done"
