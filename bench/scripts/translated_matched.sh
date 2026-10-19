#!/usr/bin/env bash
# Translated searches at matched sensitivity: the 34-query blastx / tblastn
# subsets of run_all.sh's translated section, searched on ONE machine by
# SABER 1.1.0 (full dynamic programming over six frames), every relevant
# SABER 1.2 preset (including exhaustive) and DIAMOND, so speed can be
# compared at similar sensitivity.
#
#   bench/scripts/translated_matched.sh        (same variables as run_all.sh)
#
# Needs the translated section's subset files (nxsub_queries.fna,
# pnsub_queries.fa) and scop40.dmnd in $OUT/translated.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V110=${V110:-/home/claude/work/target175_opt_v110/release/saber}
V12=${V12:-/home/claude/work/target175/release/saber}
DMD=${DMD:-/home/claude/tools/d228/diamond}
DATA=${DATA:-/home/claude/bench/data}
DER=${DER:-/home/claude/bench/derived}
OUT=${OUT:-/home/claude/bench/final}
THREADS=${THREADS:-2}
REPS=${REPS:-3}
PY="python3 -I"
TRS=$OUT/translated
SEC=$OUT/translated_matched
mkdir -p "$SEC"
rm -f "$SEC/times.jsonl"
P=$DATA/scop40.fa
NQ=$DER/scop40.bt.fna
LAB=$DATA/scop40.labels.tsv
timed() { # LABEL REPS OUTFILE -- cmd...
  local label=$1 reps=$2 out=$3; shift 4
  $PY "$HERE/runbench.py" --label "$label" --reps "$reps" --out "$SEC/times.jsonl" --stdout "$out" -- "$@" 2>>"$SEC/runbench.err"
}
NX="-q $TRS/nxsub_queries.fna -d $P -k 100000 -E 10 -t $THREADS --quiet"
PN="-q $TRS/pnsub_queries.fa -d $NQ -k 100000 -E 10 -t $THREADS --quiet"
timed nxsub_saber_v1.1.0 1 "$SEC/nxsub_saber_v1.1.0.m8" -- timeout 7200 $V110 --mode nx $NX
timed pnsub_saber_v1.1.0 1 "$SEC/pnsub_saber_v1.1.0.m8" -- timeout 7200 $V110 --mode pn $PN
for p in default more-sensitive exhaustive; do
  timed nxsub_saber_v1.2_$p $REPS "$SEC/nxsub_saber_v1.2_$p.m8" -- $V12 --mode nx $NX --preset $p
  timed pnsub_saber_v1.2_$p $REPS "$SEC/pnsub_saber_v1.2_$p.m8" -- $V12 --mode pn $PN --preset $p
done
for m in very-sensitive ultra-sensitive; do
  timed nxsub_diamond_$m $REPS "$SEC/nxsub_diamond_$m.m8" -- $DMD blastx -q $TRS/nxsub_queries.fna -d $TRS/scop40.dmnd --threads $THREADS -k 0 -e 10 --$m --quiet -o /dev/stdout
done
for f in "$SEC"/*.m8; do
  n=$(basename "$f" .m8)
  case $n in nxsub_*) q=$TRS/nxsub_queries.fna ;; *) q=$TRS/pnsub_queries.fa ;; esac
  $PY "$HERE/eval_scop.py" "$f" "$LAB" --subset all --name "$n" --queries "$q" > "$SEC/$n.all.metrics.json"
done
