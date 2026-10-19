#!/usr/bin/env bash
# Re-measure the SABER-vs-DIAMOND timing comparisons of the A1 sections in
# one later session, so that both tools are timed on the same machine at the
# same time. (In the A1 sections the DIAMOND rows were measured 14:34-17:17
# and the SABER 1.2 rows 20:10-20:37; the machine is identified only by its
# CPU model.) Same commands as run_all.sh; labels are SECTION:LABEL of the
# record being re-checked. Run with V12 = the delivered build: its output equals
# the benchmarked build's everywhere except the long-protein set, where the
# traceback fallback fixed after the main run changes 9 of 4,712 alignments.
#
#   bench/scripts/recheck.sh        (same variables as run_all.sh)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V12=${V12:-/home/claude/work/target175/release/saber}
DMD=${DMD:-/home/claude/tools/d228/diamond}
DATA=${DATA:-/home/claude/bench/data}
OUT=${OUT:-/home/claude/bench/final}
THREADS=${THREADS:-2}
REPS=${REPS:-3}
SEC=$OUT/recheck
mkdir -p "$SEC"
rm -f "$SEC/times.jsonl"
RB="python3 -I $HERE/runbench.py"
F=$DATA/scop40.fa
Q=$DATA/mmseqs_QUERY.fasta
IX=$OUT/index
t() { # LABEL OUTFILE -- cmd...
  local label=$1 out=$2; shift 3
  $RB --label "$label" --reps "$REPS" --out "$SEC/times.jsonl" --stdout "$out" -- "$@"
}
# SCOP40 all-vs-all, general tier
t scop:diamond_very-sensitive "$SEC/scop_diamond_very-sensitive.m8" -- \
  $DMD blastp -q $F -d $OUT/scop/scop40.dmnd --threads $THREADS -k 0 -e 10 --very-sensitive --quiet -o /dev/stdout
t scop:saber_v1.2_default "$SEC/scop_saber_v1.2_default.m8" -- \
  $V12 -q $F -d $F -k 100000 -E 10 -t $THREADS --preset default --quiet --stats $SEC/scop_saber_v1.2_default.stats.json
# MMseqs2 example: in-memory index, prebuilt indexes, x4 database
t mmseqs:diamond_very-sensitive "$SEC/mmseqs_diamond_very-sensitive.m8" -- \
  $DMD blastp -q $Q -d $OUT/mmseqs/db.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/stdout
t mmseqs:saber_v1.2_default "$SEC/mmseqs_saber_v1.2_default.m8" -- \
  $V12 -q $Q -d $DATA/mmseqs_DB.fasta -k 100000 -E 1e-3 -t $THREADS --preset default --quiet
t mmseqs:saber_v1.2_default_db "$SEC/mmseqs_saber_v1.2_default_db.m8" -- \
  $V12 -q $Q --db $IX/mmseqs_DB.sdx -k 100000 -E 1e-3 -t $THREADS --quiet
t index:db_warm_diamond_very-sensitive_db400 /dev/null -- \
  $DMD blastp -q $Q -d $IX/db400.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/null
t index:db_warm_saber_v1.2_db400 /dev/null -- \
  $V12 -q $Q --db $IX/db400.sdx -k 100000 -E 1e-3 -t $THREADS --quiet
# long proteins (queries from the mmseqs section)
LQ=$OUT/mmseqs/long_queries.fa
t mmseqs:long_diamond_very-sensitive "$SEC/mmseqs_long_diamond_very-sensitive.m8" -- \
  $DMD blastp -q $LQ -d $OUT/mmseqs/db.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/stdout
t mmseqs:long_saber_v1.2_default "$SEC/mmseqs_long_saber_v1.2_default.m8" -- \
  $V12 -q $LQ -d $DATA/mmseqs_DB.fasta -k 100000 -E 1e-3 -t $THREADS --quiet --stats $SEC/mmseqs_long_saber_v1.2_default.stats.json
# outputs must equal the ones evaluated in the sections
: > "$SEC/outputs_check.txt"
for p in scop/diamond_very-sensitive scop/saber_v1.2_default mmseqs/diamond_very-sensitive mmseqs/saber_v1.2_default mmseqs/saber_v1.2_default_db mmseqs/long_diamond_very-sensitive mmseqs/long_saber_v1.2_default; do
  f="$SEC/${p/\//_}.m8"
  if cmp -s "$f" "$OUT/$p.m8"; then r=identical; else r=DIFFERENT; fi
  echo "$p: output $r" >> "$SEC/outputs_check.txt"
done
