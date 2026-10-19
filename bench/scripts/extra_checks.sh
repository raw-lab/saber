#!/usr/bin/env bash
# Checks run after the suite (results in bench/results/checks/):
#  1. two-domain queries with database E-values (`--evalue-model database`),
#     to separate the effect of the E-value model from the alignments found;
#  2. few queries against a large database: 20 long proteins (> 1,500 aa) in
#     one batch and one per batch (SABER 1.3 and 1.2.0), and three of them as
#     single-query searches of prebuilt indexes (SABER 1.3, SABER 1.2.0 with
#     its postings index, DIAMOND --very-sensitive).
#
#   V13=... V12=... DMD=... DATA=... DER=... OUT=<run_all.sh output> CHK=<dir> bench/scripts/extra_checks.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V13=${V13:-/home/claude/bench/binaries/saber-1.3.0}
V12=${V12:-/home/claude/bench/binaries/saber-1.2.0-final}
DMD=${DMD:-/home/claude/tools/d228/diamond}
DATA=${DATA:-/home/claude/bench/data}
DER=${DER:-/home/claude/bench/derived}
OUT=${OUT:-/home/claude/bench/final13}
CHK=${CHK:-/home/claude/bench/checks}
THREADS=${THREADS:-2}
PY="python3 -I"
RB="$PY $HERE/runbench.py"
LAB=$DATA/scop40.labels.tsv
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
mkdir -p "$CHK"

# 1. two-domain queries, database E-values (single-domain reference searched the same way)
$RB --label scop_saber_v1.3_default_dbmodel --reps 1 --out $CHK/two_domain_database_evalues_times.jsonl --stdout $W/scop_dbmodel.m8 -- \
  $V13 -q $DATA/scop40.fa -d $DATA/scop40.fa -k 100000 -E 10 -t $THREADS --evalue-model database --quiet
$RB --label multi_saber_v1.3_default_dbmodel --reps 1 --out $CHK/two_domain_database_evalues_times.jsonl --stdout $W/multi_dbmodel.m8 -- \
  $V13 -q $DER/scop40.multi.fa -d $DATA/scop40.fa -k 100000 -E 10 -t $THREADS --evalue-model database --quiet
$PY $HERE/eval_multi.py $W/multi_dbmodel.m8 $W/scop_dbmodel.m8 $LAB --name saber_v1.3_default_dbmodel > $CHK/two_domain_database_evalues.metrics.json

# 2. few queries against the x4 database (36 M residues)
Q=$OUT/threads/long20.fa          # written by run_all.sh threads
DB=$DER/mm.db400.fa
J=$CHK/few_queries_large_database.jsonl
$RB --label long20_saber_v1.3_batch --reps 1 --out $J --stdout /dev/null -- $V13 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
$RB --label long20_saber_v1.2_onequery --reps 1 --out $J --stdout /dev/null -- $V12 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --query-batch 1 --quiet
$RB --label long20_saber_v1.2_batch --reps 1 --out $J --stdout /dev/null -- $V12 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
$V13 --makedb $W/db400.sdx -d $DB --quiet
$V12 --makedb $W/db400.v12.sdx -d $DB --quiet
$DMD makedb --in $DB -d $W/db400 --quiet
$PY - "$Q" "$W" <<'PY'
import sys
recs, h, s = [], None, []
for l in open(sys.argv[1]):
    l = l.rstrip()
    if l.startswith('>'):
        if h:
            recs.append((h, ''.join(s)))
        h, s = l, []
    else:
        s.append(l)
if h:
    recs.append((h, ''.join(s)))
for i, (h, s) in enumerate(recs[:3]):
    open(f"{sys.argv[2]}/one{i}.fa", "w").write(f"{h}\n{s}\n")
PY
for i in 0 1 2; do
  $RB --label one${i}_saber_v1.3_db --reps 1 --out $J --stdout /dev/null -- $V13 -q $W/one$i.fa --db $W/db400.sdx -k 100000 -E 1e-3 -t $THREADS --quiet
  $RB --label one${i}_saber_v1.2_db --reps 1 --out $J --stdout /dev/null -- $V12 -q $W/one$i.fa --db $W/db400.v12.sdx -k 100000 -E 1e-3 -t $THREADS --quiet
  $RB --label one${i}_diamond_very-sensitive --reps 1 --out $J --stdout /dev/null -- $DMD blastp -q $W/one$i.fa -d $W/db400.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/null
done
$V13 -q $W/one0.fa --db $W/db400.sdx -k 100000 -E 1e-3 -t $THREADS --quiet --stats $CHK/single_query_stats.json > /dev/null
rm -f $W/db400.sdx $W/db400.v12.sdx $W/db400.dmnd

# 3. the same single queries against an index with seed postings (3.4 GB):
#    query-major schedule vs the scan schedule (same alignments)
P=$CHK/single_query_postings_index.jsonl
$RB --label makedb_saber_v1.3_db400_postings --reps 1 --out $P --stdout /dev/null -- \
  $V13 --makedb $W/db400.post.sdx -d $DB --partition-residues 4000000 --index-postings --quiet
for i in 0 1 2; do
  $RB --label one${i}_saber_v1.3_db_postings_query --reps 1 --out $P --stdout $W/q$i.m8 -- \
    $V13 -q $W/one$i.fa --db $W/db400.post.sdx -k 100000 -E 1e-3 -t $THREADS --batch-order query --quiet
  $RB --label one${i}_saber_v1.3_db_postings_scan --reps 1 --out $P --stdout $W/s$i.m8 -- \
    $V13 -q $W/one$i.fa --db $W/db400.post.sdx -k 100000 -E 1e-3 -t $THREADS --quiet
  cmp -s $W/q$i.m8 $W/s$i.m8 || echo "one$i: query-major and scan outputs differ" >&2
done
