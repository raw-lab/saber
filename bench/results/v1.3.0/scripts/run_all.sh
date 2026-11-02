#!/usr/bin/env bash
# Reproducible SABER benchmark suite: SABER 1.3.0 vs 1.2.0 vs DIAMOND.
#
#   bench/scripts/run_all.sh [section ...]
#   sections (default order): env scale scop decoy mmseqs index translated
#   lowc multi threads sweep scale_extra scale_algo oneoff (the last four
#   need the scale section's indexes or database)
#
# The 1.2.0 suite (with SABER 1.0.2 / 1.1.0) and its results are kept in
# bench/results/v1.2.0/; this suite compares the 1.3.0 presets with DIAMOND's
# modes, and keeps a few SABER 1.2.0 rows as before/after references.
#
# Configuration (environment variables, defaults are the paths used for the
# published results):
#   V13               SABER 1.3.0 binary (this release)
#   V12               SABER 1.2.0 binary (previous release)
#   DMD               DIAMOND binary (2.2.8)
#   DATA              directory with scop40.fa, scop40.labels.tsv, mmseqs_*.fasta
#   DER               derived datasets (bench/scripts/make_datasets.py DATA DER;
#                     the scale section adds DER/scop40.scale.fa)
#   OUT               results directory
#   THREADS           threads for multicore runs (default 2)
#   REPS              repetitions for timed runs (default 3)
#   KEYREPS           repetitions of the interleaved headline pairs of the
#                     scale and oneoff sections (default 5)
#
# Every timed run goes through runbench.py, which records wall/CPU time and the
# child's peak RSS (getrusage) to OUT/<section>/times.jsonl. Search outputs
# are kept for evaluation (m8 = BLAST tabular, same columns for every tool).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V13=${V13:-/home/claude/bench/binaries/saber-1.3.0}
V12=${V12:-/home/claude/bench/binaries/saber-1.2.0-final}
DMD=${DMD:-/home/claude/tools/d228/diamond}
DATA=${DATA:-/home/claude/bench/data}
DER=${DER:-/home/claude/bench/derived}
OUT=${OUT:-/home/claude/bench/final13}
THREADS=${THREADS:-2}
REPS=${REPS:-3}
KEYREPS=${KEYREPS:-5}
PY="python3 -I"
RB="$PY $HERE/runbench.py"
EVAL="$PY $HERE/eval_scop.py"
LAB=$DATA/scop40.labels.tsv
mkdir -p "$OUT"

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$OUT/run.log"; }

# timed LABEL REPS OUTFILE -- cmd...   (OUTFILE receives stdout; '-' = none)
timed() {
  local label=$1 reps=$2 out=$3; shift 4
  local so=/dev/null
  [ "$out" != "-" ] && so=$out
  $RB --label "$label" --reps "$reps" --out "$SEC/times.jsonl" --stdout "$so" -- "$@" 2>>"$SEC/runbench.err" || log "FAILED: $label (exit $?)"
}

begin() { SEC=$OUT/$1; mkdir -p "$SEC"; }

evict() { # drop a file's pages from the page cache (cold-cache runs; no root needed)
  local f
  for f in "$@"; do $PY -c "import os,sys; fd=os.open(sys.argv[1],os.O_RDONLY); os.posix_fadvise(fd,0,0,os.POSIX_FADV_DONTNEED); os.close(fd)" "$f"; done
}

subset() { # subset IN OUT EVERY : keep every N-th FASTA record
  $PY -c "
import sys
recs=[];h=None;s=[]
for l in open(sys.argv[1]):
    l=l.rstrip()
    if l.startswith('>'):
        if h: recs.append((h,''.join(s)))
        h=l;s=[]
    else: s.append(l)
if h: recs.append((h,''.join(s)))
n=int(sys.argv[3])
with open(sys.argv[2],'w') as f:
    for i,(h,q) in enumerate(recs):
        if i%n==0: f.write(h+'\\n'+q+'\\n')
" "$1" "$2" "$3"
}

dflag() { [ "$1" = default ] && echo "" || echo "--$1"; }

sec_env() {
  SEC=$OUT/env; mkdir -p $SEC
  {
    echo "date: $(date -Is)"
    echo "host: $(uname -a)"
    echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
    echo "cores: $(nproc)"
    echo "mem: $(grep MemTotal /proc/meminfo)"
    echo "transparent huge pages: $(cat /sys/kernel/mm/transparent_hugepage/enabled 2>/dev/null)"
    echo "rustc (SABER builds): $(/usr/bin/rustc --version)"
    echo "flags: cargo build --release (opt-level=3, lto=thin, codegen-units=1), no target-cpu flags"
    echo "diamond: $($DMD version 2>&1 | head -1) sha256 $(sha256sum $DMD | cut -d' ' -f1)"
    for b in $V12 $V13; do echo "$b: $($b --version 2>&1 | tail -1) sha256 $(sha256sum $b | cut -d' ' -f1)"; done
    echo "datasets:"; sha256sum $DATA/* | sed 's/^/  /'
    echo "derived:"; sha256sum $DER/* | sed 's/^/  /'
  } > $SEC/environment.txt
  cat $SEC/environment.txt
}

# --------------------------------------------------------------- at scale
# SCOP40 domains against themselves plus 80,000 reversed proteins (91,206
# sequences, 38.2 M residues): every hit to a REV_ decoy is a false positive.
sec_scale() {
  begin scale
  local Q=$DATA/scop40.fa DB=$DER/scop40.scale.fa
  [ -s $DB ] || $PY $HERE/make_scale_db.py $DATA/scop40.fa $DER/mm.db400.fa $DB
  sha256sum $DB > $SEC/dataset.sha256
  log "scale: index construction"
  timed makedb_saber_v1.3 $REPS - -- $V13 --makedb $SEC/scale.sdx -d $DB --quiet
  timed makedb_diamond $REPS - -- $DMD makedb --in $DB -d $SEC/scale --threads $THREADS --quiet
  timed makedb_saber_v1.2 1 - -- $V12 --makedb $SEC/scale.v12.sdx -d $DB --quiet
  ls -l $SEC/scale.sdx $SEC/scale.dmnd $SEC/scale.v12.sdx > $SEC/sizes.txt
  log "scale: SABER 1.3 default vs DIAMOND --very-sensitive, interleaved"
  for r in $(seq $KEYREPS); do
    timed saber_v1.3_default 1 $SEC/saber_v1.3_default.m8 -- $V13 -q $Q --db $SEC/scale.sdx -k 100000 -E 10 -t $THREADS --quiet --stats $SEC/saber_v1.3_default.stats.json
    timed diamond_very-sensitive 1 $SEC/diamond_very-sensitive.m8 -- $DMD blastp -q $Q -d $SEC/scale.dmnd --threads $THREADS -k 0 -e 10 --very-sensitive --quiet -o /dev/stdout
  done
  for p in very-fast fast sensitive more-sensitive; do
    log "scale: v1.3 $p"
    timed saber_v1.3_$p $REPS $SEC/saber_v1.3_$p.m8 -- $V13 -q $Q --db $SEC/scale.sdx -k 100000 -E 10 -t $THREADS --preset $p --quiet --stats $SEC/saber_v1.3_$p.stats.json
  done
  for m in fast default mid-sensitive sensitive more-sensitive; do
    log "scale: diamond $m"
    timed diamond_$m $REPS $SEC/diamond_$m.m8 -- $DMD blastp -q $Q -d $SEC/scale.dmnd --threads $THREADS -k 0 -e 10 $(dflag $m) --quiet -o /dev/stdout
  done
  log "scale: diamond ultra-sensitive (1 rep)"
  timed diamond_ultra-sensitive 1 $SEC/diamond_ultra-sensitive.m8 -- $DMD blastp -q $Q -d $SEC/scale.dmnd --threads $THREADS -k 0 -e 10 --ultra-sensitive --quiet -o /dev/stdout
  log "scale: v1.2.0 default (1 rep)"
  timed saber_v1.2_default 1 $SEC/saber_v1.2_default.m8 -- $V12 -q $Q --db $SEC/scale.v12.sdx -k 100000 -E 10 -t $THREADS --quiet
  rm -f $SEC/scale.v12.sdx
  # practical settings: E <= 1e-3, 25 targets per query
  timed practical_saber_v1.3_default $REPS - -- $V13 -q $Q --db $SEC/scale.sdx -k 25 -E 1e-3 -t $THREADS --quiet
  timed practical_diamond_very-sensitive $REPS - -- $DMD blastp -q $Q -d $SEC/scale.dmnd --threads $THREADS -k 25 -e 1e-3 --very-sensitive --quiet -o /dev/null
  scale_eval
}

# Checks on the scale section's comparison: DIAMOND --very-sensitive with one
# seed-index chunk (-c1, DIAMOND's speed/memory knob) and SABER 1.3 default
# with the 1.2 E-value model (what the per-subject E-values contribute).
sec_scale_extra() {
  begin scale
  local Q=$DATA/scop40.fa
  [ -s $SEC/scale.sdx ] && [ -s $SEC/scale.dmnd ] || { log "scale_extra: run the scale section first"; return 1; }
  timed diamond_very-sensitive_c1 $REPS $SEC/diamond_very-sensitive_c1.m8 -- $DMD blastp -q $Q -d $SEC/scale.dmnd --threads $THREADS -k 0 -e 10 --very-sensitive -c1 --quiet -o /dev/stdout
  timed saber_v1.3_default_dbmodel 1 $SEC/saber_v1.3_default_dbmodel.m8 -- $V13 -q $Q --db $SEC/scale.sdx -k 100000 -E 10 -t $THREADS --evalue-model database --quiet
  scale_eval
}

# DIAMOND chooses its seed-search algorithm automatically (double-indexed or
# query-indexed, --algo 0/1). Check that the forced alternative is not faster
# than DIAMOND's choice: --algo 1 alternating with SABER 1.3 default, outputs
# evaluated like the scale section's.
sec_scale_algo() {
  begin scale_algo
  local Q=$DATA/scop40.fa IDX=$OUT/scale
  [ -s $IDX/scale.sdx ] && [ -s $IDX/scale.dmnd ] || { log "scale_algo: run the scale section first"; return 1; }
  for r in $(seq $REPS); do
    timed saber_v1.3_default 1 - -- $V13 -q $Q --db $IDX/scale.sdx -k 100000 -E 10 -t $THREADS --quiet
    timed diamond_very-sensitive_algo1 1 $SEC/diamond_very-sensitive_algo1.m8 -- $DMD blastp -q $Q -d $IDX/scale.dmnd --threads $THREADS -k 0 -e 10 --very-sensitive --algo 1 --quiet -o /dev/stdout
  done
  for sub in tune test all; do
    $EVAL $SEC/diamond_very-sensitive_algo1.m8 $LAB --subset $sub --name diamond_very-sensitive_algo1 --decoy-prefix REV_ --decoys 80000 > $SEC/diamond_very-sensitive_algo1.$sub.metrics.json
  done
}

# One-off search at scale, index construction included: SABER builds its
# index in memory (-d FASTA, no preparation step); DIAMOND runs makedb and
# then blastp (one timed shell command). Alternating, KEYREPS times each.
sec_oneoff() {
  begin oneoff
  local Q=$DATA/scop40.fa DB=$DER/scop40.scale.fa
  [ -s $DB ] || { log "oneoff: run the scale section first"; return 1; }
  log "oneoff: SABER 1.3 default on the FASTA vs DIAMOND makedb + blastp --very-sensitive, interleaved"
  for r in $(seq $KEYREPS); do
    timed saber_v1.3_default_fasta 1 $SEC/saber_v1.3_default_fasta.m8 -- $V13 -q $Q -d $DB -k 100000 -E 10 -t $THREADS --quiet --stats $SEC/saber_v1.3_default_fasta.stats.json
    timed diamond_very-sensitive_makedb_blastp 1 - -- bash -c "$DMD makedb --in $DB -d $SEC/oneoff --threads $THREADS --quiet && $DMD blastp -q $Q -d $SEC/oneoff.dmnd --threads $THREADS -k 0 -e 10 --very-sensitive --quiet -o /dev/null"
  done
  if cmp -s $SEC/saber_v1.3_default_fasta.m8 $OUT/scale/saber_v1.3_default.m8; then echo "in-memory index output == prebuilt index output" > $SEC/outputs_check.txt; else echo "in-memory and prebuilt index outputs DIFFER" > $SEC/outputs_check.txt; fi
  rm -f $SEC/oneoff.dmnd
}

scale_eval() {
  SEC=$OUT/scale
  for f in $SEC/*.m8; do
    n=$(basename $f .m8)
    for sub in tune test all; do
      $EVAL $f $LAB --subset $sub --name $n --decoy-prefix REV_ --decoys 80000 --roc $SEC/$n.$sub.roc.tsv > $SEC/$n.$sub.metrics.json
    done
  done
}

# --------------------------------------------------------------- SCOP40
sec_scop() {
  begin scop
  local Q=$DATA/scop40.fa
  log "scop: index construction"
  timed diamond_makedb $REPS - -- $DMD makedb --in $Q -d $SEC/scop40 --threads $THREADS
  timed saber_v1.3_makedb_default $REPS - -- $V13 --makedb $SEC/scop40_default.sdx -d $Q --quiet
  for p in very-fast fast default sensitive more-sensitive legacy; do
    log "scop: v1.3 $p"
    timed saber_v1.3_$p $REPS $SEC/saber_v1.3_$p.m8 -- $V13 -q $Q -d $Q -k 100000 -E 10 -t $THREADS --preset $p --quiet --stats $SEC/saber_v1.3_$p.stats.json
  done
  log "scop: v1.3 default with a persistent index (--db, warm)"
  timed saber_v1.3_default_db_warm $REPS $SEC/saber_v1.3_default_db.m8 -- $V13 -q $Q --db $SEC/scop40_default.sdx -k 100000 -E 10 -t $THREADS --quiet
  for m in fast default mid-sensitive sensitive more-sensitive very-sensitive ultra-sensitive; do
    log "scop: diamond $m"
    timed diamond_$m $REPS $SEC/diamond_$m.m8 -- $DMD blastp -q $Q -d $SEC/scop40.dmnd --threads $THREADS -k 0 -e 10 $(dflag $m) --quiet -o /dev/stdout
  done
  log "scop: v1.2.0"
  for p in default more-sensitive; do
    timed saber_v1.2_$p $REPS $SEC/saber_v1.2_$p.m8 -- $V12 -q $Q -d $Q -k 100000 -E 10 -t $THREADS --preset $p --quiet
  done
  # practical settings (E <= 1e-3, 25 targets) for the main configurations
  for p in fast default; do
    timed practical_saber_v1.3_$p $REPS - -- $V13 -q $Q -d $Q -k 25 -E 1e-3 -t $THREADS --preset $p --quiet
  done
  for m in default very-sensitive ultra-sensitive; do
    timed practical_diamond_$m $REPS - -- $DMD blastp -q $Q -d $SEC/scop40.dmnd --threads $THREADS -k 25 -e 1e-3 $(dflag $m) --quiet -o /dev/null
  done
  log "scop: v1.3 exhaustive (1 rep)"
  timed saber_v1.3_exhaustive 1 $SEC/saber_v1.3_exhaustive.m8 -- $V13 -q $Q -d $Q -k 100000 -E 10 -t $THREADS --preset exhaustive --quiet --stats $SEC/saber_v1.3_exhaustive.stats.json
  scop_eval
}

scop_eval() {
  SEC=$OUT/scop
  for f in $SEC/*.m8; do
    n=$(basename $f .m8)
    for sub in tune test all; do
      $EVAL $f $LAB --subset $sub --name $n --roc $SEC/$n.$sub.roc.tsv > $SEC/$n.$sub.metrics.json
    done
  done
}

# --------------------------------------------------------------- shuffled decoys
sec_decoy() {
  begin decoy
  local Q=$DATA/scop40.fa
  for d in shuf wshuf10; do
    local DB=$DER/scop40.$d.fa
    $DMD makedb --in $DB -d $SEC/$d --quiet
    for p in fast default more-sensitive; do
      timed saber_v1.3_${p}_$d 1 $SEC/saber_v1.3_${p}_$d.m8 -- $V13 -q $Q -d $DB -k 100000 -E 10 -t $THREADS --preset $p --quiet
    done
    timed saber_v1.3_default_dbmodel_$d 1 $SEC/saber_v1.3_default_dbmodel_$d.m8 -- $V13 -q $Q -d $DB -k 100000 -E 10 -t $THREADS --evalue-model database --quiet
    timed saber_v1.2_default_$d 1 $SEC/saber_v1.2_default_$d.m8 -- $V12 -q $Q -d $DB -k 100000 -E 10 -t $THREADS --quiet
    for m in default very-sensitive ultra-sensitive; do
      timed diamond_${m}_$d 1 $SEC/diamond_${m}_$d.m8 -- $DMD blastp -q $Q -d $SEC/$d.dmnd --threads $THREADS -k 0 -e 10 $(dflag $m) --quiet -o /dev/stdout
    done
  done
}

# --------------------------------------------------------------- MMseqs2 example
sec_mmseqs() {
  begin mmseqs
  local Q=$DATA/mmseqs_QUERY.fasta DB=$DATA/mmseqs_DB.fasta
  $DMD makedb --in $DB -d $SEC/db --quiet
  log "mmseqs: exhaustive reference"
  timed saber_v1.3_exhaustive 1 $SEC/saber_v1.3_exhaustive.m8 -- $V13 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --preset exhaustive --quiet
  for p in very-fast fast default sensitive more-sensitive; do
    timed saber_v1.3_$p $REPS $SEC/saber_v1.3_$p.m8 -- $V13 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --preset $p --quiet --stats $SEC/saber_v1.3_$p.stats.json
  done
  for m in fast default mid-sensitive sensitive very-sensitive ultra-sensitive; do
    timed diamond_$m $REPS $SEC/diamond_$m.m8 -- $DMD blastp -q $Q -d $SEC/db.dmnd --threads $THREADS -k 0 -e 1e-3 $(dflag $m) --quiet -o /dev/stdout
  done
  timed saber_v1.2_default $REPS $SEC/saber_v1.2_default.m8 -- $V12 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
  log "mmseqs: database-size scaling"
  for d in db25 db50 db400; do
    local D=$DER/mm.$d.fa
    $DMD makedb --in $D -d $SEC/$d --quiet
    timed scale_${d}_saber_v1.3_default 1 - -- $V13 -q $Q -d $D -k 100000 -E 1e-3 -t $THREADS --quiet
    timed scale_${d}_saber_v1.2_default 1 - -- $V12 -q $Q -d $D -k 100000 -E 1e-3 -t $THREADS --quiet
    timed scale_${d}_diamond_very-sensitive 1 - -- $DMD blastp -q $Q -d $SEC/$d.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/null
    timed scale_${d}_diamond_default 1 - -- $DMD blastp -q $Q -d $SEC/$d.dmnd --threads $THREADS -k 0 -e 1e-3 --quiet -o /dev/null
  done
  log "mmseqs: query-count scaling"
  for n in 50 100 250; do
    local QQ=$DER/mm.q$n.fa
    timed scaleq_${n}_saber_v1.3_default 1 - -- $V13 -q $QQ -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
    timed scaleq_${n}_diamond_very-sensitive 1 - -- $DMD blastp -q $QQ -d $SEC/db.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/null
  done
  log "mmseqs: long proteins (every 5th of the 494 sequences > 1,500 aa)"
  subset $DER/mm.long.fa $SEC/long_queries.fa 5
  timed long_saber_v1.3_default 1 $SEC/long_saber_v1.3_default.m8 -- $V13 -q $SEC/long_queries.fa -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet --stats $SEC/long_saber_v1.3_default.stats.json
  timed long_saber_v1.2_default 1 $SEC/long_saber_v1.2_default.m8 -- $V12 -q $SEC/long_queries.fa -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
  timed long_diamond_very-sensitive 1 $SEC/long_diamond_very-sensitive.m8 -- $DMD blastp -q $SEC/long_queries.fa -d $SEC/db.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/stdout
}

# --------------------------------------------------------------- index
sec_index() {
  begin index
  rm -f $SEC/sizes.txt
  local Q=$DATA/mmseqs_QUERY.fasta
  for d in mmseqs_DB:$DATA/mmseqs_DB.fasta db400:$DER/mm.db400.fa; do
    local n=${d%%:*} f=${d#*:}
    timed makedb_saber_v1.3_$n $REPS - -- $V13 --makedb $SEC/$n.sdx -d $f --quiet
    timed makedb_diamond_$n $REPS - -- $DMD makedb --in $f -d $SEC/$n --threads $THREADS --quiet
    timed makedb_saber_v1.2_$n 1 - -- $V12 --makedb $SEC/$n.v12.sdx -d $f --quiet
    ls -l $SEC/$n.sdx $SEC/$n.dmnd $SEC/$n.v12.sdx >> $SEC/sizes.txt || true
    rm -f $SEC/$n.v12.sdx
    timed db_warm_saber_v1.3_$n $REPS - -- $V13 -q $Q --db $SEC/$n.sdx -k 100000 -E 1e-3 -t $THREADS --quiet
    timed db_warm_diamond_very-sensitive_$n $REPS - -- $DMD blastp -q $Q -d $SEC/$n.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/null
    for r in 1 2 3; do
      evict $SEC/$n.sdx
      timed db_cold_saber_v1.3_$n 1 - -- $V13 -q $Q --db $SEC/$n.sdx -k 100000 -E 1e-3 -t $THREADS --quiet
      evict $SEC/$n.dmnd
      timed db_cold_diamond_very-sensitive_$n 1 - -- $DMD blastp -q $Q -d $SEC/$n.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/null
    done
    # an index with seed postings (several partitions): the scan schedule
    # and both index schedules must report the same alignments
    timed makedb_saber_v1.3_${n}_postings_p4M 1 - -- $V13 --makedb $SEC/$n.p4M.sdx -d $f --partition-residues 4000000 --index-postings --quiet
    ls -l $SEC/$n.p4M.sdx >> $SEC/sizes.txt
    timed db_order_scan_$n 1 $SEC/order_scan_$n.m8 -- $V13 -q $Q --db $SEC/$n.p4M.sdx -k 100000 -E 1e-3 -t $THREADS --batch-order scan --quiet
    timed db_order_query_$n 1 $SEC/order_query_$n.m8 -- $V13 -q $Q --db $SEC/$n.p4M.sdx -k 100000 -E 1e-3 -t $THREADS --batch-order query --quiet
    timed db_order_partition_$n 1 $SEC/order_partition_$n.m8 -- $V13 -q $Q --db $SEC/$n.p4M.sdx -k 100000 -E 1e-3 -t $THREADS --batch-order partition --quiet
    for o in query partition; do
      if cmp -s $SEC/order_scan_$n.m8 $SEC/order_${o}_$n.m8; then echo "$n: scan == $o-major output" >> $SEC/sizes.txt; else echo "$n: scan and $o-major outputs DIFFER" >> $SEC/sizes.txt; fi
    done
    rm -f $SEC/$n.p4M.sdx
  done
}

# --------------------------------------------------------------- translated
sec_translated() {
  begin translated
  local NQ=$DER/scop40.bt.fna P=$DATA/scop40.fa
  $DMD makedb --in $P -d $SEC/scop40 --quiet
  for p in fast default more-sensitive; do
    timed nx_saber_v1.3_$p 1 $SEC/nx_saber_v1.3_$p.m8 -- $V13 --mode nx -q $NQ -d $P -k 100000 -E 10 -t $THREADS --preset $p --quiet --stats $SEC/nx_saber_v1.3_$p.stats.json
  done
  timed nx_saber_v1.2_default 1 $SEC/nx_saber_v1.2_default.m8 -- $V12 --mode nx -q $NQ -d $P -k 100000 -E 10 -t $THREADS --quiet
  for m in default very-sensitive ultra-sensitive; do
    timed nx_diamond_$m 1 $SEC/nx_diamond_$m.m8 -- $DMD blastx -q $NQ -d $SEC/scop40.dmnd --threads $THREADS -k 0 -e 10 $(dflag $m) --quiet -o /dev/stdout
  done
  # pn (tblastn): a fixed 1-in-37 subset of SCOP40 proteins; DIAMOND has no tblastn mode.
  subset $P $SEC/pn_queries.fa 37
  for p in default more-sensitive; do
    timed pn_saber_v1.3_$p 1 $SEC/pn_saber_v1.3_$p.m8 -- $V13 --mode pn -q $SEC/pn_queries.fa -d $NQ -k 100000 -E 10 -t $THREADS --preset $p --quiet
  done
  timed pn_saber_v1.2_default 1 $SEC/pn_saber_v1.2_default.m8 -- $V12 --mode pn -q $SEC/pn_queries.fa -d $NQ -k 100000 -E 10 -t $THREADS --quiet
  for f in $SEC/*.m8; do
    n=$(basename $f .m8)
    local qopt=""
    case $n in
      pn_*) qopt="--queries $SEC/pn_queries.fa" ;;
    esac
    for sub in tune test all; do $EVAL $f $LAB --subset $sub --name $n $qopt > $SEC/$n.$sub.metrics.json; done
  done
}

# --------------------------------------------------------------- low complexity
sec_lowc() {
  begin lowc
  local Q=$DER/lowc.query.fa DB=$DER/lowc.db.fa
  $DMD makedb --in $DB -d $SEC/lowc --quiet
  timed lowc_saber_v1.3_default $REPS $SEC/lowc_saber_v1.3_default.m8 -- $V13 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet --stats $SEC/lowc_saber_v1.3_default.stats.json
  timed lowc_saber_v1.3_default_nomask 1 $SEC/lowc_saber_v1.3_default_nomask.m8 -- $V13 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --masking none --comp-bias off --quiet --stats $SEC/lowc_saber_v1.3_default_nomask.stats.json
  # adversarial: no masking, no frequency cap, every repeat hit enumerated --
  # memory must stay bounded (scan schedule: cuts by chunks, subjects and
  # queries, stops at subject boundaries, repeat-rich pairs scanned on their
  # own in slices). Every 10th query (10 queries) with the legacy preset (one
  # seed shape, as in the 1.2.0 suite); the first 4 of them with the 16-shape
  # default preset (all 10 enumerate tens of billions of seed hits and did not
  # finish within 30 minutes in development, at 0.69 GB).
  subset $Q $SEC/lowc10.fa 10
  head -8 $SEC/lowc10.fa > $SEC/lowc4.fa
  timed lowc4_saber_v1.3_default_nomask_nocap 1 $SEC/lowc4_saber_v1.3_default_nomask_nocap.m8 -- timeout 1200 $V13 -q $SEC/lowc4.fa -d $DB -k 100000 -E 1e-3 -t $THREADS --masking none --comp-bias off --seed-max-occ 1000000000 --quiet --stats $SEC/lowc4_saber_v1.3_default_nomask_nocap.stats.json
  timed lowc10_saber_v1.3_legacy_nomask_nocap 1 $SEC/lowc10_saber_v1.3_legacy_nomask_nocap.m8 -- timeout 1200 $V13 -q $SEC/lowc10.fa -d $DB -k 100000 -E 1e-3 -t $THREADS --preset legacy --masking none --comp-bias off --seed-max-occ 1000000000 --quiet --stats $SEC/lowc10_saber_v1.3_legacy_nomask_nocap.stats.json
  timed lowc_diamond_default $REPS $SEC/lowc_diamond_default.m8 -- $DMD blastp -q $Q -d $SEC/lowc.dmnd --threads $THREADS -k 0 -e 1e-3 --quiet -o /dev/stdout
  timed lowc_diamond_very-sensitive $REPS $SEC/lowc_diamond_very-sensitive.m8 -- $DMD blastp -q $Q -d $SEC/lowc.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/stdout
  # real queries against a database polluted with low-complexity sequences
  local MQ=$DATA/mmseqs_QUERY.fasta
  timed mmq_saber_v1.3_default 1 $SEC/mmq_saber_v1.3_default.m8 -- $V13 -q $MQ -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
  timed mmq_diamond_very-sensitive 1 $SEC/mmq_diamond_very-sensitive.m8 -- $DMD blastp -q $MQ -d $SEC/lowc.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/stdout
}

# --------------------------------------------------------------- two-domain queries, fragments
sec_multi() {
  begin multi
  local P=$DATA/scop40.fa
  $DMD makedb --in $P -d $SEC/scop40 --quiet
  for set in multi frag; do
    local Q=$DER/scop40.$set.fa
    for p in default more-sensitive; do
      timed ${set}_saber_v1.3_$p 1 $SEC/${set}_saber_v1.3_$p.m8 -- $V13 -q $Q -d $P -k 100000 -E 10 -t $THREADS --preset $p --quiet
    done
    for m in very-sensitive ultra-sensitive; do
      timed ${set}_diamond_$m 1 $SEC/${set}_diamond_$m.m8 -- $DMD blastp -q $Q -d $SEC/scop40.dmnd --threads $THREADS -k 0 -e 10 --$m --quiet -o /dev/stdout
    done
  done
  for a in saber_v1.3_default saber_v1.3_more-sensitive diamond_very-sensitive diamond_ultra-sensitive; do
    $PY $HERE/eval_multi.py $SEC/multi_$a.m8 $OUT/scop/$a.m8 $LAB --name $a > $SEC/multi_$a.metrics.json
    for sub in tune test all; do $EVAL $SEC/frag_$a.m8 $LAB --subset $sub --name frag_$a > $SEC/frag_$a.$sub.metrics.json; done
  done
}

# --------------------------------------------------------------- threads
sec_threads() {
  begin threads
  local Q=$DATA/scop40.fa
  $DMD makedb --in $Q -d $SEC/scop40 --quiet
  subset $DER/mm.long.fa $SEC/long20.fa 25
  for t in 1 2; do
    timed t${t}_saber_v1.3_default $REPS - -- $V13 -q $Q -d $Q -k 100000 -E 10 -t $t --quiet
    timed t${t}_diamond_very-sensitive $REPS - -- $DMD blastp -q $Q -d $SEC/scop40.dmnd --threads $t -k 0 -e 10 --very-sensitive --quiet -o /dev/null
    # one query per batch: the database is split across the threads
    timed t${t}_saber_v1.3_onequery $REPS - -- $V13 -q $SEC/long20.fa -d $DER/mm.db400.fa -k 100000 -E 1e-3 -t $t --query-batch 1 --quiet
  done
}

# --------------------------------------------------------------- preset selection
# Configurations considered for the 1.3 presets, on the tune splits of both
# benchmarks (the test splits are never used for choosing). Needs the scale
# section's index. Common arguments: Murphy-10 seeds, window filter with exact
# identities weighted 2, scoring band +-48.
DS="11110111,11100100111,110010101011,11010001001011,10101100001101,110100100100011,1010010100010011,1100101000001011,11100000100010101,11000100010010011,11010000001000111,110001001000010011,1010001000100001011,1100010100000010011,1100100000101000011,1101000010000001011"
DS8="11110111,11100100111,110010101011,11010001001011,10101100001101,110100100100011,1010010100010011,1100101000001011"
DS4="11110111,11100100111,110010101011,11010001001011"
W8="100111011,1011100011,1010011011,11110101,1101000111,1111100001,10101111,1010010111"
sec_sweep() {
  begin sweep
  local Q=$DATA/scop40.fa IDX=$OUT/scale/scale.sdx
  [ -s $IDX ] || { log "sweep: run the scale section first"; return 1; }
  local base="--seed-alphabet murphy10 --hit-filter-exact-weight 2 --score-band 48"
  while read -r name shapes T ff pf; do
    [ -z "$name" ] && continue
    local args="$base --seed-shape ${!shapes} --hit-filter-ids $T --hsp-floor-fraction $ff --score-prefilter-band $pf"
    timed sc_$name 1 $SEC/sc_$name.m8 -- $V13 -q $Q --db $IDX -k 100000 -E 10 -t $THREADS --quiet $args
    $EVAL $SEC/sc_$name.m8 $LAB --subset tune --name sc_$name --decoy-prefix REV_ --decoys 80000 > $SEC/sc_$name.tune.metrics.json
    timed av_$name 1 $SEC/av_$name.m8 -- $V13 -q $Q -d $Q -k 100000 -E 10 -t $THREADS --quiet $args
    $EVAL $SEC/av_$name.m8 $LAB --subset tune --name av_$name > $SEC/av_$name.tune.metrics.json
    rm -f $SEC/sc_$name.m8 $SEC/av_$name.m8
  done <<'CONFIGS'
ds4_T28_f70 DS4 28 0.7 0
ds8_T28_f70 DS8 28 0.7 0
ds8_T27_f70 DS8 27 0.7 0
ds16_T28_f70 DS 28 0.7 0
ds16_T27_f70 DS 27 0.7 0
ds16_T27_f65_pf16 DS 27 0.65 16
ds16_T27_f60 DS 27 0.6 0
ds16_T26_f60_pf16 DS 26 0.6 16
ds16_T25_f55_pf16 DS 25 0.55 16
w8_T27_f65_pf16 W8 27 0.65 16
CONFIGS
  # DS16 followed by W8 (24 shapes)
  local args="$base --seed-shape $DS,$W8 --hit-filter-ids 27 --hsp-floor-fraction 0.65 --score-prefilter-band 16"
  timed sc_ds16w8_T27_f65_pf16 1 $SEC/sc_ds16w8_T27_f65_pf16.m8 -- $V13 -q $Q --db $IDX -k 100000 -E 10 -t $THREADS --quiet $args
  $EVAL $SEC/sc_ds16w8_T27_f65_pf16.m8 $LAB --subset tune --name sc_ds16w8_T27_f65_pf16 --decoy-prefix REV_ --decoys 80000 > $SEC/sc_ds16w8_T27_f65_pf16.tune.metrics.json
  timed av_ds16w8_T27_f65_pf16 1 $SEC/av_ds16w8_T27_f65_pf16.m8 -- $V13 -q $Q -d $Q -k 100000 -E 10 -t $THREADS --quiet $args
  $EVAL $SEC/av_ds16w8_T27_f65_pf16.m8 $LAB --subset tune --name av_ds16w8_T27_f65_pf16 > $SEC/av_ds16w8_T27_f65_pf16.tune.metrics.json
  rm -f $SEC/sc_ds16w8_T27_f65_pf16.m8 $SEC/av_ds16w8_T27_f65_pf16.m8
}

sections=("$@")
[ ${#sections[@]} -eq 0 ] && sections=(env scale scop decoy mmseqs index translated lowc multi threads sweep scale_extra scale_algo oneoff)
for s in "${sections[@]}"; do
  log "=== section $s ==="
  sec_$s
done
log "done"
