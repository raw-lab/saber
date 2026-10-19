#!/usr/bin/env bash
# Reproducible SABER benchmark suite: SABER 1.0.2 vs 1.1.0 vs 1.2.0 vs DIAMOND.
#
#   bench/scripts/run_all.sh [section ...]
#   sections (default order): env scop decoy mmseqs index extra translated lowc
#   multi threads sweep; 'rerun' re-measures SCOP40 runs on demand. 'extra'
#   needs the index section's mmseqs_DB.sdx.
#
# Configuration (environment variables, defaults are the paths used for the
# published results):
#   V102, V110, V12   SABER binaries (all built with the same rustc)
#   DMD               DIAMOND binary (2.2.8)
#   DATA              directory with scop40.fa, scop40.labels.tsv, mmseqs_*.fasta
#   DER               derived datasets (bench/scripts/make_datasets.py DATA DER)
#   OUT               results directory
#   THREADS           threads for multicore runs (default 2)
#   REPS              repetitions for timed runs (default 3)
#
# Every timed run goes through runbench.py, which records wall/CPU time and the
# child's peak RSS (getrusage) to OUT/<section>/times.jsonl. Search outputs
# are kept for evaluation (m8 = BLAST tabular, same columns for every tool).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V102=${V102:-/home/claude/work/target175_baseline_v102/release/saber}
V110=${V110:-/home/claude/work/target175_opt_v110/release/saber}
V12=${V12:-/home/claude/work/target175/release/saber}
DMD=${DMD:-/home/claude/tools/d228/diamond}
DATA=${DATA:-/home/claude/bench/data}
DER=${DER:-/home/claude/bench/derived}
OUT=${OUT:-/home/claude/bench/final}
THREADS=${THREADS:-2}
REPS=${REPS:-3}
PY="python3 -I"
RB="$PY $HERE/runbench.py"
EVAL="$PY $HERE/eval_scop.py"
LAB=$DATA/scop40.labels.tsv
mkdir -p "$OUT"

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$OUT/run.log"; }

# timed LABEL REPS OUTFILE -- cmd...   (OUTFILE receives stdout; '-' = none)
# With REUSE=1, baseline runs (labels naming diamond, v1.0.2 or v1.1.0) whose
# successful records (and output file) already exist from a previous
# invocation on the same machine (same CPU model and cache) are kept instead
# of re-measured (the previous records are moved to times.prev.jsonl when the
# section starts); every other run (SABER 1.2) is always re-measured.
timed() {
  local label=$1 reps=$2 out=$3; shift 4
  local so=/dev/null
  [ "$out" != "-" ] && so=$out
  if [ "${REUSE:-0}" = 1 ] && [[ "$label" == *diamond* || "$label" == *v1.0.2* || "$label" == *v1.1.0* ]] && [ -f "$SEC/times.prev.jsonl" ] \
     && { [ "$out" = "-" ] || [ -s "$out" ]; }; then
    # only records measured on this machine (runbench `cpu` field) are reused
    if $PY -c "
import json,sys
sys.path.insert(0, sys.argv[4]); import runbench
cpu=runbench.cpu_id()
recs=[json.loads(l) for l in open(sys.argv[1]) if l.strip()]
ok=[r for r in recs if r['label']==sys.argv[2] and r.get('exit')==0 and r.get('cpu')==cpu]
sys.exit(0 if len(ok)>=int(sys.argv[3]) else 1)" "$SEC/times.prev.jsonl" "$label" "$reps" "$HERE"; then
      $PY -c "
import json,sys
sys.path.insert(0, sys.argv[3]); import runbench
cpu=runbench.cpu_id()
seen=set()
for l in open(sys.argv[1]):
    if l.strip():
        r=json.loads(l)
        k=(r['label'], r.get('rep'), r['timestamp'])
        if r['label']==sys.argv[2] and r.get('cpu')==cpu and k not in seen:
            seen.add(k); print(l.strip())" "$SEC/times.prev.jsonl" "$label" "$HERE" >> "$SEC/times.jsonl"
      log "reused: $label"
      return 0
    fi
  fi
  $RB --label "$label" --reps "$reps" --out "$SEC/times.jsonl" --stdout "$so" -- "$@" 2>>"$SEC/runbench.err" || log "FAILED: $label (exit $?)"
}

# start a section directory (keeps previous timings for REUSE=1)
begin() {
  SEC=$OUT/$1; mkdir -p "$SEC"
  if [ "${REUSE:-0}" = 1 ] && [ -f "$SEC/times.jsonl" ]; then
    cat "$SEC/times.jsonl" >> "$SEC/times.prev.jsonl"; rm -f "$SEC/times.jsonl"
    # keep one copy of each record (reused records would otherwise accumulate)
    $PY -c "
import json,sys
seen=set(); out=[]
for l in open(sys.argv[1]):
    if l.strip():
        r=json.loads(l); k=(r['label'], r.get('rep'), r['timestamp'])
        if k not in seen: seen.add(k); out.append(l.strip())
open(sys.argv[1],'w').write(''.join(x+'\\n' for x in out))" "$SEC/times.prev.jsonl"
  fi
}

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


sec_env() {
  SEC=$OUT/env; mkdir -p $SEC
  {
    echo "date: $(date -Is)"
    echo "host: $(uname -a)"
    echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
    echo "cores: $(nproc)"
    echo "mem: $(grep MemTotal /proc/meminfo)"
    echo "rustc (SABER builds): $(/usr/bin/rustc --version)"
    echo "flags: cargo build --release (opt-level=3, lto=thin, codegen-units=1), no target-cpu flags"
    echo "diamond: $($DMD version 2>&1 | head -1)"
    for b in $V102 $V110 $V12; do echo "$b: $($b --version 2>&1 | tail -1) sha256 $(sha256sum $b | cut -c1-16)"; done
    echo "datasets:"; sha256sum $DATA/* | sed 's/^/  /'
    echo "derived:"; sha256sum $DER/* | sed 's/^/  /'
  } > $SEC/environment.txt
  cat $SEC/environment.txt
}

sec_scop() {
  begin scop
  local Q=$DATA/scop40.fa
  log "scop: diamond makedb"
  timed diamond_makedb $REPS - -- $DMD makedb --in $Q -d $SEC/scop40 --threads $THREADS
  log "scop: saber makedb"
  timed saber_v1.2_makedb_default $REPS - -- $V12 --makedb $SEC/scop40_default.sdx -d $Q --quiet
  for p in very-fast fast default sensitive more-sensitive legacy; do
    log "scop: v1.2 $p"
    timed saber_v1.2_$p $REPS $SEC/saber_v1.2_$p.m8 -- $V12 -q $Q -d $Q -k 100000 -E 10 -t $THREADS --preset $p --quiet --stats $SEC/saber_v1.2_$p.stats.json
  done
  log "scop: v1.2 default with persistent index (--db, warm)"
  timed saber_v1.2_default_db_warm $REPS $SEC/saber_v1.2_default_db.m8 -- $V12 -q $Q --db $SEC/scop40_default.sdx -k 100000 -E 10 -t $THREADS --quiet
  for m in fast default mid-sensitive sensitive more-sensitive very-sensitive ultra-sensitive; do
    log "scop: diamond $m"
    local flag="--$m"; [ $m = default ] && flag=""
    timed diamond_$m $REPS $SEC/diamond_$m.m8 -- $DMD blastp -q $Q -d $SEC/scop40.dmnd --threads $THREADS -k 0 -e 10 $flag --quiet -o /dev/stdout
  done
  log "scop: v1.1.0"
  timed saber_v1.1.0 $REPS $SEC/saber_v1.1.0.m8 -- $V110 -q $Q -d $Q -k 100000 -E 10 -t $THREADS --quiet
  log "scop: v1.0.2"
  timed saber_v1.0.2 $REPS $SEC/saber_v1.0.2.m8 -- $V102 -q $Q -d $Q -k 100000 -E 10 -t $THREADS --quiet
  # practical settings (E <= 1e-3, 25 targets) for the main configurations
  for p in fast default; do
    timed practical_saber_v1.2_$p $REPS - -- $V12 -q $Q -d $Q -k 25 -E 1e-3 -t $THREADS --preset $p --quiet
  done
  timed practical_saber_v1.1.0 1 - -- $V110 -q $Q -d $Q -k 25 -E 1e-3 -t $THREADS --quiet
  for m in default very-sensitive ultra-sensitive; do
    local flag="--$m"; [ $m = default ] && flag=""
    timed practical_diamond_$m $REPS - -- $DMD blastp -q $Q -d $SEC/scop40.dmnd --threads $THREADS -k 25 -e 1e-3 $flag --quiet -o /dev/null
  done
  log "scop: v1.2 exhaustive (1 rep)"
  timed saber_v1.2_exhaustive 1 $SEC/saber_v1.2_exhaustive.m8 -- $V12 -q $Q -d $Q -k 100000 -E 10 -t $THREADS --preset exhaustive --quiet --stats $SEC/saber_v1.2_exhaustive.stats.json
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

sec_decoy() {
  begin decoy
  local Q=$DATA/scop40.fa
  for d in shuf wshuf10; do
    local DB=$DER/scop40.$d.fa
    $DMD makedb --in $DB -d $SEC/$d --quiet
    timed saber_v1.1.0_$d 1 $SEC/saber_v1.1.0_$d.m8 -- $V110 -q $Q -d $DB -k 100000 -E 10 -t $THREADS --quiet
    for p in fast default more-sensitive; do
      timed saber_v1.2_${p}_$d 1 $SEC/saber_v1.2_${p}_$d.m8 -- $V12 -q $Q -d $DB -k 100000 -E 10 -t $THREADS --preset $p --quiet
    done
    for m in default very-sensitive ultra-sensitive; do
      local flag="--$m"; [ $m = default ] && flag=""
      timed diamond_${m}_$d 1 $SEC/diamond_${m}_$d.m8 -- $DMD blastp -q $Q -d $SEC/$d.dmnd --threads $THREADS -k 0 -e 10 $flag --quiet -o /dev/stdout
    done
  done
}

sec_mmseqs() {
  begin mmseqs
  local Q=$DATA/mmseqs_QUERY.fasta DB=$DATA/mmseqs_DB.fasta
  $DMD makedb --in $DB -d $SEC/db --quiet
  log "mmseqs: exhaustive reference"
  timed saber_v1.2_exhaustive 1 $SEC/saber_v1.2_exhaustive.m8 -- $V12 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --preset exhaustive --quiet
  for p in very-fast fast default sensitive; do
    timed saber_v1.2_$p $REPS $SEC/saber_v1.2_$p.m8 -- $V12 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --preset $p --quiet --stats $SEC/saber_v1.2_$p.stats.json
  done
  for m in fast default mid-sensitive sensitive very-sensitive; do
    local flag="--$m"; [ $m = default ] && flag=""
    timed diamond_$m $REPS $SEC/diamond_$m.m8 -- $DMD blastp -q $Q -d $SEC/db.dmnd --threads $THREADS -k 0 -e 1e-3 $flag --quiet -o /dev/stdout
  done
  timed saber_v1.1.0 $REPS $SEC/saber_v1.1.0.m8 -- $V110 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
  timed saber_v1.0.2 1 $SEC/saber_v1.0.2.m8 -- $V102 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
  log "mmseqs: database-size scaling"
  for d in db25 db50 db400; do
    local D=$DER/mm.$d.fa
    $DMD makedb --in $D -d $SEC/$d --quiet
    timed scale_${d}_saber_v1.1.0 1 - -- $V110 -q $Q -d $D -k 100000 -E 1e-3 -t $THREADS --quiet
    timed scale_${d}_saber_v1.2_default 1 - -- $V12 -q $Q -d $D -k 100000 -E 1e-3 -t $THREADS --quiet
    timed scale_${d}_diamond_very-sensitive 1 - -- $DMD blastp -q $Q -d $SEC/$d.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/null
    timed scale_${d}_diamond_default 1 - -- $DMD blastp -q $Q -d $SEC/$d.dmnd --threads $THREADS -k 0 -e 1e-3 --quiet -o /dev/null
  done
  log "mmseqs: query-count scaling"
  for n in 50 100 250; do
    local QQ=$DER/mm.q$n.fa
    timed scaleq_${n}_saber_v1.1.0 1 - -- $V110 -q $QQ -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
    timed scaleq_${n}_saber_v1.2_default 1 - -- $V12 -q $QQ -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
    timed scaleq_${n}_diamond_very-sensitive 1 - -- $DMD blastp -q $QQ -d $SEC/db.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/null
  done
  log "mmseqs: long proteins (every 5th of the 494 sequences > 1,500 aa)"
  subset $DER/mm.long.fa $SEC/long_queries.fa 5
  timed long_saber_v1.1.0 1 $SEC/long_saber_v1.1.0.m8 -- $V110 -q $SEC/long_queries.fa -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
  timed long_saber_v1.2_default 1 $SEC/long_saber_v1.2_default.m8 -- $V12 -q $SEC/long_queries.fa -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet --stats $SEC/long_saber_v1.2_default.stats.json
  timed long_diamond_very-sensitive 1 $SEC/long_diamond_very-sensitive.m8 -- $DMD blastp -q $SEC/long_queries.fa -d $SEC/db.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/stdout
}

sec_index() {
  begin index
  rm -f $SEC/sizes.txt
  local Q=$DATA/mmseqs_QUERY.fasta
  for d in mmseqs_DB:$DATA/mmseqs_DB.fasta db400:$DER/mm.db400.fa; do
    local n=${d%%:*} f=${d#*:}
    timed makedb_saber_v1.2_$n $REPS - -- $V12 --makedb $SEC/$n.sdx -d $f --quiet
    timed makedb_diamond_$n $REPS - -- $DMD makedb --in $f -d $SEC/$n --threads $THREADS --quiet
    timed makedb_saber_v1.1.0_$n 1 - -- $V110 --makedb $SEC/$n.v110.sdx -d $f --quiet
    [ -f $SEC/$n.v110.sdx ] || $V110 --makedb $SEC/$n.v110.sdx -d $f --quiet >/dev/null 2>&1 || true
    ls -l $SEC/$n.sdx $SEC/$n.dmnd $SEC/$n.v110.sdx >> $SEC/sizes.txt || true
    timed db_warm_saber_v1.2_$n $REPS - -- $V12 -q $Q --db $SEC/$n.sdx -k 100000 -E 1e-3 -t $THREADS --quiet
    timed db_warm_diamond_very-sensitive_$n $REPS - -- $DMD blastp -q $Q -d $SEC/$n.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/null
    for r in 1 2 3; do
      evict $SEC/$n.sdx
      timed db_cold_saber_v1.2_$n 1 - -- $V12 -q $Q --db $SEC/$n.sdx -k 100000 -E 1e-3 -t $THREADS --quiet
      evict $SEC/$n.dmnd
      timed db_cold_diamond_very-sensitive_$n 1 - -- $DMD blastp -q $Q -d $SEC/$n.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/null
    done
    # partition-major vs query-major on a multi-partition persistent index
    timed makedb_saber_v1.2_${n}_p16M 1 - -- $V12 --makedb $SEC/$n.p16M.sdx -d $f --partition-residues 4000000 --quiet
    timed db_order_query_$n 1 $SEC/order_query_$n.m8 -- $V12 -q $Q --db $SEC/$n.p16M.sdx -k 100000 -E 1e-3 -t $THREADS --batch-order query --quiet
    timed db_order_partition_$n 1 $SEC/order_partition_$n.m8 -- $V12 -q $Q --db $SEC/$n.p16M.sdx -k 100000 -E 1e-3 -t $THREADS --batch-order partition --quiet
    if cmp -s $SEC/order_query_$n.m8 $SEC/order_partition_$n.m8; then echo "$n: query-major == partition-major output" >> $SEC/sizes.txt; else echo "$n: query-major and partition-major outputs DIFFER" >> $SEC/sizes.txt; fi
    rm -f $SEC/$n.v110.sdx $SEC/$n.p16M.sdx
  done
}

sec_translated() {
  begin translated
  local NQ=$DER/scop40.bt.fna P=$DATA/scop40.fa
  $DMD makedb --in $P -d $SEC/scop40 --quiet
  for p in fast default more-sensitive; do
    timed nx_saber_v1.2_$p 1 $SEC/nx_saber_v1.2_$p.m8 -- $V12 --mode nx -q $NQ -d $P -k 100000 -E 10 -t $THREADS --preset $p --quiet --stats $SEC/nx_saber_v1.2_$p.stats.json
  done
  for m in default very-sensitive ultra-sensitive; do
    local flag="--$m"; [ $m = default ] && flag=""
    timed nx_diamond_$m 1 $SEC/nx_diamond_$m.m8 -- $DMD blastx -q $NQ -d $SEC/scop40.dmnd --threads $THREADS -k 0 -e 10 $flag --quiet -o /dev/stdout
  done
  # SABER 1.1.0 runs blastx and tblastn by full DP over all six frames (no
  # index; 7-14 s per query on 2 cores here): 1.1.0 is compared with 1.2.0
  # on a fixed 1-in-333 query subset (34 queries).
  subset $NQ $SEC/nxsub_queries.fna 333
  subset $P $SEC/pn_queries.fa 37
  subset $P $SEC/pnsub_queries.fa 333
  timed nxsub_saber_v1.1.0 1 $SEC/nxsub_saber_v1.1.0.m8 -- timeout 7200 $V110 --mode nx -q $SEC/nxsub_queries.fna -d $P -k 100000 -E 10 -t $THREADS --quiet
  timed nxsub_saber_v1.2_default 1 $SEC/nxsub_saber_v1.2_default.m8 -- $V12 --mode nx -q $SEC/nxsub_queries.fna -d $P -k 100000 -E 10 -t $THREADS --quiet
  timed nxsub_diamond_very-sensitive 1 $SEC/nxsub_diamond_very-sensitive.m8 -- $DMD blastx -q $SEC/nxsub_queries.fna -d $SEC/scop40.dmnd --threads $THREADS -k 0 -e 10 --very-sensitive --quiet -o /dev/stdout
  # pn (tblastn): DIAMOND has no tblastn mode.
  for p in default more-sensitive; do
    timed pn_saber_v1.2_$p 1 $SEC/pn_saber_v1.2_$p.m8 -- $V12 --mode pn -q $SEC/pn_queries.fa -d $NQ -k 100000 -E 10 -t $THREADS --preset $p --quiet
  done
  timed pnsub_saber_v1.2_default 1 $SEC/pnsub_saber_v1.2_default.m8 -- $V12 --mode pn -q $SEC/pnsub_queries.fa -d $NQ -k 100000 -E 10 -t $THREADS --quiet
  timed pnsub_saber_v1.1.0 1 $SEC/pnsub_saber_v1.1.0.m8 -- timeout 7200 $V110 --mode pn -q $SEC/pnsub_queries.fa -d $NQ -k 100000 -E 10 -t $THREADS --quiet
  for f in $SEC/*.m8; do
    n=$(basename $f .m8)
    local qopt=""
    case $n in
      nxsub_*) qopt="--queries $SEC/nxsub_queries.fna" ;;
      pnsub_*) qopt="--queries $SEC/pnsub_queries.fa" ;;
      pn_*) qopt="--queries $SEC/pn_queries.fa" ;;
    esac
    for sub in tune test all; do $EVAL $f $LAB --subset $sub --name $n $qopt > $SEC/$n.$sub.metrics.json; done
  done
}

sec_lowc() {
  begin lowc
  local Q=$DER/lowc.query.fa DB=$DER/lowc.db.fa
  $DMD makedb --in $DB -d $SEC/lowc --quiet
  timed lowc_saber_v1.2_default $REPS $SEC/lowc_saber_v1.2_default.m8 -- $V12 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet --stats $SEC/lowc_saber_v1.2_default.stats.json
  timed lowc_saber_v1.2_default_nomask 1 $SEC/lowc_saber_v1.2_default_nomask.m8 -- $V12 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --masking none --comp-bias off --quiet --stats $SEC/lowc_saber_v1.2_default_nomask.stats.json
  # adversarial: no masking, no frequency cap, every repeat hit enumerated
  # (first 10 queries) -- memory must stay bounded by the hit buffer.
  subset $Q $SEC/lowc10.fa 10
  timed lowc10_saber_v1.2_legacy_nomask_nocap 1 $SEC/lowc10_saber_v1.2_legacy_nomask_nocap.m8 -- timeout 600 $V12 -q $SEC/lowc10.fa -d $DB -k 100000 -E 1e-3 -t $THREADS --preset legacy --masking none --comp-bias off --seed-max-occ 1000000000 --quiet --stats $SEC/lowc10_saber_v1.2_legacy_nomask_nocap.stats.json
  timed lowc10_saber_v1.1.0 1 $SEC/lowc10_saber_v1.1.0.m8 -- $HERE/limited.sh 6000000 timeout 600 $V110 -q $SEC/lowc10.fa -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
  timed lowc_diamond_default $REPS $SEC/lowc_diamond_default.m8 -- $DMD blastp -q $Q -d $SEC/lowc.dmnd --threads $THREADS -k 0 -e 1e-3 --quiet -o /dev/stdout
  timed lowc_diamond_very-sensitive $REPS $SEC/lowc_diamond_very-sensitive.m8 -- $DMD blastp -q $Q -d $SEC/lowc.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/stdout
  timed lowc_saber_v1.1.0 1 $SEC/lowc_saber_v1.1.0.m8 -- $HERE/limited.sh 6000000 timeout 600 $V110 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
  timed lowc_saber_v1.0.2 1 $SEC/lowc_saber_v1.0.2.m8 -- $HERE/limited.sh 6000000 timeout 600 $V102 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
  # real queries against a database polluted with low-complexity sequences
  local MQ=$DATA/mmseqs_QUERY.fasta
  timed mmq_saber_v1.2_default 1 $SEC/mmq_saber_v1.2_default.m8 -- $V12 -q $MQ -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
  timed mmq_diamond_very-sensitive 1 $SEC/mmq_diamond_very-sensitive.m8 -- $DMD blastp -q $MQ -d $SEC/lowc.dmnd --threads $THREADS -k 0 -e 1e-3 --very-sensitive --quiet -o /dev/stdout
  timed mmq_saber_v1.1.0 1 $SEC/mmq_saber_v1.1.0.m8 -- $V110 -q $MQ -d $DB -k 100000 -E 1e-3 -t $THREADS --quiet
}

sec_multi() {
  begin multi
  local P=$DATA/scop40.fa
  $DMD makedb --in $P -d $SEC/scop40 --quiet
  for set in multi frag; do
    local Q=$DER/scop40.$set.fa
    timed ${set}_saber_v1.1.0 1 $SEC/${set}_saber_v1.1.0.m8 -- $V110 -q $Q -d $P -k 100000 -E 10 -t $THREADS --quiet
    for p in default more-sensitive; do
      timed ${set}_saber_v1.2_$p 1 $SEC/${set}_saber_v1.2_$p.m8 -- $V12 -q $Q -d $P -k 100000 -E 10 -t $THREADS --preset $p --quiet
    done
    for m in very-sensitive ultra-sensitive; do
      timed ${set}_diamond_$m 1 $SEC/${set}_diamond_$m.m8 -- $DMD blastp -q $Q -d $SEC/scop40.dmnd --threads $THREADS -k 0 -e 10 --$m --quiet -o /dev/stdout
    done
  done
  for t in saber_v1.1.0:saber_v1.1.0 saber_v1.2_default:saber_v1.2_default saber_v1.2_more-sensitive:saber_v1.2_more-sensitive diamond_very-sensitive:diamond_very-sensitive diamond_ultra-sensitive:diamond_ultra-sensitive; do
    local a=${t%%:*} b=${t#*:}
    $PY $HERE/eval_multi.py $SEC/multi_$a.m8 $OUT/scop/$b.m8 $LAB --name $a > $SEC/multi_$a.metrics.json
    for sub in tune test all; do $EVAL $SEC/frag_$a.m8 $LAB --subset $sub --name frag_$a > $SEC/frag_$a.$sub.metrics.json; done
  done
}

sec_threads() {
  begin threads
  local Q=$DATA/scop40.fa
  $DMD makedb --in $Q -d $SEC/scop40 --quiet
  subset $DER/mm.long.fa $SEC/long20.fa 25
  for t in 1 2; do
    timed t${t}_saber_v1.2_default $REPS - -- $V12 -q $Q -d $Q -k 100000 -E 10 -t $t --quiet
    timed t${t}_diamond_very-sensitive $REPS - -- $DMD blastp -q $Q -d $SEC/scop40.dmnd --threads $t -k 0 -e 10 --very-sensitive --quiet -o /dev/null
    timed t${t}_saber_v1.1.0 1 - -- $V110 -q $Q -d $Q -k 100000 -E 10 -t $t --quiet
    # one expensive query: intra-query parallelism (database split across threads)
    timed t${t}_saber_v1.2_intraquery $REPS - -- $V12 -q $SEC/long20.fa -d $DER/mm.db400.fa -k 100000 -E 1e-3 -t $t --query-batch 1 --quiet
  done
}

sec_extra() {
  # Ablations referenced in docs/LIMITATIONS.md (appended to the mmseqs
  # section's timings without moving its existing records).
  SEC=$OUT/mmseqs
  local Q=$DATA/mmseqs_QUERY.fasta DB=$DATA/mmseqs_DB.fasta
  timed saber_v1.2_default_nofloor 1 $SEC/saber_v1.2_default_nofloor.m8 -- $V12 -q $Q -d $DB -k 100000 -E 1e-3 -t $THREADS --preset default --hsp-floor-fraction 0 --quiet
  timed saber_v1.2_default_db 3 $SEC/saber_v1.2_default_db.m8 -- $V12 -q $Q --db $OUT/index/mmseqs_DB.sdx -k 100000 -E 1e-3 -t $THREADS --quiet
}

sec_rerun() {
  # Re-measure SCOP40 runs whose first measurement overlapped other work on
  # the machine (recorded in scop/rerun.jsonl; summarize.py prefers it).
  SEC=$OUT/scop
  local Q=$DATA/scop40.fa
  for p in default sensitive more-sensitive; do
    $RB --label saber_v1.2_$p --reps $REPS --out $SEC/rerun.jsonl --stdout $SEC/saber_v1.2_$p.rerun.m8 -- \
      $V12 -q $Q -d $Q -k 100000 -E 10 -t $THREADS --preset $p --quiet 2>>$SEC/runbench.err
    cmp -s $SEC/saber_v1.2_$p.m8 $SEC/saber_v1.2_$p.rerun.m8 && echo "saber_v1.2_$p: rerun output identical" >> $SEC/rerun_check.txt \
      || echo "saber_v1.2_$p: rerun output DIFFERS" >> $SEC/rerun_check.txt
    rm -f $SEC/saber_v1.2_$p.rerun.m8
  done
}

sec_sweep() {
  # Seed / filter configurations considered for the presets, evaluated on the
  # SCOP40 tune split only (the test split is never used for choosing).
  begin sweep
  local Q=$DATA/scop40.fa
  while read -r name args; do
    [ -z "$name" ] && continue
    timed $name 1 $SEC/$name.m8 -- $V12 -q $Q -d $Q -k 100000 -E 10 -t $THREADS --quiet $args
    $EVAL $SEC/$name.m8 $LAB --subset tune --name $name > $SEC/$name.tune.metrics.json
    rm -f $SEC/$name.m8
  done <<'CONFIGS'
full_1111_hf0 --seed-alphabet full --seed-shape 1111 --hit-filter-ids 0
full_1111_hf8 --seed-alphabet full --seed-shape 1111 --hit-filter-ids 8
full_1111_hf10 --seed-alphabet full --seed-shape 1111 --hit-filter-ids 10
full_11111_hf0 --seed-alphabet full --seed-shape 11111 --hit-filter-ids 0
full_110111_hf0 --seed-alphabet full --seed-shape 110111 --hit-filter-ids 0
m10_11111_hf12 --seed-alphabet murphy10 --seed-shape 11111 --hit-filter-ids 12
m10_110111_hf11 --seed-alphabet murphy10 --seed-shape 110111 --hit-filter-ids 11
m10_110111_hf12 --seed-alphabet murphy10 --seed-shape 110111 --hit-filter-ids 12
m10_110111_hf13 --seed-alphabet murphy10 --seed-shape 110111 --hit-filter-ids 13
m10_111111_hf11 --seed-alphabet murphy10 --seed-shape 111111 --hit-filter-ids 11
m10_1111011_hf11 --seed-alphabet murphy10 --seed-shape 1111011 --hit-filter-ids 11
m10_1111011_hf12 --seed-alphabet murphy10 --seed-shape 1111011 --hit-filter-ids 12
m10_1111011_hf13 --seed-alphabet murphy10 --seed-shape 1111011 --hit-filter-ids 13
m10_11101111_hf12 --seed-alphabet murphy10 --seed-shape 11101111 --hit-filter-ids 12
m10_1111011+1101111_hf12 --seed-alphabet murphy10 --seed-shape 1111011,1101111 --hit-filter-ids 12
m10_1101011+111011_hf11 --seed-alphabet murphy10 --seed-shape 1101011,111011 --hit-filter-ids 11 --index-max-candidates 8000
m10_1101011+111011_hf12 --seed-alphabet murphy10 --seed-shape 1101011,111011 --hit-filter-ids 12 --index-max-candidates 8000
full_1111_twohit40 --seed-alphabet full --seed-shape 1111 --index-min-hits 2 --two-hit-window 40
full_1111_nbT13_twohit --neighborhood-t 13
CONFIGS
}

sections=("$@")
[ ${#sections[@]} -eq 0 ] && sections=(env scop decoy mmseqs index extra translated lowc multi threads sweep)
for s in "${sections[@]}"; do
  log "=== section $s ==="
  sec_$s
done
log "done"
