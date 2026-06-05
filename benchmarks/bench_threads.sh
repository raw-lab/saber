#!/bin/bash
# bench_threads.sh: scaling test across thread counts for BLAST, SWORD,
# DIAMOND, SABER. BLAT does not support multi-threading.
#
# Output: /tmp/bench_threads.tsv with columns:
#   tool, threads, time_s, peak_mem_mb, hits
#
# Run with several thread counts; on a typical 16-core machine, try
# THREADS="1 2 4 8 16". On a workstation with 32+ cores, try 1 2 4 8 16 32.

set -u
cd "$(dirname "$0")"
SWORD=${SWORD:-/home/claude/work/sword/build/bin/sword}
DIAMOND=${DIAMOND:-/tmp/diamond/build/diamond}
SABER=${SABER:-/home/claude/work/saber/target/release/saber}

Q=${Q:-/tmp/bio_queries100.fa}
DB_FA=${DB_FA:-/tmp/bio_db.fa}
BLAST_DB=${BLAST_DB:-/tmp/bio_blast_db}
DIAMOND_DB=${DIAMOND_DB:-/tmp/bio_diamond_db.dmnd}

THREADS=${THREADS:-"1 2 4 8 16"}
OUT=${OUT:-/tmp/bench_threads.tsv}
echo -e "tool\tthreads\ttime_s\tpeak_mem_mb\thits" > "$OUT"

# Build the SABER persistent index once; the scaling loop reuses it via --db.
SABER_IDX=${SABER_IDX:-/tmp/saber_threads.sdx}
"$SABER" --makedb "$SABER_IDX" -d "$DB_FA" --mode pp --quiet 2>/dev/null

bench() {
  local label="$1"; local nthr="$2"; local outf="$3"; shift 3
  local cmd="$*"
  local tmpf=$(mktemp)
  /usr/bin/time -v bash -c "$cmd" > /dev/null 2> "$tmpf"
  local secs=$(grep "Elapsed" "$tmpf" | awk -F': ' '{print $NF}' | awk -F: '{
    if (NF==2) { print $1*60 + $2 } else if (NF==3) { print $1*3600 + $2*60 + $3 }
  }')
  local memkb=$(grep "Maximum resident" "$tmpf" | awk '{print $NF}')
  local memmb=$(echo "scale=1; $memkb / 1024" | bc)
  local hits=$(awk '{print $1"\t"$2}' "$outf" 2>/dev/null | sort -u | wc -l)
  echo -e "$label\t$nthr\t$secs\t$memmb\t$hits" | tee -a "$OUT"
  rm -f "$tmpf"
}

for T in $THREADS; do
  echo ""
  echo "=== threads=$T ==="

  bench BLAST   "$T" /tmp/o_blast_t.tsv \
    "/usr/bin/blastp -query $Q -db $BLAST_DB -outfmt 6 -evalue 1e-5 -num_threads $T -max_target_seqs 10 -out /tmp/o_blast_t.tsv"

  bench SWORD   "$T" /tmp/o_sword_t.tsv \
    "$SWORD -i $Q -j $DB_FA --threads $T -v 1e-5 -f bm8 -o /tmp/o_sword_t.tsv"

  bench DIAMOND-sensitive "$T" /tmp/o_diamond_t.tsv \
    "$DIAMOND blastp -q $Q -d $DIAMOND_DB -o /tmp/o_diamond_t.tsv --sensitive -e 1e-5 -k 10 -p $T --quiet"

  bench SABER-default "$T" /tmp/o_saber_d_t.tsv \
    "$SABER -q $Q --db $SABER_IDX --mode pp --outfmt m8 --threads $T -E 1e-5 --quiet -o /tmp/o_saber_d_t.tsv"

  bench SABER-fast "$T" /tmp/o_saber_f_t.tsv \
    "$SABER -q $Q --db $SABER_IDX --mode pp --outfmt m8 --threads $T -E 1e-5 --hsp-threshold 50 --quiet -o /tmp/o_saber_f_t.tsv"
done

echo ""
echo "Wrote $OUT — run benchmarks/bio/plot_threads.py to visualize scaling."
