#!/bin/bash
# bench_5way.sh: BLAST vs SWORD vs DIAMOND vs BLAT vs SABER on real biology.
#
# Outputs /tmp/bench_results.tsv with columns:
#   tool, mode, time_s, peak_mem_kb, hits, disk_db_mb

set -u
cd "$(dirname "$0")"
SWORD=/home/claude/work/sword/build/bin/sword
DIAMOND=/tmp/diamond/build/diamond
BLAT=/tmp/blat/bin/blat
SABER=/home/claude/work/saber/target/release/saber

Q=/tmp/bio_queries100.fa
DB_FA=/tmp/bio_db.fa
BLAST_DB=/tmp/bio_blast_db
DIAMOND_DB=/tmp/bio_diamond_db.dmnd

OUT=${OUT:-/tmp/bench_results.tsv}
echo "Writing results to $OUT (override with OUT=/path/to/file env var)"
echo -e "tool\tmode\ttime_s\tpeak_mem_mb\thits\tdisk_mb" > "$OUT"

bench() {
  local label="$1"; shift
  local mode="$1"; shift
  local outf="$1"; shift
  local diskmb="$1"; shift
  local cmd="$*"

  # Run with /usr/bin/time -v, parse the relevant fields.
  local tmpf=$(mktemp)
  /usr/bin/time -v bash -c "$cmd" > /dev/null 2> "$tmpf"
  local secs=$(grep "Elapsed" "$tmpf" | awk -F': ' '{print $NF}' | awk -F: '{
    if (NF==2) { print $1*60 + $2 } else if (NF==3) { print $1*3600 + $2*60 + $3 }
  }')
  local memkb=$(grep "Maximum resident" "$tmpf" | awk '{print $NF}')
  local memmb=$(echo "scale=1; $memkb / 1024" | bc)
  local hits=$(awk '{print $1"\t"$2}' "$outf" 2>/dev/null | sort -u | wc -l)
  echo -e "$label\t$mode\t$secs\t$memmb\t$hits\t$diskmb" | tee -a "$OUT"
  rm -f "$tmpf"
}

echo "Running 5-way head-to-head bench..."
echo ""

# DB sizes
BLAST_DB_MB=$(du -shc /tmp/bio_blast_db.* 2>/dev/null | tail -1 | awk '{print $1}' | sed 's/M//')
DIAMOND_DB_MB=$(du -sh "$DIAMOND_DB" 2>/dev/null | awk '{print $1}' | sed 's/M//')
SOURCE_FA_MB=$(du -sh "$DB_FA" 2>/dev/null | awk '{print $1}' | sed 's/M//')

# Each tool with max_aligns ~ 10
bench BLAST    blastp     /tmp/o_blast.tsv   "$BLAST_DB_MB"   \
  "/usr/bin/blastp -query $Q -db $BLAST_DB -outfmt 6 -evalue 1e-5 -num_threads 1 -max_target_seqs 10 -out /tmp/o_blast.tsv"

bench SWORD    sw         /tmp/o_sword.tsv   "$SOURCE_FA_MB"  \
  "$SWORD -i $Q -j $DB_FA --threads 1 -v 1e-5 -f bm8 -o /tmp/o_sword.tsv"

bench DIAMOND  default    /tmp/o_diamond_d.tsv "$DIAMOND_DB_MB" \
  "$DIAMOND blastp -q $Q -d $DIAMOND_DB -o /tmp/o_diamond_d.tsv -e 1e-5 -k 10 -p 1 --quiet"

bench DIAMOND  sensitive  /tmp/o_diamond_s.tsv "$DIAMOND_DB_MB" \
  "$DIAMOND blastp -q $Q -d $DIAMOND_DB -o /tmp/o_diamond_s.tsv --sensitive -e 1e-5 -k 10 -p 1 --quiet"

bench DIAMOND  more-sens  /tmp/o_diamond_m.tsv "$DIAMOND_DB_MB" \
  "$DIAMOND blastp -q $Q -d $DIAMOND_DB -o /tmp/o_diamond_m.tsv --more-sensitive -e 1e-5 -k 10 -p 1 --quiet"

# BLAT: no E-value option; we post-filter to E≤1e-5 for fair hit count
bench BLAT     prot-prot  /tmp/o_blat_f.tsv  "$SOURCE_FA_MB"  \
  "$BLAT -t=prot -q=prot -out=blast8 $DB_FA $Q /tmp/o_blat.tsv && awk '\$11 <= 1e-5' /tmp/o_blat.tsv > /tmp/o_blat_f.tsv"

# Build a persistent SABER index once, then measure the memory-mapped search
# path (route B). This is the recommended production workflow: build the index
# one time, reuse it across many queries with minimal peak memory.
SABER_IDX=/tmp/saber_bench.sdx
"$SABER" --makedb "$SABER_IDX" -d "$DB_FA" --mode pp --quiet 2>/dev/null
SABER_IDX_MB=$(du -m "$SABER_IDX" 2>/dev/null | cut -f1)

bench SABER    default    /tmp/o_saber_d.tsv "$SABER_IDX_MB"  \
  "$SABER -q $Q --db $SABER_IDX --mode pp --outfmt m8 --threads 1 -E 1e-5 --quiet -o /tmp/o_saber_d.tsv"

bench SABER    fast       /tmp/o_saber_f.tsv "$SABER_IDX_MB"  \
  "$SABER -q $Q --db $SABER_IDX --mode pp --outfmt m8 --threads 1 -E 1e-5 --hsp-threshold 50 --quiet -o /tmp/o_saber_f.tsv"

bench SABER    sensitive  /tmp/o_saber_s.tsv "$SABER_IDX_MB"  \
  "$SABER -q $Q -d $DB_FA --mode pp --outfmt m8 --threads 1 -E 1e-5 --sensitive --quiet -o /tmp/o_saber_s.tsv"

echo ""
echo "Results written to $OUT"
cat "$OUT"
