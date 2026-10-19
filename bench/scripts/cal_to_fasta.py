#!/usr/bin/env python3
"""Convert reseek's SCOP40 .cal (sequence + C-alpha) file to FASTA + label TSV.

Input : scop40.cal.gz  (">domid/sccs" header, then one "aa\tx\ty\tz" line per residue)
Output: scop40.fa      (">domid sccs")
        scop40.labels.tsv (domid, class, fold, superfamily, family, length)
"""
import gzip, sys, collections
src, out_fa, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
recs = []
cur_id = None; cur_lab = None; seq = []
with gzip.open(src, "rt") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if cur_id is not None:
                recs.append((cur_id, cur_lab, "".join(seq)))
            hdr = line[1:]
            cur_id, cur_lab = hdr.split("/", 1)
            seq = []
        else:
            seq.append(line.split("\t", 1)[0])
    if cur_id is not None:
        recs.append((cur_id, cur_lab, "".join(seq)))
with open(out_fa, "w") as fa, open(out_tsv, "w") as tsv:
    tsv.write("domain\tclass\tfold\tsuperfamily\tfamily\tlength\n")
    for d, lab, s in recs:
        fa.write(f">{d} {lab}\n")
        for i in range(0, len(s), 80):
            fa.write(s[i:i+80] + "\n")
        p = lab.split(".")
        tsv.write(f"{d}\t{p[0]}\t{'.'.join(p[:2])}\t{'.'.join(p[:3])}\t{lab}\t{len(s)}\n")
c = collections.Counter
print("domains", len(recs))
print("classes", len(c(r[1].split('.')[0] for r in recs)))
print("folds", len(c('.'.join(r[1].split('.')[:2]) for r in recs)))
print("superfamilies", len(c('.'.join(r[1].split('.')[:3]) for r in recs)))
print("families", len(c(r[1] for r in recs)))
L = sorted(len(r[2]) for r in recs)
print("len min/median/mean/max", L[0], L[len(L)//2], sum(L)/len(L), L[-1], "total", sum(L))
aa = c("".join(r[2] for r in recs))
print("nonstd residues", {k: v for k, v in aa.items() if k not in "ACDEFGHIKLMNPQRSTVWY"})
