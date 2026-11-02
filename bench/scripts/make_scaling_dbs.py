#!/usr/bin/env python3
"""Databases of increasing size for the scaling benchmark.

    make_scaling_dbs.py SCOP40_FASTA MM_DB_FASTA MM_DB400_FASTA OUTDIR [--k 1,2,4,8,16]

OUTDIR/scale_k{K}.fa = every SCOP40 domain followed by K x 20,000 reversed
proteins REV_0 .. REV_{20000K-1}. The decoy sources are, in order, the 80,000
sequences of MM_DB400_FASTA (`mm.db400.fa` from make_datasets.py: the MMseqs2
example database and three copies at 15 % substitution) and then further
copies of MM_DB_FASTA at 15 % substitution (copy 4, 5, ...; fixed seed), so
scale_k4.fa is the at-scale benchmark database of make_scale_db.py and every
database contains the decoys of the smaller ones. Each original sequence
gives K similar decoys (one per copy).
"""
import argparse
import os
import random

AA = "ACDEFGHIKLMNPQRSTVWY"


def read_fasta(path):
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scop40")
    ap.add_argument("mmdb")
    ap.add_argument("db400")
    ap.add_argument("outdir")
    ap.add_argument("--k", default="1,2,4,8,16")
    a = ap.parse_args()
    ks = sorted({int(x) for x in a.k.split(",")})
    orig = [s for _, s in read_fasta(a.mmdb)]
    n0 = len(orig)
    need = max(ks) * n0
    sources = [s for _, s in read_fasta(a.db400)][:need]
    if len(sources) != min(need, 4 * n0):
        raise SystemExit(f"{a.db400}: expected {4 * n0} sequences")
    r = random.Random(20261019)
    c = 4
    while len(sources) < need:
        for s in orig:
            sources.append("".join(ch if r.random() >= 0.15 else r.choice(AA) for ch in s))
        c += 1
    sources = sources[:need]
    scop = [(name.split()[0], seq) for name, seq in read_fasta(a.scop40)]
    os.makedirs(a.outdir, exist_ok=True)
    for k in ks:
        path = os.path.join(a.outdir, f"scale_k{k}.fa")
        res = 0
        with open(path, "w") as out:
            for name, seq in scop:
                out.write(f">{name}\n{seq}\n")
                res += len(seq)
            for i in range(k * n0):
                out.write(f">REV_{i}\n{sources[i][::-1]}\n")
                res += len(sources[i])
        print(f"{path}\t{len(scop) + k * n0} sequences\t{res} residues\tdecoys {k * n0}")


if __name__ == "__main__":
    main()
