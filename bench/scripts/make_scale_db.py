#!/usr/bin/env python3
"""Build the at-scale accuracy benchmark database.

    make_scale_db.py SCOP40_FASTA MM_DB400_FASTA OUT_FASTA [--decoys N]

OUT = every SCOP40 domain (header reduced to its identifier, sequence on one
line) followed by N reversed protein sequences named REV_0 .. REV_{N-1}:
REV_i is the i-th sequence of MM_DB400_FASTA (bench/scripts/make_datasets.py
writes it as `mm.db400.fa`: the MMseqs2 example database and three copies at
15 % substitution) written backwards. Reversed proteins keep the length and
composition of real proteins but have no homologs among the SCOP40 domains,
so every hit to a REV_ sequence is a false positive (eval_scop.py
--decoy-prefix REV_ --decoys N). Note that the four copies of each original
sequence give four similar decoys.
"""
import argparse


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
    ap.add_argument("db400")
    ap.add_argument("out")
    ap.add_argument("--decoys", type=int, default=80000)
    a = ap.parse_args()
    n = 0
    with open(a.out, "w") as out:
        for name, seq in read_fasta(a.scop40):
            out.write(f">{name.split()[0]}\n{seq}\n")
        for i, (_, seq) in enumerate(read_fasta(a.db400)):
            if i >= a.decoys:
                break
            out.write(f">REV_{i}\n{seq[::-1]}\n")
            n += 1
    if n < a.decoys:
        raise SystemExit(f"only {n} sequences in {a.db400}, {a.decoys} decoys requested")


if __name__ == "__main__":
    main()
