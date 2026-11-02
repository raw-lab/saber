#!/usr/bin/env python3
"""At-scale database with decoys of domain length.

    make_lenmatched_db.py SCOP40_FASTA MM_DB400_FASTA OUT_FASTA [--sources N]

OUT = every SCOP40 domain followed by the first N sequences of
MM_DB400_FASTA written backwards (as in make_scale_db.py) and cut into
consecutive pieces whose lengths are drawn from the SCOP40 domain lengths
(fixed seed; a last piece shorter than 30 residues is dropped), named
REV_0, REV_1, ... True subjects (SCOP40 domains) and false subjects
(reversed pieces) then have the same length distribution, so an E-value
model that depends on subject length gains nothing from a length difference
between them.
"""
import argparse
import random


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
    ap.add_argument("--sources", type=int, default=80000)
    a = ap.parse_args()
    scop = [(n.split()[0], s) for n, s in read_fasta(a.scop40)]
    lens = [len(s) for _, s in scop]
    r = random.Random(7)
    n = res = 0
    with open(a.out, "w") as out:
        for name, s in scop:
            out.write(f">{name}\n{s}\n")
        for i, (_, s) in enumerate(read_fasta(a.db400)):
            if i >= a.sources:
                break
            rs = s[::-1]
            j = 0
            while j < len(rs):
                piece = rs[j:j + r.choice(lens)]
                j += len(piece)
                if len(piece) < 30:
                    break
                out.write(f">REV_{n}\n{piece}\n")
                n += 1
                res += len(piece)
    print(f"{a.out}\t{len(scop)} domains\t{n} decoys\t{res} decoy residues")


if __name__ == "__main__":
    main()
