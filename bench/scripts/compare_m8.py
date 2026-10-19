#!/usr/bin/env python3
"""Compare two BLAST-tabular (m8) outputs of the same search line by line and
say how they differ: same query-subject pairs? same E-values and bit scores?
which columns change?

  python3 -I compare_m8.py A.m8 B.m8 [LABEL]
"""
import collections, sys

NAMES = "query subject pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()


def rows(path):
    return [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]


def main():
    a, b = rows(sys.argv[1]), rows(sys.argv[2])
    label = sys.argv[3] if len(sys.argv) > 3 else f"{sys.argv[1]} vs {sys.argv[2]}"
    if a == b:
        print(f"{label}: identical ({len(a)} lines)")
        return
    pa = {(r[0], r[1]): r for r in a}
    pb = {(r[0], r[1]): r for r in b}
    common = set(pa) & set(pb)
    diff = [k for k in common if pa[k] != pb[k]]
    same_score = sum(1 for k in diff if pa[k][10:12] == pb[k][10:12])
    higher_b = sum(1 for k in diff if float(pb[k][11]) > float(pa[k][11]))
    lower_b = sum(1 for k in diff if float(pb[k][11]) < float(pa[k][11]))
    cols = collections.Counter(NAMES[i] for k in diff for i, (x, y) in enumerate(zip(pa[k], pb[k])) if x != y)
    order_only = len(diff) == 0 and set(pa) == set(pb)
    print(f"{label}: {len(a)} vs {len(b)} lines; pairs only in first {len(set(pa) - set(pb))}, only in second "
          f"{len(set(pb) - set(pa))}; {len(diff)} common pairs differ ({same_score} with the same E-value and bit score, "
          f"{higher_b} higher / {lower_b} lower bit score in the second)"
          + (f"; columns that change: {', '.join(f'{c} {n}' for c, n in cols.most_common())}" if cols else "")
          + ("; same lines in another order" if order_only else ""))


if __name__ == "__main__":
    main()
