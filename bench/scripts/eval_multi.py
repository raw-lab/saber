#!/usr/bin/env python3
"""Two-domain query evaluation ("A+B" queries from make_datasets.py).

For every query A+B and each constituent domain X in {A, B}:
  TP_X = reported subjects (E <= cutoff) in X's superfamily, excluding A and B
         themselves; FP = subjects whose fold differs from both A's and B's.
Compares against the same tool's single-domain searches (all-vs-all SCOP40
output) to measure how much of each domain's evidence survives when it is
embedded in a longer, two-domain query:

  retention_weaker = sum_q TP_weak(A+B) / sum_q TP_weak(alone)

where "weak" is the constituent with fewer TPs when searched alone (the one
most at risk of being crowded out). Also reports FP counts at the cutoff.
"""
import argparse, collections, json

def load_labels(path):
    lab = {}
    with open(path) as f:
        next(f)
        for line in f:
            d, cls, fold, sf, fam, ln = line.rstrip("\n").split("\t")
            lab[d] = (fold, sf, fam)
    return lab

def hits(path, cutoff):
    h = collections.defaultdict(set)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            try:
                ev = float(p[10])
            except ValueError:
                continue
            if ev <= cutoff:
                h[p[0].split()[0]].add(p[1].split()[0])
    return h

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("multi_m8")
    ap.add_argument("single_m8", help="the same tool's all-vs-all SCOP40 output")
    ap.add_argument("labels")
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--name", default="")
    a = ap.parse_args()
    lab = load_labels(a.labels)
    multi = hits(a.multi_m8, a.evalue)
    single = hits(a.single_m8, a.evalue)
    queries = set()
    with open(a.multi_m8) as f:
        pass
    # queries come from the dataset naming "A+B"; collect from both files' keys
    for q in list(multi):
        queries.add(q)
    tot_weak_m = tot_weak_s = tot_strong_m = tot_strong_s = 0
    fp_multi = 0
    nq = 0
    for q in sorted(queries):
        if "+" not in q:
            continue
        A, B = q.split("+", 1)
        if A not in lab or B not in lab:
            continue
        nq += 1
        def tps(subjects, X):
            return {s for s in subjects if s in lab and s not in (A, B) and lab[s][1] == lab[X][1]}
        sa, sb = tps(single.get(A, ()), A), tps(single.get(B, ()), B)
        ma, mb = tps(multi.get(q, ()), A), tps(multi.get(q, ()), B)
        weak, strong = ((sa, ma), (sb, mb)) if len(sa) <= len(sb) else ((sb, mb), (sa, ma))
        tot_weak_s += len(weak[0]); tot_weak_m += len(weak[1] & weak[0])
        tot_strong_s += len(strong[0]); tot_strong_m += len(strong[1] & strong[0])
        for s in multi.get(q, ()):
            if s in lab and lab[s][0] != lab[A][0] and lab[s][0] != lab[B][0]:
                fp_multi += 1
    out = dict(
        name=a.name, evalue=a.evalue, queries_with_hits=nq,
        retention_weaker=(tot_weak_m / tot_weak_s) if tot_weak_s else None,
        retention_stronger=(tot_strong_m / tot_strong_s) if tot_strong_s else None,
        weaker_tp_alone=tot_weak_s, weaker_tp_in_multidomain=tot_weak_m,
        stronger_tp_alone=tot_strong_s, stronger_tp_in_multidomain=tot_strong_m,
        fp_multidomain=fp_multi,
    )
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
