#!/usr/bin/env python3
"""1.4 vs 1.3 outputs per database size, and true positives lost as the
database grows (pairs reported at K=4 whose E-value, scaled to the larger
database, would still be <= 10, but which are missing at the larger size).
  compare_scaling_outputs.py RESULTS LABELS   (after run_all.sh scaling)"""
import filecmp, json, os, sys

R, LAB = sys.argv[1], sys.argv[2]
S = os.path.join(R, "scaling")
lab = {}
with open(LAB) as f:
    next(f)
    for line in f:
        d, cls, fold, sf, fam, ln = line.rstrip("\n").split("\t")[:6]
        lab[d] = (fold, sf)
sizes = {}
for line in open(os.path.join(S, "dbsizes.tsv")):
    k, n, res, d = line.split()
    sizes[int(k)] = (int(n), int(res))
ks = sorted(sizes)


def pairs(path, bits=None):
    """(query, subject) -> best E-value; with `bits` (a dict), also -> bit
    score of that best hit."""
    out = {}
    for line in open(path):
        p = line.split("\t")
        q, s, e = p[0], p[1], float(p[10])
        if (q, s) not in out or e < out[(q, s)]:
            out[(q, s)] = e
            if bits is not None:
                bits[(q, s)] = float(p[11])
    return out


def tp(k):
    q, s = k
    return q != s and q in lab and s in lab and lab[q][1] == lab[s][1]


def fp(k):
    q, s = k
    return s.startswith("REV_") or (q in lab and s in lab and lab[q][0] != lab[s][0])


print("1.4 default vs 1.3.0 default (same index):")
for k in ks:
    a, b = os.path.join(S, f"k{k}.saber_v1.4_default.m8"), os.path.join(S, f"k{k}.saber_v1.3_default.m8")
    if not (os.path.exists(a) and os.path.exists(b)):
        continue
    same = filecmp.cmp(a, b, shallow=False)
    ba, bb = {}, {}
    pa, pb = pairs(a, ba), pairs(b, bb)
    gained = [x for x in pa if x not in pb]
    lost = [x for x in pb if x not in pa]
    eg = [x for x in pa if x in pb and pa[x] != pb[x]]
    print(f"  K={k:<2} ({sizes[k][1]/1e6:.1f} M): {'IDENTICAL files' if same else 'different'}; pairs 1.4 {len(pa)} / 1.3 {len(pb)}; "
          f"gained {len(gained)} (TP {sum(map(tp, gained))}, FP {sum(map(fp, gained))}), lost {len(lost)} (TP {sum(map(tp, lost))}, FP {sum(map(fp, lost))}); "
          f"common pairs with a different E-value {len(eg)}")
    for x in sorted(eg):
        print(f"      {x[0]} -> {x[1]} ({'TP' if tp(x) else ('FP' if fp(x) else 'ignored')}): best hit 1.3 {bb[x]:.1f} bits, E {pb[x]:.3g}; "
              f"1.4 {ba[x]:.1f} bits, E {pa[x]:.3g}")

print("\nTrue positives lost from K=4 to K=16 although their scaled E-value stays <= 10:")
for name, ratio_kind in (("saber_v1.3_default", "seq"), ("saber_v1.4_default", "seq"), ("diamond_very-sensitive", "res"),
                         ("saber_v1.4_more-sensitive", "seq"), ("diamond_ultra-sensitive", "res")):
    for k0, k1 in ((4, 16), (4, 8), (1, 16)):
        a, b = os.path.join(S, f"k{k0}.{name}.m8"), os.path.join(S, f"k{k1}.{name}.m8")
        if not (os.path.exists(a) and os.path.exists(b)) or k0 not in sizes or k1 not in sizes:
            continue
        ratio = sizes[k1][0] / sizes[k0][0] if ratio_kind == "seq" else sizes[k1][1] / sizes[k0][1]
        pa, pb = pairs(a), pairs(b)
        tpa = {x: e for x, e in pa.items() if tp(x)}
        expect = [x for x, e in tpa.items() if e * ratio <= 10]
        lost = [x for x in expect if x not in pb]
        strong = [x for x in lost if tpa[x] * ratio <= 1e-3]
        print(f"  {name:28s} K={k0}->{k1} (E x {ratio:.2f}): TPs at K={k0} {len(tpa)}, expected to stay {len(expect)}, lost {len(lost)} "
              f"(of which E <= 1e-3 after scaling: {len(strong)}); TPs at K={k1} {sum(map(tp, pb))}")
