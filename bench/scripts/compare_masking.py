#!/usr/bin/env python3
"""Compare two tantan-style probability files (tantan -p -f 1 layout).
Reports correlation, mean |diff| and agreement of the >= 0.5 mask calls."""
import sys, math
def load(p):
    out = []
    with open(p) as f:
        for line in f:
            if line.startswith(">"):
                continue
            out.append(float(line))
    return out
a = load(sys.argv[1]); b = load(sys.argv[2])
assert len(a) == len(b), (len(a), len(b))
n = len(a)
ma = [x >= 0.5 for x in a]; mb = [x >= 0.5 for x in b]
tp = sum(1 for x, y in zip(ma, mb) if x and y)
only_a = sum(1 for x, y in zip(ma, mb) if x and not y)
only_b = sum(1 for x, y in zip(ma, mb) if y and not x)
agree = sum(1 for x, y in zip(ma, mb) if x == y) / n
mad = sum(abs(x - y) for x, y in zip(a, b)) / n
mx = sum(a)/n; my = sum(b)/n
cov = sum((x-mx)*(y-my) for x, y in zip(a, b)); va = sum((x-mx)**2 for x in a); vb = sum((y-my)**2 for y in b)
r = cov / math.sqrt(va*vb) if va > 0 and vb > 0 else float('nan')
print(f"positions={n} masked_ref={sum(ma)} ({sum(ma)/n:.3%}) masked_test={sum(mb)} ({sum(mb)/n:.3%})")
print(f"call agreement={agree:.5f} both={tp} ref_only={only_a} test_only={only_b} "
      f"jaccard={tp/max(1,tp+only_a+only_b):.4f}")
print(f"pearson r={r:.5f} mean|diff|={mad:.5f}")
