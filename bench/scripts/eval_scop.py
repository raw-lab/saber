#!/usr/bin/env python3
"""Score an all-vs-all SCOP40 search (BLAST tabular / m8) against SCOP labels.

Labels come from SCOPe 2.08 40%-identity domains (reseek test_data/scop40.cal.gz).
Pair classes (Foldseek / MMseqs2 convention):
  self            query == subject                      -> ignored
  TP_family       same family                           -> true positive
  TP_superfamily  same superfamily, different family    -> true positive
  ignore_fold     same fold, different superfamily      -> ignored (ambiguous)
  FP              different fold                        -> false positive

Metrics (all over the chosen query subset; --queries restricts it to the
queries of a subset run):
  * sens_to_1st_fp_{family,superfamily}: Foldseek/MMseqs2 "AUC1"-style per-query
    sensitivity counted down the E-value-sorted hit list until the first FP,
    averaged over queries that have >=1 possible TP at that level.
  * pooled ROC: hits from all queries pooled and sorted by E-value; at each
    threshold we know TP, FP (counts of pairs) ->
       recall    = TP / (all possible same-superfamily pairs)
       precision = TP / (TP + FP)                (ignored pairs excluded)
       FDP       = FP / (TP + FP)
       FPR       = FP / (all different-fold pairs)
       EPQ       = FP / n_queries               (errors per query)
  * recall at fixed EPQ (0.01, 0.1, 1) and at fixed FDP (1%, 5%, 10%)
  * calibration: EPQ at E-value cut-offs vs the cut-off itself.

With --decoy-prefix P (and --decoys N), database sequences whose ID starts
with P are known non-homologues (e.g. reversed background sequences added to
test accuracy on a large database): every hit to them is a false positive and
they add N negative pairs per query.

Output: JSON summary on stdout (and optional --roc TSV with the curve).
"""
import argparse, collections, json, math, sys, hashlib

def load_labels(path):
    lab = {}
    with open(path) as f:
        next(f)
        for line in f:
            d, cls, fold, sf, fam, ln = line.rstrip("\n").split("\t")
            lab[d] = (fold, sf, fam)
    return lab

def split_of(sf):
    """Deterministic tune/test split by superfamily (MD5 parity)."""
    h = int(hashlib.md5(sf.encode()).hexdigest(), 16)
    return "tune" if h % 2 == 0 else "test"

DECOY = None  # ID prefix of known non-homologous database sequences


def parse_hits(path, lab, qset):
    hits = collections.defaultdict(list)
    n_lines = 0
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            q = p[0].split()[0]; s = p[1].split()[0]
            # strip translated-frame suffixes like "|frame+1" if present
            q = q.split("|")[0]; s = s.split("|")[0]
            if q not in qset or (s not in lab and not (DECOY and s.startswith(DECOY))):
                continue
            try:
                ev = float(p[10]); bits = float(p[11])
            except ValueError:
                continue
            hits[q].append((ev, -bits, s))
            n_lines += 1
    # one entry per (q,s): keep the best HSP
    for q in hits:
        best = {}
        for ev, nb, s in hits[q]:
            if s not in best or (ev, nb) < best[s]:
                best[s] = (ev, nb)
        hits[q] = sorted((v[0], v[1], s) for s, v in best.items())
    return hits, n_lines

def classify(q, s, lab):
    if q == s:
        return "self"
    if DECOY and s.startswith(DECOY):
        return "FP"
    fq, sq, mq = lab[q]; fs, ss, ms = lab[s]
    if mq == ms:
        return "TP_family"
    if sq == ss:
        return "TP_superfamily"
    if fq == fs:
        return "ignore_fold"
    return "FP"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("m8")
    ap.add_argument("labels")
    ap.add_argument("--subset", choices=["all", "tune", "test"], default="all")
    ap.add_argument("--roc", help="write pooled ROC curve TSV here")
    ap.add_argument("--name", default="")
    ap.add_argument("--queries", help="FASTA of the queries that were searched (subset runs): "
                    "metrics are computed over these queries only")
    ap.add_argument("--decoy-prefix", help="ID prefix of known non-homologous database sequences (every hit = FP)")
    ap.add_argument("--decoys", type=int, default=0, help="number of such sequences in the database")
    a = ap.parse_args()
    global DECOY
    DECOY = a.decoy_prefix

    lab = load_labels(a.labels)
    doms = list(lab)
    if a.subset == "all":
        qset = set(doms)
    else:
        qset = {d for d in doms if split_of(lab[d][1]) == a.subset}
    if a.queries:
        searched = set()
        with open(a.queries) as f:
            for line in f:
                if line.startswith(">"):
                    searched.add(line[1:].split()[0].split("|")[0])
        qset &= searched

    # possible TP counts and negatives
    by_fam = collections.Counter(v[2] for v in lab.values())
    by_sf = collections.Counter(v[1] for v in lab.values())
    by_fold = collections.Counter(v[0] for v in lab.values())
    n_db = len(lab)
    poss_fam = {q: by_fam[lab[q][2]] - 1 for q in qset}
    poss_sf = {q: by_sf[lab[q][1]] - 1 for q in qset}           # family + superfamily TPs
    poss_sfonly = {q: poss_sf[q] - poss_fam[q] for q in qset}
    negatives = {q: n_db - by_fold[lab[q][0]] + a.decoys for q in qset}
    total_pos = sum(poss_sf.values())
    total_neg = sum(negatives.values())
    nq = len(qset)

    hits, n_lines = parse_hits(a.m8, lab, qset)

    # ---- per-query sensitivity up to the first FP ----
    s_fam, s_sf, s_sfonly = [], [], []
    for q in qset:
        tp_f = tp_s = 0
        for ev, nb, s in hits.get(q, []):
            c = classify(q, s, lab)
            if c == "FP":
                break
            if c == "TP_family":
                tp_f += 1
            elif c == "TP_superfamily":
                tp_s += 1
        if poss_fam[q] > 0:
            s_fam.append(tp_f / poss_fam[q])
        if poss_sf[q] > 0:
            s_sf.append((tp_f + tp_s) / poss_sf[q])
        if poss_sfonly[q] > 0:
            s_sfonly.append(tp_s / poss_sfonly[q])

    # ---- pooled ROC ----
    pooled = []
    for q, lst in hits.items():
        for ev, nb, s in lst:
            c = classify(q, s, lab)
            if c in ("self", "ignore_fold"):
                continue
            pooled.append((ev, nb, 1 if c.startswith("TP") else 0))
    pooled.sort()
    tp = fp = 0
    curve = []  # (evalue, tp, fp)
    i = 0
    while i < len(pooled):
        ev = pooled[i][0]
        j = i
        while j < len(pooled) and pooled[j][0] == ev:
            if pooled[j][2]:
                tp += 1
            else:
                fp += 1
            j += 1
        curve.append((ev, tp, fp))
        i = j

    def recall_at_epq(epq):
        lim = epq * nq
        best = 0
        for ev, t, f in curve:
            if f <= lim:
                best = t
            else:
                break
        return best / total_pos if total_pos else 0.0

    def recall_at_fdp(fdp):
        best = 0
        for ev, t, f in curve:
            if t + f > 0 and f / (t + f) <= fdp:
                best = max(best, t)
        return best / total_pos if total_pos else 0.0

    def fdp_at_recall(r):
        # FDP at the most stringent E-value threshold whose recall reaches r
        need = r * total_pos
        for ev, t, f in curve:
            if t >= need:
                return f / (t + f) if t + f else None
        return None

    def at_evalue(cut):
        t = f = 0
        for ev, tt, ff in curve:
            if ev <= cut:
                t, f = tt, ff
            else:
                break
        return dict(evalue=cut, tp=t, fp=f,
                    recall=t / total_pos if total_pos else 0.0,
                    precision=(t / (t + f)) if t + f else None,
                    fdp=(f / (t + f)) if t + f else None,
                    fpr=f / total_neg if total_neg else None,
                    epq=f / nq if nq else None)

    out = dict(
        name=a.name, m8=a.m8, subset=a.subset, n_queries=nq, n_db=n_db + a.decoys, decoys=a.decoys,
        hit_lines=n_lines, possible_tp_pairs=total_pos, negative_pairs=total_neg,
        sens_to_1st_fp_family=sum(s_fam) / len(s_fam) if s_fam else 0.0,
        sens_to_1st_fp_superfamily=sum(s_sf) / len(s_sf) if s_sf else 0.0,
        sens_to_1st_fp_superfamily_only=sum(s_sfonly) / len(s_sfonly) if s_sfonly else 0.0,
        recall_at_epq={str(e): recall_at_epq(e) for e in (0.001, 0.01, 0.1, 1.0)},
        recall_at_fdp={str(d): recall_at_fdp(d) for d in (0.01, 0.05, 0.10)},
        fdp_at_recall={str(r): fdp_at_recall(r) for r in (0.02, 0.04, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12)},
        at_evalue=[at_evalue(c) for c in (1e-10, 1e-5, 1e-3, 1e-2, 1e-1, 1.0, 10.0)],
    )
    if a.roc:
        with open(a.roc, "w") as f:
            f.write("evalue\ttp\tfp\trecall\tepq\tfdp\n")
            for ev, t, fpp in curve:
                f.write(f"{ev:.3g}\t{t}\t{fpp}\t{t/total_pos:.6f}\t{fpp/nq:.6f}\t{(fpp/(t+fpp)) if t+fpp else 0:.6f}\n")
    json.dump(out, sys.stdout, indent=1)
    print()

if __name__ == "__main__":
    main()
