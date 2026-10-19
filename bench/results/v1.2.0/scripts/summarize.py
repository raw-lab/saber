#!/usr/bin/env python3
"""Summarise benchmark results: machine-readable summary.json, markdown tables
(summary.md) and plots (PNG) from the directories written by run_all.sh.

  python3 -I bench/scripts/summarize.py RESULTS_DIR
"""
import collections, json, math, os, statistics, sys

R = sys.argv[1]
LAB = None

def jl(path):
    out = []
    if os.path.exists(path):
        for line in open(path):
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out

def times(sec, extra=()):
    """label -> list of runs (duplicate records from REUSE=1 removed);
    later files (re-runs) override earlier ones."""
    d = collections.OrderedDict()
    for p in [os.path.join(R, sec, "times.jsonl")] + [os.path.join(R, sec, e) for e in extra]:
        recs = jl(p)
        seen = set()
        uniq = []
        for r in recs:
            key = (r["label"], r["rep"], r.get("timestamp"), r["wall_s"])
            if key not in seen:
                seen.add(key)
                uniq.append(r)
        for r in uniq:
            d.setdefault(r["label"], [])
        for l in list(d):
            rs = [r for r in uniq if r["label"] == l]
            if rs:
                d[l] = rs
    return d

def med(runs, key):
    v = [r[key] for r in runs if r.get("exit", 0) == 0]
    return statistics.median(v) if v else None

def rng(runs, key):
    v = [r[key] for r in runs if r.get("exit", 0) == 0]
    return (min(v), max(v)) if v else (None, None)

def metrics(sec, name, sub):
    p = os.path.join(R, sec, f"{name}.{sub}.metrics.json")
    return json.load(open(p)) if os.path.exists(p) and os.path.getsize(p) > 0 else None

def fmt_t(runs):
    if not runs:
        return "–"
    m = med(runs, "wall_s")
    if m is None:
        ex = runs[0].get("exit")
        return f"failed (exit {ex})"
    lo, hi = rng(runs, "wall_s")
    return f"{m:.1f}" + (f" ({lo:.1f}–{hi:.1f})" if len(runs) > 1 else "")

def fmt_mb(runs):
    m = med(runs, "maxrss_kb") if runs else None
    if m is None and runs:
        m = max(r["maxrss_kb"] for r in runs)
    return f"{m / 1024:.0f}" if m else "–"

def f3(x, n=3):
    return "–" if x is None else f"{x:.{n}f}"

summary = {}
md = []
HOSTS = {}

def host_line(sec, T):
    """Which machine(s) the section's runs come from (`host` label when the
    records carry one, else the runbench `cpu` field)."""
    hs = sorted({(r.get("host", ""), r.get("cpu", "unknown")) for runs in T.values() for r in runs})
    names = [f"host {h} ({c})" if h else c for h, c in hs]
    HOSTS[sec] = names
    if len(names) == 1:
        return f"Measured on: {names[0]}.\n"
    return "**Measured on more than one machine — compare rows only within one host:** " + "; ".join(names) + ".\n"

# ---------------------------------------------------------------- SCOP40
T = times("scop", extra=("rerun.jsonl",))
order = (["saber_v1.0.2", "saber_v1.1.0"]
         + [f"saber_v1.2_{p}" for p in ("very-fast", "fast", "default", "sensitive", "more-sensitive", "legacy", "exhaustive")]
         + [f"diamond_{m}" for m in ("fast", "default", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive")])
rows = []
for n in order:
    if n not in T and not metrics("scop", n, "test"):
        continue
    row = dict(name=n, runs=len(T.get(n, [])), wall_s=med(T.get(n, []), "wall_s"), wall_range=rng(T.get(n, []), "wall_s"),
               cpu_s=med(T.get(n, []), "cpu_s"), rss_mb=(med(T.get(n, []), "maxrss_kb") or 0) / 1024)
    for sub in ("tune", "test"):
        m = metrics("scop", n, sub)
        if m:
            row[sub] = dict(
                sens_sf=m["sens_to_1st_fp_superfamily"], sens_fam=m["sens_to_1st_fp_family"],
                recall_epq001=m["recall_at_epq"]["0.01"], recall_epq01=m["recall_at_epq"]["0.1"],
                recall_fdp05=m["recall_at_fdp"]["0.05"], fdp_at_recall=m.get("fdp_at_recall", {}),
                fp_e1e3=m["at_evalue"][2]["fp"], tp_e1e3=m["at_evalue"][2]["tp"],
                epq_e1=m["at_evalue"][5]["epq"], fdp_e1e3=m["at_evalue"][2]["fdp"])
    rows.append(row)
summary["scop"] = rows
md.append("## SCOP40 all-vs-all (11,206 domains, 2 threads, E ≤ 10, no hit limit)\n")
md.append(host_line("scop", T))
md.append("Sensitivity columns are on the **held-out test split**; tune-split sensitivity in brackets. Time = wall seconds, median of runs (range); RSS = peak resident memory.\n")
md.append("| Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP (superfamily) | Sens. (family) | Recall @ EPQ 0.01 | Recall @ FDP 5 % | TP at E ≤ 1e-3 | FP at E ≤ 1e-3 |")
md.append("|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    te, tu = r.get("test", {}), r.get("tune", {})
    md.append(f"| {r['name']} | {fmt_t(T.get(r['name'], []))} | {f3(r['cpu_s'], 1)} | {r['rss_mb']:.0f} | "
              f"{f3(te.get('sens_sf'))} [{f3(tu.get('sens_sf'))}] | {f3(te.get('sens_fam'))} | {f3(te.get('recall_epq001'), 4)} | "
              f"{f3(te.get('recall_fdp05'), 4)} | {te.get('tp_e1e3', '–')} | {te.get('fp_e1e3', '–')} |")
md.append("")
md.append("TP / FP at E ≤ 1e-3 = same-superfamily / different-fold pairs reported by the test-split queries at that cut-off. E-values are not interchangeable between tools: the same cut-off yields more true and more false positives from DIAMOND than from the SABER preset of similar sensitivity to the first FP.\n")

# FDP at matched recall
md.append("### False-discovery proportion at matched recall (test split)\n")
FDP_GRID = ("0.06", "0.08", "0.09", "0.1", "0.11")
md.append("Pooled over all test-split queries: hits sorted by E-value, FDP = FP / (TP + FP) at the most stringent cut-off that reaches the recall (TP / all same-superfamily pairs). 'not reached' = the tool never reaches that recall at E ≤ 10.\n")
md.append("| Tool / preset | " + " | ".join(f"FDP @ recall {x}" for x in FDP_GRID) + " |")
md.append("|---" * (len(FDP_GRID) + 1) + "|")
for r in rows:
    fr = r.get("test", {}).get("fdp_at_recall", {})
    md.append(f"| {r['name']} | " + " | ".join(f3(fr.get(x), 4) if fr.get(x) is not None else "not reached" for x in FDP_GRID) + " |")
md.append("")

# practical settings
prac = [(n, v) for n, v in T.items() if n.startswith("practical_")]
if prac:
    md.append("### Practical settings (E ≤ 1e-3, 25 targets per query), SCOP40 all-vs-all\n")
    md.append("| Run | Time s | RSS MB |")
    md.append("|---|---|---|")
    for n, v in prac:
        md.append(f"| {n.replace('practical_', '')} | {fmt_t(v)} | {fmt_mb(v)} |")
    md.append("")
    summary["scop_practical"] = {n: dict(wall_s=med(v, "wall_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024) for n, v in prac}

# index build
idx = [(n, v) for n, v in T.items() if "makedb" in n or "_db_" in n]
if idx:
    md.append("### SCOP40 index build and persistent-index search\n")
    md.append("| Run | Time s | RSS MB |")
    md.append("|---|---|---|")
    for n, v in idx:
        md.append(f"| {n} | {fmt_t(v)} | {fmt_mb(v)} |")
    md.append("")

# ---------------------------------------------------------------- decoys
LABELS = {}
_lp = os.environ.get("SCOP_LABELS", "/home/claude/bench/data/scop40.labels.tsv")
if os.path.exists(_lp):
    for _line in open(_lp):
        _p = _line.rstrip("\n").split("\t")
        if _p[0] != "domain" and len(_p) >= 6:
            LABELS[_p[0]] = _p[2]  # fold

def hits_per_query(path, cuts=(1e-3, 1e-2, 0.1, 1.0, 10.0), nq=11206):
    """False hits per query against a decoy database. Hits to decoys derived
    from a sequence of the query's own fold are excluded (a windowed shuffle of
    a homolog keeps much of its similarity, so it is not a clean negative)."""
    c = [0] * len(cuts)
    excluded = 0
    if not os.path.exists(path):
        return None, None
    for line in open(path):
        p = line.split("\t")
        if len(p) < 12:
            continue
        try:
            e = float(p[10])
        except ValueError:
            continue
        q = p[0].split()[0]
        src = p[1].split()[0].rsplit("_", 1)[0]
        if q in LABELS and src in LABELS and LABELS[q] == LABELS[src]:
            excluded += 1
            continue
        for i, cut in enumerate(cuts):
            if e <= cut:
                c[i] += 1
    return [x / nq for x in c], excluded

dec = {}
D = times("decoy")
for n in sorted(D):
    h, ex = hits_per_query(os.path.join(R, "decoy", n + ".m8"))
    dec[n] = dict(hits_per_query=h, excluded_same_fold_decoy_hits=ex, wall_s=med(D[n], "wall_s"))
if dec:
    summary["decoy"] = dec
    md.append("## Null behaviour: real SCOP40 queries against shuffled databases\n")
    md.append(host_line("decoy", D))
    md.append("Real queries against composition-preserving decoys (`shuf`: whole-sequence shuffle; `wshuf10`: shuffle inside 10-residue windows). Hits to decoys derived from a sequence of the query's own fold are excluded (a windowed shuffle of a homolog keeps part of its similarity); every remaining hit is a false positive. A calibrated E-value gives at most about E false hits per query (EPQ ≤ E).\n")
    md.append("| Run | EPQ at E ≤ 1e-3 | E ≤ 0.01 | E ≤ 0.1 | E ≤ 1 | E ≤ 10 |")
    md.append("|---|---|---|---|---|---|")
    for n, v in dec.items():
        h = v["hits_per_query"]
        if h:
            md.append(f"| {n} | " + " | ".join(f"{x:.4f}" for x in h) + " |")
    md.append("")

# ---------------------------------------------------------------- MMseqs
def pairs(path, cut=1e-3):
    s = set()
    if not os.path.exists(path):
        return None
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < 12:
            continue
        try:
            if float(p[10]) <= cut:
                s.add((p[0].split()[0], p[1].split()[0]))
        except ValueError:
            pass
    return s

STRONG = 1e-10  # reference pairs far below the reporting cut-off

M = times("mmseqs")
ref = pairs(os.path.join(R, "mmseqs", "saber_v1.2_exhaustive.m8"))
ref_strong = pairs(os.path.join(R, "mmseqs", "saber_v1.2_exhaustive.m8"), STRONG)
mm = []
for n in (["saber_v1.0.2", "saber_v1.1.0"] + [f"saber_v1.2_{p}" for p in ("very-fast", "fast", "default", "default_db", "default_nofloor", "sensitive", "exhaustive")]
          + [f"diamond_{m}" for m in ("fast", "default", "mid-sensitive", "sensitive", "very-sensitive")]):
    if n not in M:
        continue
    p = pairs(os.path.join(R, "mmseqs", n + ".m8"))
    rec = (len(p & ref) / len(ref)) if (p is not None and ref) else None
    rec_strong = (len(p & ref_strong) / len(ref_strong)) if (p is not None and ref_strong) else None
    mm.append(dict(name=n, wall_s=med(M[n], "wall_s"), rss_mb=(med(M[n], "maxrss_kb") or 0) / 1024,
                   pairs=len(p) if p is not None else None, recall_vs_exhaustive=rec,
                   recall_strong_pairs=rec_strong,
                   extra_pairs=len(p - ref) if (p is not None and ref) else None))
if mm:
    summary["mmseqs"] = dict(reference_pairs=len(ref) if ref else None,
                             strong_reference_pairs=len(ref_strong) if ref_strong else None, runs=mm)
    md.append(f"## MMseqs2 example workload (500 queries × 20,000 UniProt sequences, E ≤ 1e-3)\n")
    md.append(host_line("mmseqs", M))
    md.append(f"Recall is measured against SABER's own exhaustive Smith–Waterman search (`--preset exhaustive`: {len(ref) if ref else '–'} query–subject pairs at E ≤ 1e-3). That reference shares SABER's scoring, composition correction and E-value statistics, so near the cut-off it favours SABER; 'Recall, E ≤ 1e-10 pairs' uses only the {len(ref_strong) if ref_strong else '–'} reference pairs with E ≤ 1e-10, which a tool that finds them reports at E ≤ 1e-3 whatever its statistics. 'Extra' = pairs a tool reports that the reference does not (mostly borderline E-values under different statistics). 74 of the 500 queries are also database members; self-hits count for every tool.\n")
    md.append("| Tool / preset | Time s | RSS MB | Pairs | Recall of exhaustive SW | Recall, E ≤ 1e-10 pairs | Extra pairs |")
    md.append("|---|---|---|---|---|---|---|")
    for r in mm:
        md.append(f"| {r['name']} | {fmt_t(M[r['name']])} | {r['rss_mb']:.0f} | {r['pairs']} | {f3(r['recall_vs_exhaustive'], 4)} | {f3(r['recall_strong_pairs'], 4)} | {r['extra_pairs']} |")
    md.append("")
    sc = [(n, v) for n, v in M.items() if n.startswith("scale")]
    if sc:
        md.append("### Scaling (single run each)\n")
        md.append("| Run | Time s | RSS MB |")
        md.append("|---|---|---|")
        for n, v in sc:
            md.append(f"| {n} | {fmt_t(v)} | {fmt_mb(v)} |")
        md.append("")
        summary["mmseqs_scaling"] = {n: dict(wall_s=med(v, "wall_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024) for n, v in sc}
    lg = [(n, v) for n, v in M.items() if n.startswith("long_")]
    if lg:
        lq = os.path.join(R, "mmseqs", "long_queries.fa")
        nlq = sum(1 for l in open(lq) if l.startswith(">")) if os.path.exists(lq) else "?"
        md.append(f"### Long proteins ({nlq} queries > 1,500 aa against the 20,000-sequence database)\n")
        md.append("The queries are every fifth database sequence longer than 1,500 residues, so each also finds itself (self-hits count for every tool). Single runs. The SABER 1.2 row is the benchmarked build, whose traceback memory cap lowered the score of 9 of these alignments; the delivered build restores them (sections 'Re-check' and 'Before / after').\n")
        md.append("| Run | Time s | RSS MB | Pairs at E ≤ 1e-3 |")
        md.append("|---|---|---|---|")
        for n, v in lg:
            p = pairs(os.path.join(R, "mmseqs", n + ".m8"))
            md.append(f"| {n.replace('long_', '')} | {fmt_t(v)} | {fmt_mb(v)} | {len(p) if p is not None else '–'} |")
        md.append("")

# ---------------------------------------------------------------- index
I = times("index")
if I:
    md.append("## Index construction, size and persistent-index search (MMseqs database and ×4 database)\n")
    md.append(host_line("index", I))
    md.append("| Run | Time s | RSS MB |")
    md.append("|---|---|---|")
    for n, v in I.items():
        md.append(f"| {n} | {fmt_t(v)} | {fmt_mb(v)} |")
    sz = os.path.join(R, "index", "sizes.txt")
    if os.path.exists(sz):
        md.append("\nFile sizes / checks:\n```\n" + open(sz).read().strip() + "\n```")
    md.append("")
    summary["index"] = {n: dict(wall_s=med(v, "wall_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024, runs=len(v)) for n, v in I.items()}

# ---------------------------------------------------------------- translated
TR = times("translated")
if TR:
    md.append("## Translated searches\n")
    md.append(host_line("translated", TR))
    md.append("blastx (`nx`): back-translated SCOP40 domains (random synonymous codons, random flanks and strand) against SCOP40 proteins; tblastn (`pn`): a fixed 1-in-37 subset of SCOP40 proteins (303 queries) against the back-translated set. `nxsub_*` / `pnsub_*`: a fixed 1-in-333 query subset (34 queries), used to compare with SABER 1.1.0, which searched translated data by full dynamic programming (its sensitivity numbers on 34 queries are indicative only). Metrics are over the queries of each run (both splits): superfamily-level sensitivity to the first FP, pooled recall at 0.01 false positives per query, and false positives at E ≤ 1e-3.\n")
    md.append("| Run | Queries | Time s | RSS MB | Sens. to 1st FP (SF) | Recall @ EPQ 0.01 | FP at E ≤ 1e-3 |")
    md.append("|---|---|---|---|---|---|---|")
    trs = {}
    tr_order = sorted(TR, key=lambda n: (n.split("_")[0] != "nx", n.split("_")[0] != "nxsub", n.split("_")[0] != "pn", n))
    for n in tr_order:
        v = TR[n]
        m = metrics("translated", n, "all")
        trs[n] = dict(wall_s=med(v, "wall_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024,
                      n_queries=m["n_queries"] if m else None,
                      sens_sf=m["sens_to_1st_fp_superfamily"] if m else None,
                      recall_epq001=m["recall_at_epq"]["0.01"] if m else None,
                      fp_e1e3=m["at_evalue"][2]["fp"] if m else None)
        md.append(f"| {n} | {trs[n]['n_queries'] or '–'} | {fmt_t(v)} | {fmt_mb(v)} | {f3(trs[n]['sens_sf'])} | {f3(trs[n]['recall_epq001'], 4)} | {trs[n]['fp_e1e3'] if trs[n]['fp_e1e3'] is not None else '–'} |")
    md.append("")
    summary["translated"] = trs

# ---------------------------------------------------------------- translated, matched sensitivity
TM = times("translated_matched")
if TM:
    md.append("### Translated subsets at matched sensitivity (follow-up run)\n")
    md.append(host_line("translated_matched", TM))
    md.append("The 34-query blastx / tblastn subsets above, searched again in one session on one machine after the main run (`bench/scripts/translated_matched.sh`): SABER 1.1.0 (full dynamic programming over six frames), SABER 1.2 `default`, `more-sensitive` and `exhaustive`, and DIAMOND (blastx only), so speed can be compared at similar sensitivity. Metrics as above (34 queries, indicative only); SABER 1.1.0 single runs, the others median of 3.\n")
    md.append("| Run | Time s | RSS MB | Sens. to 1st FP (SF) | Recall @ EPQ 0.01 | FP at E ≤ 1e-3 |")
    md.append("|---|---|---|---|---|---|")
    tm = {}
    for n in [f"{m}_{t}" for m in ("nxsub", "pnsub") for t in ("saber_v1.1.0", "saber_v1.2_default", "saber_v1.2_more-sensitive",
                                                                "saber_v1.2_exhaustive", "diamond_very-sensitive", "diamond_ultra-sensitive")]:
        if n not in TM:
            continue
        v = TM[n]
        m = metrics("translated_matched", n, "all")
        tm[n] = dict(wall_s=med(v, "wall_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024,
                     sens_sf=m["sens_to_1st_fp_superfamily"] if m else None,
                     recall_epq001=m["recall_at_epq"]["0.01"] if m else None,
                     fp_e1e3=m["at_evalue"][2]["fp"] if m else None)
        md.append(f"| {n} | {fmt_t(v)} | {fmt_mb(v)} | {f3(tm[n]['sens_sf'])} | {f3(tm[n]['recall_epq001'], 4)} | {tm[n]['fp_e1e3'] if tm[n]['fp_e1e3'] is not None else '–'} |")
    md.append("")
    summary["translated_matched"] = tm

# ---------------------------------------------------------------- re-check of A1 timings
RC = times("recheck")
if RC:
    md.append("## Re-check of the A1 SABER-vs-DIAMOND timings (follow-up run)\n")
    md.append(host_line("recheck", RC))
    md.append("In the A1 sections the DIAMOND rows were measured between 14:34 and 17:17 and the SABER 1.2 rows between 20:10 and 20:37. These comparisons were re-measured with the same commands, both tools in one session on one machine of A1's CPU model (`bench/scripts/recheck.sh`, SABER = the delivered build), next to the medians of the original sections. Absolute times differ between the sessions (another virtual machine instance or host load); compare the two tools within a column. The delivered build's output equals the benchmarked build's except on the long-protein set, where its traceback fallback changes 9 of 4,712 alignments (all to higher scores).\n")
    md.append("| Section: run | Time s (re-check) | RSS MB | Time s (original section) | RSS MB |")
    md.append("|---|---|---|---|---|")
    rc = {}
    for n, v in RC.items():
        sec, lab = n.split(":", 1) if ":" in n else ("scop", n)
        o = times(sec).get(lab, [])
        rc[n] = dict(wall_s=med(v, "wall_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024, section_wall_s=med(o, "wall_s"))
        md.append(f"| {sec}: {lab} | {fmt_t(v)} | {fmt_mb(v)} | {fmt_t(o)} | {fmt_mb(o)} |")
    oc = os.path.join(R, "recheck", "outputs_check.txt")
    if os.path.exists(oc):
        md.append("")
        md.append("Outputs: " + "; ".join(l.strip() for l in open(oc) if l.strip()) + ".")
    md.append("")
    summary["recheck"] = rc

# ---------------------------------------------------------------- development build comparison
DV = times("devcheck")
if DV:
    md.append("## Before / after the performance changes (follow-up run)\n")
    md.append(host_line("devcheck", DV))
    md.append("A development build from before the index-construction, sparse-offset, 16-bit traceback, filter-first and X-drop changes (`_old`, sha256 63442e7d…, full hash in `env/final_binary.txt`) against the benchmarked build (`_new`) on one machine in one session; `--stats` stage times are summed over threads.\n")
    md.append("| Run | Time s | RSS MB | traceback stage s | seed stage s |")
    md.append("|---|---|---|---|---|")
    dv = {}
    for n, v in DV.items():
        st = os.path.join(R, "devcheck", n + ".stats.json")
        c = json.load(open(st)).get("counters", {}) if os.path.exists(st) else {}
        tb, sd = c.get("t_traceback_ns"), c.get("t_seed_ns")
        dv[n] = dict(wall_s=med(v, "wall_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024, t_traceback_s=tb / 1e9 if tb else None, t_seed_s=sd / 1e9 if sd else None)
        md.append(f"| {n} | {fmt_t(v)} | {fmt_mb(v)} | {f3(dv[n]['t_traceback_s'], 2)} | {f3(dv[n]['t_seed_s'], 2)} |")
    oc = os.path.join(R, "devcheck", "outputs_check.txt")
    if os.path.exists(oc):
        md.append("")
        md.append("Outputs of the development build compared with the released outputs (first = development build):\n")
        for l in open(oc):
            if l.strip():
                md.append(f"* {l.strip()}")
    md.append("")
    summary["devcheck"] = dv

# ---------------------------------------------------------------- low complexity
L = times("lowc")
if L:
    md.append("## Adversarial low-complexity / repeat inputs\n")
    md.append(host_line("lowc", L))
    md.append("100 homopolymer / tandem-repeat queries (some with an embedded real domain) against the MMseqs database plus 2,000 low-complexity sequences. `mmq_*`: the 500 real MMseqs queries against the same polluted database — hits to the artificial `lowcdb*` sequences are false positives. SABER 1.0.2/1.1.0 ran under a 6 GB address-space limit; the adversarial `lowc10` runs and the 1.x runs under a 10-minute time limit (exit 124 = stopped by the time limit; exit −6 = aborted after a failed memory allocation, e.g. `memory allocation of 4294967296 bytes failed`).\n")
    md.append("| Run | Time s | RSS MB | Exit | Hits (E ≤ 1e-3) | Hits to artificial low-complexity sequences |")
    md.append("|---|---|---|---|---|---|")
    lc = {}
    for n, v in L.items():
        path = os.path.join(R, "lowc", n + ".m8")
        tot = art = 0
        if os.path.exists(path):
            for line in open(path):
                p = line.split("\t")
                if len(p) >= 12:
                    try:
                        if float(p[10]) <= 1e-3:
                            tot += 1
                            art += p[1].startswith("lowcdb")
                    except ValueError:
                        pass
        ex = v[0].get("exit")
        lc[n] = dict(wall_s=med(v, "wall_s") if ex == 0 else v[0]["wall_s"], rss_mb=max(r["maxrss_kb"] for r in v) / 1024, exit=ex,
                     hits=tot if ex == 0 else None, artificial_hits=art if ex == 0 else None)
        hits_s = f"{tot} | {art}" if ex == 0 else "– (aborted) | –"
        md.append(f"| {n} | {lc[n]['wall_s']:.1f} | {lc[n]['rss_mb']:.0f} | {ex} | {hits_s} |")
    md.append("")
    summary["lowc"] = lc
    for n in ("lowc_saber_v1.2_default", "lowc_saber_v1.2_default_nomask", "lowc10_saber_v1.2_legacy_nomask_nocap"):
        p = os.path.join(R, "lowc", n + ".stats.json")
        if os.path.exists(p):
            s = json.load(open(p))
            c = s["counters"]
            summary.setdefault("lowc_counters", {})[n] = {k: c[k] for k in ("hits", "subject_range_passes", "streamed_subjects", "seeds_frequent_skipped", "query_residues_masked", "regions_dropped_by_cap")}

# ---------------------------------------------------------------- multi / fragments
MU = times("multi")
if MU:
    md.append("## Two-domain queries and short fragments\n")
    md.append(host_line("multi", MU))
    md.append("Retention = fraction of a domain's true positives from its single-domain search (each tool's own, E ≤ 1e-3) that are still found when the domain is part of a two-domain query; the TP columns give the absolute numbers found in the two-domain searches. Single runs.\n")
    md.append("| Tool / preset | Retention of the weaker domain's TPs | Retention (stronger) | Weaker-domain TPs (two-domain / alone) | Stronger-domain TPs (two-domain / alone) | FP in two-domain searches (E ≤ 1e-3) | Fragments: sens. to 1st FP (SF, test) | Fragment search time s |")
    md.append("|---|---|---|---|---|---|---|---|")
    mu = {}
    for t in ("saber_v1.1.0", "saber_v1.2_default", "saber_v1.2_more-sensitive", "diamond_very-sensitive", "diamond_ultra-sensitive"):
        p = os.path.join(R, "multi", f"multi_{t}.metrics.json")
        m = json.load(open(p)) if os.path.exists(p) and os.path.getsize(p) else {}
        fr = metrics("multi", f"frag_{t}", "test")
        mu[t] = dict(retention_weaker=m.get("retention_weaker"), retention_stronger=m.get("retention_stronger"), fp=m.get("fp_multidomain"),
                     weaker_tp=(m.get("weaker_tp_in_multidomain"), m.get("weaker_tp_alone")),
                     stronger_tp=(m.get("stronger_tp_in_multidomain"), m.get("stronger_tp_alone")),
                     frag_sens_sf=fr["sens_to_1st_fp_superfamily"] if fr else None, frag_wall_s=med(MU.get(f"frag_{t}", []), "wall_s"))
        md.append(f"| {t} | {f3(mu[t]['retention_weaker'])} | {f3(mu[t]['retention_stronger'])} | {mu[t]['weaker_tp'][0]} / {mu[t]['weaker_tp'][1]} | {mu[t]['stronger_tp'][0]} / {mu[t]['stronger_tp'][1]} | {mu[t]['fp'] if mu[t]['fp'] is not None else '–'} | {f3(mu[t]['frag_sens_sf'])} | {f3(mu[t]['frag_wall_s'], 1)} |")
    md.append("")
    summary["multi"] = mu

# ---------------------------------------------------------------- threads
TH = times("threads")
if TH:
    md.append("## Thread scaling (2-core machine)\n")
    md.append(host_line("threads", TH))
    md.append("| Run | Time s | CPU s | RSS MB |")
    md.append("|---|---|---|---|")
    for n, v in TH.items():
        md.append(f"| {n} | {fmt_t(v)} | {f3(med(v, 'cpu_s'), 1)} | {fmt_mb(v)} |")
    md.append("")
    summary["threads"] = {n: dict(wall_s=med(v, "wall_s"), cpu_s=med(v, "cpu_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024) for n, v in TH.items()}

# ---------------------------------------------------------------- sweep (tune split)
SW = times("sweep")
if SW:
    md.append("## Preset selection sweep (SCOP40 **tune** split only)\n")
    md.append(host_line("sweep", SW))
    md.append("Single runs, 2 threads. These numbers chose the presets; the test-split results above are the held-out evaluation.\n")
    md.append("| Configuration | Time s | Sens. to 1st FP (SF, tune) | Sens. (family) | Recall @ EPQ 0.01 | FP at E ≤ 1e-3 |")
    md.append("|---|---|---|---|---|---|")
    sw = {}
    for n, v in SW.items():
        m = metrics("sweep", n, "tune")
        if not m:
            continue
        sw[n] = dict(wall_s=med(v, "wall_s"), sens_sf=m["sens_to_1st_fp_superfamily"], sens_fam=m["sens_to_1st_fp_family"],
                     recall_epq001=m["recall_at_epq"]["0.01"], fp_e1e3=m["at_evalue"][2]["fp"])
    for n, r in sorted(sw.items(), key=lambda kv: kv[1]["wall_s"] or 0):
        md.append(f"| {n} | {f3(r['wall_s'], 1)} | {f3(r['sens_sf'])} | {f3(r['sens_fam'])} | {f3(r['recall_epq001'], 4)} | {r['fp_e1e3']} |")
    md.append("")
    summary["sweep"] = sw

# ---------------------------------------------------------------- plots
# Categorical slots (validated set, light surface): blue, orange, aqua, yellow;
# baselines in neutral grey. Text is neutral ink; series are also labelled
# directly, so identity never depends on colour alone.
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    BLUE, ORANGE, AQUA, YELLOW, GREY, INK, INK2 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a8984", "#0b0b0b", "#52514e"
    plt.rcParams.update({"axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                         "axes.titlecolor": INK, "grid.color": "#d9d8d4", "grid.linewidth": 0.6, "font.size": 9})

    def family(n):
        return BLUE if n.startswith("saber_v1.2") else (ORANGE if n.startswith("diamond") else GREY)

    def short(n):
        return n.replace("saber_v1.2_", "").replace("saber_", "SABER ").replace("diamond_", "D ")

    if rows:
        fig, ax = plt.subplots(figsize=(7.5, 5))
        # label offsets (points) where neighbours would collide
        off = {"saber_v1.2_default": (-6, 4, "right"), "saber_v1.2_sensitive": (6, 2, "left"),
               "diamond_sensitive": (0, -13, "center"), "diamond_more-sensitive": (6, -2, "left"),
               "diamond_very-sensitive": (6, -8, "left"), "saber_v1.2_legacy": (6, -8, "left")}
        for r in rows:
            if r.get("wall_s") and r.get("test"):
                x, y = r["wall_s"], r["test"]["sens_sf"]
                ax.scatter(x, y, color=family(r["name"]), s=48, zorder=3, edgecolors="white", linewidths=1.5)
                dx, dy, ha = off.get(r["name"], (6, 2, "left"))
                ax.annotate(short(r["name"]), (x, y), fontsize=7.5, xytext=(dx, dy), textcoords="offset points", ha=ha, color=INK)
        ax.set_xscale("log")
        ax.set_xlabel("wall time, s (log scale; SCOP40 all-vs-all, 2 threads, machine A1)")
        ax.set_ylabel("sensitivity to the first FP, superfamily (test split)")
        ax.set_title("Runtime vs sensitivity")
        for lab, c in (("SABER 1.2", BLUE), ("DIAMOND 2.2.8", ORANGE), ("SABER 1.0.2 / 1.1.0", GREY)):
            ax.scatter([], [], color=c, s=48, label=lab)
        ax.legend(fontsize=8, frameon=False, loc="lower right")
        ax.grid(True)
        fig.tight_layout()
        fig.savefig(os.path.join(R, "scop_time_vs_sensitivity.png"), dpi=130)
        plt.close(fig)

        # Pooled ROC: hue = sensitivity tier, line style = tool.
        fig, ax = plt.subplots(figsize=(7.5, 5))
        series = (("saber_v1.1.0", GREY, ":"), ("saber_v1.2_fast", BLUE, "-"), ("diamond_sensitive", BLUE, "--"),
                  ("saber_v1.2_default", ORANGE, "-"), ("diamond_very-sensitive", ORANGE, "--"),
                  ("saber_v1.2_more-sensitive", AQUA, "-"), ("diamond_ultra-sensitive", AQUA, "--"))
        for n, c, ls in series:
            p = os.path.join(R, "scop", f"{n}.test.roc.tsv")
            if not os.path.exists(p):
                continue
            xs, ys = [], []
            for line in list(open(p))[1:]:
                ev, tp, fp, rec, epq, fdp = line.split("\t")
                if float(epq) > 0:
                    xs.append(float(epq)); ys.append(float(rec))
            if not xs:
                continue
            ax.plot(xs, ys, color=c, ls=ls, lw=2, label=short(n).replace("D ", "DIAMOND "))
            ax.annotate(short(n), (xs[-1], ys[-1]), fontsize=7.5, xytext=(4, 0), textcoords="offset points", va="center", color=INK)
        ax.set_xscale("log")
        ax.set_xlim(1e-3, 30)
        ax.set_xlabel("false positives per query (EPQ, pooled over test-split queries)")
        ax.set_ylabel("recall (same-superfamily pairs)")
        ax.set_title("Pooled ROC, SCOP40 test split (solid: SABER 1.2, dashed: DIAMOND; colour = matched tier)")
        ax.legend(fontsize=7.5, frameon=False, loc="upper left")
        ax.grid(True)
        fig.tight_layout()
        fig.savefig(os.path.join(R, "scop_roc_test.png"), dpi=130)
        plt.close(fig)

    if dec:
        cuts = [1e-3, 1e-2, 0.1, 1.0, 10.0]
        panels = (("SABER", (("saber_v1.1.0", GREY, ":", "1.1.0"), ("saber_v1.2_fast", BLUE, "-", "1.2 fast"),
                             ("saber_v1.2_default", ORANGE, "-", "1.2 default"), ("saber_v1.2_more-sensitive", AQUA, "-", "1.2 more-sensitive"))),
                  ("DIAMOND 2.2.8", (("diamond_default", YELLOW, "--", "default"), ("diamond_very-sensitive", ORANGE, "--", "very-sensitive"),
                                     ("diamond_ultra-sensitive", AQUA, "--", "ultra-sensitive"))))
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
        for ax, (title, ser) in zip(axes, panels):
            for n, c, ls, lab in ser:
                v = dec.get(n + "_shuf")
                if not v or not v["hits_per_query"]:
                    continue
                ys = [max(x, 1e-5) for x in v["hits_per_query"]]
                ax.plot(cuts, ys, marker="o", ms=6, lw=2, color=c, ls=ls, label=lab, markeredgecolor="white", markeredgewidth=1.2)
                nudge = {"1.2 more-sensitive": 10, "1.1.0": 3, "1.2 default": -6, "1.2 fast": -7}.get(lab, 0)
                ax.annotate(lab, (cuts[-1], ys[-1]), fontsize=7.5, xytext=(6, nudge), textcoords="offset points", va="center", color=INK)
            ax.plot(cuts, cuts, color=INK2, ls=(0, (1, 2)), lw=1, label="EPQ = E (calibrated limit)")
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlim(5e-4, 80)
            ax.set_xlabel("E-value cut-off")
            ax.set_title(title)
            ax.legend(fontsize=7.5, frameon=False, loc="upper left")
            ax.grid(True)
        axes[0].set_ylabel("false hits per query (whole-sequence shuffles; 1e-5 = none)")
        fig.suptitle("Real SCOP40 queries against shuffled decoys", color=INK)
        fig.tight_layout()
        fig.savefig(os.path.join(R, "decoy_calibration.png"), dpi=130)
        plt.close(fig)
except Exception as e:  # plotting is optional
    print("plotting skipped:", e, file=sys.stderr)

summary["hosts"] = HOSTS
json.dump(summary, open(os.path.join(R, "summary.json"), "w"), indent=1)
open(os.path.join(R, "summary.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md))
