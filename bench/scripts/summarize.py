#!/usr/bin/env python3
"""Summarise benchmark results: machine-readable summary.json, markdown tables
(summary.md) and plots (PNG) from the directories written by run_all.sh
(SABER 1.4 suite; the summaries of the 1.2.0 and 1.3.0 suites are kept in
results/v1.2.0/ and results/v1.3.0/).

  python3 -I bench/scripts/summarize.py RESULTS_DIR
"""
import collections, json, os, statistics, sys

R = sys.argv[1]
NEW = "saber_v1.4"
OLD = "saber_v1.3"
NV, OV = "1.4", "1.3.0"
PRESETS = ("very-fast", "fast", "default", "sensitive", "more-sensitive")
DMODES = ("fast", "default", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive")


def jl(path):
    out = []
    if os.path.exists(path):
        for line in open(path):
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def times(sec):
    """label -> list of runs, in file order."""
    d = collections.OrderedDict()
    for r in jl(os.path.join(R, sec, "times.jsonl")):
        d.setdefault(r["label"], []).append(r)
    return d


def ok(runs):
    return [r for r in runs if r.get("exit", 0) == 0]


def med(runs, key):
    v = [r[key] for r in ok(runs)]
    return statistics.median(v) if v else None


def rng(runs, key):
    v = [r[key] for r in ok(runs)]
    return (min(v), max(v)) if v else (None, None)


def metrics(sec, name, sub):
    p = os.path.join(R, sec, f"{name}.{sub}.metrics.json")
    return json.load(open(p)) if os.path.exists(p) and os.path.getsize(p) > 0 else None


def fmt_t(runs):
    if not runs:
        return "–"
    m = med(runs, "wall_s")
    if m is None:
        return f"failed (exit {runs[0].get('exit')})"
    lo, hi = rng(runs, "wall_s")
    return f"{m:.1f}" + (f" ({lo:.1f}–{hi:.1f})" if len(ok(runs)) > 1 else "")


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
# Machines named A, B, ... in the order of their first timing record (a
# suite whose virtual machine was replaced mid-run has more than one).
_first = {}
for _dirpath, _dirs, _files in os.walk(R):
    if "times.jsonl" in _files:
        for _r in jl(os.path.join(_dirpath, "times.jsonl")):
            _c, _t = _r.get("cpu", "unknown"), _r.get("timestamp", "")
            if _c not in _first or _t < _first[_c]:
                _first[_c] = _t
MACHINE = {c: chr(ord("A") + i) for i, c in enumerate(sorted(_first, key=lambda c: _first[c]))}


def machines(runs):
    return "".join(sorted({MACHINE.get(r.get("cpu", "unknown"), "?") for r in runs}))


def host_line(sec, T):
    cpus = sorted({r.get("cpu", "unknown") for runs in T.values() for r in runs}, key=lambda c: MACHINE.get(c, "?"))
    HOSTS[sec] = cpus
    if len(cpus) == 1:
        return f"Measured on machine {MACHINE.get(cpus[0], '?')}: {cpus[0]}.\n"
    return ("**Measured on more than one machine — compare rows only within one machine** (marked after the name): "
            + "; ".join(f"{MACHINE.get(c, '?')} = {c}" for c in cpus) + ".\n")


def mlabel(sec, name, runs):
    """`name`, with its machine when the section ran on more than one."""
    return f"{name} [{machines(runs)}]" if len(HOSTS.get(sec, [])) > 1 and runs else name


def acc(m):
    if not m:
        return {}
    return dict(sens_sf=m["sens_to_1st_fp_superfamily"], sens_fam=m["sens_to_1st_fp_family"],
                recall_epq001=m["recall_at_epq"]["0.01"], recall_epq01=m["recall_at_epq"]["0.1"],
                recall_fdp05=m["recall_at_fdp"]["0.05"], fdp_at_recall=m.get("fdp_at_recall", {}),
                tp_e1e3=m["at_evalue"][2]["tp"], fp_e1e3=m["at_evalue"][2]["fp"],
                epq_e1=m["at_evalue"][5]["epq"])


def tool_rows(sec, T, order):
    rows = []
    for n in order:
        if n not in T and not metrics(sec, n, "test"):
            continue
        # timings of a persistent-index search are recorded as `<name>_warm`
        v = T.get(n) or T.get(n + "_warm", [])
        rows.append(dict(name=n, runs=len(ok(v)), wall_s=med(v, "wall_s"), wall_range=rng(v, "wall_s"),
                         cpu_s=med(v, "cpu_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024,
                         tune=acc(metrics(sec, n, "tune")), test=acc(metrics(sec, n, "test"))))
    return rows


def acc_table(rows, T, extra_note="", sec=None):
    md.append("| Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP (superfamily) | Recall @ EPQ 0.01 | Recall @ EPQ 0.1 | Recall @ FDP 5 % | TP at E ≤ 1e-3 | FP at E ≤ 1e-3 |")
    md.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        te, tu = r["test"], r["tune"]
        runs = T.get(r['name']) or T.get(r['name'] + '_warm', [])
        md.append(f"| {mlabel(sec, r['name'], runs) if sec else r['name']} | {fmt_t(runs)} | {f3(r['cpu_s'], 1)} | {r['rss_mb']:.0f} | "
                  f"{f3(te.get('sens_sf'), 4)} [{f3(tu.get('sens_sf'), 4)}] | {f3(te.get('recall_epq001'), 4)} [{f3(tu.get('recall_epq001'), 4)}] | "
                  f"{f3(te.get('recall_epq01'), 4)} | {f3(te.get('recall_fdp05'), 4)} | {te.get('tp_e1e3', '–')} | {te.get('fp_e1e3', '–')} |")
    md.append("")
    if extra_note:
        md.append(extra_note + "\n")


def fdp_table(rows, grid):
    md.append("| Tool / preset | " + " | ".join(f"FDP @ recall {x}" for x in grid) + " |")
    md.append("|---" * (len(grid) + 1) + "|")
    for r in rows:
        fr = r["test"].get("fdp_at_recall", {})
        md.append(f"| {r['name']} | " + " | ".join(f3(fr.get(x), 4) if fr.get(x) is not None else "not reached" for x in grid) + " |")
    md.append("")


def simple_table(title, items):
    if not items:
        return
    if title:
        md.append(title + "\n")
    md.append("| Run | Time s | CPU s | RSS MB |")
    md.append("|---|---|---|---|")
    for n, v in items:
        md.append(f"| {n} | {fmt_t(v)} | {f3(med(v, 'cpu_s'), 1)} | {fmt_mb(v)} |")
    md.append("")


ACC_NOTE = ("Accuracy columns are on the **held-out test split** (tune split in brackets): sensitivity = fraction of a query's "
            "same-superfamily domains ranked above its first false positive, averaged over queries; recall @ EPQ x = pooled "
            "fraction of same-superfamily pairs found before x false positives per query; TP / FP = same-superfamily / "
            "different-fold (or decoy) pairs reported at E ≤ 1e-3. Time = wall seconds, median (range) of the runs; RSS = peak "
            "resident memory (median).")

# ---------------------------------------------------------------- at scale
ST = times("scale")
scale_rows = []
if ST:
    scale_rows = tool_rows("scale", ST, [f"{NEW}_{p}" for p in PRESETS] + [f"diamond_{m}" for m in DMODES]
                           + ["diamond_very-sensitive_c1", f"{NEW}_default_dbmodel", f"{OLD}_default"])
    summary["scale"] = scale_rows
    md.append("## At scale: SCOP40 against SCOP40 + 80,000 reversed proteins (91,206 sequences, 38.2 M residues)\n")
    md.append(host_line("scale", ST))
    singles = [r["name"] for r in scale_rows if r["runs"] == 1]
    md.append("11,206 SCOP40 queries, 2 threads, E ≤ 10, no hit limit; every hit to a reversed protein (`REV_*`) or to a domain of "
              "another fold is a false positive. Prebuilt indexes for both tools (`saber --makedb`, `diamond makedb`; built once, "
              f"timed below). SABER {NV} `default` and DIAMOND `--very-sensitive` were run alternately "
              f"({len(ok(ST.get(NEW + '_default', [])))} runs each) so that machine drift affects both alike; the other rows are "
              "median of 3" + (" (single runs: " + ", ".join(f"`{n}`" for n in singles) + ")" if singles else "") + ". "
              + ACC_NOTE + "\n")
    acc_table(scale_rows, ST, "`diamond_very-sensitive_c1` = DIAMOND `--very-sensitive -c1` (one seed-index chunk, DIAMOND's "
              f"speed/memory option); `{NEW}_default_dbmodel` = SABER {NV} `default` with `--evalue-model database` "
              "(database-wide E-values, as in 1.2: same alignments and bit scores, different E-values and therefore a different "
              "ranking across subjects). Both are run by `run_all.sh scale_extra`.", sec="scale")
    md.append("### False-discovery proportion at matched recall (test split, at scale)\n")
    fdp_table(scale_rows, ("0.06", "0.07", "0.08", "0.09", "0.1"))
    pair = [n for n in (f"{NEW}_default", "diamond_very-sensitive") if n in ST]
    if len(pair) == 2:
        a, b = ST[pair[0]], ST[pair[1]]
        wa, wb = [r["wall_s"] for r in ok(a)], [r["wall_s"] for r in ok(b)]
        ca, cb = [r["cpu_s"] for r in ok(a)], [r["cpu_s"] for r in ok(b)]
        paired = [x / y for x, y in zip(wa, wb)]
        summary["scale_pair"] = dict(saber_wall=wa, diamond_wall=wb, saber_cpu=ca, diamond_cpu=cb, ratio_wall=paired)
        md.append(f"Alternating runs, SABER {NV} `default` / DIAMOND `--very-sensitive` wall time: "
                  + ", ".join(f"{x:.1f} / {y:.1f} s" for x, y in zip(wa, wb))
                  + f" (ratio {min(paired):.2f}–{max(paired):.2f}); CPU time "
                  + ", ".join(f"{x:.0f} / {y:.0f} s" for x, y in zip(ca, cb)) + ".\n")
    simple_table("### Index construction (scale database)", [(n, v) for n, v in ST.items() if n.startswith("makedb")])
    sz = os.path.join(R, "scale", "sizes.txt")
    if os.path.exists(sz):
        md.append("```\n" + open(sz).read().strip() + "\n```\n")
    simple_table("### Practical settings at scale (E ≤ 1e-3, 25 targets per query)", [(n, v) for n, v in ST.items() if n.startswith("practical_")])

SA = times("scale_algo")
if SA:
    algo_rows = tool_rows("scale_algo", SA, [f"{NEW}_default", "diamond_very-sensitive_algo1"])
    md.append("### DIAMOND's seed-search algorithm at scale\n")
    md.append(host_line("scale_algo", SA))
    md.append("DIAMOND chooses between its double-indexed and query-indexed seed search automatically; here `--very-sensitive "
              f"--algo 1` (query-indexed) is forced, run alternately with SABER {NV} `default` (whose accuracy is in the table "
              "above; same output). " + ACC_NOTE + "\n")
    acc_table(algo_rows, SA, sec="scale_algo")
    a, b = ok(SA.get(f"{NEW}_default", [])), ok(SA.get("diamond_very-sensitive_algo1", []))
    if a and b:
        md.append("Alternating runs, SABER / DIAMOND `--algo 1` wall time: "
                  + ", ".join(f"{x['wall_s']:.1f} / {y['wall_s']:.1f} s" for x, y in zip(a, b)) + ".\n")
    summary["scale_algo"] = algo_rows

OO = times("oneoff")
if OO:
    md.append("### One-off search at scale, index construction included\n")
    md.append(host_line("oneoff", OO))
    md.append(f"No preparation step: SABER {NV} `default` searches the FASTA database directly (`-d`, index built in memory, "
              "included in its time); DIAMOND runs `makedb` and then `blastp --very-sensitive` (one timed command). "
              "Alternating runs.\n")
    simple_table("", list(OO.items()))
    a, b = ok(OO.get(f"{NEW}_default_fasta", [])), ok(OO.get("diamond_very-sensitive_makedb_blastp", []))
    if a and b:
        md.append("Alternating runs, SABER / DIAMOND wall time: "
                  + ", ".join(f"{x['wall_s']:.1f} / {y['wall_s']:.1f} s" for x, y in zip(a, b)) + ".\n")
    oc = os.path.join(R, "oneoff", "outputs_check.txt")
    if os.path.exists(oc):
        md.append("Output check: " + open(oc).read().strip() + ".\n")
    summary["oneoff"] = {n: dict(wall_s=med(v, "wall_s"), wall_range=rng(v, "wall_s"), cpu_s=med(v, "cpu_s"),
                                 rss_mb=(med(v, "maxrss_kb") or 0) / 1024, runs=len(ok(v))) for n, v in OO.items()}

# ---------------------------------------------------------------- SCOP40 all-vs-all
T = times("scop")
rows = []
if T:
    rows = tool_rows("scop", T, [f"{NEW}_{p}" for p in PRESETS + ("legacy", "exhaustive")] + [f"diamond_{m}" for m in DMODES]
                     + [f"{OLD}_default", f"{OLD}_more-sensitive", f"{NEW}_default_db"])
    summary["scop"] = rows
    md.append("## SCOP40 all-against-all (11,206 domains, 2 threads, E ≤ 10, no hit limit)\n")
    md.append(host_line("scop", T))
    md.append("SABER builds its index in memory on every run (`-d`; `default_db` = prebuilt `--db`), DIAMOND uses a prebuilt "
              "`.dmnd`. " + ACC_NOTE + "\n")
    acc_table(rows, T)
    md.append("### False-discovery proportion at matched recall (test split)\n")
    md.append("Pooled over all test-split queries: hits sorted by E-value, FDP = FP / (TP + FP) at the most stringent cut-off that reaches the recall. 'not reached' = the tool never reaches that recall at E ≤ 10.\n")
    fdp_table(rows, ("0.06", "0.08", "0.09", "0.1", "0.11"))
    simple_table("### Practical settings (E ≤ 1e-3, 25 targets per query), SCOP40 all-vs-all", [(n, v) for n, v in T.items() if n.startswith("practical_")])
    simple_table("### SCOP40 index build and persistent-index search", [(n, v) for n, v in T.items() if "makedb" in n or n.endswith("_db_warm")])

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
    md.append("Real queries against composition-preserving decoys (`shuf`: whole-sequence shuffle; `wshuf10`: shuffle inside "
              "10-residue windows). Hits to decoys derived from a sequence of the query's own fold are excluded; every remaining "
              "hit is a false positive. A calibrated E-value gives at most about E false hits per query (EPQ ≤ E). "
              f"`dbmodel` = SABER {NV} with `--evalue-model database` (database-wide E-values, as in 1.2).\n")
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


STRONG = 1e-10
M = times("mmseqs")
if M:
    ref = pairs(os.path.join(R, "mmseqs", f"{NEW}_exhaustive.m8"))
    ref_strong = pairs(os.path.join(R, "mmseqs", f"{NEW}_exhaustive.m8"), STRONG)
    mm = []
    for n in [f"{NEW}_{p}" for p in PRESETS + ("exhaustive",)] + [f"diamond_{m}" for m in DMODES] + [f"{OLD}_default"]:
        if n not in M:
            continue
        p = pairs(os.path.join(R, "mmseqs", n + ".m8"))
        mm.append(dict(name=n, wall_s=med(M[n], "wall_s"), rss_mb=(med(M[n], "maxrss_kb") or 0) / 1024,
                       pairs=len(p) if p is not None else None,
                       recall_vs_exhaustive=(len(p & ref) / len(ref)) if (p is not None and ref) else None,
                       recall_strong_pairs=(len(p & ref_strong) / len(ref_strong)) if (p is not None and ref_strong) else None,
                       extra_pairs=len(p - ref) if (p is not None and ref) else None))
    summary["mmseqs"] = dict(reference_pairs=len(ref) if ref else None, strong_reference_pairs=len(ref_strong) if ref_strong else None, runs=mm)
    md.append("## MMseqs2 example workload (500 queries × 20,000 UniProt sequences, E ≤ 1e-3)\n")
    md.append(host_line("mmseqs", M))
    md.append(f"Recall is measured against SABER {NV}'s exhaustive Smith–Waterman search (`--preset exhaustive`: {len(ref) if ref else '–'} "
              "query–subject pairs at E ≤ 1e-3). That reference shares SABER's scoring, composition correction and E-value "
              f"statistics, so near the cut-off it favours SABER; 'Recall, E ≤ 1e-10 pairs' uses only the {len(ref_strong) if ref_strong else '–'} "
              "reference pairs with E ≤ 1e-10, which a tool that finds them reports at E ≤ 1e-3 whatever its statistics. "
              "'Extra' = pairs a tool reports that the reference does not. SABER builds its index in memory on every run; "
              "DIAMOND uses a prebuilt `.dmnd`. 74 of the 500 queries are also database members.\n")
    md.append("| Tool / preset | Time s | RSS MB | Pairs | Recall of exhaustive SW | Recall, E ≤ 1e-10 pairs | Extra pairs |")
    md.append("|---|---|---|---|---|---|---|")
    for r in mm:
        md.append(f"| {r['name']} | {fmt_t(M[r['name']])} | {r['rss_mb']:.0f} | {r['pairs']} | {f3(r['recall_vs_exhaustive'], 4)} | {f3(r['recall_strong_pairs'], 4)} | {r['extra_pairs']} |")
    md.append("")
    simple_table("### Database and query-count scaling (single runs)", [(n, v) for n, v in M.items() if n.startswith("scale")])
    lg = [(n, v) for n, v in M.items() if n.startswith("long_")]
    if lg:
        lq = os.path.join(R, "mmseqs", "long_queries.fa")
        nlq = sum(1 for l in open(lq) if l.startswith(">")) if os.path.exists(lq) else "?"
        md.append(f"### Long proteins ({nlq} queries > 1,500 aa against the 20,000-sequence database)\n")
        md.append("Every fifth database sequence longer than 1,500 residues, so each also finds itself. Single runs.\n")
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
    md.append("| Run | Time s | CPU s | RSS MB |")
    md.append("|---|---|---|---|")
    for n, v in I.items():
        md.append(f"| {n} | {fmt_t(v)} | {f3(med(v, 'cpu_s'), 1)} | {fmt_mb(v)} |")
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
    md.append("blastx (`nx`): back-translated SCOP40 domains (random synonymous codons, random flanks and strand) against SCOP40 "
              "proteins; tblastn (`pn`): a fixed 1-in-37 subset of SCOP40 proteins (303 queries) against the back-translated set "
              "(DIAMOND has no tblastn mode). Metrics over the queries of each run (both splits). Single runs.\n")
    md.append("| Run | Queries | Time s | RSS MB | Sens. to 1st FP (SF) | Recall @ EPQ 0.01 | FP at E ≤ 1e-3 |")
    md.append("|---|---|---|---|---|---|---|")
    trs = {}
    for n in sorted(TR, key=lambda n: (n.split("_")[0] != "nx", n)):
        v = TR[n]
        m = metrics("translated", n, "all")
        trs[n] = dict(wall_s=med(v, "wall_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024, n_queries=m["n_queries"] if m else None,
                      sens_sf=m["sens_to_1st_fp_superfamily"] if m else None, recall_epq001=m["recall_at_epq"]["0.01"] if m else None,
                      fp_e1e3=m["at_evalue"][2]["fp"] if m else None)
        md.append(f"| {n} | {trs[n]['n_queries'] or '–'} | {fmt_t(v)} | {fmt_mb(v)} | {f3(trs[n]['sens_sf'])} | {f3(trs[n]['recall_epq001'], 4)} | {trs[n]['fp_e1e3'] if trs[n]['fp_e1e3'] is not None else '–'} |")
    md.append("")
    summary["translated"] = trs

# ---------------------------------------------------------------- low complexity
L = times("lowc")
if L:
    md.append("## Adversarial low-complexity / repeat inputs\n")
    md.append(host_line("lowc", L))
    md.append("100 homopolymer / tandem-repeat queries (some with an embedded real domain) against the MMseqs database plus 2,000 "
              "low-complexity sequences. `lowc10_*` / `lowc4_*`: every 10th query (10 queries) / the first 4 of them, with masking, "
              "composition correction and the frequency cap all off (every repeat hit enumerated; 20-minute time limit, exit 124 = "
              "stopped by it). `mmq_*`: the 500 real "
              "MMseqs queries against the same polluted database — hits to the artificial `lowcdb*` sequences are false positives.\n")
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
        hits_s = f"{tot} | {art}" if ex == 0 else "– | –"
        md.append(f"| {n} | {lc[n]['wall_s']:.1f} | {lc[n]['rss_mb']:.0f} | {ex} | {hits_s} |")
    md.append("")
    summary["lowc"] = lc
    for n in L:
        p = os.path.join(R, "lowc", n + ".stats.json")
        if os.path.exists(p):
            c = json.load(open(p))["counters"]
            summary.setdefault("lowc_counters", {})[n] = {k: c.get(k) for k in ("hits", "hits_filtered", "scan_splits", "scan_pairs_streamed", "subject_range_passes", "streamed_subjects", "seeds_frequent_skipped", "query_residues_masked")}

# ---------------------------------------------------------------- multi / fragments
MU = times("multi")
if MU:
    md.append("## Two-domain queries and short fragments\n")
    md.append(host_line("multi", MU))
    md.append("Retention = fraction of a domain's true positives from its single-domain search (each tool's own, E ≤ 1e-3) that are "
              "still found when the domain is part of a two-domain query. Fragments: 40-residue pieces of SCOP40 domains. Single runs.\n")
    md.append("| Tool / preset | Retention of the weaker domain's TPs | Retention (stronger) | Weaker-domain TPs (two-domain / alone) | Stronger-domain TPs (two-domain / alone) | FP in two-domain searches (E ≤ 1e-3) | Fragments: sens. to 1st FP (SF, test) | Fragment search time s |")
    md.append("|---|---|---|---|---|---|---|---|")
    mu = {}
    for t in (f"{NEW}_default", f"{NEW}_more-sensitive", "diamond_very-sensitive", "diamond_ultra-sensitive"):
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
    md.append("`onequery`: 20 long proteins searched one query per batch (`--query-batch 1`) against the ×4 database (index built in memory); "
              "every batch streams the whole database, in chunks shared by the threads (few-query searches: `docs/LIMITATIONS.md`).\n")
    md.append("| Run | Time s | CPU s | RSS MB |")
    md.append("|---|---|---|---|")
    for n, v in TH.items():
        md.append(f"| {n} | {fmt_t(v)} | {f3(med(v, 'cpu_s'), 1)} | {fmt_mb(v)} |")
    md.append("")
    summary["threads"] = {n: dict(wall_s=med(v, "wall_s"), cpu_s=med(v, "cpu_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024) for n, v in TH.items()}

# ---------------------------------------------------------------- sweep (tune splits)
SW = times("sweep")
if SW:
    md.append("## Preset selection (tune splits only)\n")
    md.append(host_line("sweep", SW))
    md.append("Configurations considered for the 1.3 presets (unchanged in 1.4), on the tune splits of the at-scale benchmark (`sc`) and of SCOP40 "
              "all-against-all (`av`). Single runs, 2 threads; the presets' held-out results are in the sections above. "
              "Names: seed set (DS4/DS8/DS16 = the first 4/8/16 of DIAMOND's very-sensitive weight-7 seeds over Murphy-10; "
              "W8 = eight weight-6 seeds; DS16W8 = both), window filter threshold T, region threshold fraction f, two-stage "
              "scoring prefilter band pf.\n")
    md.append("| Configuration | Scale: time s | Scale: sens. (tune) | Scale: recall @ EPQ 0.01 (tune) | All-vs-all: time s | All-vs-all: sens. (tune) | All-vs-all: recall @ EPQ 0.01 (tune) |")
    md.append("|---|---|---|---|---|---|---|")
    sw = {}
    names = sorted({n[3:] for n in SW if n[:3] in ("sc_", "av_")})
    for n in names:
        a, b = metrics("sweep", f"sc_{n}", "tune"), metrics("sweep", f"av_{n}", "tune")
        sw[n] = dict(scale_wall_s=med(SW.get(f"sc_{n}", []), "wall_s"), scale=acc(a), allvsall_wall_s=med(SW.get(f"av_{n}", []), "wall_s"), allvsall=acc(b))
    for n, r in sorted(sw.items(), key=lambda kv: kv[1]["scale_wall_s"] or 0):
        md.append(f"| {n} | {f3(r['scale_wall_s'], 1)} | {f3(r['scale'].get('sens_sf'), 4)} | {f3(r['scale'].get('recall_epq001'), 4)} | "
                  f"{f3(r['allvsall_wall_s'], 1)} | {f3(r['allvsall'].get('sens_sf'), 4)} | {f3(r['allvsall'].get('recall_epq001'), 4)} |")
    md.append("")
    summary["sweep"] = sw

# ---------------------------------------------------------------- scaling with database size
SC = times("scaling")
scaling = {}
if SC:
    sizes = {}
    dbs = os.path.join(R, "scaling", "dbsizes.tsv")
    if os.path.exists(dbs):
        for line in open(dbs):
            k, n, r, d = line.split()
            sizes[int(k)] = dict(sequences=int(n), residues=int(r), decoys=int(d))
    ks = sorted({int(n.split("_")[0][1:]) for n in SC})
    SCONF = [f"{NEW}_very-fast", f"{NEW}_default", f"{NEW}_more-sensitive", f"{OLD}_default",
             "diamond_default", "diamond_very-sensitive", "diamond_ultra-sensitive"]
    md.append("## Scaling with database size\n")
    md.append(host_line("scaling", SC))
    md.append("SCOP40 (11,206 domains) against SCOP40 + K × 20,000 reversed proteins (`make_scaling_dbs.py`; K = 4 is the at-scale "
              "database above). The decoy sources are the MMseqs2 example database and copies of it at 15 % substitution, so each "
              f"original protein gives K similar decoys. Prebuilt indexes; SABER {NV} `default` and DIAMOND `--very-sensitive` run "
              "alternately (median of the runs, range in brackets), the other rows single runs. Accuracy on the held-out test split "
              "(tune split in brackets); every hit to a decoy or another fold is a false positive.\n")
    md.append("| K | Database | Tool / preset | Time s | CPU s | RSS MB | Sens. to 1st FP | Recall @ EPQ 0.01 | TP / FP at E ≤ 1e-3 |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    for k in ks:
        z = sizes.get(k, {})
        db = f"{z.get('sequences', 0):,} seq., {z.get('residues', 0) / 1e6:.1f} M res." if z else "–"
        scaling[k] = dict(size=z, runs={})
        for c in SCONF:
            v = SC.get(f"k{k}_{c}", [])
            if not v:
                continue
            # a failed (e.g. timed-out) run's output is incomplete: no accuracy
            te, tu = (acc(metrics("scaling", f"k{k}.{c}", "test")), acc(metrics("scaling", f"k{k}.{c}", "tune"))) if ok(v) else ({}, {})
            scaling[k]["runs"][c] = dict(wall_s=med(v, "wall_s"), wall_range=rng(v, "wall_s"), cpu_s=med(v, "cpu_s"),
                                         rss_mb=(med(v, "maxrss_kb") or 0) / 1024, runs=len(ok(v)), exit=v[0].get("exit"), test=te, tune=tu)
            md.append(f"| {k} | {db} | {c} | {fmt_t(v)} | {f3(med(v, 'cpu_s'), 1)} | {fmt_mb(v)} | "
                      f"{f3(te.get('sens_sf'), 4)} [{f3(tu.get('sens_sf'), 4)}] | {f3(te.get('recall_epq001'), 4)} [{f3(tu.get('recall_epq001'), 4)}] | "
                      f"{te.get('tp_e1e3', '–')} / {te.get('fp_e1e3', '–')} |")
            db = ""
    md.append("")
    md.append("### Index construction and the one-off search (no preparation step) by database size\n")
    md.append(f"One-off: SABER {NV} `default` searching the FASTA (index built in memory, included) and DIAMOND `makedb` + `blastp "
              "--very-sensitive` as one timed command; single runs. The prebuilt-index times are the medians above.\n")
    md.append("| K | SABER makedb s / MB | DIAMOND makedb s / MB | SABER one-off s / MB | DIAMOND makedb + blastp s / MB | SABER prebuilt s | DIAMOND prebuilt s |")
    md.append("|---|---|---|---|---|---|---|")
    for k in ks:
        g = lambda n: SC.get(f"k{k}_{n}", [])
        cells = []
        for n in (f"makedb_saber_v{NV}", "makedb_diamond", f"{NEW}_default_fasta", "diamond_very-sensitive_makedb_blastp"):
            v = g(n)
            cells.append(f"{fmt_t(v)} / {fmt_mb(v)}" if v else "–")
        md.append(f"| {k} | " + " | ".join(cells) + f" | {f3(med(g(NEW + '_default'), 'wall_s'), 1)} | {f3(med(g('diamond_very-sensitive'), 'wall_s'), 1)} |")
        for n in (f"makedb_saber_v{NV}", "makedb_diamond", f"{NEW}_default_fasta", "diamond_very-sensitive_makedb_blastp"):
            v = g(n)
            if v:
                scaling[k]["runs"][n] = dict(wall_s=med(v, "wall_s"), cpu_s=med(v, "cpu_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024, runs=len(ok(v)))
    md.append("")
    oc = os.path.join(R, "scaling", "outputs_check.txt")
    if os.path.exists(oc):
        md.append("Output checks: " + "; ".join(l.strip() for l in open(oc) if l.strip()) + ".\n")
    summary["scaling"] = scaling

# ---------------------------------------------------------------- length-matched decoys
LM = times("lenmatch")
if LM:
    lm_rows = tool_rows("lenmatch", LM, [f"{NEW}_default", f"{NEW}_default_dbmodel", f"{NEW}_more-sensitive",
                                         "diamond_very-sensitive", "diamond_ultra-sensitive"])
    summary["lenmatch"] = lm_rows
    dsz = os.path.join(R, "lenmatch", "datasets.txt")
    md.append("## Length-matched decoys\n")
    md.append(host_line("lenmatch", LM))
    md.append("The at-scale database's decoys are whole proteins, longer than the SCOP40 domains, and per-subject E-values make a "
              "hit in a long subject less significant — which favours SABER there. Here the same 80,000 reversed proteins are cut into "
              "pieces with the lengths of SCOP40 domains (`make_lenmatched_db.py`), so true and false subjects have the same length "
              f"distribution. `{NEW}_default_dbmodel` = database-wide E-values. " + ACC_NOTE + "\n")
    if os.path.exists(dsz):
        md.append("```\n" + open(dsz).read().strip() + "\n```\n")
    acc_table(lm_rows, LM, sec="lenmatch")
    fdp_table(lm_rows, ("0.06", "0.07", "0.08", "0.09", "0.1"))

# ---------------------------------------------------------------- scaling with the number of queries
QS = times("qscale")
qscale = {}
if QS:
    md.append("## Scaling with the number of queries\n")
    md.append(host_line("qscale", QS))
    md.append(f"1 to 11,206 SCOP40 domains (every N-th) against the at-scale database (38.2 M residues), prebuilt indexes; median of "
              f"the runs. SABER {NV} `default`, SABER {OV} `default`, DIAMOND `--very-sensitive`.\n")
    md.append("| Queries | SABER " + NV + " s | SABER " + OV + " s | DIAMOND s | SABER " + NV + " MB | DIAMOND MB |")
    md.append("|---|---|---|---|---|---|")
    for n in (1, 10, 100, 1000, 11206):
        a, b, c = QS.get(f"q{n}_{NEW}_default", []), QS.get(f"q{n}_{OLD}_default", []), QS.get(f"q{n}_diamond_very-sensitive", [])
        if not (a or b or c):
            continue
        qscale[n] = {k: dict(wall_s=med(v, "wall_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024) for k, v in (("new", a), ("old", b), ("diamond", c)) if v}
        md.append(f"| {n:,} | {fmt_t(a)} | {fmt_t(b)} | {fmt_t(c)} | {fmt_mb(a)} | {fmt_mb(c)} |")
    md.append("")
    summary["qscale"] = qscale

# ---------------------------------------------------------------- cross-machine check
CM = times("checks_machine")
if CM:
    md.append("## Cross-machine check: the headline pairs on a second machine\n")
    md.append(host_line("checks_machine", CM))
    md.append("`cross_machine_check.sh`: SABER " + NV + " `default` and DIAMOND `--very-sensitive` alternately on the machine that ran "
              "the later sections, at scale (prebuilt indexes of the scale section; 5 runs each) and on the K = 16 scaling database "
              "(indexes built again; 3 runs each). Medians; the first machine's numbers are those of the scale and scaling sections.\n")
    md.append("| Database | Machine | SABER s | DIAMOND s | ratio | SABER MB | DIAMOND MB |")
    md.append("|---|---|---|---|---|---|---|")
    cm = {}
    def _row(db, mach, a, b):
        if not (a and b and med(a, "wall_s") and med(b, "wall_s")):
            return
        cm[f"{db} {mach}"] = dict(saber_s=med(a, "wall_s"), diamond_s=med(b, "wall_s"), saber_mb=med(a, "maxrss_kb") / 1024, diamond_mb=med(b, "maxrss_kb") / 1024)
        md.append(f"| {db} | {mach} | {fmt_t(a)} | {fmt_t(b)} | {med(a, 'wall_s') / med(b, 'wall_s'):.2f} | {fmt_mb(a)} | {fmt_mb(b)} |")
    mB = machines(next(iter(CM.values())))
    pick = lambda runs, m: [r for r in runs if MACHINE.get(r.get("cpu", "unknown")) == m]
    for m in sorted(set(MACHINE.values())):
        if m == mB:
            continue
        _row("at scale (38.2 M residues)", m, pick(ST.get(f"{NEW}_default", []), m), pick(ST.get("diamond_very-sensitive", []), m))
    _row("at scale (38.2 M residues)", mB, CM.get(f"scale_{NEW}_default", []), CM.get("scale_diamond_very-sensitive", []))
    for m in sorted(set(MACHINE.values())):
        if m == mB:
            continue
        _row("K = 16 (146.8 M residues)", m, pick(SC.get(f"k16_{NEW}_default", []), m), pick(SC.get("k16_diamond_very-sensitive", []), m))
    _row("K = 16 (146.8 M residues)", mB, CM.get(f"k16_{NEW}_default", []), CM.get("k16_diamond_very-sensitive", []))
    md.append("")
    oc = os.path.join(R, "checks_machine", "outputs_check.txt")
    if os.path.exists(oc):
        md.append("Outputs: " + "; ".join(l.strip() for l in open(oc) if l.strip()) + ".\n")
    summary["cross_machine"] = cm

# ---------------------------------------------------------------- one-off check (SABER 1.4, 1.3.0, DIAMOND)
OC = times("checks_oneoff")
if OC:
    md.append("## One-off search at scale: SABER " + NV + ", SABER " + OV + " and DIAMOND, alternating\n")
    md.append(host_line("checks_oneoff", OC).replace("Measured on machine ", "Measured on a new virtual machine of type "))
    md.append("`oneoff_check.sh`: no preparation step, at scale (38.2 M residues). SABER " + NV + " and SABER " + OV + " `default` "
              "search the FASTA database (index built in memory, included in the time); DIAMOND runs `makedb` and then "
              "`blastp --very-sensitive` as one timed command. The three were run alternately. Machines are named by CPU type: "
              "this check ran on a new virtual machine (the second replacement of the session), so a machine of a type seen "
              "earlier is not necessarily the same host.\n")
    md.append("| Run | Time s | CPU s | RSS MB |")
    md.append("|---|---|---|---|")
    oo = {}
    for n in (f"{NEW}_default_fasta", f"{OLD}_default_fasta", "diamond_very-sensitive_makedb_blastp"):
        v = OC.get(n, [])
        if not v:
            continue
        oo[n] = dict(wall_s=med(v, "wall_s"), cpu_s=med(v, "cpu_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024, runs=len(ok(v)))
        md.append(f"| {n} | {fmt_t(v)} | {f3(med(v, 'cpu_s'), 1)} | {fmt_mb(v)} |")
    md.append("")
    a, b, c = (OC.get(n, []) for n in (f"{NEW}_default_fasta", f"{OLD}_default_fasta", "diamond_very-sensitive_makedb_blastp"))
    if a and c and len(a) == len(c):
        pairs = ", ".join(f"{x['wall_s']:.1f} / {y['wall_s']:.1f} / {z['wall_s']:.1f} s" if b and len(b) == len(a) else f"{x['wall_s']:.1f} / {z['wall_s']:.1f} s"
                          for x, y, z in zip(a, b if b and len(b) == len(a) else a, c))
        md.append(("Alternating runs, SABER " + NV + " / SABER " + OV + " / DIAMOND" if b and len(b) == len(a) else "Alternating runs, SABER / DIAMOND")
                  + f" wall time: {pairs}.\n")
    oc = os.path.join(R, "checks_oneoff", "outputs_check.txt")
    if os.path.exists(oc):
        md.append("Outputs: " + "; ".join(l.strip() for l in open(oc) if l.strip()) + ".\n")
    summary["oneoff_check"] = oo

# ---------------------------------------------------------------- extra checks
XC = times("checks_extra")
if XC:
    def _counters(path):
        if not os.path.exists(path):
            return {}
        d = json.load(open(path))
        return d.get("counters", d)

    xc = {"counters": {}, "cap": {}}
    md.append("## Additional checks: SABER " + OV + "'s counters, the frequency cap, " + OV + " vs " + NV + " outputs\n")
    md.append(host_line("checks_extra", XC).replace("Measured on machine ", "Measured on a new virtual machine of type "))
    md.append("`extra_checks_14.sh`. Counters from `--stats` do not depend on the machine: candidate regions (`regions`) and "
              "query seeds skipped by the frequency cap (`seeds_frequent_skipped`) of SABER " + OV + " (run here on indexes built "
              "by " + NV + ") and of SABER " + NV + " (the scaling section's runs).\n")
    md.append("| Database | SABER " + OV + " regions | SABER " + NV + " regions | ratio | SABER " + OV + " seeds skipped | SABER " + NV + " seeds skipped |")
    md.append("|---|---|---|---|---|---|")
    for k in ks if SC else (8, 16):
        a = _counters(os.path.join(R, "checks_extra", f"k{k}.{OLD}_default.stats.json"))
        b = _counters(os.path.join(R, "scaling", f"k{k}.{NEW}_default.stats.json"))
        if not (a and b):
            continue
        xc["counters"][k] = dict(old_regions=a.get("regions"), new_regions=b.get("regions"),
                                 old_skipped=a.get("seeds_frequent_skipped"), new_skipped=b.get("seeds_frequent_skipped"))
        size = f"K = {k} ({scaling[k]['size']['residues'] / 1e6:.1f} M residues)" if k in scaling and "size" in scaling[k] and "residues" in scaling[k]["size"] else f"K = {k}"
        md.append(f"| {size} | {a.get('regions', 0):,} | {b.get('regions', 0):,} | {b.get('regions', 0) / max(1, a.get('regions', 0)):.2f} | "
                  f"{a.get('seeds_frequent_skipped', 0):,} | {b.get('seeds_frequent_skipped', 0):,} |")
    md.append("")
    md.append("The frequency cap at scale (38.2 M residues): SABER " + NV + " `default` with the automatic cap and without it "
              "(`--seed-max-occ 1000000000`), alternately; accuracy on both splits (the capped run's accuracy is the scale section's, "
              "its output is identical).\n")
    md.append("| Run | Time s | CPU s | RSS MB | Sens. test [tune] | Recall @ EPQ 0.01 test [tune] |")
    md.append("|---|---|---|---|---|---|")
    for n, msec, mname in ((f"scale_{NEW}_default", "scale", f"{NEW}_default"),
                           (f"scale_{NEW}_default_nocap", "checks_extra", f"scale_{NEW}_default_nocap")):
        v = XC.get(n, [])
        if not v:
            continue
        te, tu = acc(metrics(msec, mname, "test")), acc(metrics(msec, mname, "tune"))
        xc["cap"][n] = dict(wall_s=med(v, "wall_s"), cpu_s=med(v, "cpu_s"), rss_mb=(med(v, "maxrss_kb") or 0) / 1024, runs=len(ok(v)),
                            test=te, tune=tu)
        sens = f"{te['sens_sf']:.4f} [{tu['sens_sf']:.4f}]" if te and tu else "–"
        rec = f"{te['recall_epq001']:.4f} [{tu['recall_epq001']:.4f}]" if te and tu else "–"
        md.append(f"| {n} | {fmt_t(v)} | {f3(med(v, 'cpu_s'), 1)} | {fmt_mb(v)} | {sens} | {rec} |")
    md.append("")
    oc = os.path.join(R, "checks_extra", "outputs_check.txt")
    if os.path.exists(oc):
        md.append("Outputs: " + "; ".join(l.strip() for l in open(oc) if l.strip()) + ".\n")
    summary["extra_checks"] = xc

# ---------------------------------------------------------------- plots
# Categorical slots (validated set, light surface): blue = SABER 1.4,
# orange = DIAMOND, grey = SABER 1.3.0. Text is neutral ink; points are also
# labelled directly, so identity never depends on colour alone.
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, NullFormatter

    BLUE, ORANGE, AQUA, YELLOW, GREY, INK, INK2 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a8984", "#0b0b0b", "#52514e"
    plt.rcParams.update({"axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                         "axes.titlecolor": INK, "grid.color": "#d9d8d4", "grid.linewidth": 0.6, "font.size": 9})

    def family(n):
        return BLUE if n.startswith(NEW) else (ORANGE if n.startswith("diamond") else GREY)

    def short(n):
        return n.replace(NEW + "_", "").replace(OLD + "_", OV[:3] + " ").replace("diamond_", "D ")

    def log_ticks(ax, lo, hi):
        ticks = [t for t in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000) if lo <= t <= hi]
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(NullFormatter())

    # label offsets (points) for neighbours that would collide
    OFF = {"diamond_sensitive": (6, -10), "diamond_more-sensitive": (6, 3)}

    def tradeoff(rows_, fname, xlabel, title):
        pts = [r for r in rows_ if r.get("wall_s") and r["test"] and "exhaustive" not in r["name"]
               and not r["name"].endswith(("_db", "_c1", "_dbmodel"))]
        if not pts:
            return
        fig, ax = plt.subplots(figsize=(7.5, 5))
        for r in pts:
            x, y = r["wall_s"], r["test"]["sens_sf"]
            ax.scatter(x, y, color=family(r["name"]), s=52, zorder=3, edgecolors="white", linewidths=1.5)
            dx, dy = OFF.get(r["name"], (6, 2))
            ax.annotate(short(r["name"]), (x, y), fontsize=7.5, xytext=(dx, dy), textcoords="offset points", ha="left", color=INK)
        ax.set_xscale("log")
        xs = [r["wall_s"] for r in pts]
        ax.set_xlim(min(xs) / 1.3, max(xs) * 1.6)
        log_ticks(ax, min(xs) / 1.3, max(xs) * 1.6)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("sensitivity to the first FP, superfamily (test split)")
        ax.set_title(title)
        for lab, c in ((f"SABER {NV}", BLUE), ("DIAMOND 2.2.8", ORANGE), (f"SABER {OV}", GREY)):
            ax.scatter([], [], color=c, s=52, label=lab)
        ax.legend(fontsize=8, frameon=False, loc="lower right")
        ax.grid(True)
        fig.tight_layout()
        fig.savefig(os.path.join(R, fname), dpi=130)
        plt.close(fig)

    def roc(sec, fname, title, series):
        fig, ax = plt.subplots(figsize=(7.5, 5))
        drawn = 0
        for n, c, ls in series:
            p = os.path.join(R, sec, f"{n}.test.roc.tsv")
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
            drawn += 1
        if not drawn:
            plt.close(fig)
            return
        ax.set_xscale("log")
        ax.set_xlim(1e-3, max(l.get_xdata().max() for l in ax.get_lines()) * 1.5)
        ax.set_xlabel("false positives per query (EPQ, pooled over test-split queries)")
        ax.set_ylabel("recall (same-superfamily pairs)")
        ax.set_title(title)
        ax.legend(fontsize=7.5, frameon=False, loc="upper left")
        ax.grid(True)
        fig.tight_layout()
        fig.savefig(os.path.join(R, fname), dpi=130)
        plt.close(fig)

    tiers = ((f"{NEW}_fast", BLUE, "-"), ("diamond_mid-sensitive", BLUE, "--"),
             (f"{NEW}_default", ORANGE, "-"), ("diamond_very-sensitive", ORANGE, "--"),
             (f"{NEW}_more-sensitive", AQUA, "-"), ("diamond_ultra-sensitive", AQUA, "--"))
    if scale_rows:
        tradeoff(scale_rows, "scale_time_vs_sensitivity.png", "wall time, s (log scale; 11,206 queries vs 91,206 sequences, 2 threads)",
                 "At scale: runtime vs sensitivity")
        roc("scale", "scale_roc_test.png", f"Pooled ROC at scale, test split (solid: SABER {NV}, dashed: DIAMOND)", tiers)
    if rows:
        tradeoff(rows, "scop_time_vs_sensitivity.png", "wall time, s (log scale; SCOP40 all-vs-all, 2 threads)", "SCOP40 all-vs-all: runtime vs sensitivity")
        roc("scop", "scop_roc_test.png", f"Pooled ROC, SCOP40 test split (solid: SABER {NV}, dashed: DIAMOND)", tiers)
    if dec:
        cuts = [1e-3, 1e-2, 0.1, 1.0, 10.0]
        panels = (("SABER", ((f"{OLD}_default", GREY, ":", f"{OV} default"), (f"{NEW}_fast", BLUE, "-", f"{NV} fast"),
                             (f"{NEW}_default", ORANGE, "-", f"{NV} default"), (f"{NEW}_more-sensitive", AQUA, "-", f"{NV} more-sensitive"))),
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
    # ------------------------------------------------ scaling figures
    # Colour = tool everywhere in these figures (blue SABER 1.4, orange
    # DIAMOND, grey SABER 1.3.0), line style and marker repeat it, and the
    # lines are labelled directly in neutral ink.
    def end_labels(ax, items, gap=0.06):
        """Label line ends (x, y, text in data coordinates) to the right of
        the last point, spread vertically so the labels never overlap."""
        to_ax = ax.transData + ax.transAxes.inverted()
        pts = sorted(((tuple(to_ax.transform((x, y))), t) for x, y, t in items), key=lambda p: p[0][1])
        ys = [p[0][1] for p in pts]
        for i in range(1, len(ys)):
            ys[i] = max(ys[i], ys[i - 1] + gap)
        if ys and ys[-1] > 0.97:
            d = ys[-1] - 0.97
            ys = [y - d for y in ys]
            for i in range(len(ys) - 2, -1, -1):
                ys[i] = min(ys[i], ys[i + 1] - gap)
        for (xy, t), y in zip(pts, ys):
            ax.annotate(t, xy=xy, xycoords="axes fraction", xytext=(xy[0] + 0.025, y), textcoords="axes fraction",
                        fontsize=7.5, color=INK, va="center", ha="left")

    def size_axis(ax, xs, room=2.4):
        ax.set_xscale("log")
        ax.set_xlim(xs[0] / 1.25, xs[-1] * room)
        ax.xaxis.set_major_locator(FixedLocator(xs))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}"))
        ax.xaxis.set_minor_formatter(NullFormatter())

    def fix_limits(ax, ylog=False, zero=False):
        ax.autoscale_view()
        lo, hi = ax.get_ylim()
        if ylog:
            ax.set_ylim(lo / 1.15, hi * 1.15)
            ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else f"{v:g}"))
            ax.yaxis.set_minor_formatter(NullFormatter())
        else:
            ax.set_ylim(0 if zero else lo - 0.06 * (hi - lo), hi + 0.06 * (hi - lo))

    MAIN = ((f"{NEW}_default", BLUE, "-", "o", f"SABER {NV} default"),
            ("diamond_very-sensitive", ORANGE, "--", "s", "DIAMOND very-sensitive"),
            (f"{OLD}_default", GREY, ":", "^", f"SABER {OV} default"))
    if scaling:
        ks = sorted(k for k in scaling if scaling[k]["size"].get("residues"))
        xs = [scaling[k]["size"]["residues"] / 1e6 for k in ks]
        if ks:
            # 1. the matched comparison: time, memory and accuracy against database size
            fig, axes = plt.subplots(4, 2, figsize=(11, 15))
            panels = ((axes[0][0], "wall time, s (log)", "wall_s", True),
                      (axes[0][1], "peak resident memory, MB (log)", "rss_mb", True),
                      (axes[2][0], "sensitivity to the first false positive (test split)", "sens_sf", False),
                      (axes[2][1], "recall at 0.01 false positives per query (test split)", "recall_epq001", False),
                      (axes[3][0], "true positives at E ≤ 1e-3 (test split)", "tp_e1e3", False),
                      (axes[3][1], "false positives at E ≤ 1e-3 (test split)", "fp_e1e3", False))
            # relative to DIAMOND --very-sensitive at the same size (1 = parity)
            for ax, key, ylab in ((axes[1][0], "wall_s", "wall time relative to DIAMOND --very-sensitive"),
                                  (axes[1][1], "rss_mb", "peak memory relative to DIAMOND --very-sensitive")):
                ends = []
                rser = [(f"{NEW}_default", "diamond_very-sensitive", BLUE, "-", "o", True, f"SABER {NV} default, prebuilt indexes", f"SABER {NV}"),
                        (f"{OLD}_default", "diamond_very-sensitive", GREY, ":", "^", True, f"SABER {OV} default, prebuilt indexes", f"SABER {OV}")]
                if key == "wall_s":
                    rser.append((f"{NEW}_default_fasta", "diamond_very-sensitive_makedb_blastp", BLUE, "--", "o", False,
                                 f"SABER {NV} on the FASTA vs DIAMOND makedb + blastp", f"SABER {NV}, no preparation"))
                for a_, b_, c, ls, mk, fill, lab, short_ in rser:
                    pts = []
                    for x, k in zip(xs, ks):
                        ra, rb = scaling[k]["runs"].get(a_), scaling[k]["runs"].get(b_)
                        if ra and rb and ra.get(key) and rb.get(key):
                            pts.append((x, ra[key] / rb[key]))
                    if not pts:
                        continue
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=c, ls=ls, lw=2, marker=mk, ms=7,
                            markerfacecolor=c if fill else "white", markeredgecolor="white" if fill else c, markeredgewidth=1.3, label=lab, zorder=3)
                    ends.append((pts[-1][0], pts[-1][1], short_))
                ax.axhline(1.0, color=INK2, lw=1, ls=(0, (4, 3)), zorder=1)
                size_axis(ax, xs)
                fix_limits(ax)
                lo, hi = ax.get_ylim()
                ax.set_ylim(min(lo, 0.0 if key == "rss_mb" else 0.7), max(hi, 1.1))
                ax.annotate("DIAMOND = 1", xy=(0.02, 1.0), xycoords=("axes fraction", "data"), xytext=(0, 3), textcoords="offset points",
                            fontsize=7.5, color=INK2, va="bottom")
                end_labels(ax, ends)
                ax.set_xlabel("database residues, millions (log)")
                ax.set_ylabel(ylab)
                ax.legend(fontsize=7, frameon=False, loc="lower left" if key == "wall_s" else "upper right")
                ax.grid(True)
            for ax, ylab, key, logy in panels:
                ends = []
                for n, c, ls, mk, lab in MAIN:
                    pts = []
                    for x, k in zip(xs, ks):
                        r = scaling[k]["runs"].get(n)
                        if not r or not r.get("runs"):
                            continue
                        y = r.get(key) if key in ("wall_s", "rss_mb") else (r.get("test") or {}).get(key)
                        if y is not None:
                            pts.append((x, y, r))
                    if not pts:
                        continue
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=c, ls=ls, lw=2, marker=mk, ms=7,
                            markeredgecolor="white", markeredgewidth=1.2, label=lab, zorder=3)
                    if key == "wall_s":
                        for x, y, r in pts:
                            lo, hi = r.get("wall_range", (None, None))
                            if r.get("runs", 0) > 1 and lo is not None and hi > lo:
                                ax.plot([x, x], [lo, hi], color=c, lw=1, zorder=2)
                    ends.append((pts[-1][0], pts[-1][1], lab.replace(" default", "").replace(" very-sensitive", " v.-sens.")))
                size_axis(ax, xs)
                if logy:
                    ax.set_yscale("log")
                fix_limits(ax, ylog=logy, zero=key == "fp_e1e3")
                end_labels(ax, ends)
                ax.set_xlabel("database residues, millions (log)")
                ax.set_ylabel(ylab)
                ax.grid(True)
            h, l = axes[0][0].get_legend_handles_labels()
            fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8.5)
            fig.suptitle("SCOP40 (11,206 queries) against SCOP40 + K × 20,000 reversed proteins, K = " + ", ".join(map(str, ks))
                         + "\nprebuilt indexes, 2 threads; time: median and range of the alternating runs", color=INK, fontsize=10)
            fig.tight_layout(rect=(0, 0.03, 1, 0.975))
            fig.savefig(os.path.join(R, "scaling_database_size.png"), dpi=130)
            plt.close(fig)

            # 2. speed / sensitivity trade-off at every size (small multiples, shared axes)
            # colour = tool, marker = tier (fastest / matched / most sensitive)
            TIER_MK = ("v", "o", "D")
            TIERS = (("SABER", BLUE, "-", "o", ((f"{NEW}_very-fast", "very-fast"), (f"{NEW}_default", "default"), (f"{NEW}_more-sensitive", "more-sens."))),
                     ("DIAMOND", ORANGE, "--", "s", (("diamond_default", "default"), ("diamond_very-sensitive", "very-sens."), ("diamond_ultra-sensitive", "ultra"))),
                     (f"SABER {OV}", GREY, ":", "^", ((f"{OLD}_default", "1.3.0"),)))
            ncol = 3
            nrow = (len(ks) + 1 + ncol - 1) // ncol
            fig, axes = plt.subplots(nrow, ncol, figsize=(12, 4.1 * nrow), sharex=True, sharey=True, squeeze=False)
            flat = [a for row in axes for a in row]
            allx = [r["wall_s"] for k in ks for r in scaling[k]["runs"].values() if r.get("wall_s") and r.get("test")]
            for ax, k in zip(flat, ks):
                for tool, c, ls, mk, confs in TIERS:
                    pts = []
                    for ti, (n, short_) in enumerate(confs):
                        r = scaling[k]["runs"].get(n)
                        if r and r.get("wall_s") and r.get("test") and r.get("runs"):
                            pts.append((r["wall_s"], r["test"]["sens_sf"], TIER_MK[ti] if len(confs) > 1 else mk))
                    if not pts:
                        continue
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=c, ls=ls, lw=1.6, zorder=2)
                    for x, y, m in pts:
                        ax.plot([x], [y], color=c, marker=m, ms=8, markeredgecolor="white", markeredgewidth=1.2, ls="none", zorder=3)
                z = scaling[k]["size"]
                ax.set_title(f"K = {k}: {z['sequences']:,} sequences, {z['residues'] / 1e6:.1f} M residues", fontsize=9)
                ax.set_xscale("log")
                ax.grid(True)
            for ax in flat[len(ks):]:
                ax.axis("off")
            if allx:
                lo, hi = min(allx) / 1.4, max(allx) * 2.2
                flat[0].set_xlim(lo, hi)
                ticks = [t for t in (2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000) if lo <= t <= hi]
                for ax in flat[:len(ks)]:
                    ax.xaxis.set_major_locator(FixedLocator(ticks))
                    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
                    ax.xaxis.set_minor_formatter(NullFormatter())
            for i, ax in enumerate(flat[:len(ks)]):
                if i % ncol == 0:
                    ax.set_ylabel("sensitivity to the first FP (test split)")
                if i >= len(ks) - ncol:
                    ax.set_xlabel("wall time, s (log)")
                    ax.xaxis.set_tick_params(labelbottom=True)
            from matplotlib.lines import Line2D
            h = [Line2D([], [], color=BLUE, ls="-", lw=1.6, label=f"SABER {NV}"),
                 Line2D([], [], color=ORANGE, ls="--", lw=1.6, label="DIAMOND 2.2.8"),
                 Line2D([], [], color=GREY, ls="none", marker="^", ms=8, label=f"SABER {OV} default"),
                 Line2D([], [], color=INK2, ls="none", marker="v", ms=8, label="fastest: very-fast / DIAMOND default"),
                 Line2D([], [], color=INK2, ls="none", marker="o", ms=8, label="matched: default / --very-sensitive"),
                 Line2D([], [], color=INK2, ls="none", marker="D", ms=7, label="most sensitive: more-sensitive / --ultra-sensitive")]
            if len(ks) < len(flat):
                flat[len(ks)].legend(handles=h, loc="center", frameon=False, fontsize=8.5,
                                     title="2 threads, prebuilt indexes;\nup and to the left is better", title_fontsize=8.5)
            else:
                fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=8.5)
            fig.suptitle("Speed and sensitivity at every database size", color=INK)
            fig.tight_layout()
            fig.savefig(os.path.join(R, "scaling_tradeoff.png"), dpi=130)
            plt.close(fig)

            # 3. no preparation step: one-off searches and index construction
            fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
            ax = axes[0]
            ends = []
            for n, c, ls, mk, fill, lab, short_ in (
                    (f"{NEW}_default", BLUE, "-", "o", True, f"SABER {NV} default, prebuilt index", "SABER, prebuilt"),
                    (f"{NEW}_default_fasta", BLUE, "--", "o", False, f"SABER {NV} default on the FASTA (index built in memory)", "SABER, FASTA"),
                    ("diamond_very-sensitive", ORANGE, "-", "s", True, "DIAMOND very-sensitive, prebuilt .dmnd", "DIAMOND, prebuilt"),
                    ("diamond_very-sensitive_makedb_blastp", ORANGE, "--", "s", False, "DIAMOND makedb + blastp very-sensitive", "DIAMOND, makedb + blastp")):
                pts = [(x, scaling[k]["runs"][n]["wall_s"]) for x, k in zip(xs, ks) if scaling[k]["runs"].get(n, {}).get("wall_s")]
                if not pts:
                    continue
                ax.plot([p[0] for p in pts], [p[1] for p in pts], color=c, ls=ls, lw=2, marker=mk, ms=7,
                        markerfacecolor=c if fill else "white", markeredgecolor=c if not fill else "white", markeredgewidth=1.4, label=lab, zorder=3)
                ends.append((pts[-1][0], pts[-1][1], short_))
            size_axis(ax, xs, room=4.2)
            ax.set_yscale("log")
            fix_limits(ax, ylog=True)
            end_labels(ax, ends)
            ax.set_xlabel("database residues, millions (log)")
            ax.set_ylabel("wall time, s (log)")
            ax.set_title("Search, with and without index construction", fontsize=9.5)
            ax.legend(fontsize=7.5, frameon=False, loc="lower right")
            ax.grid(True)
            ax = axes[1]
            ends = []
            for n, c, ls, mk, lab in ((f"makedb_saber_v{NV}", BLUE, "-", "o", f"SABER {NV} --makedb"), ("makedb_diamond", ORANGE, "--", "s", "DIAMOND makedb")):
                pts = [(x, scaling[k]["runs"][n]["wall_s"]) for x, k in zip(xs, ks) if scaling[k]["runs"].get(n, {}).get("wall_s")]
                if not pts:
                    continue
                ax.plot([p[0] for p in pts], [p[1] for p in pts], color=c, ls=ls, lw=2, marker=mk, ms=7,
                        markeredgecolor="white", markeredgewidth=1.2, label=lab, zorder=3)
                ends.append((pts[-1][0], pts[-1][1], lab.split(" ")[0]))
            size_axis(ax, xs)
            ax.set_yscale("log")
            fix_limits(ax, ylog=True)
            end_labels(ax, ends)
            ax.set_xlabel("database residues, millions (log)")
            ax.set_ylabel("wall time, s (log)")
            ax.set_title("Index construction alone (single runs)", fontsize=9.5)
            ax.legend(fontsize=7.5, frameon=False, loc="lower right")
            ax.grid(True)
            fig.suptitle("No preparation step vs prebuilt indexes (SCOP40 queries, 2 threads)", color=INK)
            fig.tight_layout()
            fig.savefig(os.path.join(R, "scaling_no_preparation.png"), dpi=130)
            plt.close(fig)
    if qscale:
        ns = sorted(qscale)
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
        for ax, key, ylab, logy in ((axes[0], "wall_s", "wall time, s (log)", True), (axes[1], "rss_mb", "peak resident memory, MB (log)", True)):
            ends = []
            for tk, c, ls, mk, lab in (("new", BLUE, "-", "o", f"SABER {NV} default"), ("diamond", ORANGE, "--", "s", "DIAMOND very-sensitive"),
                                       ("old", GREY, ":", "^", f"SABER {OV} default")):
                pts = [(n, qscale[n][tk][key]) for n in ns if qscale[n].get(tk, {}).get(key)]
                if not pts:
                    continue
                ax.plot([p[0] for p in pts], [p[1] for p in pts], color=c, ls=ls, lw=2, marker=mk, ms=7,
                        markeredgecolor="white", markeredgewidth=1.2, label=lab, zorder=3)
                ends.append((pts[-1][0], pts[-1][1], lab.replace(" default", "").replace(" very-sensitive", "")))
            ax.set_xscale("log")
            ax.set_xlim(ns[0] / 1.6, ns[-1] * 4.5)
            ax.xaxis.set_major_locator(FixedLocator(ns))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
            ax.xaxis.set_minor_formatter(NullFormatter())
            if logy:
                ax.set_yscale("log")
            fix_limits(ax, ylog=logy)
            end_labels(ax, ends)
            ax.set_xlabel("number of query sequences (log)")
            ax.set_ylabel(ylab)
            ax.legend(fontsize=7.5, frameon=False, loc="upper left" if key == "wall_s" else "center left")
            ax.grid(True)
        fig.suptitle("Scaling with the number of queries (38.2 M-residue database, prebuilt indexes, 2 threads)", color=INK)
        fig.tight_layout()
        fig.savefig(os.path.join(R, "scaling_queries.png"), dpi=130)
        plt.close(fig)
except Exception as e:  # plotting is optional
    print("plotting skipped:", e, file=sys.stderr)

summary["hosts"] = HOSTS
json.dump(summary, open(os.path.join(R, "summary.json"), "w"), indent=1)
open(os.path.join(R, "summary.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md))
