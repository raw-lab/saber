#!/usr/bin/env python3
"""Reproduce the verified v1.0.2/v1.1.0 scoring/reporting defects on tiny
inputs (bench/defects/*.fa) and show the v1.2 result next to an independent
reference where one applies.

  python3 -I bench/scripts/defects.py V102 V110 V12 bench/defects > defects.json
"""
import json, subprocess, sys, os

CASES = {
    "amb":  "ambiguous residues (X/B/Z) in the subject",
    "amb2": "120-residue pair identical except for an X at every fifth subject position "
            "(SABER 1.x default search drops it in its ungapped filter; 1.1.0 --no-index finds it at E = 5e-73)",
    "gap":  "alignment with one 5-residue and one 3-residue gap",
    "long": "3,000-residue identical pair (v1.x traceback window cap = 2,048)",
    "sat":  "4,000-residue poly-W identical pair (exact score 44,000, above the 16-bit range). "
            "SABER 1.2's default search reports no hit by design: tantan masks the homopolymer for "
            "seeding and the composition-bias correction removes the score of a pure poly-W "
            "alignment; the extra run with both disabled shows the exact score",
}
# Additional runs: (version, label suffix, extra arguments).
EXTRA = {
    "amb":  [("v1.2.0", " --comp-bias off", ["--comp-bias", "off"])],
    "amb2": [("v1.1.0", " --no-index", ["--no-index"])],
    "gap":  [("v1.2.0", " --comp-bias off", ["--comp-bias", "off"])],
    "sat":  [("v1.2.0", " --masking none --comp-bias off", ["--masking", "none", "--comp-bias", "off"])],
}

def run(binary, q, d, extra=()):
    p = subprocess.run([binary, "-q", q, "-d", d, "--quiet", "-t", "1", "--outfmt", "json", "-E", "10", *extra],
                       capture_output=True, text=True)
    out = p.stdout
    # keep only the JSON array; any other stdout line is counted as noise
    i, j = out.find("["), out.rfind("]")
    try:
        hits = json.loads(out[i:j + 1]) if i >= 0 else []
    except json.JSONDecodeError:
        hits = []
    stray = [l for l in out.splitlines() if l.startswith("[") and not l.strip().startswith("[{") and l.strip() not in ("[", "]")]
    return hits, stray

def cigar_lengths(c):
    q = s = n = 0
    for ch in c:
        if ch.isdigit():
            n = n * 10 + int(ch)
        else:
            if ch in "M=X":
                q += n; s += n
            elif ch == "I":
                q += n
            elif ch == "D":
                s += n
            n = 0
    return q, s

def main():
    bins = dict(zip(["v1.0.2", "v1.1.0", "v1.2.0"], sys.argv[1:4]))
    d = sys.argv[4]
    res = {}
    for case, desc in CASES.items():
        q, db = os.path.join(d, f"{case}_q.fa"), os.path.join(d, f"{case}_db.fa")
        qlen = sum(len(l.strip()) for l in open(q) if not l.startswith(">"))
        slen = sum(len(l.strip()) for l in open(db) if not l.startswith(">"))
        res[case] = {"description": desc, "query_length": qlen, "subject_length": slen}
        runs = [(name, b, ()) for name, b in bins.items()]
        runs += [(ver + suffix, bins[ver], args) for ver, suffix, args in EXTRA.get(case, [])]
        for name, b, extra in runs:
            hits, stray = run(b, q, db, extra)
            if not hits:
                res[case][name] = {"hits": 0}
                continue
            h = hits[0]
            qr, sr = h.get("query_range"), h.get("subject_range")
            cq, cs = cigar_lengths(h.get("cigar", ""))
            res[case][name] = {
                "score": h.get("score"), "evalue": h.get("evalue"), "bitscore": round(h.get("bitscore", 0), 1),
                "align_length": h.get("align_length"), "query_range": qr, "subject_range": sr,
                "cigar": h.get("cigar") if len(h.get("cigar", "")) < 40 else h.get("cigar")[:37] + "...",
                "cigar_consumes": [cq, cs],
                "cigar_consistent": bool(qr and sr and cq == qr[1] - qr[0] + 1 and cs == sr[1] - sr[0] + 1),
                "gaps_field": h.get("gaps"), "gap_opens": h.get("gap_opens"),
                "query_cover": h.get("query_cover", h.get("query_coverage")),
                "stdout_noise_lines": len(stray),
            }
    json.dump(res, sys.stdout, indent=1)
    print()

if __name__ == "__main__":
    main()
