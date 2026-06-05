#!/usr/bin/env python3
"""Generate benchmark plots.

Reads /tmp/bench_results.tsv (output of bench_5way.sh) by default; falls
back to the bundled results_single_thread.tsv next to this script when
the temp file is absent. That way you can re-render plots from the
captured results without rerunning the bench.

Writes four PNGs into benchmarks/bio/plots/:
  - speed.png         (wall-clock time, log scale)
  - memory.png        (peak RSS)
  - disk.png          (DB index size)
  - sensitivity.png   (hits vs BLAST)

And a composite figure benchmarks/bio/plots/summary.png.
"""

import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
OUT = HERE / "plots"
OUT.mkdir(exist_ok=True)

# ---- Load data ----
# Prefer /tmp/bench_results.tsv (fresh bench output); fall back to bundled file.
tsv_paths = [Path("/tmp/bench_results.tsv"), HERE / "results_single_thread.tsv"]
tsv = next((p for p in tsv_paths if p.exists()), None)
if tsv is None:
    raise SystemExit(
        f"No bench results found. Looked in: {[str(p) for p in tsv_paths]}.\n"
        f"Run bench_5way.sh first, or ensure results_single_thread.tsv is present."
    )
print(f"Reading results from {tsv}")

rows = []
with open(tsv) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        rows.append({
            "tool": r["tool"],
            "mode": r["mode"],
            "label": f'{r["tool"]} ({r["mode"]})' if r["tool"] not in ("BLAST",) else r["tool"],
            "time_s": float(r["time_s"]),
            "peak_mem_mb": float(r["peak_mem_mb"]),
            "hits": int(r["hits"]),
            "disk_mb": float(r["disk_mb"]),
        })

# Friendly labels and an ordered presentation matching how people talk about these tools.
ORDER = [
    ("BLAST", "blastp"),
    ("SWORD", "sw"),
    ("DIAMOND", "default"),
    ("DIAMOND", "sensitive"),
    ("DIAMOND", "more-sens"),
    ("BLAT", "prot-prot"),
    ("SABER", "default"),
    ("SABER", "fast"),
]
LABELS = {
    ("BLAST", "blastp"):       "BLAST",
    ("SWORD", "sw"):           "SWORD",
    ("DIAMOND", "default"):    "DIAMOND\n(default)",
    ("DIAMOND", "sensitive"):  "DIAMOND\n--sensitive",
    ("DIAMOND", "more-sens"):  "DIAMOND\n--more-sensitive",
    ("BLAT", "prot-prot"):     "BLAT",
    ("SABER", "default"):      "SABER v1.0.0\n(default)",
    ("SABER", "fast"):         "SABER v1.0.0\n(fast)",
}
COLORS = {
    "BLAST":   "#3b82f6",   # blue — reference
    "SWORD":   "#a855f7",   # purple
    "DIAMOND": "#f59e0b",   # amber
    "BLAT":    "#10b981",   # green
    "SABER":   "#ef4444",   # red — us
}

ordered = []
for key in ORDER:
    for r in rows:
        if (r["tool"], r["mode"]) == key:
            ordered.append(r)
            break

labels = [LABELS[(r["tool"], r["mode"])] for r in ordered]
colors = [COLORS[r["tool"]] for r in ordered]

# Use a clean style with reasonable defaults.
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": "--",
    "figure.dpi": 120,
})

# ---------------------------------------------------------------------------
# 1. Speed (log scale because BLAT is 50× faster than BLAST)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 5.5))
times = [r["time_s"] for r in ordered]
bars = ax.bar(range(len(ordered)), times, color=colors, edgecolor="black", linewidth=0.6)
ax.set_yscale("log")
ax.set_ylabel("Wall-clock time (s, log scale)")
ax.set_xticks(range(len(ordered)))
ax.set_xticklabels(labels, rotation=0, fontsize=9.5)
ax.set_title("Real-biology benchmark: 100 UniProt queries × 20K UniProt subjects, E ≤ 1e-5, single-thread",
             fontsize=11, pad=12)
for bar, t in zip(bars, times):
    ax.text(bar.get_x() + bar.get_width()/2, t * 1.08, f"{t:.2f}s",
            ha="center", va="bottom", fontsize=10, fontweight="bold")
# horizontal reference at BLAST time
blast_t = next(r["time_s"] for r in ordered if r["tool"] == "BLAST")
ax.axhline(blast_t, color="#3b82f6", linestyle=":", alpha=0.5, linewidth=1)
ax.text(len(ordered)-0.4, blast_t * 1.15, "BLAST baseline", color="#3b82f6",
        fontsize=9, ha="right", style="italic")
plt.tight_layout()
plt.savefig(OUT / "speed.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Wrote {OUT / 'speed.png'}")

# ---------------------------------------------------------------------------
# 2. Peak memory
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 5.5))
mem = [r["peak_mem_mb"] for r in ordered]
bars = ax.bar(range(len(ordered)), mem, color=colors, edgecolor="black", linewidth=0.6)
ax.set_ylabel("Peak resident memory (MB)")
ax.set_xticks(range(len(ordered)))
ax.set_xticklabels(labels, rotation=0, fontsize=9.5)
ax.set_title("Peak memory usage — same workload as speed plot",
             fontsize=11, pad=12)
for bar, m in zip(bars, mem):
    ax.text(bar.get_x() + bar.get_width()/2, m + max(mem)*0.015, f"{m:.0f} MB",
            ha="center", va="bottom", fontsize=10, fontweight="bold")
plt.tight_layout()
plt.savefig(OUT / "memory.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Wrote {OUT / 'memory.png'}")

# ---------------------------------------------------------------------------
# 3. Disk (DB index size)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 5.5))
disk = [r["disk_mb"] for r in ordered]
bars = ax.bar(range(len(ordered)), disk, color=colors, edgecolor="black", linewidth=0.6)
ax.set_ylabel("Disk footprint of DB (MB)")
ax.set_xticks(range(len(ordered)))
ax.set_xticklabels(labels, rotation=0, fontsize=9.5)
ax.set_title("On-disk DB size — BLAST/DIAMOND pre-build; SWORD/BLAT/SABER work from raw FASTA",
             fontsize=11, pad=12)
for bar, d in zip(bars, disk):
    ax.text(bar.get_x() + bar.get_width()/2, d + max(disk)*0.015, f"{d:.0f} MB",
            ha="center", va="bottom", fontsize=10, fontweight="bold")
plt.tight_layout()
plt.savefig(OUT / "disk.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Wrote {OUT / 'disk.png'}")

# ---------------------------------------------------------------------------
# 4. Sensitivity (hits found vs BLAST baseline)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 5.5))
hits = [r["hits"] for r in ordered]
bars = ax.bar(range(len(ordered)), hits, color=colors, edgecolor="black", linewidth=0.6)
ax.set_ylabel("Unique (query, subject) hits at E ≤ 1e-5")
ax.set_xticks(range(len(ordered)))
ax.set_xticklabels(labels, rotation=0, fontsize=9.5)
ax.set_title("Recall: number of unique query–subject pairs reported (max_aligns=10)",
             fontsize=11, pad=12)
for bar, h in zip(bars, hits):
    ax.text(bar.get_x() + bar.get_width()/2, h + max(hits)*0.012, str(h),
            ha="center", va="bottom", fontsize=10, fontweight="bold")
# BLAST baseline reference
blast_h = next(r["hits"] for r in ordered if r["tool"] == "BLAST")
ax.axhline(blast_h, color="#3b82f6", linestyle=":", alpha=0.5, linewidth=1)
ax.text(len(ordered)-0.4, blast_h + max(hits)*0.03, "BLAST baseline",
        color="#3b82f6", fontsize=9, ha="right", style="italic")
plt.tight_layout()
plt.savefig(OUT / "sensitivity.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Wrote {OUT / 'sensitivity.png'}")

# ---------------------------------------------------------------------------
# 5. Composite summary — 2×2 panel
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(16, 11))

# Speed
ax = axes[0, 0]
ax.bar(range(len(ordered)), times, color=colors, edgecolor="black", linewidth=0.6)
ax.set_yscale("log")
ax.set_ylabel("Time (s, log)")
ax.set_xticks(range(len(ordered)))
ax.set_xticklabels(labels, rotation=0, fontsize=8.5)
ax.set_title("Wall-clock time")
for i, t in enumerate(times):
    ax.text(i, t * 1.08, f"{t:.1f}s", ha="center", va="bottom", fontsize=9, fontweight="bold")
ax.axhline(blast_t, color="#3b82f6", linestyle=":", alpha=0.5, linewidth=1)

# Memory
ax = axes[0, 1]
ax.bar(range(len(ordered)), mem, color=colors, edgecolor="black", linewidth=0.6)
ax.set_ylabel("Peak memory (MB)")
ax.set_xticks(range(len(ordered)))
ax.set_xticklabels(labels, rotation=0, fontsize=8.5)
ax.set_title("Peak resident memory")
for i, m in enumerate(mem):
    ax.text(i, m + max(mem)*0.015, f"{m:.0f}", ha="center", va="bottom",
            fontsize=9, fontweight="bold")

# Disk
ax = axes[1, 0]
ax.bar(range(len(ordered)), disk, color=colors, edgecolor="black", linewidth=0.6)
ax.set_ylabel("DB on-disk (MB)")
ax.set_xticks(range(len(ordered)))
ax.set_xticklabels(labels, rotation=0, fontsize=8.5)
ax.set_title("DB disk footprint")
for i, d in enumerate(disk):
    ax.text(i, d + max(disk)*0.015, f"{d:.0f}", ha="center", va="bottom",
            fontsize=9, fontweight="bold")

# Hits
ax = axes[1, 1]
ax.bar(range(len(ordered)), hits, color=colors, edgecolor="black", linewidth=0.6)
ax.set_ylabel("Hits (E ≤ 1e-5)")
ax.set_xticks(range(len(ordered)))
ax.set_xticklabels(labels, rotation=0, fontsize=8.5)
ax.set_title("Recall (unique pairs found)")
for i, h in enumerate(hits):
    ax.text(i, h + max(hits)*0.012, str(h), ha="center", va="bottom",
            fontsize=9, fontweight="bold")
ax.axhline(blast_h, color="#3b82f6", linestyle=":", alpha=0.5, linewidth=1)

fig.suptitle("Real-biology benchmark: 100 UniProt queries × 20K UniProt subjects, single-thread\n"
             "SABER v1.0.0 vs BLAST / SWORD / DIAMOND",
             fontsize=13, fontweight="bold", y=0.995)
plt.tight_layout()
plt.savefig(OUT / "summary.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Wrote {OUT / 'summary.png'}")

print("\nAll plots done.")
