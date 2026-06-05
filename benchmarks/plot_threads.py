#!/usr/bin/env python3
"""Plot multi-thread scaling from /tmp/bench_threads.tsv.

Writes benchmarks/bio/plots/threads_speed.png and threads_speedup.png.
"""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
OUT = HERE / "plots"
OUT.mkdir(exist_ok=True)

# ---- Load ----
data = defaultdict(list)  # tool -> [(threads, time_s, hits), ...]
with open("/tmp/bench_threads.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        data[r["tool"]].append((int(r["threads"]), float(r["time_s"]), int(r["hits"])))

# Sort each tool's series by thread count.
for k in data:
    data[k].sort(key=lambda t: t[0])

COLORS = {
    "BLAST":   "#3b82f6",
    "SWORD":   "#a855f7",
    "DIAMOND-sensitive": "#f59e0b",
    "SABER-default": "#ef4444",
    "SABER-fast":    "#dc2626",
}

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.3, "grid.linestyle": "--",
    "figure.dpi": 120,
})

# ---- Plot 1: raw time vs threads ----
fig, ax = plt.subplots(figsize=(10, 6))
for tool, series in data.items():
    ts = [t for t, _, _ in series]
    times = [s for _, s, _ in series]
    color = COLORS.get(tool, "#888")
    ax.plot(ts, times, marker="o", markersize=8, linewidth=2.2,
            label=tool, color=color)
    # value labels
    for t, s in zip(ts, times):
        ax.annotate(f"{s:.1f}s", (t, s), textcoords="offset points",
                    xytext=(0, 8), fontsize=9, ha="center", color=color)

ax.set_xlabel("Threads")
ax.set_ylabel("Wall-clock time (s)")
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_title("Multi-thread scaling: 100 UniProt queries × 20K UniProt subjects",
             pad=12)
ax.legend(loc="upper right", framealpha=0.95)
all_ts = sorted({t for series in data.values() for t, _, _ in series})
ax.set_xticks(all_ts)
ax.set_xticklabels([str(t) for t in all_ts])
plt.tight_layout()
plt.savefig(OUT / "threads_speed.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Wrote {OUT / 'threads_speed.png'}")

# ---- Plot 2: speedup vs 1-thread baseline ----
fig, ax = plt.subplots(figsize=(10, 6))
for tool, series in data.items():
    # find the t=1 baseline
    base = next((s for t, s, _ in series if t == 1), None)
    if base is None:
        continue
    ts = [t for t, _, _ in series]
    speedup = [base / s for _, s, _ in series]
    color = COLORS.get(tool, "#888")
    ax.plot(ts, speedup, marker="o", markersize=8, linewidth=2.2,
            label=tool, color=color)

# Ideal-scaling reference
max_t = max(t for series in data.values() for t, _, _ in series) if data else 16
xs = sorted({t for series in data.values() for t, _, _ in series})
if xs:
    ax.plot(xs, xs, "--", color="#666", alpha=0.6, label="ideal (linear)")

ax.set_xlabel("Threads")
ax.set_ylabel("Speedup vs 1-thread")
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=2)
ax.set_title("Multi-thread scaling efficiency", pad=12)
ax.legend(loc="upper left", framealpha=0.95)
ax.set_xticks(xs)
ax.set_xticklabels([str(t) for t in xs])
ax.set_yticks(xs)
ax.set_yticklabels([str(t) + "×" for t in xs])
plt.tight_layout()
plt.savefig(OUT / "threads_speedup.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Wrote {OUT / 'threads_speedup.png'}")

print("\nDone.")
