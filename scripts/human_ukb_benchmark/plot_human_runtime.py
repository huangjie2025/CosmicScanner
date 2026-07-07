#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt


BASE = Path("/disk228/huangjie/PanOmics/CosmicScanner/results/human_ukb_benchmark")
OUT = BASE / "figures"
OUT.mkdir(parents=True, exist_ok=True)

cosmic_trait_rows = [
    (2, 882.804, "14:43.15", 4301.94),
    (4, 849.028, "14:09.40", 4558.31),
    (6, 930.568, "15:31.00", 4976.61),
    (8, 983.350, "16:23.90", 5563.96),
    (10, 1102.55, "18:23.21", 6293.55),
]
cosmic_seconds = cosmic_trait_rows[-1][1]
cosmic_hours = cosmic_seconds / 3600.0

rows = [
    ("GEMMA", "paper reference", 200.0, "near 200 h in paper"),
    ("EMmvGWAS", "paper reference", 40.0, "within 40 h in paper"),
    ("CosmicScanner", "measured", cosmic_hours, "1102.55 s / 18:23.21 wall"),
]

tsv = OUT / "human_10trait_runtime_comparison.tsv"
with tsv.open("w") as handle:
    handle.write("method\tsource\truntime_hours\tnote\n")
    for method, source, hours, note in rows:
        handle.write(f"{method}\t{source}\t{hours:.8g}\t{note}\n")

curve_tsv = OUT / "human_cosmicscanner_trait_scale_runtime.tsv"
with curve_tsv.open("w") as handle:
    handle.write("traits\tcomputation_seconds\tcomputation_hours\twall_time\tpeak_memory_mb\tsnps\tsamples\n")
    for traits, seconds, wall, mem in cosmic_trait_rows:
        handle.write(f"{traits}\t{seconds:.6g}\t{seconds / 3600.0:.8g}\t{wall}\t{mem:.6g}\t2026362\t5000\n")

plt.rcParams.update({
    "font.size": 12,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

methods = [r[0] for r in rows]
hours = [r[2] for r in rows]
colors = ["#d95f5f", "#2b9ca3", "#5b5bd6"]

fig, ax = plt.subplots(figsize=(7.2, 4.9), dpi=180)
bars = ax.bar(methods, hours, color=colors, width=0.62)
ax.set_yscale("log")
ax.set_ylabel("Runtime (hours, log scale)")
ax.set_title("Human 10-trait multivariate GWAS")
ax.grid(axis="y", which="both", alpha=0.25)

labels = ["~200 h", "<40 h", "0.306 h\n18.4 min"]
for bar, label in zip(bars, labels):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() * 1.18,
        label,
        ha="center",
        va="bottom",
        fontsize=11,
    )

ax.text(
    0.5,
    -0.22,
    "GEMMA/EMmvGWAS: EMmvGWAS paper Figure 1 text; CosmicScanner: local full PGEN run, n=5000, m=2,026,362, p=10, 10 threads.",
    transform=ax.transAxes,
    fontsize=8.5,
    color="#444444",
    ha="center",
)
fig.subplots_adjust(bottom=0.24, top=0.88)
fig.savefig(OUT / "human_10trait_runtime_comparison_loghours.png", bbox_inches="tight")
fig.savefig(OUT / "human_10trait_runtime_comparison_loghours.pdf", bbox_inches="tight")

fig2, ax2 = plt.subplots(figsize=(7.2, 4.2), dpi=180)
speedups = [200.0 / cosmic_hours, 40.0 / cosmic_hours]
bars2 = ax2.bar(["vs GEMMA\npaper ~200 h", "vs EMmvGWAS\npaper <40 h"], speedups, color=["#d95f5f", "#2b9ca3"], width=0.55)
ax2.set_ylabel("Runtime ratio to CosmicScanner")
ax2.set_title("Human 10-trait speed comparison")
ax2.grid(axis="y", alpha=0.25)
for bar, value in zip(bars2, speedups):
    ax2.text(bar.get_x() + bar.get_width() / 2, value * 1.02, f"{value:.0f}x", ha="center", va="bottom", fontsize=12)
ax2.text(
    0.5,
    -0.22,
    "Ratios use paper reference values; EMmvGWAS ratio is against the 40 h upper reference.",
    transform=ax2.transAxes,
    fontsize=8.5,
    color="#444444",
    ha="center",
)
fig2.subplots_adjust(bottom=0.24, top=0.88)
fig2.savefig(OUT / "human_10trait_speedup_vs_paper_reference.png", bbox_inches="tight")
fig2.savefig(OUT / "human_10trait_speedup_vs_paper_reference.pdf", bbox_inches="tight")

fig3, ax3 = plt.subplots(figsize=(7.4, 4.6), dpi=180)
traits = [r[0] for r in cosmic_trait_rows]
hours = [r[1] / 3600.0 for r in cosmic_trait_rows]
mem_gb = [r[3] / 1024.0 for r in cosmic_trait_rows]
ax3.plot(traits, hours, color="#5b5bd6", marker="s", linewidth=2.4, markersize=6, label="Runtime")
ax3.set_xlabel("Number of traits")
ax3.set_ylabel("Runtime (hours)")
ax3.set_xticks(traits)
ax3.grid(alpha=0.25)
ax3.set_title("CosmicScanner human full-scale runtime")
for x, y, row in zip(traits, hours, cosmic_trait_rows):
    ax3.text(x, y + 0.012, row[2].replace(":", "m ", 1) + "s", ha="center", fontsize=8.5)
ax4 = ax3.twinx()
ax4.plot(traits, mem_gb, color="#4b4b4b", marker="o", linewidth=1.8, linestyle="--", label="Peak memory")
ax4.set_ylabel("Peak memory (GB)")
ax4.set_ylim(0, max(mem_gb) * 1.25)
lines, labels = ax3.get_legend_handles_labels()
lines2, labels2 = ax4.get_legend_handles_labels()
ax3.legend(lines + lines2, labels + labels2, loc="upper left", frameon=False)
ax3.text(
    0.5,
    -0.23,
    "All points use n=5000, m=2,026,362, PGEN input, --mv-method emmv --mv-compatible emmvgwas, --no-loco, 10 threads.",
    transform=ax3.transAxes,
    fontsize=8.5,
    color="#444444",
    ha="center",
)
fig3.subplots_adjust(bottom=0.25, top=0.88, right=0.86)
fig3.savefig(OUT / "human_cosmicscanner_trait_scale_runtime.png", bbox_inches="tight")
fig3.savefig(OUT / "human_cosmicscanner_trait_scale_runtime.pdf", bbox_inches="tight")

print(tsv)
print(curve_tsv)
print(OUT / "human_10trait_runtime_comparison_loghours.png")
print(OUT / "human_10trait_speedup_vs_paper_reference.png")
print(OUT / "human_cosmicscanner_trait_scale_runtime.png")
