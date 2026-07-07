#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


OUT = Path("/disk228/huangjie/PanOmics/CosmicScanner/results/paper_plus_cosmic_runtime")
OUT.mkdir(parents=True, exist_ok=True)

traits = [2, 4, 6, 8, 10]

# Paper values are approximate visual reads from EMmvGWAS Figure 1.
# The paper text explicitly states human GEMMA ~200 h and human EMmvGWAS <40 h.
paper_rows = []
for dataset, method, values in [
    ("mouse", "GEMMA", [0.05, 0.15, 1.05, 1.60, 6.05]),
    ("mouse", "EMmvGWAS", [0.90, 1.20, 1.05, 1.00, 1.30]),
    ("human", "GEMMA", [6.0, 8.0, 20.0, 100.0, 195.0]),
    ("human", "EMmvGWAS", [25.0, 33.0, 34.0, 36.0, 36.0]),
]:
    for t, hours in zip(traits, values):
        paper_rows.append(
            {
                "dataset": dataset,
                "method": method,
                "traits": t,
                "runtime_hours": hours,
                "source": "paper_figure1_approx",
            }
        )

cosmic_mouse_seconds = {
    2: 23.00,
    4: 26.17,
    6: 28.09,
    8: 31.92,
    10: 38.36,  # same Dryad dosage trait-scale benchmark series
}
cosmic_human_seconds = {
    2: 882.804,
    4: 849.028,
    6: 930.568,
    8: 983.350,
    10: 1102.55,
}
cosmic_rows = []
for dataset, seconds_by_trait in [("mouse", cosmic_mouse_seconds), ("human", cosmic_human_seconds)]:
    for t in traits:
        cosmic_rows.append(
            {
                "dataset": dataset,
                "method": "CosmicScanner",
                "traits": t,
                "runtime_hours": seconds_by_trait[t] / 3600.0,
                "source": "local_measured",
            }
        )

df = pd.DataFrame(paper_rows + cosmic_rows)
df.to_csv(OUT / "mouse_human_paper_plus_cosmic_runtime.tsv", sep="\t", index=False)

styles = {
    "GEMMA": {"color": "#e45d5d", "marker": "o", "linestyle": "-", "linewidth": 2.2},
    "EMmvGWAS": {"color": "#22a7a7", "marker": "^", "linestyle": "--", "linewidth": 2.2},
    "CosmicScanner": {"color": "#5b5bd6", "marker": "s", "linestyle": "-.", "linewidth": 2.4},
}


def draw(path, log_scale=False):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), dpi=180, sharex=True)
    for ax, dataset in zip(axes, ["mouse", "human"]):
        sub = df[df["dataset"] == dataset]
        for method in ["GEMMA", "EMmvGWAS", "CosmicScanner"]:
            m = sub[sub["method"] == method].sort_values("traits")
            ax.plot(m["traits"], m["runtime_hours"], label=method, **styles[method])
        ax.set_title(dataset)
        ax.set_xlabel("Number of traits")
        ax.set_xticks(traits)
        ax.grid(alpha=0.25)
        if log_scale:
            ax.set_yscale("log")
            ax.set_ylabel("Runtime (hours, log scale)")
        else:
            ax.set_ylabel("Runtime (hours)")
            ax.set_ylim(bottom=0)
    axes[0].legend(loc="upper left", frameon=False)
    if log_scale:
        title = "EMmvGWAS Figure 1 plus CosmicScanner (log scale)"
        note = "GEMMA/EMmvGWAS: approximate values digitized from paper Figure 1; CosmicScanner: local measured full runs."
    else:
        title = "EMmvGWAS Figure 1 plus CosmicScanner"
        note = "Linear scale preserves the paper style; CosmicScanner is near the x-axis because it is measured in seconds/minutes."
    fig.suptitle(title, y=0.98, fontsize=16)
    fig.text(0.5, 0.02, note, ha="center", fontsize=8.5, color="#444444")
    fig.subplots_adjust(bottom=0.18, top=0.84, wspace=0.22)
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")


draw(OUT / "figure1_mouse_human_plus_cosmic_linear.png", log_scale=False)
draw(OUT / "figure1_mouse_human_plus_cosmic_log.png", log_scale=True)

print(OUT / "mouse_human_paper_plus_cosmic_runtime.tsv")
print(OUT / "figure1_mouse_human_plus_cosmic_linear.png")
print(OUT / "figure1_mouse_human_plus_cosmic_log.png")
