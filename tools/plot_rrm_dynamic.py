#!/usr/bin/env python3
"""Plot RRM dynamic GWAS effect curves, heatmaps, and method comparisons."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_dynamic(curve_path: Path, summary_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    curve = pd.read_csv(curve_path, sep="\t")
    summary = pd.read_csv(summary_path, sep="\t")
    for col in ["Rank", "Time", "Beta", "SE", "P", "JointP"]:
        if col in curve.columns:
            curve[col] = pd.to_numeric(curve[col], errors="coerce")
    for col in ["Rank", "JointP", "BestTime", "BestTimeP", "DeltaBeta", "BetaRange", "MaxAbsBeta", "EffectRatio", "PeakTime"]:
        if col in summary.columns:
            summary[col] = pd.to_numeric(summary[col], errors="coerce")
    return curve, summary


def plot_effect_curves(curve: pd.DataFrame, summary: pd.DataFrame, outdir: Path, top: int) -> None:
    sub = curve[curve["Rank"] <= top].copy()
    if sub.empty:
        return
    pattern_map = summary.set_index("SNP")["Pattern"].to_dict() if "Pattern" in summary.columns else {}

    fig, ax = plt.subplots(figsize=(8, 5))
    for snp, snp_df in sub.groupby("SNP", sort=False):
        snp_df = snp_df.sort_values("Time")
        rank = int(snp_df["Rank"].iloc[0])
        label = f"{rank}. {snp} ({pattern_map.get(snp, 'NA')})"
        ax.plot(snp_df["Time"], snp_df["Beta"], marker="o", linewidth=1.8, label=label)
        if "SE" in snp_df.columns and snp_df["SE"].notna().any():
            lo = snp_df["Beta"] - 1.96 * snp_df["SE"]
            hi = snp_df["Beta"] + 1.96 * snp_df["SE"]
            ax.fill_between(snp_df["Time"].to_numpy(), lo.to_numpy(), hi.to_numpy(), alpha=0.12)
    ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.6)
    ax.set_xlabel("Time / environmental gradient")
    ax.set_ylabel("SNP effect beta")
    ax.set_title(f"RRM Dynamic GWAS Top {top} Effect Curves")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(outdir / "rrm_dynamic_effect_curves.png", dpi=300)
    fig.savefig(outdir / "rrm_dynamic_effect_curves.pdf")
    plt.close(fig)


def plot_effect_heatmap(curve: pd.DataFrame, summary: pd.DataFrame, outdir: Path, top: int) -> None:
    sub = curve[curve["Rank"] <= top].copy()
    if sub.empty:
        return
    pivot = sub.pivot_table(index=["Rank", "SNP"], columns="Time", values="Beta", aggfunc="first").sort_index()
    data = pivot.to_numpy(dtype=float)
    if data.size == 0:
        return
    vmax = np.nanmax(np.abs(data))
    vmax = vmax if np.isfinite(vmax) and vmax > 0 else 1.0

    labels = []
    pattern_map = summary.set_index("SNP")["Pattern"].to_dict() if "Pattern" in summary.columns else {}
    for rank, snp in pivot.index:
        labels.append(f"{int(rank)} {snp}\n{pattern_map.get(snp, 'NA')}")

    fig_h = max(4.0, 0.45 * len(labels) + 1.2)
    fig, ax = plt.subplots(figsize=(8, fig_h))
    im = ax.imshow(data, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels([f"{float(x):g}" for x in pivot.columns])
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Time / environmental gradient")
    ax.set_title(f"RRM Dynamic GWAS Top {top} Effect Heatmap")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Beta")
    fig.tight_layout()
    fig.savefig(outdir / "rrm_dynamic_effect_heatmap.png", dpi=300)
    fig.savefig(outdir / "rrm_dynamic_effect_heatmap.pdf")
    plt.close(fig)


def plot_comparison(comparison_path: Path, outdir: Path) -> None:
    if not comparison_path or not comparison_path.exists():
        return
    df = pd.read_csv(comparison_path, sep="\t")
    df["RankNumeric"] = pd.to_numeric(df["Rank"], errors="coerce")
    if df.empty:
        return
    plot_df = df.copy()
    missing_rank = plot_df["RankNumeric"].max(skipna=True)
    missing_rank = 1.0 if not np.isfinite(missing_rank) else missing_rank + 5.0
    plot_df["RankPlot"] = plot_df["RankNumeric"].fillna(missing_rank)
    plot_df["Label"] = plot_df["Method"] + "\n" + plot_df["SNP"]

    fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(plot_df)), 5))
    colors = ["#4C78A8" if m == "RRM_dynamic" else "#F58518" if m == "MvGWAS" else "#54A24B" for m in plot_df["Method"]]
    ax.bar(np.arange(len(plot_df)), plot_df["RankPlot"], color=colors)
    ax.set_xticks(np.arange(len(plot_df)))
    ax.set_xticklabels(plot_df["Label"], rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Causal SNP rank (lower is better)")
    ax.set_title("Dynamic GWAS Simulation: Method Rank Comparison")
    ax.grid(True, axis="y", alpha=0.25)
    for i, row in plot_df.iterrows():
        if pd.isna(row["RankNumeric"]):
            ax.text(i, row["RankPlot"], "NA", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "rrm_dynamic_method_comparison.png", dpi=300)
    fig.savefig(outdir / "rrm_dynamic_method_comparison.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curve", type=Path, required=True, help="*.rrm.dynamic_curve.txt")
    parser.add_argument("--summary", type=Path, required=True, help="*.rrm.dynamic_summary.txt")
    parser.add_argument("--comparison", type=Path, help="comparison.tsv from simulation benchmark")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--top", type=int, default=5)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    curve, summary = read_dynamic(args.curve, args.summary)
    top = max(1, args.top)
    plot_effect_curves(curve, summary, args.outdir, top)
    plot_effect_heatmap(curve, summary, args.outdir, top)
    plot_comparison(args.comparison, args.outdir)
    summary.to_csv(args.outdir / "rrm_dynamic_plot_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
