#!/usr/bin/env python3
"""Plot repeatability GWAS benchmark and RepeatABEL concordance figures."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_benchmark(paths: list[Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        if path.exists():
            frames.append(pd.read_csv(path, sep="\t"))
    if not frames:
        raise SystemExit("No benchmark summary files were found.")
    df = pd.concat(frames, ignore_index=True)
    df["runtime_sec"] = pd.to_numeric(df["runtime_sec"], errors="coerce")
    df["max_rss_mb"] = pd.to_numeric(df["max_rss_mb"], errors="coerce")
    for col in ["n_subjects", "n_records_per_subject", "n_records", "n_snps", "n_output_snps"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def method_label(method: str) -> str:
    return {
        "exact": "CS Exact/Woodbury",
        "fast": "CS Fast approximation",
        "mean": "Subject mean LMM",
        "repeatabel": "RepeatABEL rGLS",
    }.get(method, method)


def plot_runtime_memory(df: pd.DataFrame, outdir: Path) -> None:
    ok = df[df["status"].isin(["ok", "timeout"])].copy()
    ok["method_label"] = ok["method"].map(method_label)

    fig, ax = plt.subplots(figsize=(8, 5))
    for label, sub in ok.groupby("method_label"):
        sub = sub.sort_values(["n_subjects", "n_records_per_subject", "n_snps"])
        marker = "x" if "RepeatABEL" in label else "o"
        ax.scatter(sub["n_records"], sub["runtime_sec"], label=label, marker=marker, s=55)
        for _, row in sub.iterrows():
            if row["status"] == "timeout" and math.isfinite(row["runtime_sec"]):
                ax.annotate("timeout", (row["n_records"], row["runtime_sec"]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Number of phenotype records")
    ax.set_ylabel("Runtime (seconds)")
    ax.set_yscale("log")
    ax.set_title("Repeatability GWAS Runtime Scaling")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "repeatability_runtime_scaling.png", dpi=300)
    fig.savefig(outdir / "repeatability_runtime_scaling.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for label, sub in ok.groupby("method_label"):
        sub = sub.sort_values(["n_subjects", "n_records_per_subject", "n_snps"])
        marker = "x" if "RepeatABEL" in label else "o"
        ax.scatter(sub["n_records"], sub["max_rss_mb"], label=label, marker=marker, s=55)
    ax.set_xlabel("Number of phenotype records")
    ax.set_ylabel("Peak memory (MB)")
    ax.set_yscale("log")
    ax.set_title("Repeatability GWAS Memory Scaling")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "repeatability_memory_scaling.png", dpi=300)
    fig.savefig(outdir / "repeatability_memory_scaling.pdf")
    plt.close(fig)


def read_accuracy_pair(cs_path: Path, rb_path: Path) -> pd.DataFrame:
    cs = pd.read_csv(cs_path, sep="\t")
    rb = pd.read_csv(rb_path, sep="\t")
    cs = cs.rename(columns={"BETA": "cs_beta", "SE": "cs_se", "P": "cs_p"})
    rb = rb.rename(columns={"effB": "rb_beta", "se_effB": "rb_se", "P1df": "rb_p"})
    merged = cs[["SNP", "cs_beta", "cs_se", "cs_p"]].merge(
        rb[["SNP", "rb_beta", "rb_se", "rb_p"]], on="SNP", how="inner"
    )
    for col in ["cs_beta", "cs_se", "cs_p", "rb_beta", "rb_se", "rb_p"]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce")
    merged = merged.replace([np.inf, -np.inf], np.nan).dropna()
    merged = merged[(merged["cs_p"] > 0) & (merged["rb_p"] > 0)]
    merged["cs_log10p"] = -np.log10(merged["cs_p"])
    merged["rb_log10p"] = -np.log10(merged["rb_p"])
    if merged["cs_beta"].corr(merged["rb_beta"]) < 0:
        merged["cs_beta_oriented"] = -merged["cs_beta"]
    else:
        merged["cs_beta_oriented"] = merged["cs_beta"]
    return merged


def plot_concordance(cs_path: Path, rb_path: Path, outdir: Path, label: str) -> None:
    df = read_accuracy_pair(cs_path, rb_path)
    if df.empty:
        raise SystemExit(f"No overlapping valid SNPs for {label}")

    panels = [
        ("cs_beta_oriented", "rb_beta", "Effect size beta", "repeatability_beta_concordance"),
        ("cs_se", "rb_se", "Standard error", "repeatability_se_concordance"),
        ("cs_log10p", "rb_log10p", "-log10(P)", "repeatability_log10p_concordance"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, (xcol, ycol, title, _) in zip(axes, panels):
        ax.scatter(df[xcol], df[ycol], s=12, alpha=0.65, edgecolor="none")
        lo = float(np.nanmin([df[xcol].min(), df[ycol].min()]))
        hi = float(np.nanmax([df[xcol].max(), df[ycol].max()]))
        ax.plot([lo, hi], [lo, hi], color="black", linewidth=1, alpha=0.7)
        corr = df[xcol].corr(df[ycol])
        ax.set_title(f"{title}\nr={corr:.4f}")
        ax.set_xlabel("CosmicScanner")
        ax.set_ylabel("RepeatABEL")
        ax.grid(True, alpha=0.25)
    fig.suptitle(f"CS Exact/Woodbury vs RepeatABEL ({label}, n={len(df)} SNPs)")
    fig.tight_layout()
    fig.savefig(outdir / f"repeatability_concordance_{label}.png", dpi=300)
    fig.savefig(outdir / f"repeatability_concordance_{label}.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", action="append", type=Path, default=[], help="benchmark.summary.tsv or combined TSV")
    parser.add_argument("--accuracy-cs", type=Path, help="CS exact assoc.txt for accuracy scatter")
    parser.add_argument("--accuracy-repeatabel", type=Path, help="RepeatABEL full TSV for accuracy scatter")
    parser.add_argument("--accuracy-label", default="accuracy", help="label for accuracy output filenames")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.benchmark:
        df = read_benchmark(args.benchmark)
        df.to_csv(args.outdir / "repeatability_benchmark_plot_data.tsv", sep="\t", index=False)
        plot_runtime_memory(df, args.outdir)
    if args.accuracy_cs and args.accuracy_repeatabel:
        plot_concordance(args.accuracy_cs, args.accuracy_repeatabel, args.outdir, args.accuracy_label)


if __name__ == "__main__":
    main()
