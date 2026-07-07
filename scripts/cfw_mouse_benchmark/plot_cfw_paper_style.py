#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_assoc(path, method):
    df = pd.read_csv(path, sep=("," if str(path).endswith(".csv") else "\t"))
    if method == "gemma":
        df = df.rename(columns={"rs": "SNP", "p_wald": "P"})
    elif method == "emmv":
        df = df.rename(columns={"SNP_ID": "SNP", "p": "P"})
    else:
        df = df.rename(columns={"P": "P"})
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    return df[["SNP", "P"]].dropna()


def neglogp(p):
    return -np.log10(pd.to_numeric(p, errors="coerce").clip(lower=1e-300, upper=1))


def corr_logp(a, b):
    x = neglogp(a)
    y = neglogp(b)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() <= 2:
        return np.nan, int(ok.sum())
    return float(np.corrcoef(x[ok], y[ok])[0, 1]), int(ok.sum())


def load_runtime(set_dir):
    rows = []
    for tpath in sorted(set_dir.glob("traits_*/logs/*.time")):
        traits = int(tpath.parts[-3].split("_")[-1])
        method = tpath.stem
        vals = {"seconds": np.nan, "max_rss_kb": np.nan}
        for line in tpath.read_text(encoding="utf-8").splitlines():
            parts = line.split()
            if len(parts) == 2:
                vals[parts[0]] = float(parts[1])
        rows.append({"traits": traits, "method": method, **vals})
    return pd.DataFrame(rows)


def plot_runtime(rt, out):
    mapping = {
        "gemma": ("GEMMA", "#d95f5f", "-", "o"),
        "emmvgwas": ("EMmvGWAS", "#2b9c9c", "--", "^"),
        "cosmic_emmv": ("CosmicScanner", "#5b5bd6", "-.", "s"),
    }
    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    for key, (label, color, ls, marker) in mapping.items():
        g = rt[rt["method"] == key].sort_values("traits")
        if g.empty:
            continue
        ax.plot(g["traits"], g["seconds"] / 3600.0, color=color, ls=ls, marker=marker, label=label)
    ax.set_title("mouse")
    ax.set_xlabel("Number of traits")
    ax.set_ylabel("Runtime (hours)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def plot_scatter(joined, out):
    specs = [
        ("P_gemma", "P_emmv", "EMmvGWAS", "GEMMA"),
        ("P_gemma", "P_cs_emmv", "CosmicScanner", "GEMMA"),
        ("P_emmv", "P_cs_emmv", "CosmicScanner", "EMmvGWAS"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.2))
    rows = []
    for ax, (xcol, ycol, ylab, xlab) in zip(axes, specs):
        x = neglogp(joined[xcol])
        y = neglogp(joined[ycol])
        ok = np.isfinite(x) & np.isfinite(y)
        r = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 2 else np.nan
        ax.scatter(x[ok], y[ok], s=4, c="#222222", alpha=0.45, linewidths=0)
        if ok.sum() > 2:
            coef = np.polyfit(x[ok], y[ok], 1)
            xx = np.linspace(0, max(x[ok].max(), 1), 100)
            ax.plot(xx, coef[0] * xx + coef[1], color="#d95f5f", lw=1)
        ax.set_xlabel(f"{xlab} -log10(p)")
        ax.set_ylabel(f"{ylab} -log10(p)")
        ax.set_title(f"mouse\nr = {r:.3f}", fontsize=10)
        ax.grid(alpha=0.2)
        rows.append({"comparison": f"{ylab}_vs_{xlab}", "n": int(ok.sum()), "pearson_logp": r})
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return pd.DataFrame(rows)


def add_genome_x(bim):
    b = bim.copy()
    b["CHR_NUM"] = pd.to_numeric(b["CHR"], errors="coerce")
    b["BP_NUM"] = pd.to_numeric(b["BP"], errors="coerce")
    b = b.sort_values(["CHR_NUM", "BP_NUM", "SNP"]).reset_index(drop=True)
    offsets = {}
    offset = 0
    ticks = []
    labels = []
    for chrom, g in b.groupby("CHR_NUM", sort=True):
        if not np.isfinite(chrom):
            continue
        offsets[chrom] = offset
        bp = g["BP_NUM"].to_numpy(dtype=float)
        fallback = np.arange(len(g), dtype=float)
        x = np.where(np.isfinite(bp), bp, fallback) + offset
        b.loc[g.index, "GENOME_X"] = x
        ticks.append((x.min() + x.max()) / 2.0)
        labels.append(str(int(chrom)))
        offset = x.max() + 1_000_000
    return b, ticks, labels


def plot_manhattan(bim, assoc, out):
    methods = [
        ("EMmvGWAS", assoc["emmv"]),
        ("GEMMA", assoc["gemma"]),
        ("CosmicScanner", assoc["cosmic"]),
    ]
    b, ticks, labels = add_genome_x(bim)
    threshold = -np.log10(0.05 / max(1, len(b)))
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
    colors = {0: "#222222", 1: "#1687b9"}
    for ax, (name, df) in zip(axes, methods):
        m = b.merge(df, on="SNP", how="inner")
        y = neglogp(m["P"])
        c = [colors[int(chrom) % 2] for chrom in m["CHR_NUM"]]
        ax.scatter(m["GENOME_X"], y, c=c, s=3, alpha=0.75, linewidths=0)
        ax.axhline(threshold, color="#d95f5f", lw=1)
        ax.set_title(name)
        ax.set_xlabel("Chromosome")
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, fontsize=7)
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("-log10(p)")
    fig.suptitle("Mouse", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_qq(assoc, out):
    methods = [("EMmvGWAS", assoc["emmv"]), ("GEMMA", assoc["gemma"]), ("CosmicScanner", assoc["cosmic"])]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.2), sharex=True, sharey=True)
    rows = []
    for ax, (name, df) in zip(axes, methods):
        p = pd.to_numeric(df["P"], errors="coerce").dropna().clip(lower=1e-300, upper=1).sort_values().to_numpy()
        exp = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
        obs = -np.log10(p)
        ax.scatter(exp, obs, s=3, c="#222222", alpha=0.55, linewidths=0)
        mx = max(float(exp.max()), float(obs.max()), 1.0)
        ax.plot([0, mx], [0, mx], color="#d95f5f", lw=1)
        ax.set_title(name)
        ax.set_xlabel("Expected -log10(p)")
        ax.grid(alpha=0.2)
        chisq = np.nanmedian(np.square(np.sqrt(2) * 0 + 0))
        rows.append({"method": name, "n": len(p), "max_logp": float(obs.max()), "lambda_gc": np.nan if not np.isfinite(chisq) else np.nan})
    axes[0].set_ylabel("Observed -log10(p)")
    fig.suptitle("Mouse", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--traits", type=int, default=10)
    args = ap.parse_args()

    set_dir = Path(args.set_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    res = set_dir / f"traits_{args.traits}" / "results"
    bim = pd.read_csv(set_dir / "plink.bim", sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    assoc = {
        "gemma": read_assoc(res / "gemma/output/gemma.assoc.txt", "gemma"),
        "emmv": read_assoc(res / "emmvgwas/emmv_association_with_ids.csv", "emmv"),
        "cosmic": read_assoc(res / "cosmicscanner_emmv/cosmic_emmv.assoc.mv.txt", "cosmic"),
    }
    joined = assoc["gemma"].rename(columns={"P": "P_gemma"}) \
        .merge(assoc["emmv"].rename(columns={"P": "P_emmv"}), on="SNP", how="inner") \
        .merge(assoc["cosmic"].rename(columns={"P": "P_cs_emmv"}), on="SNP", how="inner")
    joined.to_csv(out_dir / "cfw_mouse_paper_style_joined.tsv", sep="\t", index=False)

    rt = load_runtime(set_dir)
    rt.to_csv(out_dir / "cfw_mouse_paper_style_runtime.tsv", sep="\t", index=False)
    plot_runtime(rt, out_dir / "figure1_mouse_runtime_plus_cosmic.png")
    scatter_summary = plot_scatter(joined, out_dir / "figure2_mouse_scatter_plus_cosmic.png")
    plot_manhattan(bim, assoc, out_dir / "figure3_mouse_manhattan_plus_cosmic.png")
    qq_summary = plot_qq(assoc, out_dir / "figure4_mouse_qq_observed_plus_cosmic.png")

    summary = scatter_summary.drop_duplicates("comparison").reset_index(drop=True)
    summary.to_csv(out_dir / "cfw_mouse_paper_style_summary.tsv", sep="\t", index=False)
    qq_summary.to_csv(out_dir / "cfw_mouse_paper_style_qq_summary.tsv", sep="\t", index=False)

    report = [
        "# CFW Mouse Paper-Style Figures Plus CosmicScanner",
        "",
        f"Input set: `{set_dir.name}`",
        f"Traits: {args.traits}",
        f"SNPs in BIM: {len(bim)}",
        f"Joined SNPs: {len(joined)}",
        "",
        "## Figures",
        "",
        f"- Figure 1 style runtime: `{out_dir / 'figure1_mouse_runtime_plus_cosmic.png'}`",
        f"- Figure 2 style scatter: `{out_dir / 'figure2_mouse_scatter_plus_cosmic.png'}`",
        f"- Figure 3 style Manhattan: `{out_dir / 'figure3_mouse_manhattan_plus_cosmic.png'}`",
        f"- Figure 4 style observed QQ: `{out_dir / 'figure4_mouse_qq_observed_plus_cosmic.png'}`",
        "",
        "## Summary",
        "",
        "```tsv",
        summary.to_csv(sep="\t", index=False).strip(),
        "```",
        "",
        "## Limitation",
        "",
        "These figures use the validated harmonized hard-call/imputed CFW set, not the original Dryad dosage `geno.txt`/`map.txt` files. They reproduce the EMmvGWAS paper figure types and add CosmicScanner, but are not a byte-for-byte reproduction of the manuscript figures.",
    ]
    (out_dir / "cfw_mouse_paper_style_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
