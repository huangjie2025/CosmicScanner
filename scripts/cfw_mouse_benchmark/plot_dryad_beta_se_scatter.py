#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path("/disk228/huangjie/PanOmics/CosmicScanner/results/cfw_mouse_benchmark")
TRAIT_DIR = BASE / "scaled_harmonized/dryad_dosage/traits_10/results"
OUT = BASE / "paper_style_dryad_dosage"

GEMMA = TRAIT_DIR / "gemma_bimbam/output/gemma.assoc.txt"
EMMV = TRAIT_DIR / "emmvgwas/emmv_association_with_ids.csv"
COSMIC = TRAIT_DIR / "cosmicscanner_emmv_bgen/cosmic_emmv_bgen.assoc.mv.txt"

TRAITS = [
    "TA",
    "EDL",
    "bw1",
    "PPIweight",
    "Tail length",
    "tibia",
    "Testis weight",
    "gastroc",
    "Plantaris",
    "Soleus",
]


def pearson(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() <= 2:
        return np.nan, int(ok.sum())
    return float(np.corrcoef(x[ok], y[ok])[0, 1]), int(ok.sum())


def draw_identity_scatter(ax, x, y, xlab, ylab, title, value_label):
    ok = np.isfinite(x) & np.isfinite(y)
    r, n = pearson(x, y)
    ax.scatter(x[ok], y[ok], s=1.2, c="#222222", alpha=0.18, linewidths=0)
    if ok.sum() > 2:
        lo = float(min(np.nanmin(x[ok]), np.nanmin(y[ok])))
        hi = float(max(np.nanmax(x[ok]), np.nanmax(y[ok])))
        pad = 0.03 * max(hi - lo, 1e-12)
        lo -= pad
        hi += pad
        ax.plot([lo, hi], [lo, hi], color="#d95f5f", lw=1.2)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(f"{title}\nr = {r:.3f}", fontsize=10)
    ax.grid(alpha=0.2)
    return {"comparison": value_label, "n": n, "pearson": r}


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    gemma = pd.read_csv(GEMMA, sep="\t")
    emmv = pd.read_csv(EMMV)
    cosmic = pd.read_csv(COSMIC, sep="\t")

    gemma = gemma.rename(columns={"rs": "SNP"})
    emmv = emmv.rename(columns={"SNP_ID": "SNP"})

    gemma_cols = ["SNP"] + [f"beta_{i}" for i in range(1, 11)] + [f"Vbeta_{i}_{i}" for i in range(1, 11)]
    emmv_cols = ["SNP"] + TRAITS
    cosmic_cols = ["SNP"] + [f"BETA_{i}" for i in range(1, 11)] + [f"SE_{i}" for i in range(1, 11)]

    joined = (
        gemma[gemma_cols]
        .merge(emmv[emmv_cols], on="SNP", how="inner")
        .merge(cosmic[cosmic_cols], on="SNP", how="inner")
    )

    beta_gemma = []
    beta_emmv = []
    beta_cosmic = []
    se_gemma = []
    se_cosmic = []
    per_trait_rows = []

    for i, trait in enumerate(TRAITS, start=1):
        bg = pd.to_numeric(joined[f"beta_{i}"], errors="coerce").to_numpy(float)
        be = pd.to_numeric(joined[trait], errors="coerce").to_numpy(float)
        bc = pd.to_numeric(joined[f"BETA_{i}"], errors="coerce").to_numpy(float)
        sg = np.sqrt(pd.to_numeric(joined[f"Vbeta_{i}_{i}"], errors="coerce").clip(lower=0).to_numpy(float))
        sc = pd.to_numeric(joined[f"SE_{i}"], errors="coerce").to_numpy(float)

        beta_gemma.append(bg)
        beta_emmv.append(be)
        beta_cosmic.append(bc)
        se_gemma.append(sg)
        se_cosmic.append(sc)

        for label, x, y in [
            ("beta_EMmvGWAS_vs_GEMMA", bg, be),
            ("beta_CosmicScanner_vs_GEMMA", bg, bc),
            ("beta_CosmicScanner_vs_EMmvGWAS", be, bc),
            ("se_CosmicScanner_vs_GEMMA", sg, sc),
        ]:
            r, n = pearson(x, y)
            per_trait_rows.append({"trait": trait, "comparison": label, "n": n, "pearson": r})

    beta_gemma = np.concatenate(beta_gemma)
    beta_emmv = np.concatenate(beta_emmv)
    beta_cosmic = np.concatenate(beta_cosmic)
    se_gemma = np.concatenate(se_gemma)
    se_cosmic = np.concatenate(se_cosmic)

    summary_rows = []

    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.2))
    summary_rows.append(draw_identity_scatter(
        axes[0], beta_gemma, beta_emmv,
        "GEMMA beta", "EMmvGWAS beta", "mouse", "beta_EMmvGWAS_vs_GEMMA_all_traits"
    ))
    summary_rows.append(draw_identity_scatter(
        axes[1], beta_gemma, beta_cosmic,
        "GEMMA beta", "CosmicScanner beta", "mouse", "beta_CosmicScanner_vs_GEMMA_all_traits"
    ))
    summary_rows.append(draw_identity_scatter(
        axes[2], beta_emmv, beta_cosmic,
        "EMmvGWAS beta", "CosmicScanner beta", "mouse", "beta_CosmicScanner_vs_EMmvGWAS_all_traits"
    ))
    fig.tight_layout()
    fig.savefig(OUT / "figure2b_mouse_beta_scatter_dryad_dosage_plus_cosmic.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    summary_rows.append(draw_identity_scatter(
        ax, se_gemma, se_cosmic,
        "GEMMA SE", "CosmicScanner SE", "mouse", "se_CosmicScanner_vs_GEMMA_all_traits"
    ))
    fig.tight_layout()
    fig.savefig(OUT / "figure2c_mouse_se_scatter_dryad_dosage_gemma_vs_cosmic.png", dpi=300)
    plt.close(fig)

    pd.DataFrame(summary_rows).to_csv(
        OUT / "dryad_dosage_beta_se_scatter_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(per_trait_rows).to_csv(
        OUT / "dryad_dosage_beta_se_scatter_by_trait.tsv", sep="\t", index=False
    )

    print(OUT / "figure2b_mouse_beta_scatter_dryad_dosage_plus_cosmic.png")
    print(OUT / "figure2c_mouse_se_scatter_dryad_dosage_gemma_vs_cosmic.png")
    print(OUT / "dryad_dosage_beta_se_scatter_summary.tsv")


if __name__ == "__main__":
    main()
