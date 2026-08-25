#!/usr/bin/env python3
"""Publication-quality figure generation for CosmicScanner manuscript."""
import csv, math, os, sys
from pathlib import Path
import numpy as np

ROOT = Path("/disk228/huangjie/PanOmics/CosmicScanner")
RESULT_DIR = ROOT / "tests" / "results"
FIG_DIR = ROOT / "docs" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

os.environ["MPLBACKEND"] = "Agg"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.ticker as ticker

# ── Style configuration ──────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "figure.dpi": 180,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
})

CS_BLUE = "#2166ac"
CS_TEAL = "#35978f"
CS_ORANGE = "#d6604d"
CS_PURPLE = "#762a83"
CS_GRAY = "#666666"
REF_COLORS = {"CosmicScanner": CS_BLUE, "GEMMA": CS_ORANGE, "EMmvGWAS": CS_PURPLE, "RepeatABEL": CS_ORANGE}
CB_PALETTE = ["#2166ac", "#35978f", "#d6604d", "#762a83", "#999999"]


def read_tsv(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_csv(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ═══════════════════════════════════════════════════════════════════
# Figure 1: Architecture overview
# ═══════════════════════════════════════════════════════════════════
def fig1():
    fig, ax = plt.subplots(1, 1, figsize=(8.5, 4.5))
    ax.set_xlim(0, 10); ax.set_ylim(0, 5.5)
    ax.axis("off")

    C_IN = "#d4e6f1"
    C_CORE = "#d5f5e3"
    C_MODEL = "#fef9e7"
    C_OUT = "#fadbd8"
    C_INN = "#fdebd0"
    BORDER = "#2c3e50"

    def rbox(x, y, w, h, color, title, sub=", ", fs=8, sfs=6):
        r = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                           facecolor=color, edgecolor=BORDER, linewidth=0.8)
        ax.add_patch(r)
        ax.text(x+w/2, y+h*0.65, title, ha="center", va="center", fontsize=fs, fontweight="bold")
        ax.text(x+w/2, y+h*0.25, sub, ha="center", va="center", fontsize=sfs, color="#555555")

    # Row 0: Genotype inputs
    rbox(0.15, 4.7, 3.0, 0.55, C_IN, "Genotype Input",
         "PLINK 1.9  |  PLINK 2.0 (.pgen)  |  BGEN")

    # Row 1: Core infrastructure
    rbox(3.8, 3.9, 2.5, 0.55, C_CORE, "Shared Core Library",
         "Memory-mapped I/O  |  OpenMP  |  AVX2/AVX-512")
    rbox(6.8, 3.9, 2.5, 0.55, C_CORE, "GRM Builder",
         "VanRaden  |  6 paths  |  LOCO")

    # Row 2: Engine components
    rbox(0.15, 2.9, 1.7, 0.5, C_CORE, "LMM Engine",
         "Eigen-decomp  |  Rotation  |  Block scan")
    rbox(2.0, 2.9, 1.5, 0.5, C_CORE, "Null Model",
         "IDUL / AIDUL / Brent")
    rbox(3.7, 2.9, 1.5, 0.5, C_CORE, "VCE Engine",
         "AI-REML  |  EM-REML  |  Grid")
    rbox(5.5, 2.9, 1.5, 0.5, C_CORE, "Biobank",
         "REGENIE 2-step  |  Ridge")
    rbox(7.3, 2.9, 2.0, 0.5, C_CORE, "Input Adapter",
         "CSV/TSV  |  Wide/Long  |  Auto-convert")

    # Row 3: Model classes
    models = ["Single-trait\nLMM", "Multivariate\nLMM", "Repeatability\n(Woodbury\u2217)",
              "RRM\n(BLUP-Kron\u2217)", "Maternal", "ssGWAS",
              "Threshold", "Count PSR", "Ordinal PSR", "Survival", "GxE"]
    for i, name in enumerate(models):
        x = 0.08 + i * 0.88
        y = 1.55 if i < 6 else 0.35
        rbox(x, y, 0.82, 0.5, C_MODEL, name, "", fs=6.5 if "\n" in name else 7)

    # Row 4: Output
    rbox(7.5, 0.95, 2.0, 0.5, C_OUT, "Output",
         ".assoc.txt  |  Dynamic curves  |  Summary")

    # Innovation callouts
    for xy, txt, clr in [
        ((3.4, 2.15), "Woodbury\nReduction", CS_ORANGE),
        ((4.45, 2.15), "BLUP-Kronecker\nRRM", CS_TEAL),
    ]:
        ax.annotate("", xy=(xy[0]+0.1, 2.9), xytext=(xy[0]+0.1, 2.6),
                    arrowprops=dict(arrowstyle="->", color=clr, lw=1.2))
        ax.text(xy[0]+0.1, xy[1], txt, ha="center", va="center",
                fontsize=6.5, fontweight="bold", color=clr,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="#fff9f0", edgecolor=clr, linewidth=0.8))

    # Arrows
    for xs in [0.5, 4.5, 7.5]:
        ax.annotate("", xy=(xs, 4.7), xytext=(xs, 4.25),
                    arrowprops=dict(arrowstyle="->", color=BORDER, lw=0.8))

    ax.set_title("Figure 1. CosmicScanner Architecture and Model Coverage",
                 fontsize=10, fontweight="bold", pad=6, loc="left")
    fig.savefig(FIG_DIR / "figure1_architecture.png")
    fig.savefig(FIG_DIR / "figure1_architecture.pdf")
    plt.close()
    print("✓ Figure 1")


# ═══════════════════════════════════════════════════════════════════
# Figure 2: Woodbury reduction + repeatability performance
# ═══════════════════════════════════════════════════════════════════
def fig2():
    fig, axes = plt.subplots(1, 3, figsize=(8.5, 3.0),
                             gridspec_kw={"width_ratios": [1, 1.1, 1.3]})

    # A: Direct approach
    ax = axes[0]; ax.axis("off"); ax.set_xlim(0, 4); ax.set_ylim(0, 3)
    r = FancyBboxPatch((0.2, 0.8), 1.8, 1.4, boxstyle="round,pad=0.1",
                       facecolor="#fce4ec", edgecolor="#c62828", linewidth=1.2)
    ax.add_patch(r)
    ax.text(1.1, 2.0, "V", ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(1.1, 1.5, r"$n_{\text{records}} \times n_{\text{records}}$",
            ha="center", va="center", fontsize=7.5)
    ax.text(1.1, 0.95, r"$O(n_{\text{records}}^3)$", ha="center", va="center",
            fontsize=8, fontweight="bold", color="#c62828")
    ax.set_title("A. Direct record-level", fontsize=8, fontweight="bold", pad=4)

    # B: Woodbury reduction
    ax = axes[1]; ax.axis("off"); ax.set_xlim(0, 4.5); ax.set_ylim(0, 3)
    rr = FancyBboxPatch((0.15, 0.4), 2.4, 0.65, boxstyle="round,pad=0.08",
                        facecolor="#e8f5e9", edgecolor="#2e7d32", linewidth=1.2)
    ax.add_patch(rr)
    ax.text(1.35, 0.9, r"$\mathbf{M} = \mathbf{K}^{-1} + \mathbf{Z}'\mathbf{R}^{-1}\mathbf{Z}$",
            ha="center", va="center", fontsize=7.5)
    ax.text(1.35, 0.55, r"$n_{\text{subj}} \times n_{\text{subj}}$", ha="center", va="center", fontsize=7)
    ax.text(1.35, 2.5, r"$\mathbf{V}^{-1} = \mathbf{R}^{-1} - \mathbf{R}^{-1}\mathbf{Z}\mathbf{M}^{-1}\mathbf{Z}'\mathbf{R}^{-1}$",
            ha="center", va="center", fontsize=7.5, color="#2e7d32", fontweight="bold")
    ax.annotate("", xy=(1.35, 2.2), xytext=(1.35, 1.05),
                arrowprops=dict(arrowstyle="->", color="#2e7d32", lw=1.5))
    ax.set_title("B. Woodbury subject-level", fontsize=8, fontweight="bold", pad=4)

    # C: Performance bars (using real benchmark data)
    ax = axes[2]
    data = read_tsv(RESULT_DIR / "repeatability_perf_bench" / "benchmark.summary.tsv")
    cs_exact = {f"{r['n_subjects']}r{r['n_records_per_subject']}": r
                for r in data if r["method"] == "exact"}
    labels = ["100/2", "100/5", "200/2", "200/5", "500/2", "500/5"]
    cs_times = [float(cs_exact.get(f"{l.replace('/','r')}", {}).get("runtime_sec", 0)) for l in labels]
    x = np.arange(len(labels))
    ax.bar(x, cs_times, width=0.55, color=CS_BLUE, edgecolor="white", linewidth=0.3)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7)
    ax.set_xlabel("N$_{\\text{subj}}$ / Records per subject")
    ax.set_ylabel("Runtime (s)")
    for i, v in enumerate(cs_times):
        ax.text(i, v + 0.3, f"{v:.1f}s", ha="center", fontsize=6.5, fontweight="bold")
    ax.set_title("C. Woodbury runtime", fontsize=8, fontweight="bold", pad=4)

    fig.suptitle("Figure 2. Woodbury Subject-Level Reduction for Repeatability GWAS",
                 fontsize=10, fontweight="bold", y=1.04)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figure2_woodbury.png")
    fig.savefig(FIG_DIR / "figure2_woodbury.pdf")
    plt.close()
    print("✓ Figure 2")


# ═══════════════════════════════════════════════════════════════════
# Figure 3: BLUP-based Kronecker RRM
# ═══════════════════════════════════════════════════════════════════
def fig3():
    fig, axes = plt.subplots(1, 3, figsize=(8.5, 2.8),
                             gridspec_kw={"width_ratios": [1, 1.3, 1]})

    # A: Observation-level
    ax = axes[0]; ax.axis("off"); ax.set_xlim(0, 4); ax.set_ylim(0, 3)
    r = FancyBboxPatch((0.2, 0.8), 1.8, 1.3, boxstyle="round,pad=0.1",
                       facecolor="#fce4ec", edgecolor="#c62828", linewidth=1.2)
    ax.add_patch(r)
    ax.text(1.1, 1.8, r"$\mathbf{V} = \mathbf{Z}(\mathbf{K}\otimes\mathbf{G})\mathbf{Z}' + \mathbf{R}$",
            ha="center", va="center", fontsize=6.5)
    ax.text(1.1, 1.15, r"$n_{\text{obs}} \times n_{\text{obs}}$", ha="center", va="center", fontsize=7.5)
    ax.text(1.1, 2.5, "19,392  observations", ha="center", va="center",
            fontsize=7.5, color="#c62828", fontweight="bold")
    ax.set_title("A. Observation-level", fontsize=8, fontweight="bold", pad=4)

    # B: BLUP-based Kronecker
    ax = axes[1]; ax.axis("off"); ax.set_xlim(0, 5); ax.set_ylim(0, 3)
    b1 = FancyBboxPatch((0.15, 1.7), 2.2, 0.5, boxstyle="round,pad=0.08",
                        facecolor="#d4e6f1", edgecolor="#1b4f72", linewidth=1)
    ax.add_patch(b1); ax.text(1.25, 1.95, "Step 1: Individual BLUP", ha="center", fontsize=7, fontweight="bold")
    ax.text(1.25, 1.75, r"$\hat{\mathbf{a}}_i = \text{BLUP}(\boldsymbol{\phi}(t_{ij}), y_{ij})$", ha="center", fontsize=6.5)

    b2 = FancyBboxPatch((0.15, 0.4), 2.2, 0.5, boxstyle="round,pad=0.08",
                        facecolor="#e8f5e9", edgecolor="#2e7d32", linewidth=1)
    ax.add_patch(b2); ax.text(1.25, 0.65, "Step 2: Kronecker mvLMM", ha="center", fontsize=7, fontweight="bold")
    ax.text(1.25, 0.45, r"$\text{vec}(\hat{\mathbf{A}}) \sim N(\mathbf{0}, \mathbf{K}_0\otimes\mathbf{G} + \mathbf{P}_0\otimes\mathbf{I})$", ha="center", fontsize=5.5)

    ax.annotate("", xy=(1.25, 1.7), xytext=(1.25, 0.9),
                arrowprops=dict(arrowstyle="->", color="#2e7d32", lw=1.5))

    ax.text(2.8, 1.1, "1,212 subjects\n213x faster\n93x less memory",
            ha="center", va="center", fontsize=8, fontweight="bold", color="#2e7d32",
            bbox=dict(boxstyle="round", facecolor="#e8f5e9", edgecolor="#2e7d32", linewidth=0.8))

    ax.set_title("B. BLUP-based Kronecker", fontsize=8, fontweight="bold", pad=4)

    # C: Memory comparison
    ax = axes[2]; ax.axis("off")
    ax.barh([0], [195], height=0.35, color=CS_BLUE, edgecolor="white")
    ax.barh([1], [17700], height=0.35, color=CS_ORANGE, edgecolor="white")
    ax.set_yticks([0, 1]); ax.set_yticklabels(["CS", "Observation"], fontsize=7)
    ax.set_xlim(0, 20000)
    ax.text(195+200, 0, "195 MB", va="center", fontsize=7.5, color=CS_BLUE, fontweight="bold")
    ax.text(17700-800, 1, "17,700 MB", va="center", fontsize=7.5, color=CS_ORANGE, fontweight="bold")
    ax.text(5000, -0.8, "93x reduction", ha="center", fontsize=8, fontweight="bold", color="#2e7d32")
    ax.set_title("C. Memory comparison", fontsize=8, fontweight="bold", pad=4)

    fig.suptitle("Figure 3. BLUP-Based Kronecker RRM for Longitudinal GWAS",
                 fontsize=10, fontweight="bold", y=1.04)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figure3_blup_rrm.png")
    fig.savefig(FIG_DIR / "figure3_blup_rrm.pdf")
    plt.close()
    print("✓ Figure 3")


# ═══════════════════════════════════════════════════════════════════
# Figure 4: Accuracy scatter plots (4 panels)
# ═══════════════════════════════════════════════════════════════════
def fig4():
    fig, axes = plt.subplots(2, 2, figsize=(7, 6.5))

    # ── Panel A: Single-trait vs GEMMA ──
    ax = axes[0, 0]
    data = read_tsv(RESULT_DIR / "gemma_accuracy" / "gemma_accuracy.compare.tsv")
    cs_p = [float(r["log10p_cosmic"]) for r in data]
    gm_p = [float(r["log10p_gemma"]) for r in data]
    ax.scatter(gm_p, cs_p, s=2, alpha=0.35, color=CS_BLUE, edgecolors="none")
    lims = [min(gm_p + cs_p) - 0.1, max(gm_p + cs_p) + 0.1]
    ax.plot(lims, lims, "k-", lw=0.6, alpha=0.4)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("GEMMA $-$log$_{10}$(P)"); ax.set_ylabel("CosmicScanner $-$log$_{10}$(P)")
    r_val = np.corrcoef(cs_p, gm_p)[0, 1]
    ax.text(0.05, 0.92, f"r = {r_val:.4f}", transform=ax.transAxes, fontsize=8,
            fontweight="bold", bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.8))
    ax.set_title("A. Single-trait LMM (vs GEMMA)", fontsize=8.5, fontweight="bold")

    # ── Panel B: Multivariate vs EMmvGWAS ──
    ax = axes[0, 1]
    try:
        mv_dir = RESULT_DIR / "multivariate_system_benchmark_n240_m1200_t2"
        cs_fast = read_tsv(mv_dir / "cs_fast.assoc.mv.txt")
        emmv = read_csv(mv_dir / "emmvgwas" / "emmv_association_with_ids.csv")
        cs_p_map = {r["SNP"]: -math.log10(max(float(r["P"]), 1e-300)) for r in cs_fast}
        emmv_p_map = {r["SNP_ID"]: -math.log10(max(float(r["p"]), 1e-300)) for r in emmv}
        common = sorted(set(cs_p_map) & set(emmv_p_map))
        if len(common) > 10:
            xs_v = [emmv_p_map[s] for s in common]
            ys_v = [cs_p_map[s] for s in common]
            ax.scatter(xs_v, ys_v, s=2, alpha=0.35, color=CS_TEAL, edgecolors="none")
            lims = [min(xs_v + ys_v) - 0.1, max(xs_v + ys_v) + 0.1]
            ax.plot(lims, lims, "k-", lw=0.6, alpha=0.4)
            ax.set_xlim(lims); ax.set_ylim(lims)
            r_v = np.corrcoef(ys_v, xs_v)[0, 1]
            ax.text(0.05, 0.92, f"r = {r_v:.4f}", transform=ax.transAxes, fontsize=8,
                    fontweight="bold", bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.8))
    except Exception as e:
        ax.text(0.5, 0.5, f"Data not available", ha="center", va="center", fontsize=8)
    ax.set_xlabel("EMmvGWAS $-$log$_{10}$(P)"); ax.set_ylabel("CosmicScanner $-$log$_{10}$(P)")
    ax.set_title("B. Multivariate LMM (vs EMmvGWAS)", fontsize=8.5, fontweight="bold")

    # ── Panel C: Repeatability benchmark ──
    ax = axes[1, 0]
    rp_data = read_tsv(RESULT_DIR / "repeatability_perf_bench" / "benchmark.summary.tsv")
    configs = [(100, 2), (100, 5), (200, 2), (200, 5), (500, 2), (500, 5)]
    labels = ["100/2", "100/5", "200/2", "200/5", "500/2", "500/5"]
    exact_times = []
    for n, r in configs:
        match = [x for x in rp_data if int(x["n_subjects"]) == n
                 and int(x["n_records_per_subject"]) == r and x["method"] == "exact"]
        if match:
            exact_times.append(float(match[-1]["runtime_sec"]))
    if exact_times:
        x = np.arange(len(labels))
        ax.bar(x, exact_times, width=0.55, color=CS_BLUE, edgecolor="white")
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7)
        ax.set_xlabel("N$_{\\text{subj}}$ / Records")
        ax.set_ylabel("Runtime (s)")
        for i, v in enumerate(exact_times):
            ax.text(i, v + 0.3, f"{v:.1f}s", ha="center", fontsize=6.5, va="bottom")
        ax.set_title("C. Repeatability runtime", fontsize=8.5, fontweight="bold")

    # ── Panel D: RRM dynamic curves ──
    ax = axes[1, 1]
    try:
        curve_data = read_tsv(RESULT_DIR / "rrm_dynamic_simulation"
                              / "rrm_dynamic_sim.rrm.dynamic_curve.txt")
        causal_snps = ["rs_dyn_causal", "rs_stable_causal", "rs_early_causal",
                       "rs_late_causal", "rs_transient_causal"]
        colors_d = {"rs_dyn_causal": "#d6604d", "rs_stable_causal": "#2166ac",
                    "rs_early_causal": "#35978f", "rs_late_causal": "#d95f02",
                    "rs_transient_causal": "#762a83"}
        for snp in causal_snps:
            snp_rows = [r for r in curve_data if r["SNP"] == snp]
            if snp_rows:
                times = [float(r["Time"]) for r in snp_rows]
                betas = [float(r["Beta"]) for r in snp_rows]
                ax.plot(times, betas, "o-", color=colors_d.get(snp, "#666"),
                        linewidth=1.5, markersize=3, label=snp.split("_causal")[0].replace("_", ""))
        ax.axhline(0, color="gray", linestyle="--", lw=0.5)
        ax.set_xlabel("Time"); ax.set_ylabel("SNP Effect ($\\beta$)")
        ax.legend(fontsize=6, ncol=2, framealpha=0.8)
        ax.set_title("D. RRM dynamic effects", fontsize=8.5, fontweight="bold")
    except Exception as e:
        ax.text(0.5, 0.5, "Data not available", ha="center", va="center", fontsize=8)

    fig.suptitle("Figure 4. Validation and Performance of CosmicScanner",
                 fontsize=10, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figure4_accuracy.png")
    fig.savefig(FIG_DIR / "figure4_accuracy.pdf")
    plt.close()
    print("✓ Figure 4")


# ═══════════════════════════════════════════════════════════════════
# Figure 5: Multivariate benchmark
# ═══════════════════════════════════════════════════════════════════
def fig5():
    mv_dir = RESULT_DIR / "multivariate_system_benchmark_n240_m1200_t2"
    data = read_tsv(mv_dir / "method_metrics.tsv")
    methods = [r["method"] for r in data]
    seconds = [float(r["seconds"]) for r in data]
    rss = [float(r["max_rss_mb"]) for r in data]
    corr = [float(r["log10p_corr_vs_gemma"]) if r["log10p_corr_vs_gemma"] != "NA" else float("nan") for r in data]
    emmv_corr = [float(r["log10p_corr_vs_emmvgwas"]) if r["log10p_corr_vs_emmvgwas"] != "NA" else float("nan") for r in data]
    top10 = [int(r["causal_in_top10"]) for r in data]
    top50 = [int(r["causal_in_top50"]) for r in data]

    fig, axes = plt.subplots(2, 2, figsize=(7, 5.5))

    # A: Runtime
    ax = axes[0, 0]
    colors = [CS_BLUE, CS_TEAL, CS_ORANGE, CS_PURPLE, CS_GRAY]
    bars = ax.bar(methods, seconds, color=colors[:len(methods)], edgecolor="white", linewidth=0.3)
    ax.set_ylabel("Runtime (s, log scale)")
    ax.set_yscale("log")
    for bar, v in zip(bars, seconds):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()*1.05,
                f"{v:.1f}s", ha="center", va="bottom", fontsize=7, fontweight="bold")
    ax.set_title("A. Runtime comparison", fontsize=8.5, fontweight="bold")

    # B: Memory
    ax = axes[0, 1]
    ax.bar(methods, rss, color=colors[:len(methods)], edgecolor="white", linewidth=0.3)
    ax.set_ylabel("Peak memory (MB)")
    ax.set_xticklabels(methods, fontsize=6.5, rotation=20, ha="right")
    for i, v in enumerate(rss):
        ax.text(i, v*1.02, f"{v:.0f}", ha="center", fontsize=7)
    ax.set_title("B. Memory usage", fontsize=8.5, fontweight="bold")

    # C: Accuracy
    ax = axes[1, 0]
    x = np.arange(len(methods))
    w = 0.3
    valid_c = [(i, c) for i, c in enumerate(corr) if not math.isnan(c)]
    if valid_c:
        ax.bar([v[0] for v in valid_c], [v[1] for v in valid_c], width=w,
               color=CS_TEAL, edgecolor="white", label="vs GEMMA")
    valid_e = [(i, c) for i, c in enumerate(emmv_corr) if not math.isnan(c)]
    if valid_e:
        ax.bar([v[0]+w for v in valid_e], [v[1] for v in valid_e], width=w,
               color=CS_ORANGE, edgecolor="white", label="vs EMmvGWAS")
    ax.set_xticks(x + w/2); ax.set_xticklabels(methods, fontsize=6.5, rotation=20, ha="right")
    ax.set_ylabel("P correlation")
    ax.legend(fontsize=7, framealpha=0.8)
    ax.set_title("C. Accuracy concordance", fontsize=8.5, fontweight="bold")

    # D: Causal recovery
    ax = axes[1, 1]
    x = np.arange(len(methods))
    w2 = 0.3
    ax.bar(x - w2/2, top10, w2, color=CS_ORANGE, edgecolor="white", label="Top 10")
    ax.bar(x + w2/2, top50, w2, color=CS_TEAL, edgecolor="white", label="Top 50")
    ax.set_xticks(x); ax.set_xticklabels(methods, fontsize=6.5, rotation=20, ha="right")
    ax.set_ylabel("Causal SNPs (of 4)")
    ax.set_ylim(0, 5.5)
    ax.legend(fontsize=7, framealpha=0.8)
    ax.set_title("D. Causal SNP recovery", fontsize=8.5, fontweight="bold")

    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Figure 5. Multivariate GWAS System Benchmark (N=240, M=1,200, 2 traits)",
                 fontsize=10, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figure5_multivariate_benchmark.png")
    fig.savefig(FIG_DIR / "figure5_multivariate_benchmark.pdf")
    plt.close()
    print("✓ Figure 5")


# ═══════════════════════════════════════════════════════════════════
# Figure 6: Repeatability scaling
# ═══════════════════════════════════════════════════════════════════
def fig6():
    data = read_tsv(RESULT_DIR / "repeatability_perf_bench" / "benchmark.summary.tsv")
    cs_exact = [r for r in data if r["method"] == "exact"]

    fig, axes = plt.subplots(1, 3, figsize=(8.5, 3.0))

    # A: Runtime vs subjects (bubble)
    ax = axes[0]
    n_subj = sorted(set(int(r["n_subjects"]) for r in cs_exact))
    for r_val, marker, color in [(2, "o", CS_BLUE), (5, "s", CS_TEAL)]:
        pts = [(int(r["n_subjects"]), float(r["runtime_sec"]))
               for r in cs_exact if int(r["n_records_per_subject"]) == r_val]
        if pts:
            pts.sort()
            ax.plot([p[0] for p in pts], [p[1] for p in pts],
                    marker=marker, color=color, label=f"{r_val} rec/subj", linewidth=1.2)
    ax.set_xlabel("Number of subjects"); ax.set_ylabel("Runtime (s)")
    ax.legend(fontsize=7, framealpha=0.8)
    ax.set_title("A. Runtime vs subjects", fontsize=8.5, fontweight="bold")

    # B: Runtime vs records/subject
    ax = axes[1]
    for n in sorted(set(int(r["n_subjects"]) for r in cs_exact)):
        pts = [(int(r["n_records_per_subject"]), float(r["runtime_sec"]))
               for r in cs_exact if int(r["n_subjects"]) == n]
        if pts:
            pts.sort()
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-",
                    label=f"N={n}", linewidth=1.2)
    ax.set_xlabel("Records per subject"); ax.set_ylabel("Runtime (s)")
    ax.legend(fontsize=7, framealpha=0.8)
    ax.set_title("B. Runtime vs records", fontsize=8.5, fontweight="bold")

    # C: Memory
    ax = axes[2]
    for n in sorted(set(int(r["n_subjects"]) for r in cs_exact)):
        pts = [(int(r["n_records_per_subject"]), float(r["max_rss_mb"]))
               for r in cs_exact if int(r["n_subjects"]) == n]
        if pts:
            pts.sort()
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-",
                    label=f"N={n}", linewidth=1.2)
    ax.set_xlabel("Records per subject"); ax.set_ylabel("Peak memory (MB)")
    ax.legend(fontsize=7, framealpha=0.8)
    ax.set_title("C. Memory scaling", fontsize=8.5, fontweight="bold")

    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Figure 6. Repeatability GWAS Scaling Behavior",
                 fontsize=10, fontweight="bold", y=1.04)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figure6_repeatability_scaling.png")
    fig.savefig(FIG_DIR / "figure6_repeatability_scaling.pdf")
    plt.close()
    print("✓ Figure 6")


# ═══════════════════════════════════════════════════════════════════
# Figure 7: RRM dynamic effect curves
# ═══════════════════════════════════════════════════════════════════
def fig7():
    curve_path = RESULT_DIR / "rrm_dynamic_simulation" / "rrm_dynamic_sim.rrm.dynamic_curve.txt"
    curve_data = read_tsv(curve_path)

    causal_info = {
        "rs_dyn_causal": ("Sign-switch", "#d6604d"),
        "rs_stable_causal": ("Stable", "#2166ac"),
        "rs_early_causal": ("Early-specific", "#35978f"),
        "rs_late_causal": ("Late-specific", "#d95f02"),
        "rs_transient_causal": ("Transient", "#762a83"),
    }

    fig, axes = plt.subplots(2, 3, figsize=(8, 5))
    axes_flat = axes.flat

    for idx, (snp, (label, color)) in enumerate(causal_info.items()):
        ax = axes_flat[idx]
        snp_rows = [r for r in curve_data if r["SNP"] == snp]
        if not snp_rows:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            continue
        times = [float(r["Time"]) for r in snp_rows]
        betas = [float(r["Beta"]) for r in snp_rows]
        ses = [float(r["SE"]) for r in snp_rows]

        ax.fill_between(times,
                        [b - 1.96*s for b, s in zip(betas, ses)],
                        [b + 1.96*s for b, s in zip(betas, ses)],
                        alpha=0.15, color=color)
        ax.plot(times, betas, "o-", color=color, linewidth=2, markersize=5)
        ax.axhline(0, color="gray", linestyle="--", lw=0.6)
        ax.set_xlabel("Time", fontsize=7.5)
        ax.set_ylabel("Effect ($\\beta$)", fontsize=7.5)
        ax.set_title(label, fontsize=8, fontweight="bold", color=color)
        ax.spines[["top", "right"]].set_visible(False)

    axes_flat[5].axis("off")

    fig.suptitle("Figure 7. RRM Dynamic Effect Curves for Five QTL Trajectory Types",
                 fontsize=10, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figure7_rrm_dynamic.png")
    fig.savefig(FIG_DIR / "figure7_rrm_dynamic.pdf")
    plt.close()
    print("✓ Figure 7")


if __name__ == "__main__":
    print("Generating publication-quality figures...")
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    fig7()
    print(f"\nAll figures saved to {FIG_DIR}")
