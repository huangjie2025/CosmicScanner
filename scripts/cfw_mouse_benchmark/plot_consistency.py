#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def load_assoc(path):
    p = Path(path)
    if not p.exists():
        return pd.DataFrame(columns=["SNP", "P"])
    sep = "," if p.suffix == ".csv" else "\t"
    df = pd.read_csv(p, sep=sep)
    snp_col = next((c for c in ["SNP", "SNP_ID", "rs", "snp"] if c in df.columns), df.columns[0])
    p_col = next((c for c in ["P", "p_wald", "p_lrt", "p_score", "p"] if c in df.columns), None)
    if p_col is None:
        p_col = next(c for c in df.columns if c.lower().startswith("p"))
    return df[[snp_col, p_col]].rename(columns={snp_col: "SNP", p_col: "P"})

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gemma")
    ap.add_argument("--emmv")
    ap.add_argument("--cosmic")
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()
    base = load_assoc(args.gemma) if args.gemma else pd.DataFrame()
    comparisons = [("EMmvGWAS", load_assoc(args.emmv) if args.emmv else pd.DataFrame()), ("CosmicScanner", load_assoc(args.cosmic) if args.cosmic else pd.DataFrame())]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    rows = []
    for ax, (name, df) in zip(axes, comparisons):
        if base.empty or df.empty:
            ax.text(0.5, 0.5, "missing results", ha="center", va="center")
            ax.set_axis_off()
            rows.append({"comparison": f"{name} vs GEMMA", "n": 0, "pearson_r": np.nan})
            continue
        m = base.merge(df, on="SNP", suffixes=("_x", "_y"))
        x = -np.log10(pd.to_numeric(m["P_x"], errors="coerce").clip(lower=1e-300))
        y = -np.log10(pd.to_numeric(m["P_y"], errors="coerce").clip(lower=1e-300))
        ok = np.isfinite(x) & np.isfinite(y)
        r = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 2 else np.nan
        ax.scatter(x[ok], y[ok], s=4, alpha=0.35)
        ax.set_xlabel("GEMMA -log10(p)")
        ax.set_ylabel(f"{name} -log10(p)")
        ax.set_title(f"r={r:.4f}")
        rows.append({"comparison": f"{name} vs GEMMA", "n": int(ok.sum()), "pearson_r": r})
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    pd.DataFrame(rows).to_csv(args.summary, sep="\t", index=False)

if __name__ == "__main__":
    main()
