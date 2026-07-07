#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plot_consistency import load_assoc

ap = argparse.ArgumentParser()
ap.add_argument("--bim", required=True)
ap.add_argument("--gemma")
ap.add_argument("--emmv")
ap.add_argument("--cosmic")
ap.add_argument("--out", required=True)
args = ap.parse_args()
bim = pd.read_csv(args.bim, sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
methods = [("GEMMA", args.gemma), ("EMmvGWAS", args.emmv), ("CosmicScanner", args.cosmic)]
fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
threshold = -np.log10(0.05 / max(1, len(bim)))
for ax, (name, path) in zip(axes, methods):
    df = load_assoc(path) if path else pd.DataFrame()
    if df.empty:
        ax.text(0.5, 0.5, "missing results", ha="center", va="center")
        continue
    m = bim.merge(df, on="SNP", how="inner")
    y = -np.log10(pd.to_numeric(m["P"], errors="coerce").clip(lower=1e-300))
    x = np.arange(len(m))
    colors = ["#2f6fbb" if int(c) % 2 else "#111111" for c in m["CHR"]]
    ax.scatter(x, y, c=colors, s=3, alpha=0.7)
    ax.axhline(threshold, color="red", lw=1)
    ax.set_ylabel(name)
axes[-1].set_xlabel("SNP order")
plt.tight_layout()
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(args.out, dpi=200)
