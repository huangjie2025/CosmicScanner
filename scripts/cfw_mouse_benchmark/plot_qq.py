#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plot_consistency import load_assoc

ap = argparse.ArgumentParser()
ap.add_argument("--inputs", nargs="*", default=[])
ap.add_argument("--labels", nargs="*", default=[])
ap.add_argument("--out", required=True)
ap.add_argument("--summary", required=True)
args = ap.parse_args()
fig, ax = plt.subplots(figsize=(5, 5))
rows = []
for label, path in zip(args.labels, args.inputs):
    df = load_assoc(path)
    p = pd.to_numeric(df["P"], errors="coerce").dropna().clip(lower=1e-300, upper=1).sort_values().to_numpy()
    if len(p) == 0:
        rows.append({"method": label, "n": 0, "lambda_gc": np.nan})
        continue
    exp = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
    obs = -np.log10(p)
    ax.plot(exp, obs, marker=".", lw=0, ms=3, label=label)
    rows.append({"method": label, "n": len(p), "lambda_gc": np.nan})
mx = max(ax.get_xlim()[1], ax.get_ylim()[1], 1)
ax.plot([0, mx], [0, mx], color="black", lw=1)
ax.set_xlabel("Expected -log10(p)")
ax.set_ylabel("Observed -log10(p)")
ax.legend(frameon=False)
plt.tight_layout()
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(args.out, dpi=200)
pd.DataFrame(rows).to_csv(args.summary, sep="\t", index=False)
