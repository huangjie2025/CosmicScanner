#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--runtime-tsv", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()
df = pd.read_csv(args.runtime_tsv, sep="\t") if Path(args.runtime_tsv).exists() else pd.DataFrame()
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
plt.figure(figsize=(6, 4))
if not df.empty:
    df["runtime_hours"] = pd.to_numeric(df["seconds"], errors="coerce") / 3600.0
    for method, g in df.groupby("method"):
        g = g.sort_values("traits")
        plt.plot(g["traits"], g["runtime_hours"], marker="o", label=method)
plt.xlabel("Number of traits")
plt.ylabel("Runtime hours")
plt.xticks([2, 4, 6, 8, 10])
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(args.out, dpi=200)
