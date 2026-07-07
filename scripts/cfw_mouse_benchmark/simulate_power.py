#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--template-pheno", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--pve", default="0.005,0.01,0.015")
ap.add_argument("--ratios", default="0.2,0.4")
ap.add_argument("--seed", type=int, default=20260627)
args = ap.parse_args()
rng = np.random.default_rng(args.seed)
ph = pd.read_csv(args.template_pheno, sep="\t")
traits = list(ph.columns[2:12])
Path(args.out).mkdir(parents=True, exist_ok=True)
rows = []
for pve in map(float, args.pve.split(",")):
    for ratio in map(float, args.ratios.split(",")):
        for sign in ["+", "-"]:
            y = ph.copy()
            effect = np.sqrt(pve) * (1 if sign == "+" else -1)
            n_causal_traits = max(1, int(round(len(traits) * ratio)))
            for t in traits[:n_causal_traits]:
                y[t] = pd.to_numeric(y[t], errors="coerce") + rng.normal(0, effect, len(y))
            name = f"alt_pve{pve:g}_ratio{ratio:g}_sign{sign.replace('+','pos').replace('-','neg')}.pheno.tsv"
            y.to_csv(Path(args.out) / name, sep="\t", index=False, na_rep="NA")
            rows.append({"file": name, "pve": pve, "ratio": ratio, "sign": sign, "causal_traits": n_causal_traits})
pd.DataFrame(rows).to_csv(Path(args.out) / "power_design.tsv", sep="\t", index=False)
