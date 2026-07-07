#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--template-pheno", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--traits", type=int, default=10)
ap.add_argument("--replicates", type=int, default=10)
ap.add_argument("--seed", type=int, default=20260627)
args = ap.parse_args()
rng = np.random.default_rng(args.seed)
ph = pd.read_csv(args.template_pheno, sep="\t")
ids = ph[["FID", "IID"]]
Y = ph.iloc[:, 2:2 + args.traits].apply(pd.to_numeric, errors="coerce")
cov = np.cov(Y.fillna(Y.mean()).to_numpy().T)
Path(args.out).mkdir(parents=True, exist_ok=True)
for r in range(args.replicates):
    sim = rng.multivariate_normal(np.zeros(args.traits), cov, size=len(ids))
    out = ids.copy()
    for i in range(args.traits):
        out[f"null_trait_{i+1}"] = sim[:, i]
    out.to_csv(Path(args.out) / f"null_{r+1}.pheno.tsv", sep="\t", index=False)
