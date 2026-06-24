#!/bin/bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

ROOT="/disk228/huangjie/PanOmics/CosmicScanner"
BIN="${ROOT}/build/cosmicscanner"
DATA_DIR="${ROOT}/tests/data/repeatability_exact_vs_fast_sim"
RESULT_DIR="${ROOT}/tests/results/repeatability_exact_vs_fast_sim"

SIM_N=${SIM_N:-80}
SIM_M=${SIM_M:-240}
SIM_REPS=${SIM_REPS:-3}
SIM_ALPHA=${SIM_ALPHA:-1e-3}

mkdir -p "${DATA_DIR}" "${RESULT_DIR}"

python3 - <<PY
from pathlib import Path
import math
import numpy as np

root = Path("${ROOT}")
data_dir = Path("${DATA_DIR}")
result_dir = Path("${RESULT_DIR}")
n = int("${SIM_N}")
m = int("${SIM_M}")
reps = int("${SIM_REPS}")
rng = np.random.default_rng(20260608)

prefix = data_dir / "sim"
causal_idx = [5, 37, 96]
causal_idx = [i for i in causal_idx if i < m]
if not causal_idx:
    causal_idx = [min(5, m - 1)]

maf = rng.uniform(0.08, 0.45, size=m)
G = rng.binomial(2, maf, size=(n, m)).astype(float)
missing = rng.random((n, m)) < 0.01
G[missing] = np.nan
G_imp = G.copy()
for j in range(m):
    col = G_imp[:, j]
    mean = np.nanmean(col)
    if not np.isfinite(mean):
        mean = 2 * maf[j]
    col[np.isnan(col)] = mean
    G_imp[:, j] = col

with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for i in range(n):
        fid = iid = f"ID{i+1:04d}"
        f.write(f"{fid}\t{iid}\t0\t0\t0\t-9\n")

with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for j in range(m):
        chrom = 1 + (j % 2)
        f.write(f"{chrom}\trs{j+1:06d}\t0\t{1000 + 10*j}\tA\tG\n")

bed = bytearray(b"\x6c\x1b\x01")
bytes_per_snp = (n + 3) // 4
for j in range(m):
    block = bytearray(bytes_per_snp)
    for i in range(n):
        val = G[i, j]
        if np.isnan(val):
            code = 0b01
        elif val <= 0.5:
            code = 0b00
        elif val <= 1.5:
            code = 0b10
        else:
            code = 0b11
        block[i // 4] |= code << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix(".bed").write_bytes(bed)

np.savetxt(data_dir / "causal_snps.txt", [f"rs{i+1:06d}" for i in causal_idx], fmt="%s")

scenarios = {
    "balanced": dict(counts=lambda i: [1, 2, 3], year_scale=0.0, confounded=False),
    "unbalanced": dict(counts=lambda i: list(range(1, 2 + (i % 5))), year_scale=0.0, confounded=False),
    "unbalanced_year": dict(counts=lambda i: list(range(1, 2 + (i % 5))), year_scale=1.6, confounded=True),
}

for rep in range(reps):
    additive = np.zeros(n)
    effects = [0.65, -0.55, 0.50]
    for eff, j in zip(effects, causal_idx):
        additive += eff * (G_imp[:, j] - np.mean(G_imp[:, j]))
    polygenic = rng.normal(0, 0.45, size=n)
    pe = rng.normal(0, 0.45, size=n)
    for name, cfg in scenarios.items():
        pheno = data_dir / f"{name}.rep{rep+1}.pheno"
        qcovar = data_dir / f"{name}.rep{rep+1}.qcovar"
        mean_pheno = data_dir / f"{name}.rep{rep+1}.mean.pheno"
        rows = []
        with pheno.open("w", encoding="utf-8") as pf, qcovar.open("w", encoding="utf-8") as cf:
            for i in range(n):
                fid = iid = f"ID{i+1:04d}"
                subj_rows = []
                for t in cfg["counts"](i):
                    if cfg["confounded"]:
                        year = 1.0 if (G_imp[i, causal_idx[0]] > np.nanmedian(G_imp[:, causal_idx[0]]) and t >= 2) else 0.0
                    else:
                        year = 1.0 if t % 2 == 0 else 0.0
                    noise = rng.normal(0, 0.70)
                    y = additive[i] + polygenic[i] + pe[i] + cfg["year_scale"] * year + noise
                    pf.write(f"{fid}\t{iid}\t{t}\t{y:.8f}\n")
                    cf.write(f"{fid}\t{iid}\t{t}\t{year:.1f}\n")
                    subj_rows.append(y)
                rows.append((fid, iid, float(np.mean(subj_rows))))
        with mean_pheno.open("w", encoding="utf-8") as mf:
            for fid, iid, y in rows:
                mf.write(f"{fid}\t{iid}\t{y:.8f}\n")
PY

run_repeat_exact() {
  local scenario=$1 rep=$2 out=$3
  "${BIN}" --bfile "${DATA_DIR}/sim" --pheno "${DATA_DIR}/${scenario}.rep${rep}.pheno" \
    --tpheno 1 --mpheno 2 --qcovar "${DATA_DIR}/${scenario}.rep${rep}.qcovar" \
    --repeatability-model --repeatability-mode exact --gwas --no-loco --maf 0.01 \
    --out "${out}" >/"${out}.log" 2>&1
}

run_repeat_fast() {
  local scenario=$1 rep=$2 out=$3
  "${BIN}" --bfile "${DATA_DIR}/sim" --pheno "${DATA_DIR}/${scenario}.rep${rep}.pheno" \
    --tpheno 1 --mpheno 2 --qcovar "${DATA_DIR}/${scenario}.rep${rep}.qcovar" \
    --repeatability-model --repeatability-mode fast --gwas --no-loco --maf 0.01 \
    --out "${out}" >/"${out}.log" 2>&1
}

run_mean_lmm() {
  local scenario=$1 rep=$2 out=$3
  "${BIN}" --bfile "${DATA_DIR}/sim" --pheno "${DATA_DIR}/${scenario}.rep${rep}.mean.pheno" \
    --mpheno 1 --gwas --no-loco --maf 0.01 --out "${out}" >/"${out}.log" 2>&1
}

for rep in $(seq 1 "${SIM_REPS}"); do
  for scenario in balanced unbalanced unbalanced_year; do
    run_repeat_exact "${scenario}" "${rep}" "${RESULT_DIR}/${scenario}.rep${rep}.exact"
    run_repeat_fast "${scenario}" "${rep}" "${RESULT_DIR}/${scenario}.rep${rep}.fast"
    run_mean_lmm "${scenario}" "${rep}" "${RESULT_DIR}/${scenario}.rep${rep}.mean"
  done
done

python3 - <<PY
from pathlib import Path
import csv, math, statistics

data_dir = Path("${DATA_DIR}")
result_dir = Path("${RESULT_DIR}")
reps = int("${SIM_REPS}")
alpha = float("${SIM_ALPHA}")
causal = set((data_dir / "causal_snps.txt").read_text(encoding="utf-8").split())

def read_assoc(path):
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def metrics(path):
    rows = read_assoc(path)
    pvals = []
    null_p = []
    causal_hit = 0
    top = None
    top_p = 2.0
    for r in rows:
        try:
            p = float(r["P"])
        except Exception:
            continue
        if not math.isfinite(p) or p <= 0:
            continue
        pvals.append(p)
        snp = r["SNP"]
        if p < top_p:
            top_p = p
            top = snp
        if snp in causal:
            if p <= alpha:
                causal_hit += 1
        else:
            null_p.append(p)
    if not pvals:
        raise SystemExit(f"No valid p-values in {path}")
    false_pos = sum(1 for p in null_p if p <= alpha)
    null_rate = false_pos / max(1, len(null_p))
    lam = None
    if null_p:
        chisq = []
        for p in null_p:
            # Wilson-Hilferty approximation to chi-square(1) inverse survival.
            from statistics import NormalDist
            z = NormalDist().inv_cdf(1 - p / 2)
            chisq.append(z * z)
        lam = statistics.median(chisq) / 0.4549364231
    return dict(
        n_snps=len(pvals), causal_hits=causal_hit, false_pos=false_pos,
        null_rate=null_rate, lambda_gc=lam if lam is not None else float("nan"),
        top_snp=top or "NA", top_p=top_p,
    )

summary_path = result_dir / "simulation.summary.tsv"
with summary_path.open("w", encoding="utf-8", newline="") as f:
    fields = ["scenario", "rep", "method", "n_snps", "causal_hits", "false_pos", "null_rate", "lambda_gc", "top_snp", "top_p"]
    w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
    w.writeheader()
    for scenario in ["balanced", "unbalanced", "unbalanced_year"]:
        for rep in range(1, reps + 1):
            for method in ["exact", "fast", "mean"]:
                row = metrics(result_dir / f"{scenario}.rep{rep}.{method}.assoc.txt")
                row.update(scenario=scenario, rep=rep, method=method)
                w.writerow(row)

print(f"Wrote {summary_path}")
PY

cat "${RESULT_DIR}/simulation.summary.tsv"
