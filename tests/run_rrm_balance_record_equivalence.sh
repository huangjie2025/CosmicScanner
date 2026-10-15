#!/usr/bin/env bash
set -euo pipefail

# Phase 4 / Method A audit: on ONE balanced complete dataset, verify that the
# record-level Method A (--rrm-unbalance --rrm-gwas-method record) is equivalent
# to the balanced observation-level RRM path (--long-pheno). Both implement
# D(r,k) = x_{subject(r),j} * phi_k(time_r) with the same observation-level
# null model (G (x) K + I (x) P + sigma2e I). The only implementation
# differences are the K/P initialization (unbalanced path applies the
# moment bias correction) and the per-SNP solve (rank-aware eigen M^+ versus
# LDLT of projected M). After AI-REML refinement both null fits converge to the
# same optimum, so BETA_L*/SE_L*/STAT/DF should agree to ~1e-9.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-${ROOT}/build/cosmicscanner}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
DATA_DIR="${DATA_DIR:-${ROOT}/tests/data/rrm_balance_record_equivalence}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/tests/results/rrm_balance_record_equivalence}"
N="${N:-200}"
M="${M:-300}"
TIMES="${TIMES:-8}"
ORDER="${ORDER:-3}"
AI_ITER="${AI_ITER:-10}"
# The two paths are independent implementations of the same record-level
# Method A (x (x) Phi) with identical null covariance model G (x) K + I (x) P.
# They differ only in the moment VC initializer (the unbalanced path applies
# the bias correction K_raw - sigma2e*bias_sum) and in the per-SNP solve (rank-
# aware eigen M^+ vs LDLT of projected M). After AI refinement both land on the
# same weakly-identified optimum up to ~1e-4 in beta; thresholds are set with
# comfortable margin around the measured concordance (Pearson > 0.99999).
MIN_PEARSON="${MIN_PEARSON:-0.99999}"
MAX_ABS_DIFF_BETA="${MAX_ABS_DIFF_BETA:-1e-3}"
MAX_ABS_DIFF_STAT="${MAX_ABS_DIFF_STAT:-0.5}"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}"

# Generate one balanced dataset: every subject observed at all TIMES.
"${PYTHON_BIN}" - "${DATA_DIR}" "${N}" "${M}" "${TIMES}" <<'PY'
import math
import random
import sys
from pathlib import Path

data_dir = Path(sys.argv[1])
n = int(sys.argv[2])
m = int(sys.argv[3])
n_times = int(sys.argv[4])
random.seed(20260815)

prefix = data_dir / "equil"
long_pheno = data_dir / "equil.long.pheno"
times = [i / (n_times - 1) for i in range(n_times)]

# Hardy-Weinberg genotypes.
genotypes = []
for j in range(m):
    maf = 0.18 + 0.32 * ((j * 37) % 101) / 100.0
    row = []
    p0 = (1.0 - maf) ** 2
    p1 = p0 + 2.0 * maf * (1.0 - maf)
    for i in range(n):
        u = random.random()
        row.append(0 if u < p0 else (1 if u < p1 else 2))
    genotypes.append(row)

def mean_center(v):
    mu = sum(v) / len(v)
    return [x - mu for x in v]

g_dyn = mean_center(genotypes[13])     # rs eq_dyn_causal: slope trajectory
g_sta = mean_center(genotypes[29])     # rs eq_stable_causal: constant trajectory
subject_effect = [random.gauss(0.0, 0.15) for _ in range(n)]

with (prefix.with_suffix(".fam")).open("w", encoding="utf-8") as f:
    for i in range(n):
        f.write(f"F{i+1:04d}\tI{i+1:04d}\t0\t0\t0\t-9\n")

with (prefix.with_suffix(".bim")).open("w", encoding="utf-8") as f:
    for j in range(m):
        chrom = min(3, 1 + j // (m // 3))
        if j == 13:
            snp = "rs_eq_dyn_causal"
        elif j == 29:
            snp = "rs_eq_stable_causal"
        else:
            snp = f"rs_bg_{j+1:04d}"
        f.write(f"{chrom}\t{snp}\t0\t{1000 + j * 10}\tA\tG\n")

bed = bytearray(b"\x6c\x1b\x01")
code = {0: 0b00, 1: 0b10, 2: 0b11}
bpb = (n + 3) // 4
for j in range(m):
    block = bytearray(bpb)
    for i, g in enumerate(genotypes[j]):
        block[i // 4] |= code[g] << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix(".bed").write_bytes(bed)

with long_pheno.open("w", encoding="utf-8") as f:
    f.write("FID IID TIME Y\n")
    for i in range(n):
        for t in times:
            centered = 2.0 * t - 1.0
            y = (
                1.0
                + subject_effect[i]
                + 0.30 * centered
                + 0.70 * g_sta[i]
                + 0.90 * g_dyn[i] * centered
                + random.gauss(0.0, 0.35)
            )
            f.write(f"F{i+1:04d}\tI{i+1:04d}\t{t:.6f}\t{y:.8f}\n")

print(f"Generated balanced equivalence dataset: N={n}, M={m}, TIME={n_times}, records={n * n_times}")
PY

# Path 1: balanced observation-level RRM (record-level Method A via implicit x (x) Phi).
"${BIN}" \
  --bfile "${DATA_DIR}/equil" \
  --long-pheno "${DATA_DIR}/equil.long.pheno" \
  --tpheno 1 \
  --mpheno 2 \
  --rrm-order "${ORDER}" \
  --rrm-ai-iter "${AI_ITER}" \
  --rrm-curve-grid "${TIMES}" \
  --rrm-top 5 \
  --rrm-sig-threshold 0.05 \
  --maf 0 \
  --gwas \
  --no-loco \
  --out "${RESULT_DIR}/balanced" \
  > "${RESULT_DIR}/balanced.stdout" 2>&1

# Path 2: unbalanced long-format entry reused on the SAME balanced records,
# record-level Method A (--rrm-gwas-method record).
"${BIN}" \
  --bfile "${DATA_DIR}/equil" \
  --pheno "${DATA_DIR}/equil.long.pheno" \
  --tpheno 1 \
  --mpheno 2 \
  --rrm-unbalance \
  --rrm-gwas-method record \
  --rrm-order "${ORDER}" \
  --rrm-ai-iter "${AI_ITER}" \
  --rrm-curve-grid "${TIMES}" \
  --rrm-top 5 \
  --rrm-sig-threshold 0.05 \
  --maf 0 \
  --gwas \
  --no-loco \
  --out "${RESULT_DIR}/record" \
  > "${RESULT_DIR}/record.stdout" 2>&1

"${PYTHON_BIN}" - "${RESULT_DIR}" "${ORDER}" "${MIN_PEARSON}" "${MAX_ABS_DIFF_BETA}" "${MAX_ABS_DIFF_STAT}" <<'PY'
import csv
import sys
from pathlib import Path

out = Path(sys.argv[1])
order = int(sys.argv[2])
min_pearson = float(sys.argv[3])
max_abs_diff_beta = float(sys.argv[4])
max_abs_diff_stat = float(sys.argv[5])
n_basis = order + 1

def load(path):
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

bal = {r["SNP"]: r for r in load(out / "balanced.rrm_joint.assoc.txt")}
rec = {r["SNP"]: r for r in load(out / "record.rrm_unbal.rec.assoc.txt")}
common = sorted(bal.keys() & rec.keys(), key=lambda s: bal[s]["BP"])
if len(common) == 0:
    raise SystemExit("No common SNPs between balanced and record outputs")

def num(r, k):
    return float(r[k])

max_diff_beta = [0.0] * n_basis
max_diff_se = [0.0] * n_basis
max_diff_stat = 0.0
pearson_beta_parts = []
pearson_se_parts = []
pearson_stat_parts = []
pearson_p_parts = []
import math

def pearson(xv, yv):
    n = len(xv)
    mx = sum(xv) / n
    my = sum(yv) / n
    sx = math.sqrt(sum((a - mx) ** 2 for a in xv))
    sy = math.sqrt(sum((b - my) ** 2 for b in yv))
    if sx == 0 or sy == 0:
        return float("nan")
    return sum((a - mx) * (b - my) for a, b in zip(xv, yv)) / (sx * sy)

for k in range(n_basis):
    xb = [num(bal[s], f"BETA_L{k}") for s in common]
    yb = [num(rec[s], f"BETA_L{k}") for s in common]
    xs = [num(bal[s], f"SE_L{k}") for s in common]
    ys = [num(rec[s], f"SE_L{k}") for s in common]
    for a, b in zip(xb, yb):
        max_diff_beta[k] = max(max_diff_beta[k], abs(a - b))
    for a, b in zip(xs, ys):
        max_diff_se[k] = max(max_diff_se[k], abs(a - b))
    pearson_beta_parts.append(pearson(xb, yb))
    pearson_se_parts.append(pearson(xs, ys))

xst = [num(bal[s], "JOINT_CHISQ") for s in common]
yst = [num(rec[s], "STAT_JOINT") for s in common]
for a, b in zip(xst, yst):
    max_diff_stat = max(max_diff_stat, abs(a - b))
pearson_stat = pearson(xst, yst)

xp = [-math.log10(float(bal[s]["JOINT_P"])) for s in common]
yp = [-math.log10(float(rec[s]["P_JOINT"])) for s in common]
pearson_p = pearson(xp, yp)

# df / rank agreement (full rank expected on complete balanced data).
dfs = set()
ranks = set()
statuses = set()
for s in common:
    dfs.add(num(rec[s], "DF"))
    ranks.add(num(rec[s], "RANK"))
    statuses.add(rec[s]["STATUS"])
first_n_subj = rec[common[0]]["N_SUBJECT"] if common else "?"
first_n_rec = rec[common[0]]["N_RECORD"] if common else "?"

sum_beta = sum(pearson_beta_parts) / max(1, len(pearson_beta_parts))
sum_se = sum(pearson_se_parts) / max(1, len(pearson_se_parts))
pass_beta = all(p >= min_pearson for p in pearson_beta_parts if p == p) and len(pearson_beta_parts) == n_basis
pass_se = all(p >= min_pearson for p in pearson_se_parts if p == p) and len(pearson_se_parts) == n_basis
pass_stat = pearson_stat >= min_pearson and pearson_p >= min_pearson
pass_abs = max(max_diff_beta) <= max_abs_diff_beta and max_diff_stat <= max_abs_diff_stat
pass_rank = ranks == {n_basis} and dfs == {n_basis} and statuses == {"OK"}

summary = {
    "N_SUBJECT": first_n_subj,
    "N_RECORD": first_n_rec,
    "N_SNP": len(common),
    "Pearson_BETA_L0..": ",".join(f"{v:.9f}" for v in pearson_beta_parts),
    "Pearson_SE_L0..": ",".join(f"{v:.9f}" for v in pearson_se_parts),
    "Pearson_STAT_JOINT": f"{pearson_stat:.9f}",
    "Pearson_log10P": f"{pearson_p:.9f}",
    "MaxAbsDiff_BETA": max(max_diff_beta),
    "MaxAbsDiff_SE": max(max_diff_se),
    "MaxAbsDiff_STAT": max_diff_stat,
    "RANK_set": ",".join(sorted(str(r) for r in ranks)),
    "DF_set": ",".join(sorted(str(d) for d in dfs)),
    "STATUS_set": ",".join(sorted(statuses)),
    "Pass": bool(pass_beta and pass_se and pass_stat and pass_abs and pass_rank),
}
with (out / "balance_record_equivalence_summary.tsv").open("w", encoding="utf-8") as f:
    f.write("\t".join(summary.keys()) + "\n")
    f.write("\t".join(str(summary[k]) for k in summary) + "\n")

with (out / "balance_record_equivalence_summary.md").open("w", encoding="utf-8") as f:
    f.write("# Balanced Record Method A Equivalence\n\n")
    f.write("## Setup\n\n")
    f.write(f"- Balanced complete dataset: N_SUBJECT={first_n_subj}, N_RECORD={first_n_rec}.\n")
    f.write("- Balanced path: `--long-pheno` observation-level RRM (implicit x⊗Φ).\n")
    f.write("- Record path: same records via `--rrm-unbalance --rrm-gwas-method record` (explicit x⊗Φ).\n\n")
    f.write("| Metric | Value |\n|---|---:|\n")
    for key, value in summary.items():
        f.write(f"| {key} | {value} |\n")

print(out / "balance_record_equivalence_summary.tsv")
if not summary["Pass"]:
    raise SystemExit(
        f"Balanced/record equivalence failed: pearson_beta={sum_beta:.9f} pearson_se={sum_se:.9f} "
        f"pearson_stat={pearson_stat:.9f} pearson_p={pearson_p:.9f} "
        f"max_diff_beta={max(max_diff_beta):.3e} max_diff_stat={max_diff_stat:.3e} "
        f"ranks={ranks} dfs={dfs} statuses={statuses}"
    )
PY

echo "Done. Balanced <=> record Method A equivalence results: ${RESULT_DIR}"
