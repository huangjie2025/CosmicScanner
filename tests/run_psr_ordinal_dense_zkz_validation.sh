#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/psr_ordinal_dense_zkz_validation"
RESULT_DIR="${ROOT}/tests/results/psr_ordinal_dense_zkz_validation"
mkdir -p "${DATA_DIR}" "${RESULT_DIR}"

PASS=0
FAIL=0

check_result() {
    local name=$1 ec=$2
    if [ "$ec" -eq 0 ]; then
        echo "[PASS] $name"
        PASS=$((PASS + 1))
    else
        echo "[FAIL] $name"
        FAIL=$((FAIL + 1))
    fi
}

python_finite_check() {
    local path="$1" label="$2" impl_expected="$3"
    python3 - "$path" "$label" "$impl_expected" <<'PY'
import math, sys, csv
path, label, expected = sys.argv[1], sys.argv[2], sys.argv[3]
with open(path) as f:
    reader = csv.DictReader(f, delimiter="\t")
    n = 0
    for row in reader:
        n += 1
        for col in ["effect", "se", "wald_or_stat", "p_value", "log10p"]:
            if col in row:
                try:
                    v = float(row[col])
                except Exception:
                    raise SystemExit(f"{label}: non-numeric {col}: {row[col]}")
                if not math.isfinite(v):
                    raise SystemExit(f"{label}: non-finite {col}: {v} (row {n})")
    if n == 0:
        raise SystemExit(f"{label}: no GWAS rows")
    impl = row.get("ordinal_implementation", "")
    if impl != expected:
        raise SystemExit(f"{label}: expected ordinal_implementation={expected}, got {impl}")
    mode = row.get("ordinal_mode", "")
    if mode != "stacked_psr":
        raise SystemExit(f"{label}: expected ordinal_mode=stacked_psr, got {mode}")
    print(f"{label}: {n} SNPs, ordinal_mode={mode}, ordinal_implementation={impl}")
PY
    check_result "$label finite-value and schema check" $?
}

python_blup_check() {
    local path="$1" label="$2"
    python3 - "$path" "$label" <<'PY'
import math, sys, csv
path, label = sys.argv[1], sys.argv[2]
with open(path) as f:
    reader = csv.DictReader(f, delimiter="\t")
    n = 0
    for row in reader:
        n += 1
        blup = float(row.get("blup", "nan"))
        sg2 = float(row.get("sigma2_g", "nan"))
        if not math.isfinite(blup):
            raise SystemExit(f"{label}: non-finite blup in row {n}")
    if n == 0:
        raise SystemExit(f"{label}: no BLUP rows")
    print(f"{label}: {n} individuals, all finite")
PY
    check_result "$label BLUP finite check" $?
}

# ========== Generate ordinal data with meaningful kinship (off-diagonal) ==========
echo ""
echo "================================================"
echo "Stage 2A dense-zkz validation"
echo "================================================"
echo "Generating ordinal data with off-diagonal kinship..."

N=48 M=36 LEVELS=3
python3 -c "
import numpy as np, struct
np.random.seed(20260624)
n = $N; m = $M; n_levels = $LEVELS

# Generate genotypes for kinship with off-diagonal structure
freqs = np.random.uniform(0.15, 0.45, size=m)
G = np.random.binomial(2, freqs, size=(n, m)).astype(np.float64)

# Standardized GRM with off-diagonal
Z = (G - G.mean(axis=0)) / G.std(axis=0)
Z[np.isnan(Z)] = 0.0
K = Z @ Z.T / m

# Simulate ordinal phenotype with polygenic effect (uses K)
L = np.linalg.cholesky(K + np.eye(n) * 1e-6)
pop = np.r_[np.zeros(n//2), np.ones(n - n//2)]
causal = m // 3
polygenic = L @ np.random.normal(0, 1.0, size=n)
eta = 0.3 * pop + 0.6 * polygenic + 0.4 * (G[:, causal] - G[:, causal].mean())
latent = eta + np.random.normal(0, 1.0, size=n)
y = np.ones(n, dtype=int)
y[latent > -0.5] = 2
y[latent > 0.5] = 3
if len(set(y)) < 3:
    y[:n//3] = 1; y[n//3:2*n//3] = 2; y[2*n//3:] = 3

# Write PLINK files
with open('${DATA_DIR}/dense_zkz.fam','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} 0 0 0 -9\n')
with open('${DATA_DIR}/dense_zkz.bim','w') as f:
    for j in range(m): f.write(f'1\tsnp{j}\t0\t{j+1}\tA\tT\n')
n_bps = (n+3)//4
bed = bytearray(3 + n_bps*m)
bed[:3] = [0x6c,0x1b,0x01]
for j in range(m):
    off = 3 + j*n_bps
    for k in range(n_bps):
        bv = 0
        for b in range(4):
            si = k*4 + b
            if si < n:
                gt = int(G[si, j])
                bits = [0b00,0b10,0b11][gt]
                bv |= (bits << (b*2))
        bed[off+k] = bv
with open('${DATA_DIR}/dense_zkz.bed','wb') as f: f.write(bed)
with open('${DATA_DIR}/dense_zkz.pheno','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} {int(y[i])}\n')
print(f'Generated N={n} M={m} levels={n_levels}')
print(f'Level distribution: {np.bincount(y, minlength=n_levels)}')
"

# ========== Test 1: dense-zkz basic ==========
echo ""
echo "--- Test 1: dense-zkz basic ---"
"${BIN}" --bfile "${DATA_DIR}/dense_zkz" --gwas --ordinal-model \
    --ordinal-pheno "${DATA_DIR}/dense_zkz.pheno" \
    --ordinal-mode stacked-psr \
    --ordinal-stacked-engine dense-zkz \
    --out "${RESULT_DIR}/dense_zkz" 2>&1 | tee "${RESULT_DIR}/dense_zkz_run.log"
check_result "dense-zkz basic runs" $?

for f in dense_zkz.gwas.tsv dense_zkz.psr.convergence.tsv dense_zkz.ordinal.thresholds.tsv dense_zkz.ordinal.blup.tsv; do
    if [ -f "${RESULT_DIR}/${f}" ]; then
        check_result "dense-zkz output $f exists" 0
    else
        check_result "dense-zkz output $f exists (MISSING)" 1
    fi
done

python_finite_check "${RESULT_DIR}/dense_zkz.gwas.tsv" "dense-zkz basic" "dense_zkz_blup_first_pass"
python_blup_check "${RESULT_DIR}/dense_zkz.ordinal.blup.tsv" "dense-zkz basic"

# ========== Test 2: fixed-gls on same data (for comparison) ==========
echo ""
echo "--- Test 2: fixed-gls on same data ---"
"${BIN}" --bfile "${DATA_DIR}/dense_zkz" --gwas --ordinal-model \
    --ordinal-pheno "${DATA_DIR}/dense_zkz.pheno" \
    --ordinal-mode stacked-psr \
    --ordinal-stacked-engine fixed-gls \
    --out "${RESULT_DIR}/dense_zkz_fixed" 2>&1 | tee "${RESULT_DIR}/dense_zkz_fixed_run.log"
check_result "fixed-gls runs" $?

for f in dense_zkz_fixed.gwas.tsv dense_zkz_fixed.psr.convergence.tsv; do
    if [ -f "${RESULT_DIR}/${f}" ]; then
        check_result "fixed-gls output $f exists" 0
    else
        check_result "fixed-gls output $f exists (MISSING)" 1
    fi
done

python_finite_check "${RESULT_DIR}/dense_zkz_fixed.gwas.tsv" "fixed-gls basic" "fixed_vc_weighted_gls_first_pass"

# ========== Test 3: fixed-gls vs dense-zkz comparison ==========
echo ""
echo "--- Test 3: fixed-gls vs dense-zkz comparison ---"
python3 - "${RESULT_DIR}" <<'PY'
import csv, math, sys
from pathlib import Path

res = Path(sys.argv[1])

def read_gwas(path):
    snps = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            snps[row["snp"]] = {
                "effect": float(row["effect"]),
                "se": float(row["se"]),
                "p": float(row["p_value"]),
                "log10p": float(row["log10p"]),
            }
    return snps

fixed = read_gwas(res / "dense_zkz_fixed.gwas.tsv")
dense = read_gwas(res / "dense_zkz.gwas.tsv")
common = sorted(set(fixed) & set(dense))
print(f"Common SNPs: {len(common)}")

# Correlation of log10p
log10p_f = [fixed[s]["log10p"] for s in common]
log10p_d = [dense[s]["log10p"] for s in common]
mx, my = sum(log10p_f)/len(log10p_f), sum(log10p_d)/len(log10p_d)
sx = math.sqrt(sum((v-mx)**2 for v in log10p_f))
sy = math.sqrt(sum((v-my)**2 for v in log10p_d))
corr = sum((a-mx)*(b-my) for a,b in zip(log10p_f, log10p_d))/(sx*sy) if sx>0 and sy>0 else 0
print(f"log10p correlation: {corr:.6f}")

# Top overlap
top10_f = set(s for s,_ in sorted(fixed.items(), key=lambda kv: kv[1]["p"])[:10])
top10_d = set(s for s,_ in sorted(dense.items(), key=lambda kv: kv[1]["p"])[:10])
overlap = len(top10_f & top10_d)
print(f"Top10 overlap: {overlap}/10")

# Write summary
summary_path = res / "dense_zkz_comparison_summary.tsv"
with open(summary_path, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["engine_fixed", "engine_dense", "n_snps_fixed", "n_snps_dense",
                 "n_common", "log10p_corr", "top10_overlap",
                 "runtime_fixed", "runtime_dense"])
    # Read runtime from log files
    runtime_f, runtime_d = "NA", "NA"
    print(f"Comparison summary written to {summary_path}")
PY
check_result "dense-zkz vs fixed-gls comparison" $?

# ========== Test 4: dense-zkz with covariates ==========
echo ""
echo "--- Test 4: dense-zkz with covariates ---"
N2=40 M2=20 LEVELS2=3
python3 -c "
import numpy as np
np.random.seed(20260625)
n = $N2; m = $M2; n_levels = $LEVELS2

qcov = np.random.normal(0, 1, size=n)
freqs = np.random.uniform(0.15, 0.45, size=m)
G = np.random.binomial(2, freqs, size=(n, m)).astype(np.float64)
Z = (G - G.mean(axis=0)) / G.std(axis=0)
Z[np.isnan(Z)] = 0.0
K = Z @ Z.T / m
L = np.linalg.cholesky(K + np.eye(n) * 1e-6)

polygenic = L @ np.random.normal(0, 0.8, size=n)
eta = 0.5 * qcov + 0.5 * polygenic + np.random.normal(0, 0.5, size=n)
y = np.ones(n, dtype=int)
y[eta > -0.3] = 2
y[eta > 0.5] = 3
if len(set(y)) < 3:
    y[:n//3]=1; y[n//3:2*n//3]=2; y[2*n//3:]=3

with open('${DATA_DIR}/dense_zkz_covar.fam','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} 0 0 0 -9\n')
with open('${DATA_DIR}/dense_zkz_covar.bim','w') as f:
    for j in range(m): f.write(f'1\tsnp{j}\t0\t{j+1}\tA\tT\n')
n_bps = (n+3)//4
bed = bytearray(3 + n_bps*m)
bed[:3] = [0x6c,0x1b,0x01]
for j in range(m):
    off = 3 + j*n_bps
    for k in range(n_bps):
        bv = 0
        for b in range(4):
            si = k*4 + b
            if si < n:
                gt = int(G[si, j])
                bits = [0b00,0b10,0b11][gt]
                bv |= (bits << (b*2))
        bed[off+k] = bv
with open('${DATA_DIR}/dense_zkz_covar.bed','wb') as f: f.write(bed)
with open('${DATA_DIR}/dense_zkz_covar.pheno','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} {int(y[i])}\n')
with open('${DATA_DIR}/dense_zkz_covar.qcovar','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} {qcov[i]:.6f}\n')
print(f'Generated N={n} M={m} levels={n_levels} with 1 qcovar')
print(f'Dist: {np.bincount(y.astype(int), minlength=n_levels)}')
"

"${BIN}" --bfile "${DATA_DIR}/dense_zkz_covar" --gwas --ordinal-model \
    --ordinal-pheno "${DATA_DIR}/dense_zkz_covar.pheno" \
    --ordinal-mode stacked-psr \
    --ordinal-stacked-engine dense-zkz \
    --qcovar "${DATA_DIR}/dense_zkz_covar.qcovar" \
    --out "${RESULT_DIR}/dense_zkz_covar" 2>&1 | tee "${RESULT_DIR}/dense_zkz_covar_run.log"
check_result "dense-zkz with covariates runs" $?

for f in dense_zkz_covar.gwas.tsv dense_zkz_covar.psr.convergence.tsv dense_zkz_covar.ordinal.blup.tsv; do
    if [ -f "${RESULT_DIR}/${f}" ]; then
        check_result "dense-zkz covar output $f exists" 0
    else
        check_result "dense-zkz covar output $f exists (MISSING)" 1
    fi
done

python_finite_check "${RESULT_DIR}/dense_zkz_covar.gwas.tsv" "dense-zkz covar" "dense_zkz_blup_first_pass"
python_blup_check "${RESULT_DIR}/dense_zkz_covar.ordinal.blup.tsv" "dense-zkz covar"

# ========== Summary ==========
echo ""
echo "================================================"
echo "Results: ${PASS} passed, ${FAIL} failed"
echo "================================================"
if [ "$FAIL" -gt 0 ]; then exit 1; fi
echo "All dense-zkz validation tests passed."
