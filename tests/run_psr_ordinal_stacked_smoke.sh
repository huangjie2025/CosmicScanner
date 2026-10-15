#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${ROOT}/build/cosmicscanner"
OUT_DIR="${ROOT}/tests/results/psr_ordinal_stacked_smoke"
DATA_DIR="${ROOT}/tests/data/psr_ordinal_stacked_smoke"
mkdir -p "${OUT_DIR}" "${DATA_DIR}"

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

# ========== Test 1: stacked-psr basic smoke ==========
echo ""
echo "================================================"
echo "Test 1: stacked-psr basic smoke (no covariates)"
echo "================================================"

N=60 M=50 LEVELS=4
python3 -c "
import numpy as np, struct
np.random.seed(42)
n = $N; m = $M; n_levels = $LEVELS
pheno = np.random.choice(np.arange(n_levels), size=n, p=[0.25,0.30,0.25,0.20])
with open('${DATA_DIR}/ordinal_stacked.fam','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} 0 0 0 -9\n')
with open('${DATA_DIR}/ordinal_stacked.bim','w') as f:
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
                gt = np.random.choice([0,1,2], p=[0.25,0.50,0.25])
                bits = [0b00,0b10,0b11][gt]
                bv |= (bits << (b*2))
        bed[off+k] = bv
with open('${DATA_DIR}/ordinal_stacked.bed','wb') as f: f.write(bed)
with open('${DATA_DIR}/ordinal_stacked.pheno','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} {int(pheno[i])}\n')
print(f'Generated N={n} M={m} levels={n_levels}')
"
"${BIN}" --bfile "${DATA_DIR}/ordinal_stacked" --gwas --ordinal-model \
    --ordinal-pheno "${DATA_DIR}/ordinal_stacked.pheno" \
    --ordinal-mode stacked-psr \
    --out "${OUT_DIR}/ordinal_stacked" 2>&1 | tee "${OUT_DIR}/stacked_run.log"
check_result "stacked-psr basic runs" $?

# === File existence checks ===
for f in ordinal_stacked.gwas.tsv ordinal_stacked.psr.convergence.tsv ordinal_stacked.ordinal.thresholds.tsv; do
    if [ -f "${OUT_DIR}/${f}" ]; then
        check_result "output file $f exists" 0
    else
        check_result "output file $f exists (MISSING)" 1
        echo "Error: ${OUT_DIR}/${f} not found"; exit 1
    fi
done

# === Python finite-value check ===
COSMIC_TEST_OUT_DIR="${OUT_DIR}" python3 - <<'PY'
import os
import math, sys, csv
path = os.path.join(os.environ["COSMIC_TEST_OUT_DIR"], "ordinal_stacked.gwas.tsv")
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
                    raise SystemExit(f"non-numeric {col}: {row[col]}")
                if not math.isfinite(v):
                    raise SystemExit(f"non-finite {col}: {v} (row {n})")
    if n == 0:
        raise SystemExit("no GWAS rows")
    mode = row.get("ordinal_mode", "")
    impl = row.get("ordinal_implementation", "")
    if mode != "stacked_psr":
        raise SystemExit(f"expected ordinal_mode=stacked_psr, got {mode}")
    if impl != "fixed_vc_weighted_gls_first_pass":
        raise SystemExit(f"expected ordinal_implementation=fixed_vc_weighted_gls_first_pass, got {impl}")
    print(f"stacked ordinal finite-value check: {n} SNPs finite, mode={mode}, impl={impl}")
PY
check_result "stacked-psr no-NaN finite-value check" $?

# ========== Test 2: stacked-psr with covariates ==========
echo ""
echo "================================================"
echo "Test 2: stacked-psr with covariates"
echo "================================================"

N_Q=60 M_Q=30 LEVELS_Q=3
python3 -c "
import numpy as np
np.random.seed(7)
n = $N_Q; m = $M_Q; n_levels = $LEVELS_Q

# Continuous covariate
qcov = np.random.normal(0, 1, size=n)

# Ordinal trait with some covariate effect
logit = -1.0 + 0.5 * qcov + np.random.normal(0, 0.5, size=n)
probs = np.zeros((n, n_levels))
# Simple proportional odds: 3 thresholds
thresh = np.array([-0.5, 0.5, 1.5])
for k in range(n_levels):
    if k == 0:
        probs[:, 0] = 1.0 / (1.0 + np.exp(-(thresh[0] - logit)))
    elif k == n_levels-1:
        probs[:, k] = 1.0 - 1.0/(1.0 + np.exp(-(thresh[k-1] - logit)))
    else:
        probs[:, k] = (1.0/(1.0 + np.exp(-(thresh[k] - logit))) -
                       1.0/(1.0 + np.exp(-(thresh[k-1] - logit))))
pheno = np.array([np.random.choice(n_levels, p=probs[i]) for i in range(n)])

with open('${DATA_DIR}/ordinal_covar.fam','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} 0 0 0 -9\n')
with open('${DATA_DIR}/ordinal_covar.bim','w') as f:
    for j in range(m): f.write(f'1\tsnp{j}\t0\t{j+1}\tA\tT\n')
import struct
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
                gt = np.random.choice([0,1,2], p=[0.25,0.50,0.25])
                bits = [0b00,0b10,0b11][gt]
                bv |= (bits << (b*2))
        bed[off+k] = bv
with open('${DATA_DIR}/ordinal_covar.bed','wb') as f: f.write(bed)
with open('${DATA_DIR}/ordinal_covar.pheno','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} {int(pheno[i])}\n')
with open('${DATA_DIR}/ordinal_covar.qcovar','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} {qcov[i]:.6f}\n')
print(f'Generated N={n} M={m} levels={n_levels} with 1 qcovar')
print(f'Dist: {np.bincount(pheno.astype(int), minlength=n_levels)}')
"
"${BIN}" --bfile "${DATA_DIR}/ordinal_covar" --gwas --ordinal-model \
    --ordinal-pheno "${DATA_DIR}/ordinal_covar.pheno" \
    --ordinal-mode stacked-psr \
    --qcovar "${DATA_DIR}/ordinal_covar.qcovar" \
    --out "${OUT_DIR}/ordinal_covar" 2>&1 | tee "${OUT_DIR}/ordinal_covar_run.log"
check_result "stacked-psr with covariate runs" $?

# === Check output exists and all finite ===
COSMIC_TEST_OUT_DIR="${OUT_DIR}" python3 - <<'PY'
import os
import math, sys, csv
path = os.path.join(os.environ["COSMIC_TEST_OUT_DIR"], "ordinal_covar.gwas.tsv")
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
                    raise SystemExit(f"non-numeric {col}: {row[col]}")
                if not math.isfinite(v):
                    raise SystemExit(f"non-finite {col}: {v} (row {n})")
    if n == 0:
        raise SystemExit("no GWAS rows")
    print(f"covariate test: {n} SNPs all finite")
PY
check_result "stacked-psr covariate finite-value check" $?

# Check convergence file exists and has at least 1 iteration
if [ -f "${OUT_DIR}/ordinal_covar.psr.convergence.tsv" ]; then
    N_ITER=$(tail -n +2 "${OUT_DIR}/ordinal_covar.psr.convergence.tsv" | wc -l)
    if [ "$N_ITER" -ge 1 ]; then
        check_result "covariate convergence file has >=1 iteration" 0
    else
        check_result "covariate convergence file empty" 1
    fi
else
    check_result "covariate convergence file missing" 1
fi

# ========== Test 3: fallback ==========
echo ""
echo "================================================"
echo "Test 3: fallback error path"
echo "================================================"

# Test 3a: stacked-psr without --ordinal-allow-fallback when validation fails
# Use <3 levels -> validateInputs throws before runStackedPsr
python3 -c "
n = 20; m = 10
import numpy as np
np.random.seed(0)
with open('${DATA_DIR}/ordinal_bad.fam','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} 0 0 0 -9\n')
with open('${DATA_DIR}/ordinal_bad.bim','w') as f:
    for j in range(m): f.write(f'1\tsnp{j}\t0\t{j+1}\tA\tT\n')
import struct
n_bps = (n+3)//4
bed = bytearray(3 + n_bps*m)
bed[:3] = [0x6c,0x1b,0x01]
for j in range(m):
    off = 3 + j*n_bps
    for k in range(n_bps): bed[off+k] = 0
with open('${DATA_DIR}/ordinal_bad.bed','wb') as f: f.write(bed)
with open('${DATA_DIR}/ordinal_bad.pheno','w') as f:
    for i in range(n): f.write(f'pop{i} sample{i} 0\n')
"
set +e
"${BIN}" --bfile "${DATA_DIR}/ordinal_bad" --gwas --ordinal-model \
    --ordinal-pheno "${DATA_DIR}/ordinal_bad.pheno" \
    --ordinal-mode stacked-psr \
    --out "${OUT_DIR}/ordinal_fallback_noallow" > /dev/null 2>&1
if [ $? -ne 0 ]; then
    check_result "stacked-psr no-fallback correctly errors on bad data" 0
else
    check_result "stacked-psr no-fallback (expected non-zero exit)" 1
fi

# Test 3b: also fails at prepare stage (validateInputs), not caught by run() try/catch
# This tests that the help message mentions --ordinal-allow-fallback.
"${BIN}" --bfile "${DATA_DIR}/ordinal_bad" --gwas --ordinal-model \
    --ordinal-pheno "${DATA_DIR}/ordinal_bad.pheno" \
    --ordinal-mode stacked-psr --ordinal-allow-fallback \
    --out "${OUT_DIR}/ordinal_fallback_allow" > /dev/null 2>&1
if [ $? -ne 0 ]; then
    check_result "stacked-psr with fallback (validation fail still errors, expected)" 0
else
    check_result "stacked-psr with fallback (expected non-zero on bad data)" 1
fi
set -e

# Test 3c: proper stacked success with --ordinal-allow-fallback (no fallback needed)
"${BIN}" --bfile "${DATA_DIR}/ordinal_stacked" --gwas --ordinal-model \
    --ordinal-pheno "${DATA_DIR}/ordinal_stacked.pheno" \
    --ordinal-mode stacked-psr --ordinal-allow-fallback \
    --out "${OUT_DIR}/ordinal_fallback_notrigger" 2>&1 | tee "${OUT_DIR}/fallback_notrigger.log"
check_result "stacked-psr success with --ordinal-allow-fallback" $?
COSMIC_TEST_OUT_DIR="${OUT_DIR}" python3 - <<'PY'
import os
import csv
path = os.path.join(os.environ["COSMIC_TEST_OUT_DIR"], "ordinal_fallback_notrigger.gwas.tsv")
with open(path) as f:
    reader = csv.DictReader(f, delimiter="\t")
    n = 0
    for row in reader:
        fb = row.get("ordinal_fallback_used", "?")
        mode = row.get("ordinal_mode", "?")
        n += 1
    if fb != "FALSE":
        raise SystemExit(f"expected ordinal_fallback_used=FALSE (no fallback needed), got {fb}")
    if mode != "stacked_psr":
        raise SystemExit(f"expected ordinal_mode=stacked_psr, got {mode}")
    print(f"fallback not triggered: {n} SNPs, ordinal_fallback_used={fb}")
PY
check_result "stacked-psr fallback-not-triggered check" $?

# ========== Test 4: dense-zkz engine smoke ==========
echo ""
echo "================================================"
echo "Test 4: dense-zkz engine (no covariates)"
echo "================================================"

"${BIN}" --bfile "${DATA_DIR}/ordinal_stacked" --gwas --ordinal-model \
    --ordinal-pheno "${DATA_DIR}/ordinal_stacked.pheno" \
    --ordinal-mode stacked-psr \
    --ordinal-stacked-engine dense-zkz \
    --out "${OUT_DIR}/ordinal_stacked_densezkz" 2>&1 | tee "${OUT_DIR}/stacked_densezkz_run.log"
check_result "dense-zkz basic runs" $?

# === File existence checks ===
for f in ordinal_stacked_densezkz.gwas.tsv ordinal_stacked_densezkz.psr.convergence.tsv ordinal_stacked_densezkz.ordinal.thresholds.tsv ordinal_stacked_densezkz.ordinal.blup.tsv; do
    if [ -f "${OUT_DIR}/${f}" ]; then
        check_result "dense-zkz output file $f exists" 0
    else
        check_result "dense-zkz output file $f exists (MISSING)" 1
    fi
done

# === Python finite-value check for dense-zkz results ===
COSMIC_TEST_OUT_DIR="${OUT_DIR}" python3 - <<'PY'
import os
import math, sys, csv
path = os.path.join(os.environ["COSMIC_TEST_OUT_DIR"], "ordinal_stacked_densezkz.gwas.tsv")
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
                    raise SystemExit(f"non-numeric {col}: {row[col]}")
                if not math.isfinite(v):
                    raise SystemExit(f"non-finite {col}: {v} (row {n})")
    if n == 0:
        raise SystemExit("no GWAS rows")
    mode = row.get("ordinal_mode", "")
    impl = row.get("ordinal_implementation", "")
    if mode != "stacked_psr":
        raise SystemExit(f"expected ordinal_mode=stacked_psr, got {mode}")
    if impl != "dense_zkz_blup_first_pass":
        raise SystemExit(f"expected ordinal_implementation=dense_zkz_blup_first_pass, got {impl}")
    print(f"dense-zkz finite-value check: {n} SNPs finite, mode={mode}, impl={impl}")
PY
check_result "dense-zkz no-NaN finite-value check" $?

# === BLUP file check ===
COSMIC_TEST_OUT_DIR="${OUT_DIR}" python3 - <<'PY'
import os
import math, csv
path = os.path.join(os.environ["COSMIC_TEST_OUT_DIR"], "ordinal_stacked_densezkz.ordinal.blup.tsv")
with open(path) as f:
    reader = csv.DictReader(f, delimiter="\t")
    n = 0
    for row in reader:
        n += 1
        blup = float(row.get("blup", "nan"))
        sg2 = float(row.get("sigma2_g", "nan"))
        if not math.isfinite(blup):
            raise SystemExit(f"non-finite blup in row {n}: {blup}")
    if n == 0:
        raise SystemExit("no BLUP rows")
    print(f"dense-zkz BLUP: {n} individuals, all finite")
PY
check_result "dense-zkz BLUP finite check" $?

# === dense-zkz with covariates ===
echo ""
echo "================================================"
echo "Test 5: dense-zkz with covariates"
echo "================================================"
"${BIN}" --bfile "${DATA_DIR}/ordinal_covar" --gwas --ordinal-model \
    --ordinal-pheno "${DATA_DIR}/ordinal_covar.pheno" \
    --ordinal-mode stacked-psr \
    --qcovar "${DATA_DIR}/ordinal_covar.qcovar" \
    --ordinal-stacked-engine dense-zkz \
    --out "${OUT_DIR}/ordinal_covar_densezkz" 2>&1 | tee "${OUT_DIR}/ordinal_covar_densezkz_run.log"
check_result "dense-zkz with covariates runs" $?

COSMIC_TEST_OUT_DIR="${OUT_DIR}" python3 - <<'PY'
import os
import math, sys, csv
path = os.path.join(os.environ["COSMIC_TEST_OUT_DIR"], "ordinal_covar_densezkz.gwas.tsv")
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
                    raise SystemExit(f"non-numeric {col}: {row[col]}")
                if not math.isfinite(v):
                    raise SystemExit(f"non-finite {col}: {v} (row {n})")
    if n == 0:
        raise SystemExit("no GWAS rows")
    print(f"dense-zkz covariate test: {n} SNPs all finite")
PY
check_result "dense-zkz covariate finite-value check" $?

# ========== Summary ==========
echo ""
echo "================================================"
echo "Results: ${PASS} passed, ${FAIL} failed"
echo "================================================"
if [ "$FAIL" -gt 0 ]; then exit 1; fi
echo "All stacked ordinal PSR smoke tests passed."
