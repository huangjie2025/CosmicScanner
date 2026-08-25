#!/bin/bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="/disk228/huangjie/PanOmics/CosmicScanner"
RESULT_DIR="${ROOT}/tests/results/repeatability_edge_cases"
DATA_DIR="${ROOT}/tests/data/repeatability_mouse_clean"
BIN="${ROOT}/build/cosmicscanner"
SUBSET_PREFIX="${DATA_DIR}/mouse_repeatability_subset"

mkdir -p "${RESULT_DIR}"

echo "=== Testing Repeatability Edge Cases ==="

# Test 1: Invalid phenotype file (empty)
echo "Test 1: Empty phenotype file"
EMPTY_PHENO="${RESULT_DIR}/empty.pheno"
touch "${EMPTY_PHENO}"
if ! "${BIN}" --bfile "${SUBSET_PREFIX}" --pheno "${EMPTY_PHENO}" --repeatability-model --gwas --no-loco --out "${RESULT_DIR}/test_empty" 2>&1; then
    echo "  PASS: Empty phenotype file correctly rejected"
else
    echo "  FAIL: Empty phenotype file should be rejected"
    exit 1
fi

# Test 2: Phenotype with all NA values
echo "Test 2: Phenotype with all NA values"
NA_PHENO="${RESULT_DIR}/all_na.pheno"
python3 - <<PY
with open("${NA_PHENO}", "w") as f:
    for i in range(10):
        f.write(f"FAM{i}\tIND{i}\tNA\n")
PY
if ! "${BIN}" --bfile "${SUBSET_PREFIX}" --pheno "${NA_PHENO}" --repeatability-model --gwas --no-loco --out "${RESULT_DIR}/test_all_na" 2>&1; then
    echo "  PASS: All-NA phenotype file correctly rejected"
else
    echo "  FAIL: All-NA phenotype file should be rejected"
    exit 1
fi

# Test 3: Single observation per subject (minimum valid case)
echo "Test 3: Single observation per subject"
SINGLE_OBS_PHENO="${RESULT_DIR}/single_obs.pheno"
python3 - <<PY
from pathlib import Path
src_fam = Path("${SUBSET_PREFIX}").with_suffix(".fam")
with src_fam.open("r") as f:
    lines = [line.strip().split()[:2] for line in f if line.strip()]
with open("${SINGLE_OBS_PHENO}", "w") as f:
    for fid, iid in lines[:10]:
        f.write(f"{fid}\t{iid}\t1.0\n")
PY
"${BIN}" --bfile "${SUBSET_PREFIX}" --pheno "${SINGLE_OBS_PHENO}" --repeatability-model --gwas --no-loco --out "${RESULT_DIR}/test_single_obs"
if [ -f "${RESULT_DIR}/test_single_obs.assoc.txt" ]; then
    echo "  PASS: Single observation per subject works"
else
    echo "  FAIL: Single observation per subject should work"
    exit 1
fi

# Test 4: MAF filter boundary
echo "Test 4: MAF filter boundary values"
"${BIN}" --bfile "${SUBSET_PREFIX}" --pheno "${SINGLE_OBS_PHENO}" --repeatability-model --gwas --no-loco --maf 0.0 --out "${RESULT_DIR}/test_maf_0"
"${BIN}" --bfile "${SUBSET_PREFIX}" --pheno "${SINGLE_OBS_PHENO}" --repeatability-model --gwas --no-loco --maf 0.5 --out "${RESULT_DIR}/test_maf_05"
echo "  PASS: MAF filter boundary tests completed"

# Test 5: Verify output files exist and have expected format
echo "Test 5: Output file validation"
if [ -f "${RESULT_DIR}/test_single_obs.repeatability.null.txt" ]; then
    if grep -q "n_records" "${RESULT_DIR}/test_single_obs.repeatability.null.txt" && \
       grep -q "n_subjects" "${RESULT_DIR}/test_single_obs.repeatability.null.txt" && \
       grep -q "h2" "${RESULT_DIR}/test_single_obs.repeatability.null.txt"; then
        echo "  PASS: Null summary file has expected format"
    else
        echo "  FAIL: Null summary file missing expected fields"
        exit 1
    fi
else
    echo "  FAIL: Null summary file missing"
    exit 1
fi

# Test 6: Test with covariates
echo "Test 6: Repeatability with covariates"
COVAR_FILE="${RESULT_DIR}/covar.txt"
python3 - <<PY
from pathlib import Path
src_fam = Path("${SUBSET_PREFIX}").with_suffix(".fam")
with src_fam.open("r") as f:
    lines = [line.strip().split()[:2] for line in f if line.strip()]
with open("${COVAR_FILE}", "w") as f:
    f.write("FID\tIID\tage\tsex\n")
    for idx, (fid, iid) in enumerate(lines[:10]):
        f.write(f"{fid}\t{iid}\t{20 + idx}\t{idx % 2}\n")
PY
"${BIN}" --bfile "${SUBSET_PREFIX}" --pheno "${SINGLE_OBS_PHENO}" --covar "${COVAR_FILE}" --repeatability-model --gwas --no-loco --out "${RESULT_DIR}/test_with_covar"
if [ -f "${RESULT_DIR}/test_with_covar.assoc.txt" ]; then
    echo "  PASS: Repeatability with covariates works"
else
    echo "  FAIL: Repeatability with covariates failed"
    exit 1
fi

# Test 7: Test cache invalidation (corrupted cache file)
echo "Test 7: Cache invalidation (corrupted cache)"
CORRUPTED_CACHE="${DATA_DIR}/corrupted_cache"
mkdir -p "${CORRUPTED_CACHE}"
echo "invalid binary data" > "${CORRUPTED_CACHE}/grm.global.f32grm.bin"
python3 - <<PY
with open("${CORRUPTED_CACHE}/grm.global.meta.txt", "w") as f:
    f.write("scale 100.0\nchroms 1 2\n")
PY
"${BIN}" --bfile "${SUBSET_PREFIX}" --pheno "${SINGLE_OBS_PHENO}" --repeatability-model --repeatability-grm-cache "${CORRUPTED_CACHE}/grm" --gwas --no-loco --out "${RESULT_DIR}/test_corrupted_cache"
if [ -f "${RESULT_DIR}/test_corrupted_cache.assoc.txt" ]; then
    echo "  PASS: Corrupted cache handled correctly (recomputed)"
else
    echo "  FAIL: Corrupted cache should be handled gracefully"
    exit 1
fi

# Test 8: Test with quantitative covariates
echo "Test 8: Repeatability with quantitative covariates"
QCOVAR_FILE="${RESULT_DIR}/qcovar.txt"
python3 - <<PY
from pathlib import Path
src_fam = Path("${SUBSET_PREFIX}").with_suffix(".fam")
with src_fam.open("r") as f:
    lines = [line.strip().split()[:2] for line in f if line.strip()]
with open("${QCOVAR_FILE}", "w") as f:
    f.write("FID\tIID\tpc1\tpc2\n")
    for idx, (fid, iid) in enumerate(lines[:10]):
        f.write(f"{fid}\t{iid}\t{0.1*idx}\t{-0.05*idx}\n")
PY
"${BIN}" --bfile "${SUBSET_PREFIX}" --pheno "${SINGLE_OBS_PHENO}" --qcovar "${QCOVAR_FILE}" --repeatability-model --gwas --no-loco --out "${RESULT_DIR}/test_with_qcovar"
if [ -f "${RESULT_DIR}/test_with_qcovar.assoc.txt" ]; then
    echo "  PASS: Repeatability with quantitative covariates works"
else
    echo "  FAIL: Repeatability with quantitative covariates failed"
    exit 1
fi

echo ""
echo "=== All edge case tests passed! ==="
