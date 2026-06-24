#!/bin/bash
# P2c: Large-scale (>1000 samples) smoke test for all models
# Usage: bash tests/run_large_scale_smoke.sh [--rebuild]

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="${ROOT}/build"
RESULTS_DIR="${ROOT}/tests/results/large_smoke"
mkdir -p "$RESULTS_DIR"

if [[ "${1:-}" == "--rebuild" ]]; then
    echo "=== Rebuilding ==="
    cmake --build "$BUILD_DIR" -j"$(nproc)" 2>&1 | tail -3
fi

SCANNER="${BUILD_DIR}/cosmicscanner"
if [[ ! -x "$SCANNER" ]]; then
    echo "ERROR: cosmicscanner not found. Build first."
    exit 1
fi

PASS=0
FAIL=0
TOTAL=0

run_test() {
    local name="$1" desc="$2" cmd="$3"
    TOTAL=$((TOTAL + 1))
    echo ""
    echo "========================================================================"
    echo "  [$TOTAL] $name — $desc"
    echo "========================================================================"
    local start=$(date +%s.%N)
    if eval "$cmd" > "${RESULTS_DIR}/${name}.log" 2>&1; then
        local end=$(date +%s.%N)
        local elapsed=$(echo "$end - $start" | bc)
        echo "  PASS (${elapsed}s)" | tee -a "${RESULTS_DIR}/${name}.log"
        PASS=$((PASS + 1))
    else
        local end=$(date +%s.%N)
        local elapsed=$(echo "$end - $start" | bc)
        tail -20 "${RESULTS_DIR}/${name}.log"
        echo "  FAIL (${elapsed}s)" | tee -a "${RESULTS_DIR}/${name}.log"
        FAIL=$((FAIL + 1))
    fi
}

# Datasets
N1000_PREFIX="${ROOT}/tests/data/gwas_speed_n1000_m10000/speed"
N2000_PREFIX="${ROOT}/tests/data/gwas_speed_n2000_m20000/speed"
N5000_PREFIX="${ROOT}/tests/data/gwas_speed_n5000_m50000/speed"
PIG1000="${ROOT}/examples/pig_repro_test/s1000/pig_s1000"
PIG1686="${ROOT}/examples/pig_repro_test/s1686/pig_s1686"

# =====================================================================
# 1. Standard LMM — n1000 (baseline)
# =====================================================================
run_test "lmm_n1000" "Standard LMM, 1000 samples, 10000 SNPs" \
    "$SCANNER --bfile $N1000_PREFIX --pheno ${N1000_PREFIX}.pheno --gwas --no-loco --out ${RESULTS_DIR}/lmm_n1000"

# =====================================================================
# 2. Standard LMM — n2000
# =====================================================================
run_test "lmm_n2000" "Standard LMM, 2000 samples, 20000 SNPs" \
    "$SCANNER --bfile $N2000_PREFIX --pheno ${N2000_PREFIX}.pheno --gwas --no-loco --out ${RESULTS_DIR}/lmm_n2000"

# =====================================================================
# 3. Standard LMM — n5000 (5000 x 50000)
# =====================================================================
run_test "lmm_n5000" "Standard LMM, 5000 samples, 50000 SNPs" \
    "$SCANNER --bfile $N5000_PREFIX --pheno ${N5000_PREFIX}.pheno --gwas --no-loco --out ${RESULTS_DIR}/lmm_n5000"

# =====================================================================
# 4. ThresholdModel — n1000 (binary pheno from FAM col 6)
# =====================================================================
awk '{print $1, $2, ($6>0.5) ? 1 : 0}' "${N1000_PREFIX}.fam" > "${RESULTS_DIR}/threshold_n1000.pheno"
run_test "threshold_n1000" "Threshold probit, 1000 samples, 10000 SNPs" \
    "$SCANNER --bfile $N1000_PREFIX --pheno ${RESULTS_DIR}/threshold_n1000.pheno --gwas --threshold-model --link probit --no-loco --out ${RESULTS_DIR}/threshold_n1000"

# =====================================================================
# 5. CountModel — n1000
# =====================================================================
run_test "count_n1000" "Count Poisson, 1000 samples, 10000 SNPs" \
    "$SCANNER --bfile $N1000_PREFIX --pheno ${N1000_PREFIX}.pheno --gwas --count-model --family poisson --no-loco --out ${RESULTS_DIR}/count_n1000"

# =====================================================================
# 6. SurvivalModel — n1000 (TIME=1, EVENT from FAM col 6)
# =====================================================================
awk '{print $1, $2, 1, ($6>0) ? 1 : 0}' "${N1000_PREFIX}.fam" > "${RESULTS_DIR}/survival_n1000.pheno"
run_test "survival_n1000" "Cox survival, 1000 samples, 10000 SNPs" \
    "$SCANNER --bfile $N1000_PREFIX --pheno ${RESULTS_DIR}/survival_n1000.pheno --gwas --survival-model --no-loco --out ${RESULTS_DIR}/survival_n1000"

# =====================================================================
# 7. MaternalModel — n1000 (dummy dam ID = 0 for all)
# =====================================================================
awk '{print $1, $2, 0, $6}' "${N1000_PREFIX}.fam" > "${RESULTS_DIR}/maternal_n1000.pheno"
run_test "maternal_n1000" "Maternal model, 1000 samples, 10000 SNPs" \
    "$SCANNER --bfile $N1000_PREFIX --pheno ${RESULTS_DIR}/maternal_n1000.pheno --gwas --maternal-model --no-loco --out ${RESULTS_DIR}/maternal_n1000"

# =====================================================================
# 8. SingleStepModel — n1000 (dummy pedigree)
# =====================================================================
awk '{print $2, 0, 0}' "${N1000_PREFIX}.fam" > "${RESULTS_DIR}/singlestep_n1000.ped"
awk '{print $2}' "${N1000_PREFIX}.fam" > "${RESULTS_DIR}/singlestep_n1000.genolist"
run_test "singlestep_n1000" "ssGBLUP, 1000 samples, 10000 SNPs" \
    "$SCANNER --bfile $N1000_PREFIX --pheno ${N1000_PREFIX}.pheno --gwas --single-step --ped ${RESULTS_DIR}/singlestep_n1000.ped --genotyped-list ${RESULTS_DIR}/singlestep_n1000.genolist --no-loco --out ${RESULTS_DIR}/singlestep_n1000"

# =====================================================================
# 9. GxEModel — n1000 (env-covar = pheno column)
# =====================================================================
run_test "gxe_n1000" "GxE reaction norm, 1000 samples, 10000 SNPs" \
    "$SCANNER --bfile $N1000_PREFIX --pheno ${N1000_PREFIX}.pheno --gwas --gxe-model --env-covar 3 --no-loco --out ${RESULTS_DIR}/gxe_n1000"

# =====================================================================
# 10. RepeatabilityModel — pig s1000 (multi-record)
# =====================================================================
run_test "repeatability_pig1000" "Repeatability, 1000 pigs, 3697 records, 78 SNPs" \
    "$SCANNER --bfile $PIG1000 --pheno ${PIG1000}_pheno.txt --gwas --repeatability-model --no-loco --out ${RESULTS_DIR}/repeatability_pig1000"

# =====================================================================
# 11. RepeatabilityModel — pig s1686 (larger multi-record)
# =====================================================================
run_test "repeatability_pig1686" "Repeatability, 1686 pigs, 4143 records, 78 SNPs" \
    "$SCANNER --bfile $PIG1686 --pheno ${PIG1686}_pheno.txt --gwas --repeatability-model --no-loco --out ${RESULTS_DIR}/repeatability_pig1686"

# =====================================================================
# Summary
# =====================================================================
echo ""
echo "========================================================================"
echo "  P2c Large-Scale Smoke Test Summary"
echo "========================================================================"
echo "  Total: $TOTAL  Pass: $PASS  Fail: $FAIL"
if [[ "$FAIL" -eq 0 ]]; then
    echo "  ALL PASSED"
else
    echo "  FAILURES DETECTED — check logs in ${RESULTS_DIR}/"
fi
echo "========================================================================"
exit "$FAIL"
