#!/usr/bin/env bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULT_DIR="${ROOT}/tests/results/survival_smoke"
DATA_DIR="${ROOT}/tests/data/gemma_accuracy"
DATA_PREFIX="${DATA_DIR}/gemma_accuracy"
OUT_PREFIX="${RESULT_DIR}/survival_smoke"
BIN="${ROOT}/build/cosmicscanner"

mkdir -p "${RESULT_DIR}"

echo "=== SurvivalModel Smoke Test ==="

# Generate survival phenotype: FID IID TIME EVENT
# TIME = 1 (constant), EVENT = 1 if FAM pheno > 0, else 0
awk '{print $1, $2, 1, ($6 > 0 ? 1 : 0)}' "${DATA_PREFIX}.fam" > "${RESULT_DIR}/survival.pheno"

"${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${RESULT_DIR}/survival.pheno" \
  --gwas \
  --survival-model \
  --no-loco \
  --out "${OUT_PREFIX}"

# Validate null model summary
NULL_FILE="${OUT_PREFIX}_survival_null.txt"
if [ ! -f "${NULL_FILE}" ]; then
  echo "FAIL: missing ${NULL_FILE}"
  exit 1
fi
HEADER=$(head -1 "${NULL_FILE}")
if [ "${HEADER}" != "sigma2_g	sigma2_e" ]; then
  echo "FAIL: unexpected header in ${NULL_FILE}: ${HEADER}"
  exit 1
fi

# Validate assoc output
ASSOC_FILE="${OUT_PREFIX}.assoc.txt"
if [ ! -f "${ASSOC_FILE}" ]; then
  echo "FAIL: missing ${ASSOC_FILE}"
  exit 1
fi
LINE_COUNT=$(wc -l < "${ASSOC_FILE}")
if [ "${LINE_COUNT}" -lt 2 ]; then
  echo "FAIL: ${ASSOC_FILE} has only ${LINE_COUNT} lines"
  exit 1
fi

echo "PASS: SurvivalModel smoke test"
exit 0
