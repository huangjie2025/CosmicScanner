#!/usr/bin/env bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULT_DIR="${ROOT}/tests/results/threshold_smoke"
DATA_DIR="${ROOT}/tests/data/gemma_accuracy"
DATA_PREFIX="${DATA_DIR}/gemma_accuracy"
OUT_PREFIX="${RESULT_DIR}/threshold_smoke"
BIN="${ROOT}/build/cosmicscanner"

mkdir -p "${RESULT_DIR}"

echo "=== ThresholdModel Smoke Test ==="

# Generate binary phenotype: 1 if FAM pheno > 0.5, else 0
awk '{print $1, $2, ($6 > 0.5 ? 1 : 0)}' "${DATA_PREFIX}.fam" > "${RESULT_DIR}/binary.pheno"

"${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${RESULT_DIR}/binary.pheno" \
  --gwas \
  --threshold-model \
  --link probit \
  --no-loco \
  --out "${OUT_PREFIX}"

# Validate null model summary
NULL_FILE="${OUT_PREFIX}_threshold_null.txt"
if [ ! -f "${NULL_FILE}" ]; then
  echo "FAIL: missing ${NULL_FILE}"
  exit 1
fi
HEADER=$(head -1 "${NULL_FILE}")
if [ "${HEADER}" != "sigma2_g	sigma2_e	h2_liability" ]; then
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

echo "PASS: ThresholdModel smoke test"
exit 0
