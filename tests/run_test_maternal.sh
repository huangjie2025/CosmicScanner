#!/usr/bin/env bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULT_DIR="${ROOT}/tests/results/maternal_smoke"
DATA_DIR="${ROOT}/tests/data/gemma_accuracy"
DATA_PREFIX="${DATA_DIR}/gemma_accuracy"
OUT_PREFIX="${RESULT_DIR}/maternal_smoke"
BIN="${ROOT}/build/cosmicscanner"

mkdir -p "${RESULT_DIR}"

echo "=== MaternalModel Smoke Test ==="

# Generate pheno with 4 columns: FID IID DAM_ID PHENO
# Use "0" as dummy dam ID, pheno from FAM column 6
awk '{print $1, $2, 0, $6}' "${DATA_PREFIX}.fam" > "${RESULT_DIR}/maternal.pheno"

"${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${RESULT_DIR}/maternal.pheno" \
  --gwas \
  --maternal-model \
  --no-loco \
  --out "${OUT_PREFIX}"

# Validate null model summary
NULL_FILE="${OUT_PREFIX}_maternal_null.txt"
if [ ! -f "${NULL_FILE}" ]; then
  echo "FAIL: missing ${NULL_FILE}"
  exit 1
fi
HEADER=$(head -1 "${NULL_FILE}")
if [ "${HEADER}" != "sigma2_ad	sigma2_am	sigma_adm	sigma2_c	sigma2_e" ]; then
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

echo "PASS: MaternalModel smoke test"
exit 0
