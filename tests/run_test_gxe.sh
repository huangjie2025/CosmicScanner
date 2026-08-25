#!/usr/bin/env bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULT_DIR="${ROOT}/tests/results/gxe_smoke"
DATA_DIR="${ROOT}/tests/data/gemma_accuracy"
DATA_PREFIX="${DATA_DIR}/gemma_accuracy"
OUT_PREFIX="${RESULT_DIR}/gxe_smoke"
BIN="${ROOT}/build/cosmicscanner"

mkdir -p "${RESULT_DIR}"

echo "=== GxEModel Smoke Test ==="

# Use the phenotype column itself as the environmental covariate (col 3 in pheno file)
"${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${DATA_PREFIX}.pheno" \
  --gwas \
  --gxe-model \
  --env-covar 3 \
  --no-loco \
  --out "${OUT_PREFIX}"

# Validate null model summary
NULL_FILE="${OUT_PREFIX}_gxe_null.txt"
if [ ! -f "${NULL_FILE}" ]; then
  echo "FAIL: missing ${NULL_FILE}"
  exit 1
fi
HEADER=$(head -1 "${NULL_FILE}")
if [ "${HEADER}" != "sigma2_g0	sigma2_g1	sigma2_e" ]; then
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

echo "PASS: GxEModel smoke test"
exit 0
