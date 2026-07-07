#!/usr/bin/env bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULT_DIR="${ROOT}/tests/results/singlestep_smoke"
DATA_DIR="${ROOT}/tests/data/gemma_accuracy"
DATA_PREFIX="${DATA_DIR}/gemma_accuracy"
OUT_PREFIX="${RESULT_DIR}/singlestep_smoke"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"

mkdir -p "${RESULT_DIR}"

echo "=== SingleStepModel Smoke Test ==="

# Create a dummy pedigree from FAM (all unrelated)
awk '{print $1":"$2, 0, 0}' "${DATA_PREFIX}.fam" > "${RESULT_DIR}/pedigree.txt"

# Create genotyped list (all samples have genotypes)
awk '{print $1":"$2}' "${DATA_PREFIX}.fam" > "${RESULT_DIR}/genotyped_list.txt"

"${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${DATA_PREFIX}.pheno" \
  --gwas \
  --single-step \
  --ped "${RESULT_DIR}/pedigree.txt" \
  --genotyped-list "${RESULT_DIR}/genotyped_list.txt" \
  --no-loco \
  --out "${OUT_PREFIX}"

# Validate null model summary
NULL_FILE="${OUT_PREFIX}_ssgwas_null.txt"
if [ ! -f "${NULL_FILE}" ]; then
  echo "FAIL: missing ${NULL_FILE}"
  exit 1
fi
HEADER=$(head -1 "${NULL_FILE}")
if [ "${HEADER}" != "sigma2_g	sigma2_e	h2	a22_correction	h_tau	h_omega	hinv_symmetry_max_abs" ]; then
  echo "FAIL: unexpected header in ${NULL_FILE}: ${HEADER}"
  exit 1
fi
if ! awk 'NR==2 {exit !($4 == 1 && $7 < 1e-8)}' "${NULL_FILE}"; then
  echo "FAIL: A22 correction diagnostics failed in ${NULL_FILE}"
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

echo "PASS: SingleStepModel smoke test"
exit 0
