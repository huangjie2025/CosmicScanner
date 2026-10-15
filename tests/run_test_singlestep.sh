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

# First two samples are founders; the third is their offspring. Remaining samples are founders.
awk 'NR==1 {p1=$1":"$2} NR==2 {p2=$1":"$2} NR==3 {print $1":"$2, p1, p2; next} {print $1":"$2, 0, 0}' \
  "${DATA_PREFIX}.fam" > "${RESULT_DIR}/pedigree.txt"

# Create genotyped list (all samples have genotypes)
awk '{print $1":"$2}' "${DATA_PREFIX}.fam" > "${RESULT_DIR}/genotyped_list.txt"

"${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${DATA_PREFIX}.pheno" \
  --gwas \
  --single-step \
  --ped "${RESULT_DIR}/pedigree.txt" \
  --genotyped-list "${RESULT_DIR}/genotyped_list.txt" \
  --write-agh \
  --no-loco \
  --out "${OUT_PREFIX}"

# Validate null model summary
NULL_FILE="${OUT_PREFIX}_ssgwas_null.txt"
if [ ! -f "${NULL_FILE}" ]; then
  echo "FAIL: missing ${NULL_FILE}"
  exit 1
fi

# Validate the Henderson A-inverse block for two founders and their offspring.
python3 - "${OUT_PREFIX}.Ainv.bin" <<'PY'
import struct
import sys

with open(sys.argv[1], "rb") as handle:
    values = struct.unpack("<" + "f" * (160 * 161 // 2), handle.read())

def lower(i, j):
    if i < j:
        i, j = j, i
    return values[i * (i + 1) // 2 + j]

expected = (
    (1.5, 0.5, -1.0),
    (0.5, 1.5, -1.0),
    (-1.0, -1.0, 2.0),
)
for i in range(3):
    for j in range(3):
        if abs(lower(i, j) - expected[i][j]) > 1e-5:
            raise SystemExit(f"A-inverse mismatch at ({i},{j}): {lower(i,j)} != {expected[i][j]}")
PY
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
