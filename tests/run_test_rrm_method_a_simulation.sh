#!/usr/bin/env bash
set -euo pipefail

# Phase 4 / Method A record-level RRM-GWAS simulation wrapper.
# Runs scenarios A (null calibration) ... G (rank-deficient) and asserts the
# joint summary passes all acceptance criteria.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
BIN="${BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${DATA_DIR:-${ROOT}/tests/data/rrm_method_a_simulation}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/tests/results/rrm_method_a_simulation}"
SCENARIOS="${SCENARIOS:-A,B,C,D,E,F,G}"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

"${PYTHON_BIN}" "${ROOT}/tests/run_rrm_snp_trajectory_simulation.py" \
  --bin "${BIN}" \
  --data-dir "${DATA_DIR}" \
  --result-dir "${RESULT_DIR}" \
  --scenarios "${SCENARIOS}"

echo "Done. Method A record-level RRM-GWAS simulation results: ${RESULT_DIR}"
