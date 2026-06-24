#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  echo "Build first: make -C ${ROOT}/build cosmicscanner -j4" >&2
  exit 2
fi

run_step() {
  local name="$1"
  shift
  echo
  echo "[repeatability-release] START ${name}"
  local start end
  start=$(date +%s)
  "$@"
  end=$(date +%s)
  echo "[repeatability-release] PASS ${name} ($((end - start))s)"
}

run_step "repeatability smoke exact/cache/LOCO/fast" bash "${ROOT}/tests/run_test_repeatability.sh"
run_step "repeatability edge cases" bash "${ROOT}/tests/run_test_repeatability_edge_cases.sh"

SIM_N="${SIM_N:-50}" \
SIM_M="${SIM_M:-120}" \
SIM_REPS="${SIM_REPS:-1}" \
run_step "repeatability exact-vs-fast simulation smoke" bash "${ROOT}/tests/run_repeatability_exact_vs_fast_simulation.sh"

echo
echo "[repeatability-release] REQUIRED VALIDATION PASSED"
echo "[repeatability-release] Optional: bash ${ROOT}/tests/run_repeatability_format_consistency.sh"
echo "[repeatability-release] Optional: bash ${ROOT}/tests/run_repeatability_performance_benchmark.sh"
