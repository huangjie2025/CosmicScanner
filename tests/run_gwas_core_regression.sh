#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"

RUN_BUILD=0
RUN_PGEN=0
RUN_MATRIX=1
RUN_GEMMA_ACC=1
RUN_EB_SMOKE=1
RUN_REPEATABILITY_EB=1
RUN_MV_EB=1
RUN_MV_SMOKE=1
RUN_MV_VALIDATION=1
RUN_GEMMA_COMPAT_SMOKE=0
RUN_GEMMA_REFERENCE=0
RUN_RRM=0
RUN_EB_VALIDATION=0
RUN_EB_LOCO=0
RUN_EB_CALIBRATION=0
RUN_EB_REPORT=0
RUN_EB_FORMATS=0
RUN_EB_PERF=0
RUN_EB_SUITE=0
RUN_EB_GRID=0
RUN_PSR=0

usage() {
  cat <<'EOF'
Usage: tests/run_gwas_core_regression.sh [options]

Options:
  --build            Build cosmicscanner before running tests
  --with-pgen        Also run PLINK vs PGEN consistency test
  --with-rrm         Also run RRM dynamic GWAS smoke + simulation validation
  --with-eb-validation Also run EB GWAS null/power validation matrix
  --with-eb-report  Also generate standard LMM vs EB-LMM comparison report
  --with-eb-loco    Also run EB LOCO validation
  --with-eb-calibration Also run EB calibration benchmark
  --with-eb-formats Also run PLINK/PGEN/BGEN EB consistency validation
  --with-eb-perf    Also run EB-LMM performance benchmark
  --with-eb-suite   Run board-level EB GWAS suite/dashboard
  --with-eb-grid    Also run EB calibration grid benchmark
  --with-psr        Also run PSR-GWAS binary/binomial/Poisson smoke and R baseline summary
  --with-gemma-compat-smoke Also run local GEMMA/EMmvGWAS compatibility smoke without external GEMMA/R
  --with-gemma-reference Also run external GEMMA BSLMM2 reference validation when GEMMA is available
  --compat-smoke-only Run only local GEMMA/EMmvGWAS compatibility smoke
  --rrm-only         Run only RRM dynamic GWAS smoke + simulation validation
  --fast             Run EB smoke tests + multivariate smoke + GEMMA/EMmvGWAS validation
  --no-matrix        Skip multivariate validation matrix
  --no-gemma-acc     Skip single-trait GEMMA accuracy test
  -h, --help         Show this help

Default suite:
  1. Single-trait GEMMA accuracy
  2. EB-LMM single-trait smoke test
  3. Repeatability EB-LMM smoke test
  4. Multivariate EB-LMM smoke test
  5. Multivariate smoke test (emmv/gemma/gemma-full/hybrid)
  6. Multivariate GEMMA + EMmvGWAS validation
  7. Multivariate validation matrix (3 traits + covariates + missing + LOCO)

Optional:
  --with-rrm adds RRM dynamic GWAS smoke, multi-pattern dynamic QTL simulation,
  and GMAT unbalanced alignment when the reference data are available.

Use --with-pgen after reader/PGEN/BGEN-related changes.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --build) RUN_BUILD=1 ;;
    --with-pgen) RUN_PGEN=1 ;;
    --with-rrm) RUN_RRM=1 ;;
    --with-eb-validation) RUN_EB_VALIDATION=1 ;;
    --with-eb-report) RUN_EB_REPORT=1 ;;
    --with-eb-loco) RUN_EB_LOCO=1 ;;
    --with-eb-calibration) RUN_EB_CALIBRATION=1 ;;
    --with-eb-formats) RUN_EB_FORMATS=1 ;;
    --with-eb-perf) RUN_EB_PERF=1 ;;
    --with-eb-suite) RUN_EB_SUITE=1 ;;
    --with-eb-grid) RUN_EB_GRID=1 ;;
    --with-psr) RUN_PSR=1 ;;
    --with-gemma-compat-smoke) RUN_GEMMA_COMPAT_SMOKE=1 ;;
    --with-gemma-reference) RUN_GEMMA_REFERENCE=1 ;;
    --compat-smoke-only) RUN_GEMMA_ACC=0; RUN_EB_SMOKE=0; RUN_REPEATABILITY_EB=0; RUN_MV_EB=0; RUN_MATRIX=0; RUN_MV_SMOKE=0; RUN_MV_VALIDATION=0; RUN_PGEN=0; RUN_RRM=0; RUN_PSR=0; RUN_GEMMA_COMPAT_SMOKE=1 ;;
    --rrm-only) RUN_GEMMA_ACC=0; RUN_EB_SMOKE=0; RUN_REPEATABILITY_EB=0; RUN_MV_EB=0; RUN_MATRIX=0; RUN_MV_SMOKE=0; RUN_MV_VALIDATION=0; RUN_PGEN=0; RUN_RRM=1 ;;
    --fast) RUN_GEMMA_ACC=0; RUN_MATRIX=0 ;;
    --no-matrix) RUN_MATRIX=0 ;;
    --no-gemma-acc) RUN_GEMMA_ACC=0 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
  esac
  shift
done

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

run_step() {
  local name="$1"
  shift
  echo
  echo "[GWAS-core-regression] START ${name}"
  local start end
  start=$(date +%s)
  "$@"
  end=$(date +%s)
  echo "[GWAS-core-regression] PASS ${name} ($((end - start))s)"
}

run_optional_script() {
  local name="$1"
  local script="$2"
  if [[ -x "${script}" || -f "${script}" ]]; then
    run_step "${name}" bash "${script}"
  else
    echo "[GWAS-core-regression] SKIP ${name} (missing ${script})"
  fi
}

if [[ ${RUN_BUILD} -eq 1 ]]; then
  run_step "build" make -C "${ROOT}/build" cosmicscanner -j4
fi

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner binary: ${BIN}" >&2
  echo "Run with --build or set COSMICSCANNER_BIN." >&2
  exit 2
fi

if [[ ${RUN_GEMMA_ACC} -eq 1 ]]; then
  run_optional_script "single-trait GEMMA accuracy" "${ROOT}/tests/run_gemma_accuracy.sh"
fi

if [[ ${RUN_EB_SMOKE} -eq 1 ]]; then
  run_optional_script "single-trait EB-LMM smoke" "${ROOT}/tests/run_eb_lmm_smoke.sh"
fi

if [[ ${RUN_REPEATABILITY_EB} -eq 1 ]]; then
  run_optional_script "repeatability EB-LMM smoke" "${ROOT}/tests/run_repeatability_eb_smoke.sh"
fi

if [[ ${RUN_MV_EB} -eq 1 ]]; then
  run_optional_script "multivariate EB-LMM smoke" "${ROOT}/tests/run_mv_eb_smoke.sh"
fi

if [[ ${RUN_MV_SMOKE} -eq 1 ]]; then
  run_optional_script "multivariate smoke" "${ROOT}/tests/run_test_multivariate.sh"
fi

if [[ ${RUN_MV_VALIDATION} -eq 1 ]]; then
  run_optional_script "multivariate GEMMA/EMmvGWAS validation" "${ROOT}/tests/run_gemma_multivariate_validation.sh"
fi

if [[ ${RUN_GEMMA_COMPAT_SMOKE} -eq 1 ]]; then
  run_optional_script "local GEMMA/EMmvGWAS compatibility smoke" "${ROOT}/tests/run_gemma_emmvgwas_compat_smoke.sh"
fi

if [[ ${RUN_GEMMA_REFERENCE} -eq 1 ]]; then
  run_optional_script "external GEMMA BSLMM reference validation" "${ROOT}/tests/run_gemma_bslmm_reference_validation.sh"
fi

if [[ ${RUN_MATRIX} -eq 1 ]]; then
  run_optional_script "multivariate validation matrix" "${ROOT}/tests/run_multivariate_validation_matrix.sh"
fi

if [[ ${RUN_RRM} -eq 1 ]]; then
  run_step "RRM dynamic smoke" bash "${ROOT}/tests/run_test_rrm.sh"
  run_step "RRM dynamic simulation" bash "${ROOT}/tests/run_rrm_dynamic_simulation.sh"
  GMAT_DATA_DIR="${GMAT_DATA:-/tmp/opencode/gmat_data}"
  if [[ -f "${GMAT_DATA_DIR}/data/mouse_long/plink.bed" && -f "${GMAT_DATA_DIR}/data/mouse_long/phe.cs_format.txt" && -f "${GMAT_DATA_DIR}/gmat_unbal_out.0_11833.res" ]]; then
    run_step "GMAT unbalanced RRM alignment" bash "${ROOT}/tests/run_gmat_unbalanced_alignment.sh"
  else
    echo "[GWAS-core-regression] SKIP GMAT unbalanced RRM alignment (missing GMAT_DATA=${GMAT_DATA_DIR})"
  fi
fi

if [[ ${RUN_EB_VALIDATION} -eq 1 ]]; then
  run_optional_script "EB GWAS validation matrix" "${ROOT}/tests/run_eb_gwas_validation.sh"
fi

if [[ ${RUN_EB_REPORT} -eq 1 ]]; then
  run_optional_script "EB standard-vs-EB report" "${ROOT}/tests/run_eb_comparison_report.sh"
fi

if [[ ${RUN_EB_LOCO} -eq 1 ]]; then
  run_optional_script "EB GWAS LOCO validation" "${ROOT}/tests/run_eb_loco_validation.sh"
fi

if [[ ${RUN_EB_CALIBRATION} -eq 1 ]]; then
  run_optional_script "EB GWAS calibration benchmark" "${ROOT}/tests/run_eb_calibration_benchmark.sh"
fi

if [[ ${RUN_EB_FORMATS} -eq 1 ]]; then
  run_optional_script "EB PGEN/BGEN consistency" "${ROOT}/tests/run_eb_format_consistency.sh"
fi

if [[ ${RUN_EB_PERF} -eq 1 ]]; then
  run_optional_script "EB-LMM performance benchmark" "${ROOT}/tests/run_eb_performance_benchmark.sh"
fi

if [[ ${RUN_EB_SUITE} -eq 1 ]]; then
  run_optional_script "EB GWAS board suite" "${ROOT}/tests/run_eb_gwas_suite.sh"
fi

if [[ ${RUN_EB_GRID} -eq 1 ]]; then
  run_optional_script "EB calibration grid" "${ROOT}/tests/run_eb_calibration_grid.sh"
fi

if [[ ${RUN_PSR} -eq 1 ]]; then
  run_step "PSR-GWAS smoke" bash "${ROOT}/tests/run_psr_gwas_validation.sh"
  run_step "PSR-GWAS ordinal dense-zkz" bash "${ROOT}/tests/run_psr_ordinal_dense_zkz_validation.sh"
fi

if [[ ${RUN_PGEN} -eq 1 ]]; then
  run_optional_script "PGEN consistency" "${ROOT}/tests/run_pgen_consistency.sh"
fi

echo
echo "[GWAS-core-regression] ALL REQUESTED TESTS PASSED"
