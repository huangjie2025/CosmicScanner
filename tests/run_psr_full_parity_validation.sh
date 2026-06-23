#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

run_step() {
  local name="$1"
  shift
  echo
  echo "[PSR-GWAS v0.4] START ${name}"
  local start end
  start=$(date +%s)
  "$@"
  end=$(date +%s)
  echo "[PSR-GWAS v0.4] PASS ${name} ($((end - start))s)"
}

run_step "stable CI" bash "${ROOT}/tests/run_psr_stable_ci.sh"
run_step "v0.3 PSR validation" bash "${ROOT}/tests/run_psr_gwas_validation.sh"
run_step "normal special case" bash "${ROOT}/tests/run_psr_normal_special_case.sh"
run_step "ordinal experimental MVP" bash "${ROOT}/tests/run_psr_ordinal_validation.sh"
run_step "ordinal stacked-psr smoke" bash "${ROOT}/tests/run_psr_ordinal_stacked_smoke.sh"
run_step "ordinal dense-zkz validation" bash "${ROOT}/tests/run_psr_ordinal_dense_zkz_validation.sh"
run_step "ordinal VC update validation" bash "${ROOT}/tests/run_psr_ordinal_vc_update_validation.sh"
run_step "PSR-GLMM/R availability tracker" bash "${ROOT}/tests/check_psr_glmm_r_availability.sh"
run_step "formal simulation" python3 "${ROOT}/tests/run_psr_formal_simulation.py"
run_step "accuracy/speed benchmark" python3 "${ROOT}/tests/run_psr_accuracy_speed_benchmark.py"
if [[ -x "${ROOT}/tests/summarize_psr_gwas_benchmark.py" ]]; then
  run_step "benchmark summary" python3 "${ROOT}/tests/summarize_psr_gwas_benchmark.py" "${ROOT}/tests/results/psr_gwas_smoke"
fi
run_step "full ordinal design check" bash "${ROOT}/tests/run_psr_full_ordinal_design_check.sh"

cp "${ROOT}/tests/psr_gwas_parity_matrix.tsv" "${ROOT}/tests/results/psr_gwas_smoke/psr_gwas_parity_matrix.tsv"

echo
echo "[PSR-GWAS v0.4] Full parity validation passed."
echo "[PSR-GWAS v0.4] Outputs:"
echo "  ${ROOT}/tests/results/psr_gwas_smoke/psr_validation_report.md"
echo "  ${ROOT}/tests/results/psr_gwas_smoke/psr_accuracy_speed_summary.tsv"
echo "  ${ROOT}/tests/results/psr_gwas_smoke/psr_accuracy_speed_report.md"
echo "  ${ROOT}/tests/results/psr_gwas_smoke/psr_gwas_benchmark.summary.tsv"
echo "  ${ROOT}/tests/results/psr_gwas_smoke/psr_gwas_parity_matrix.tsv"
echo "  ${ROOT}/tests/results/psr_formal_simulation/psr_formal_simulation_report.md"
echo "  ${ROOT}/tests/results/psr_gwas_smoke/psr_method_availability.tsv"
