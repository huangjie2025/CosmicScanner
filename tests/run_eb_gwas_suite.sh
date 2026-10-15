#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULT_DIR="${ROOT}/tests/results/eb_gwas_suite"
RUN_PERF=1
RUN_FORMATS=1
RUN_CALIBRATION=1
RUN_GRID=1
RUN_LOCO=1

usage() {
  cat <<'EOF'
Usage: tests/run_eb_gwas_suite.sh [options]

Options:
  --quick           Run smoke + validation + standard-vs-EB report only
  --full            Run all EB checks (default)
  --no-perf         Skip EB performance benchmark
  --no-formats      Skip PGEN/BGEN EB consistency validation
  --no-calibration  Skip EB calibration benchmark
  --no-grid         Skip EB calibration grid benchmark
  --no-loco         Skip EB LOCO validation
  -h, --help        Show this help

This is the EB GWAS board-level validation entrypoint. It runs smoke,
validation, report, LOCO, calibration, format consistency, and performance
checks, then writes a single dashboard under tests/results/eb_gwas_suite/.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --quick) RUN_LOCO=0; RUN_CALIBRATION=0; RUN_GRID=0; RUN_FORMATS=0; RUN_PERF=0 ;;
    --full) RUN_LOCO=1; RUN_CALIBRATION=1; RUN_GRID=1; RUN_FORMATS=1; RUN_PERF=1 ;;
    --no-perf) RUN_PERF=0 ;;
    --no-formats) RUN_FORMATS=0 ;;
    --no-calibration) RUN_CALIBRATION=0 ;;
    --no-grid) RUN_GRID=0 ;;
    --no-loco) RUN_LOCO=0 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
  esac
  shift
done

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

mkdir -p "${RESULT_DIR}"
STATUS_TSV="${RESULT_DIR}/eb_gwas_suite.status.tsv"
printf 'step\tstatus\telapsed_sec\toutput\n' >"${STATUS_TSV}"

run_step() {
  local name="$1"
  local output="$2"
  shift 2
  local start end elapsed
  echo "[EB-GWAS-suite] START ${name}"
  start=$(date +%s)
  "$@"
  end=$(date +%s)
  elapsed=$((end - start))
  printf '%s\tPASS\t%s\t%s\n' "${name}" "${elapsed}" "${output}" >>"${STATUS_TSV}"
  echo "[EB-GWAS-suite] PASS ${name} (${elapsed}s)"
}

run_step "single-trait EB smoke" "tests/results/eb_lmm_smoke/eb_lmm_smoke.eb.summary.txt" \
  bash "${ROOT}/tests/run_eb_lmm_smoke.sh"
run_step "repeatability EB smoke" "tests/results/repeatability_eb_smoke/repeatability_eb.repeatability.eb.summary.txt" \
  bash "${ROOT}/tests/run_repeatability_eb_smoke.sh"
run_step "multivariate EB smoke" "tests/results/mv_eb_smoke/mv_eb.eb.mv.summary.txt" \
  bash "${ROOT}/tests/run_mv_eb_smoke.sh"
run_step "EB validation matrix" "tests/results/eb_gwas_validation/eb_gwas_validation.summary.tsv" \
  bash "${ROOT}/tests/run_eb_gwas_validation.sh"
run_step "standard-vs-EB report" "tests/results/eb_comparison_report/standard_vs_eb.summary.tsv" \
  bash "${ROOT}/tests/run_eb_comparison_report.sh"
if [[ ${RUN_LOCO} -eq 1 ]]; then
  run_step "EB LOCO validation" "tests/results/eb_loco_validation/loco.eb.summary.txt" \
    bash "${ROOT}/tests/run_eb_loco_validation.sh"
fi

if [[ ${RUN_CALIBRATION} -eq 1 ]]; then
  run_step "EB calibration benchmark" "tests/results/eb_calibration_benchmark/eb_calibration.aggregate.tsv" \
    bash "${ROOT}/tests/run_eb_calibration_benchmark.sh"
fi

if [[ ${RUN_GRID} -eq 1 ]]; then
  run_step "EB calibration grid" "tests/results/eb_calibration_grid/eb_calibration_grid.summary.tsv" \
    bash "${ROOT}/tests/run_eb_calibration_grid.sh"
fi

if [[ ${RUN_FORMATS} -eq 1 ]]; then
  run_step "EB PGEN/BGEN consistency" "tests/results/eb_format_consistency/eb_format_consistency.summary.tsv" \
    bash "${ROOT}/tests/run_eb_format_consistency.sh"
fi

if [[ ${RUN_PERF} -eq 1 ]]; then
  run_step "EB performance benchmark" "tests/results/eb_performance_benchmark/eb_performance.summary.tsv" \
    bash "${ROOT}/tests/run_eb_performance_benchmark.sh"
fi

python3 - "${ROOT}" "${RESULT_DIR}" <<'PY'
import csv, sys
from pathlib import Path

root = Path(sys.argv[1])
out = Path(sys.argv[2])

def read_kv(path):
    vals = {}
    if not path.exists():
        return vals
    for line in path.read_text().splitlines():
        k, _, v = line.partition('\t')
        vals[k] = v
    return vals

def read_tsv(path):
    if not path.exists():
        return []
    return list(csv.DictReader(open(path), delimiter='\t'))

status = read_tsv(out/'eb_gwas_suite.status.tsv')
validation = read_tsv(root/'tests/results/eb_gwas_validation/eb_gwas_validation.summary.tsv')
comparison = read_kv(root/'tests/results/eb_comparison_report/standard_vs_eb.summary.tsv')
calibration = read_tsv(root/'tests/results/eb_calibration_benchmark/eb_calibration.aggregate.tsv')
calibration_grid = read_tsv(root/'tests/results/eb_calibration_grid/eb_calibration_grid.summary.tsv')
formats = read_tsv(root/'tests/results/eb_format_consistency/eb_format_consistency.summary.tsv')
perf = read_tsv(root/'tests/results/eb_performance_benchmark/eb_performance.summary.tsv')

metrics = []
metrics.append(('suite_steps_passed', str(sum(1 for r in status if r.get('status') == 'PASS'))))
metrics.append(('suite_elapsed_sec', str(sum(int(float(r.get('elapsed_sec', 0))) for r in status))))
if validation:
    null_cases = [r for r in validation if r.get('scenario') == 'null']
    power_cases = [r for r in validation if r.get('scenario') == 'power']
    metrics.append(('validation_null_eb_mbonf_hits', str(sum(int(r['eb_mbonf_hits']) for r in null_cases))))
    metrics.append(('validation_power_causal_top1', str(sum(1 for r in power_cases if r.get('causal_rank') == '1'))))
if comparison:
    metrics.append(('standard_vs_eb_top10_overlap', comparison.get('top10_overlap', 'NA')))
    metrics.append(('standard_top1', comparison.get('standard_top1', 'NA')))
    metrics.append(('eb_top1', comparison.get('eb_top1', 'NA')))
if calibration:
    for r in calibration:
        prefix = 'calibration_' + r['scenario']
        metrics.append((prefix + '_eb_mbonf_reps_with_hits', r['eb_mbonf_reps_with_hits']))
        metrics.append((prefix + '_mean_effective_tests', r['mean_effective_tests']))
        metrics.append((prefix + '_power_causal_top10', r['power_causal_top10']))
if calibration_grid:
    grid_power = [r for r in calibration_grid if r['scenario'] == 'power']
    grid_null = [r for r in calibration_grid if r['scenario'] == 'null']
    metrics.append(('calibration_grid_cells', str(len(set(r['label'] for r in calibration_grid)))))
    metrics.append(('calibration_grid_power_causal_top10_total', str(sum(int(r['power_causal_top10']) for r in grid_power))))
    metrics.append(('calibration_grid_null_eb_mbonf_hit_reps_total', str(sum(int(r['eb_mbonf_reps_with_hits']) for r in grid_null))))
if formats:
    max_format_diff = max(float(r['max_abs_diff']) for r in formats)
    metrics.append(('format_consistency_max_abs_diff', f'{max_format_diff:.6g}'))
if perf:
    eb_rows = [r for r in perf if r['model'] == 'lmm-eb']
    if eb_rows:
        best = max(eb_rows, key=lambda r: float(r['rows_per_sec']))
        metrics.append(('performance_best_eb_case', best['case']))
        metrics.append(('performance_best_eb_rows_per_sec', best['rows_per_sec']))

with (out/'eb_gwas_suite.dashboard.tsv').open('w') as f:
    f.write('metric\tvalue\n')
    for k, v in metrics:
        f.write(f'{k}\t{v}\n')

with (out/'eb_gwas_suite.dashboard.md').open('w') as f:
    f.write('# EB GWAS Suite Dashboard\n\n')
    f.write('## Status\n\n')
    f.write('| Step | Status | Elapsed sec | Output |\n')
    f.write('|---|---|---:|---|\n')
    for r in status:
        f.write(f"| {r['step']} | {r['status']} | {r['elapsed_sec']} | `{r['output']}` |\n")
    f.write('\n## Key Metrics\n\n')
    f.write('| Metric | Value |\n|---|---:|\n')
    for k, v in metrics:
        f.write(f'| {k} | {v} |\n')
    if validation:
        f.write('\n## Validation Matrix\n\n')
        f.write('| Case | Model | Scenario | Effective tests | Bonf hits | EB-mBonf hits | Causal rank | Top SNP |\n')
        f.write('|---|---|---|---:|---:|---:|---:|---|\n')
        for r in validation:
            f.write(f"| {r['case']} | {r['model']} | {r['scenario']} | {float(r['effective_tests']):.4g} | {r['bonf_hits']} | {r['eb_mbonf_hits']} | {r['causal_rank']} | {r['top_snp']} |\n")
    if calibration_grid:
        f.write('\n## Calibration Grid\n\n')
        f.write('| Label | Scenario | Reps | EB-mBonf hit reps | Mean effective tests | Causal top10 reps |\n')
        f.write('|---|---|---:|---:|---:|---:|\n')
        for r in calibration_grid:
            f.write(f"| {r['label']} | {r['scenario']} | {r['reps']} | {r['eb_mbonf_reps_with_hits']} | {float(r['mean_effective_tests']):.4g} | {r['power_causal_top10']} |\n")
    f.write('\n## Reporting Note\n\n')
    f.write('EB-GWAS default routes are validated auxiliary scans. Report standard LMM/mvLMM/repeatability GWAS, ordinary Bonferroni/FDR, and EB-mBonf together; do not report EB-mBonf alone.\n')

print(f'EB GWAS suite dashboard written: {out / "eb_gwas_suite.dashboard.md"}')
PY

echo "[EB-GWAS-suite] ALL REQUESTED EB GWAS CHECKS PASSED"
