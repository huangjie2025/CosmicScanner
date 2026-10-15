#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULT_DIR="${ROOT}/tests/results/eb_calibration_grid"
REPS="${REPS:-3}"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

mkdir -p "${RESULT_DIR}"

# Compact default grid for CI/smoke. Increase REPS or override GRID_SPEC for paper-scale runs.
GRID_SPEC="${GRID_SPEC:-n96_m320_e0.60_noise0.80_pop0.04 n96_m320_e0.80_noise0.70_pop0.08 n120_m420_e0.60_noise0.80_pop0.08 n120_m420_e0.90_noise0.65_pop0.04}"

SUMMARY="${RESULT_DIR}/eb_calibration_grid.summary.tsv"
printf 'label\tscenario\treps\tstandard_bonf_reps_with_hits\teb_bonf_reps_with_hits\teb_mbonf_reps_with_hits\tmean_effective_tests\tpower_causal_top10\n' >"${SUMMARY}"

for spec in ${GRID_SPEC}; do
  n=""; m=""; effect=""; noise=""; pop=""
  IFS='_' read -ra parts <<<"${spec}"
  for part in "${parts[@]}"; do
    case "${part}" in
      noise*) noise="${part#noise}" ;;
      n*) n="${part#n}" ;;
      m*) m="${part#m}" ;;
      e*) effect="${part#e}" ;;
      pop*) pop="${part#pop}" ;;
    esac
  done
  if [[ -z "${n}" || -z "${m}" || -z "${effect}" || -z "${noise}" || -z "${pop}" ]]; then
    echo "Invalid GRID_SPEC item: ${spec}" >&2
    exit 2
  fi
  echo "[EB-calibration-grid] START ${spec}"
  CAL_LABEL="grid_${spec}" N="${n}" M="${m}" EFFECT="${effect}" NOISE="${noise}" POP_SHIFT="${pop}" REPS="${REPS}" \
    bash "${ROOT}/tests/run_eb_calibration_benchmark.sh" >/dev/null
  agg="${ROOT}/tests/results/eb_calibration_benchmark_grid_${spec}/eb_calibration.aggregate.tsv"
  python3 - "${spec}" "${agg}" "${SUMMARY}" <<'PY'
import csv, sys
label, agg_path, summary_path = sys.argv[1:]
rows = list(csv.DictReader(open(agg_path), delimiter='\t'))
with open(summary_path, 'a') as out:
    for r in rows:
        out.write('\t'.join([
            label,
            r['scenario'],
            r['reps'],
            r['standard_bonf_reps_with_hits'],
            r['eb_bonf_reps_with_hits'],
            r['eb_mbonf_reps_with_hits'],
            r['mean_effective_tests'],
            r['power_causal_top10'],
        ]) + '\n')
PY
  echo "[EB-calibration-grid] PASS ${spec}"
done

python3 - "${SUMMARY}" "${RESULT_DIR}/eb_calibration_grid.summary.md" <<'PY'
import csv, sys
tsv, md = sys.argv[1:]
rows = list(csv.DictReader(open(tsv), delimiter='\t'))
if not rows:
    raise SystemExit('empty EB calibration grid summary')
with open(md, 'w') as f:
    f.write('# EB Calibration Grid Summary\n\n')
    f.write('| Label | Scenario | Reps | Standard Bonf hit reps | EB Bonf hit reps | EB-mBonf hit reps | Mean effective tests | Causal top10 reps |\n')
    f.write('|---|---|---:|---:|---:|---:|---:|---:|\n')
    for r in rows:
        f.write(f"| {r['label']} | {r['scenario']} | {r['reps']} | {r['standard_bonf_reps_with_hits']} | {r['eb_bonf_reps_with_hits']} | {r['eb_mbonf_reps_with_hits']} | {float(r['mean_effective_tests']):.4g} | {r['power_causal_top10']} |\n")

power = [r for r in rows if r['scenario'] == 'power']
if not power or sum(int(r['power_causal_top10']) for r in power) < len(power):
    raise SystemExit('weak EB calibration grid power signal')
print(f'EB calibration grid passed: {tsv}')
PY
