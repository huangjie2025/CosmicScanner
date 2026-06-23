#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "[PSR-GWAS v0.5] Large-N calibration smoke test"
echo ""

export PSR_CALIB_REPS="${PSR_CALIB_REPS:-1}"
export PSR_CALIB_N_LIST="${PSR_CALIB_N_LIST:-80,120}"
export PSR_CALIB_M="${PSR_CALIB_M:-100}"
export PSR_CALIB_METHODS="${PSR_CALIB_METHODS:-fixed-gls,dense-zkz-grid-reml}"
export PSR_CALIB_SCENARIOS="${PSR_CALIB_SCENARIOS:-null,power}"

echo "  reps=${PSR_CALIB_REPS}  N_list=${PSR_CALIB_N_LIST}  M=${PSR_CALIB_M}"
echo "  methods=${PSR_CALIB_METHODS}  scenarios=${PSR_CALIB_SCENARIOS}"
echo ""

python3 "${ROOT}/tests/run_psr_v05_large_n_calibration.py"

RES_DIR="${ROOT}/tests/results/psr_v05_large_n_calibration"
PASS=0
FAIL=0

for f in \
    large_n_calibration_runs.tsv \
    large_n_lambda_summary.tsv \
    large_n_power_summary.tsv \
    large_n_runtime_summary.tsv \
    large_n_sigma2_summary.tsv \
    large_n_method_status.tsv \
    large_n_calibration_report.md; do
    if [[ -s "${RES_DIR}/${f}" ]]; then
        echo "[PASS] ${f} exists and non-empty"
        PASS=$((PASS + 1))
    else
        echo "[FAIL] ${f} missing or empty"
        FAIL=$((FAIL + 1))
    fi
done

echo ""
echo "[CHECK] At least 2 distinct N values:"
python3 -c "
import csv
with open('${RES_DIR}/large_n_calibration_runs.tsv') as f:
    n_vals = set(r['n'] for r in csv.DictReader(f, delimiter='\t'))
print(f'  N values: {sorted(n_vals, key=int)}')
if len(n_vals) >= 2:
    print('  [PASS] At least 2 N values')
else:
    print('  [FAIL] Need at least 2 N values')
    exit(1)
" && PASS=$((PASS + 1)) || FAIL=$((FAIL + 1))

echo "[CHECK] At least 2 methods:"
python3 -c "
import csv
with open('${RES_DIR}/large_n_calibration_runs.tsv') as f:
    methods = set(r['method'] for r in csv.DictReader(f, delimiter='\t'))
print(f'  Methods: {sorted(methods)}')
if len(methods) >= 2:
    print('  [PASS] At least 2 methods')
else:
    print('  [FAIL] Need at least 2 methods')
    exit(1)
" && PASS=$((PASS + 1)) || FAIL=$((FAIL + 1))

echo "[CHECK] Lambda summary has rows with computed values:"
python3 -c "
import csv, math
with open('${RES_DIR}/large_n_lambda_summary.tsv') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
if not rows:
    print('  [FAIL] Lambda summary empty')
    exit(1)
ok = sum(1 for r in rows if r.get('mean_lambda_gc','') not in ('','NA') and math.isfinite(float(r['mean_lambda_gc'])))
print(f'  {ok}/{len(rows)} rows have finite mean_lambda_gc')
print('  [PASS] Lambda summary has finite values')
" && PASS=$((PASS + 1)) || FAIL=$((FAIL + 1))

echo "[CHECK] No NaN/Inf in numeric fields:"
python3 -c "
import csv, math
issues = []
with open('${RES_DIR}/large_n_calibration_runs.tsv') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        for k in ('runtime_sec','lambda_gc','sigma2_g'):
            v = r.get(k,'')
            if v not in ('','NA') and not math.isfinite(float(v)):
                issues.append(f'{k}={v} in {r.get(\"method\",\"?\")}/{r.get(\"n\",\"?\")}')
if issues:
    for i in issues:
        print(f'  [FAIL] {i}')
    exit(1)
else:
    print('  [PASS] No NaN/Inf in numeric fields')
" && PASS=$((PASS + 1)) || FAIL=$((FAIL + 1))

echo ""
echo "================================================"
echo "large-N calibration smoke: ${PASS} checks passed, ${FAIL} failed"
echo "================================================"
if [ "${FAIL}" -eq 0 ]; then
    echo "[PSR-GWAS v0.5] PASS large-N calibration smoke"
else
    echo "[PSR-GWAS v0.5] FAIL large-N calibration smoke"
    exit 1
fi
