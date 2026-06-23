#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "[PSR-GWAS v0.5] Score calibration smoke test"
echo ""

export PSR_CALIB_REPS="${PSR_CALIB_REPS:-1}"
export PSR_CALIB_N_LIST="${PSR_CALIB_N_LIST:-80,120}"
export PSR_CALIB_M="${PSR_CALIB_M:-100}"
export PSR_CALIB_METHODS="${PSR_CALIB_METHODS:-fixed-gls,dense-zkz-grid-reml}"
export PSR_CALIB_TESTS="${PSR_CALIB_TESTS:-all}"

echo "  reps=${PSR_CALIB_REPS}  N_list=${PSR_CALIB_N_LIST}  M=${PSR_CALIB_M}"
echo "  methods=${PSR_CALIB_METHODS}  tests=${PSR_CALIB_TESTS}"
echo ""

python3 "${ROOT}/tests/run_psr_v05_large_n_calibration.py"

RES_DIR="${ROOT}/tests/results/psr_v05_large_n_calibration"
PASS=0
FAIL=0

echo "[CHECK] Score-specific output files:"
for f in \
    large_n_score_calibration_runs.tsv \
    large_n_score_lambda_summary.tsv \
    large_n_score_power_summary.tsv \
    large_n_score_diagnostic_report.md; do
    if [[ -s "${RES_DIR}/${f}" ]]; then
        echo "  [PASS] ${f} exists and non-empty"
        PASS=$((PASS + 1))
    else
        echo "  [FAIL] ${f} missing or empty"
        FAIL=$((FAIL + 1))
    fi
done

echo "[CHECK] Score lambda summary has psr_score + raw_score + cluster_score:"
python3 -c "
import csv
with open('${RES_DIR}/large_n_score_lambda_summary.tsv') as f:
    test_types = set(r['test_type'] for r in csv.DictReader(f, delimiter='\t'))
print(f'  Test types: {sorted(test_types)}')
expected = {'psr_score', 'raw_score', 'cluster_score'}
missing = expected - test_types
if missing:
    print(f'  [FAIL] Missing test types: {sorted(missing)}')
    exit(1)
print(f'  [PASS] psr_score, raw_score, cluster_score present')
" && PASS=$((PASS + 1)) || FAIL=$((FAIL + 1))

echo "[CHECK] Score runs file has score lambda GC columns:"
python3 -c "
import csv
with open('${RES_DIR}/large_n_score_calibration_runs.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    fields = reader.fieldnames or []
expected_cols = ['lambda_gc', 'psr_score_lambda_gc', 'raw_score_lambda_gc', 'cluster_score_lambda_gc']
missing = [c for c in expected_cols if c not in fields]
if missing:
    print(f'  [FAIL] Missing columns: {missing}')
    exit(1)
print(f'  [PASS] Wald (lambda_gc) + psr/raw/cluster score lambda GC columns present')
" && PASS=$((PASS + 1)) || FAIL=$((FAIL + 1))

echo "[CHECK] wald_psr_score_equivalent properly recorded:"
python3 -c "
import csv
with open('${RES_DIR}/large_n_score_calibration_runs.tsv') as f:
    vals = set(r.get('wald_psr_score_equivalent','') for r in csv.DictReader(f, delimiter='\t'))
print(f'  Equivalence values: {vals}')
if 'TRUE' in vals:
    print(f'  [PASS] wald_psr_score_equivalent=TRUE present')
else:
    print(f'  [FAIL] No TRUE equivalence found')
    exit(1)
" && PASS=$((PASS + 1)) || FAIL=$((FAIL + 1))

echo "[CHECK] No NaN/Inf in score numeric fields:"
python3 -c "
import csv, math
issues = []
with open('${RES_DIR}/large_n_score_calibration_runs.tsv') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        for k in ('lambda_gc','psr_score_lambda_gc','raw_score_lambda_gc','cluster_score_lambda_gc'):
            v = r.get(k,'')
            if v not in ('','NA') and not math.isfinite(float(v)):
                issues.append(f'{k}={v} in {r.get(\"method\",\"?\")}')
if issues:
    for i in issues:
        print(f'  [FAIL] {i}')
    exit(1)
else:
    print(f'  [PASS] No NaN/Inf in score lambda GC columns')
" && PASS=$((PASS + 1)) || FAIL=$((FAIL + 1))

echo ""
echo "================================================"
echo "score calibration smoke: ${PASS} checks passed, ${FAIL} failed"
echo "================================================"
if [ "${FAIL}" -eq 0 ]; then
    echo "[PSR-GWAS v0.5] PASS score calibration smoke"
else
    echo "[PSR-GWAS v0.5] FAIL score calibration smoke"
    exit 1
fi
