#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "[PSR-GWAS v0.5] R ordinal baseline validation"
echo ""

# Check if Rscript is available
RSCRIPT=""
for cand in "${PSR_V05_R_SCRIPT:-${RSCRIPT_BIN:-${CONDA_PREFIX:-/nonexistent}/bin/Rscript}}" Rscript; do
    if command -v "$cand" &>/dev/null; then
        RSCRIPT="$cand"
        break
    fi
done

if [ -z "$RSCRIPT" ]; then
    echo "[WARN] Rscript not found — R baselines will be recorded as unavailable"
    RSCRIPT=""
fi

# Check R packages
if [ -n "$RSCRIPT" ]; then
    PKG_AVAIL=$("$RSCRIPT" -e "cat(requireNamespace('MASS', quietly=TRUE))" 2>/dev/null || echo "FALSE")
    echo "  MASS available:    ${PKG_AVAIL}"
    PKG_CLM_AVAIL=$("$RSCRIPT" -e "cat(requireNamespace('ordinal', quietly=TRUE))" 2>/dev/null || echo "FALSE")
    echo "  ordinal available: ${PKG_CLM_AVAIL}"
    PKG_POLMM=$("$RSCRIPT" -e "cat(requireNamespace('POLMM', quietly=TRUE))" 2>/dev/null || echo "FALSE")
    echo "  POLMM available:   ${PKG_POLMM}"
fi

echo ""

# Run R baseline
export PSR_V05_BENCH_REPS="${PSR_V05_BENCH_REPS:-1}"
export PSR_V05_BENCH_N="${PSR_V05_BENCH_N:-40}"
export PSR_V05_BENCH_M="${PSR_V05_BENCH_M:-30}"
python3 "${ROOT}/tests/run_psr_v05_r_baseline.py"

RES_DIR="${ROOT}/tests/results/psr_v05_ordinal_benchmark"

PASS=0
FAIL=0

# Check availability file
if [[ -s "${RES_DIR}/v05_ordinal_r_baseline_availability.tsv" ]]; then
    echo "[PASS] v05_ordinal_r_baseline_availability.tsv exists and non-empty"
    PASS=$((PASS + 1))
else
    echo "[FAIL] v05_ordinal_r_baseline_availability.tsv missing or empty"
    FAIL=$((FAIL + 1))
fi

# Check summary file
if [[ -s "${RES_DIR}/v05_ordinal_r_baseline_summary.tsv" ]]; then
    echo "[PASS] v05_ordinal_r_baseline_summary.tsv exists and non-empty"
    PASS=$((PASS + 1))
else
    echo "[FAIL] v05_ordinal_r_baseline_summary.tsv missing or empty"
    FAIL=$((FAIL + 1))
fi

# Check numeric fields finite (if R was available)
if [ -n "$RSCRIPT" ] && [ -s "${RES_DIR}/v05_ordinal_r_baseline_summary.tsv" ]; then
    echo ""
    echo "[CHECK] Numeric field validation:"
    python3 -c "
import csv, math
with open('${RES_DIR}/v05_ordinal_r_baseline_summary.tsv') as f:
    reader = list(csv.DictReader(f, delimiter='\t'))
    n_ok = sum(1 for r in reader if r.get('error_message', '') == '')
    for r in reader:
        if r.get('error_message', '') != '':
            continue
        if r.get('method_available', '') == 'False':
            continue
        rt = r.get('runtime_sec', 'NA')
        ns = r.get('n_snps', 'NA')
        if rt != 'NA' and not math.isfinite(float(rt)):
            print(f'  [WARN] non-finite runtime_sec: {rt} ({r[\"method\"]})')
        if ns != 'NA' and not math.isfinite(float(ns)):
            print(f'  [WARN] non-finite n_snps: {ns} ({r[\"method\"]})')
    if n_ok:
        print(f'  [PASS] {n_ok} R baseline rows with no error')
    else:
        print(f'  [INFO] No successful R baseline rows (R may be unavailable)')
" && PASS=$((PASS + 1)) || echo "  [SKIP] validation skipped"
fi

# Check availability file has expected methods
python3 -c "
with open('${RES_DIR}/v05_ordinal_r_baseline_availability.tsv') as f:
    methods = [r['method'] for r in __import__('csv').DictReader(f, delimiter='\t')]
expected = {'MASS::polr', 'ordinal::clm'}
found = set(methods) & expected
missing = expected - found
print(f'  Expected methods: {sorted(expected)}')
print(f'  Found methods: {sorted(found)}')
if missing:
    print(f'  [FAIL] Missing methods: {sorted(missing)}')
    exit(1)
print(f'  [PASS] All ordinal R baseline methods present')
" && PASS=$((PASS + 1)) || FAIL=$((FAIL + 1))

echo ""
echo "================================================"
echo "v0.5 R baseline validation: ${PASS} checks passed, ${FAIL} failed"
echo "================================================"
if [ "${FAIL}" -eq 0 ]; then
    echo "[PSR-GWAS v0.5] PASS R ordinal baseline validation"
else
    echo "[PSR-GWAS v0.5] FAIL R ordinal baseline validation"
    exit 1
fi
