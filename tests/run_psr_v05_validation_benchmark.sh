#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"

if [[ ! -x "${BIN}" ]]; then
    echo "Missing executable: ${BIN}" >&2
    exit 2
fi

REPS="${PSR_V05_BENCH_REPS:-2}"
N="${PSR_V05_BENCH_N:-60}"
M="${PSR_V05_BENCH_M:-80}"
METHODS="${PSR_V05_BENCH_METHODS:-all}"

export PSR_V05_BENCH_REPS="$REPS"
export PSR_V05_BENCH_N="$N"
export PSR_V05_BENCH_M="$M"
export PSR_V05_BENCH_METHODS="$METHODS"

echo "[PSR-GWAS v0.5] START v0.5 ordinal validation benchmark"
echo "  reps=${REPS}  n=${N}  m=${M}  methods=${METHODS}"

python3 "${ROOT}/tests/run_psr_v05_ordinal_benchmark.py"

RES_DIR="${ROOT}/tests/results/psr_v05_ordinal_benchmark"

# Check all expected outputs
PASS=0
FAIL=0
for f in \
    v05_ordinal_benchmark_runs.tsv \
    v05_ordinal_method_summary.tsv \
    v05_ordinal_runtime_summary.tsv \
    v05_ordinal_convergence_summary.tsv \
    v05_ordinal_sigma2_summary.tsv \
    v05_ordinal_top_snp_summary.tsv \
    v05_ordinal_benchmark_report.md; do
    if [[ -s "${RES_DIR}/${f}" ]]; then
        echo "[PASS] ${f} exists and non-empty"
        PASS=$((PASS + 1))
    else
        echo "[FAIL] ${f} missing or empty"
        FAIL=$((FAIL + 1))
    fi
done

# Check method coverage
echo ""
echo "[CHECK] Method coverage in runs file:"
python3 -c "
import csv
with open('${RES_DIR}/v05_ordinal_benchmark_runs.tsv') as f:
    methods = set(r['ordinal_method'] for r in csv.DictReader(f, delimiter='\t'))
expected = {'cumulative-binary', 'fixed-gls', 'dense-zkz-fixed', 'dense-zkz-grid-reml'}
found = methods & expected
missing = expected - methods
print(f'  Found methods: {sorted(found)}')
if missing:
    print(f'  MISSING methods: {sorted(missing)}')
    print(f'  [FAIL] Missing methods')
    exit(1)
else:
    print(f'  [PASS] All 4 ordinal methods present')
" && PASS=$((PASS + 1)) || FAIL=$((FAIL + 1))

echo ""
echo "================================================"
echo "v0.5 validation benchmark: ${PASS} checks passed, ${FAIL} failed"
echo "================================================"
if [ "${FAIL}" -eq 0 ]; then
    echo "[PSR-GWAS v0.5] PASS v0.5 ordinal validation benchmark"
else
    echo "[PSR-GWAS v0.5] FAIL v0.5 ordinal validation benchmark"
    exit 1
fi
