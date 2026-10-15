#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"

SMOKE_DIR="${ROOT}/tests/data/psr_ordinal_stacked_smoke"
RES_DIR="${ROOT}/tests/results/psr_ordinal_vc_update_validation"
mkdir -p "${RES_DIR}"

PASS=0
FAIL=0

check_result() {
    local name=$1 ec=$2
    if [ "$ec" -eq 0 ]; then
        echo "[PASS] $name"
        PASS=$((PASS + 1))
    else
        echo "[FAIL] $name"
        FAIL=$((FAIL + 1))
    fi
}

# ========== Test 1: grid-reml basic with dense-zkz ==========
echo ""
echo "--- Test 1: grid-reml basic ---"
OUT1="${RES_DIR}/vc_test1"
"${BIN}" --bfile "${SMOKE_DIR}/ordinal_stacked" --gwas \
  --ordinal-pheno "${SMOKE_DIR}/ordinal_stacked.pheno" \
  --ordinal-model --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml \
  --ordinal-vc-grid-points 11 \
  --out "${OUT1}" 2>&1 | tail -5
check_result "grid-reml basic runs" $?

[ -f "${OUT1}.gwas.tsv" ] && check_result "grid-reml gwas exists" 0 || check_result "grid-reml gwas exists" 1
[ -f "${OUT1}.psr.convergence.tsv" ] && check_result "grid-reml convergence exists" 0 || check_result "grid-reml convergence exists" 1
[ -f "${OUT1}.ordinal.thresholds.tsv" ] && check_result "grid-reml thresholds exists" 0 || check_result "grid-reml thresholds exists" 1
[ -f "${OUT1}.ordinal.blup.tsv" ] && check_result "grid-reml blup exists" 0 || check_result "grid-reml blup exists" 1
[ -f "${OUT1}.ordinal.vc.tsv" ] && check_result "grid-reml vc log exists" 0 || check_result "grid-reml vc log exists" 1

# Finite check on GWAS
python3 -c "
import csv, sys
with open('${OUT1}.gwas.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    cols = ['effect','se','wald_or_stat','p_value','log10p']
    for row in r:
        for c in cols:
            v = float(row[c])
            if not (v == v and abs(v) < 1e100):
                print(f'FAIL: {c}={v}')
                sys.exit(1)
print(f'grid-reml basic: GWAS finite OK')
" && check_result "grid-reml GWAS finite check" 0 || check_result "grid-reml GWAS finite check" 1

# Finite check on BLUP
python3 -c "
import csv, sys
with open('${OUT1}.ordinal.blup.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        v = float(row['blup'])
        if not (v == v and abs(v) < 1e100):
            print(f'FAIL: blup={v}')
            sys.exit(1)
        v = float(row['sigma2_g'])
        if not (v == v and v >= 0):
            print(f'FAIL: sigma2_g={v}')
            sys.exit(1)
print(f'grid-reml basic: BLUP finite OK')
" && check_result "grid-reml BLUP finite check" 0 || check_result "grid-reml BLUP finite check" 1

# Check VC log is valid
python3 -c "
import csv, sys
with open('${OUT1}.ordinal.vc.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    rows = list(r)
    if len(rows) == 0:
        print('FAIL: empty VC log')
        sys.exit(1)
    for row in rows:
        sg = float(row['sigma2_g'])
        ll = float(row['reml_loglik'])
        if not (sg >= 0 and ll == ll):
            print(f'FAIL: sigma2_g={sg}, reml_loglik={ll}')
            sys.exit(1)
print(f'grid-reml VC log: {len(rows)} rows, sigma2_g={rows[-1][\"sigma2_g\"]}, logLik={rows[-1][\"reml_loglik\"]}')
" && check_result "grid-reml VC log check" 0 || check_result "grid-reml VC log check" 1

# ========== Test 2: fixed vs grid-reml comparison ==========
echo ""
echo "--- Test 2: fixed vs grid-reml comparison ---"
OUT_FIXED="${RES_DIR}/vc_test2_fixed"
OUT_GRID="${RES_DIR}/vc_test2_grid"

"${BIN}" --bfile "${SMOKE_DIR}/ordinal_stacked" --gwas \
  --ordinal-pheno "${SMOKE_DIR}/ordinal_stacked.pheno" \
  --ordinal-model --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method fixed \
  --ordinal-sigma2-g 0.5 \
  --out "${OUT_FIXED}" 2>&1 | tail -3
check_result "fixed dense-zkz runs" $?

"${BIN}" --bfile "${SMOKE_DIR}/ordinal_stacked" --gwas \
  --ordinal-pheno "${SMOKE_DIR}/ordinal_stacked.pheno" \
  --ordinal-model --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml \
  --ordinal-vc-grid-min 0.0 \
  --ordinal-vc-grid-max 2.0 \
  --ordinal-vc-grid-points 21 \
  --out "${OUT_GRID}" 2>&1 | tail -3
check_result "grid-reml dense-zkz runs" $?

# Compare
python3 -c "
import csv, sys
def load_gwas(path):
    with open(path) as f:
        r = csv.DictReader(f, delimiter='\t')
        return {row['snp']: row for row in r}

fixed = load_gwas('${OUT_FIXED}.gwas.tsv')
grid = load_gwas('${OUT_GRID}.gwas.tsv')

common = [s for s in fixed if s in grid]
n_common = len(common)
n_fixed = len(fixed)
n_grid = len(grid)

# sigma2_g
sg_fixed = float(list(fixed.values())[0]['sigma2_g'])
sg_grid = float(list(grid.values())[0]['sigma2_g'])

# log10p correlation
import math
fx = [math.log10(max(float(fixed[s]['p_value']), 1e-300)) for s in common]
gx = [math.log10(max(float(grid[s]['p_value']), 1e-300)) for s in common]
mx = sum(fx)/len(fx); my = sum(gx)/len(gx)
cov = sum((fx[i]-mx)*(gx[i]-my) for i in range(len(fx)))
varx = sum((fx[i]-mx)**2 for i in range(len(fx)))
vary = sum((gx[i]-my)**2 for i in range(len(gx)))
corr = cov / math.sqrt(varx * vary) if varx > 0 and vary > 0 else 0

# Top10 overlap
top_f = set(s for s,_ in sorted(fixed.items(), key=lambda kv: float(kv[1]['p_value']))[:10])
top_g = set(s for s,_ in sorted(grid.items(), key=lambda kv: float(kv[1]['p_value']))[:10])
overlap = len(top_f & top_g)

summary_path = '${RES_DIR}/vc_comparison_summary.tsv'
with open(summary_path, 'w') as sf:
    sf.write('engine_fixed\tengine_dense\tsigma2_fixed\tsigma2_grid\tn_snps_fixed\tn_snps_grid\tn_common\tlog10p_corr\ttop10_overlap\n')
    sf.write(f'fixed_vc\tgrid_reml\t{sg_fixed}\t{sg_grid}\t{n_fixed}\t{n_grid}\t{n_common}\t{corr:.6f}\t{overlap}\n')
print(f'sigma2_fixed={sg_fixed}, sigma2_grid={sg_grid}')
print(f'log10p corr={corr:.6f}, top10 overlap={overlap}/10')
print(f'Summary written to {summary_path}')
" && check_result "fixed vs grid-reml comparison" 0 || check_result "fixed vs grid-reml comparison" 1

# ========== Test 3: boundary test (near-zero K effect) ==========
echo ""
echo "--- Test 3: boundary test ---"
# Use ordinal_stacked data but with many fewer SNPs to reduce polygenic signal
OUT3="${RES_DIR}/vc_test3_boundary"
"${BIN}" --bfile "${SMOKE_DIR}/ordinal_stacked" --gwas \
  --ordinal-pheno "${SMOKE_DIR}/ordinal_stacked.pheno" \
  --ordinal-model --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml \
  --ordinal-vc-grid-min 0.0 \
  --ordinal-vc-grid-max 1.0 \
  --ordinal-vc-grid-points 11 \
  --out "${OUT3}" 2>&1 | tail -3
check_result "boundary test runs" $?

python3 -c "
import csv
with open('${OUT3}.ordinal.vc.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    rows = list(r)
    sg = float(rows[-1]['sigma2_g'])
    boundary = rows[-1]['boundary']
    print(f'sigma2_g={sg}, boundary={boundary}')
" && check_result "boundary check" 0 || check_result "boundary check" 1

# ========== Test 4: fixed-gls + grid-reml should error ==========
echo ""
echo "--- Test 4: fixed-gls + grid-reml error ---"
OUT4="${RES_DIR}/vc_test4_error"
set +e
"${BIN}" --bfile "${SMOKE_DIR}/ordinal_stacked" --gwas \
  --ordinal-pheno "${SMOKE_DIR}/ordinal_stacked.pheno" \
  --ordinal-model --ordinal-mode stacked-psr \
  --ordinal-stacked-engine fixed-gls \
  --ordinal-vc-method grid-reml \
  --out "${OUT4}" 2>&1 | tail -5
EC=${PIPESTATUS[0]}
set -e
# Expected to fail
if [ "$EC" -ne 0 ]; then
    check_result "fixed-gls + grid-reml correctly errors" 0
else
    check_result "fixed-gls + grid-reml correctly errors" 1
fi

# ========== Summary ==========
echo ""
echo "================================================"
echo "Results: ${PASS} passed, ${FAIL} failed"
echo "================================================"
if [ "${FAIL}" -eq 0 ]; then
    echo "All VC update validation tests passed."
else
    echo "WARNING: Some VC update validation tests failed."
    exit 1
fi
