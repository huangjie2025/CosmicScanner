#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${1:-${ROOT_DIR}/build/cosmicscanner}"
OUT_DIR="${COSMICSCANNER_MVLMM_REFERENCE_OUT:-${ROOT_DIR}/build/mvlmm-gemma-compat-reference}"

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

"${BIN}" --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 --multivariate --gwas --mv-method gemma-compat \
  --chr 1 --from-bp 1000 --to-bp 3000 --threads 2 --no-loco \
  --out "${OUT_DIR}/mv_gemma_compat" >"${OUT_DIR}/mv_gemma_compat.log" 2>&1

test -s "${OUT_DIR}/mv_gemma_compat.assoc.mv.txt"
test -s "${OUT_DIR}/mv_gemma_compat.mv.trace.tsv"
grep -q $'compatibility_target\tgemma_default_finite_em' "${OUT_DIR}/mv_gemma_compat.scan.meta.txt"
grep -q $'requested_method\tgemma-compat' "${OUT_DIR}/mv_gemma_compat.scan.meta.txt"
grep -q $'resolved_method\tgemma-compat' "${OUT_DIR}/mv_gemma_compat.scan.meta.txt"
grep -q $'deprecated_alias\tfalse' "${OUT_DIR}/mv_gemma_compat.scan.meta.txt"

awk 'NR == FNR { if (FNR > 1) { w[$1]=$2; p[$1]=$3; b1[$1]=$4; b2[$1]=$5; it[$1]=$6; stop[$1]=$7 } next }
     FNR > 1 {
       snp=$2;
       dw=$7-w[snp]; if (dw < 0) dw=-dw;
       dp=$8-p[snp]; if (dp < 0) dp=-dp;
       db1=$9-b1[snp]; if (db1 < 0) db1=-db1;
       db2=$11-b2[snp]; if (db2 < 0) db2=-db2;
       if (dw > 1e-5 || dp > 1e-6 || db1 > 1e-6 || db2 > 1e-6) exit 1;
       seen++;
     }
     END { if (seen != 3) exit 1 }' \
  "${ROOT_DIR}/test/gemma_compat_internal_expected.tsv" \
  "${OUT_DIR}/mv_gemma_compat.assoc.mv.txt"

awk 'NR == FNR { if (FNR > 1) { it[$1]=$6; stop[$1]=$7 } next }
     FNR > 1 {
       if ($6 != it[$2] || $8 != stop[$2]) exit 1;
       p=$10; w=-2*log(p); d=w-$0;
       seen++;
     }
     END { if (seen != 3) exit 1 }' \
  "${ROOT_DIR}/test/gemma_compat_internal_expected.tsv" \
  "${OUT_DIR}/mv_gemma_compat.mv.trace.tsv"
