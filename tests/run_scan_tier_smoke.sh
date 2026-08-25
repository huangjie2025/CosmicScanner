#!/bin/bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="/disk228/huangjie/PanOmics/CosmicScanner"
BIN="${ROOT}/build/cosmicscanner"
OUT_ROOT="${SCAN_TIER_RESULT_ROOT:-/tmp/opencode/cosmic_scan_tier_smoke}"
mkdir -p "${OUT_ROOT}"

UNI="${ROOT}/tests/data/gemma_accuracy/gemma_accuracy"
for tier in approx semi exact; do
  "${BIN}" --bfile "${UNI}" --pheno "${UNI}.pheno" --mpheno 1 \
    --gwas --no-loco --maf 0 --scan-tier "${tier}" --out "${OUT_ROOT}/uni.${tier}"
  test -s "${OUT_ROOT}/uni.${tier}.assoc.txt"
  test -s "${OUT_ROOT}/uni.${tier}.scan.meta.txt"
done

python3 - "${OUT_ROOT}/uni.approx.assoc.txt" "${OUT_ROOT}/uni.semi.assoc.txt" <<'PY'
import csv
import sys

def rows(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

approx, semi = rows(sys.argv[1]), rows(sys.argv[2])
if len(approx) != len(semi) or not approx:
    raise SystemExit("Single-trait approx/semi row mismatch")
if max(abs(float(a["BETA"]) - float(s["BETA"])) for a, s in zip(approx, semi)) > 1e-10:
    raise SystemExit("Single-trait semi changed beta under a fixed covariance shape")
if not any(abs(float(a["SE"]) - float(s["SE"])) > 1e-8 for a, s in zip(approx, semi)):
    raise SystemExit("Single-trait semi did not profile SNP-specific scale")
PY

REP="${ROOT}/tests/data/repeatability_exact_vs_fast_sim/sim"
REP_PHENO="${ROOT}/tests/data/repeatability_exact_vs_fast_sim/balanced.rep1.pheno"
for tier in approx semi; do
  "${BIN}" --bfile "${REP}" --pheno "${REP_PHENO}" --tpheno 1 --mpheno 2 \
    --repeatability-model --repeatability-mode exact --scan-tier "${tier}" \
    --no-loco --maf 0 --out "${OUT_ROOT}/rep.${tier}"
  test -s "${OUT_ROOT}/rep.${tier}.assoc.txt"
  grep -q $'scan_tier\t'"${tier}" "${OUT_ROOT}/rep.${tier}.repeatability.null.txt"
done

python3 - "${OUT_ROOT}/rep.approx.assoc.txt" "${OUT_ROOT}/rep.semi.assoc.txt" <<'PY'
import csv
import sys

def rows(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

approx, semi = rows(sys.argv[1]), rows(sys.argv[2])
if len(approx) != len(semi) or not approx:
    raise SystemExit("Repeatability approx/semi row mismatch")
if max(abs(float(a["BETA"]) - float(s["BETA"])) for a, s in zip(approx, semi)) > 1e-10:
    raise SystemExit("Repeatability semi changed beta under a fixed covariance shape")
if not any(abs(float(a["SE"]) - float(s["SE"])) > 1e-8 for a, s in zip(approx, semi)):
    raise SystemExit("Repeatability semi did not profile SNP-specific scale")
PY

"${BIN}" --bfile "${REP}" --pheno "${REP_PHENO}" --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode exact --scan-tier exact \
  --repeatability-refit-limit 10 --no-loco --maf 0 --out "${OUT_ROOT}/rep.exact"
test -s "${OUT_ROOT}/rep.exact.assoc.txt"
grep -q $'scan_tier\texact' "${OUT_ROOT}/rep.exact.repeatability.null.txt"
test "$(wc -l < "${OUT_ROOT}/rep.exact.assoc.txt")" -eq 11

MV="${ROOT}/tests/data/mv_eb_smoke/mv_eb"
for tier in approx semi; do
  "${BIN}" --bfile "${MV}" --pheno "${MV}.pheno" --mpheno-list 1,2 \
    --multivariate --scan-tier "${tier}" --no-loco --maf 0 --out "${OUT_ROOT}/mv.${tier}"
  test -s "${OUT_ROOT}/mv.${tier}.assoc.mv.txt"
  test -s "${OUT_ROOT}/mv.${tier}.scan.meta.txt"
done

"${BIN}" --bfile "${MV}" --pheno "${MV}.pheno" --mpheno-list 1,2 \
  --multivariate --scan-tier exact --no-loco --maf 0 --chr 1 --to-bp 1200 \
  --out "${OUT_ROOT}/mv.exact"
test -s "${OUT_ROOT}/mv.exact.assoc.mv.txt"
grep -q $'resolved_scan_tier\texact' "${OUT_ROOT}/mv.exact.scan.meta.txt"

if "${BIN}" --bfile "${REP}" --pheno "${REP_PHENO}" --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode fast --scan-tier semi --no-loco \
  --out "${OUT_ROOT}/invalid" >/dev/null 2>&1; then
  echo "Expected incompatible repeatability tier/engine command to fail" >&2
  exit 1
fi

echo "Scan-tier smoke tests passed."
