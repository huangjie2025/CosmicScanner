#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${COSMICSCANNER_BUILD_DIR:-${ROOT_DIR}/build}"
OUT_DIR="${COSMICSCANNER_VALIDATION_OUT:-${ROOT_DIR}/test/results}"

"${ROOT_DIR}/scripts/build.sh"
BIN="${BUILD_DIR}/cosmicscanner"
rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

run_case() {
  local name="$1"
  shift
  echo "[validate] START ${name}"
  "$@"
  echo "[validate] PASS  ${name}"
}

run_case "single-trait LMM" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/test_data.pheno" \
  --mpheno 1 --gwas --scan-tier semi --no-loco \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/single_trait"
test -s "${OUT_DIR}/single_trait.assoc.txt"
test "$(wc -l < "${OUT_DIR}/single_trait.assoc.txt")" -eq 12

run_case "multivariate LMM" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 \
  --multivariate --gwas --scan-tier semi --no-loco \
  --out "${OUT_DIR}/multivariate"
test -s "${OUT_DIR}/multivariate.assoc.mv.txt"

run_case "repeatability semi" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_repeatability.pheno" \
  --qcovar "${ROOT_DIR}/examples/demo_repeatability.qcovar" \
  --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode exact --gwas --scan-tier semi --no-loco \
  --out "${OUT_DIR}/repeatability_semi"
test -s "${OUT_DIR}/repeatability_semi.assoc.txt"
test -s "${OUT_DIR}/repeatability_semi.repeatability.null.txt"
test -s "${OUT_DIR}/repeatability_semi.scan.meta.txt"

run_case "repeatability fast comparison" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_repeatability.pheno" \
  --qcovar "${ROOT_DIR}/examples/demo_repeatability.qcovar" \
  --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode fast --no-loco --gwas \
  --out "${OUT_DIR}/repeatability_fast"
test -s "${OUT_DIR}/repeatability_fast.assoc.txt"

run_case "A/G/H matrix construction" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --make-agh \
  --ped "${ROOT_DIR}/examples/demo.ped" \
  --genotyped-list "${ROOT_DIR}/examples/demo.genotyped.txt" \
  --out "${OUT_DIR}/agh"
test -s "${OUT_DIR}/agh.Ainv.bin"
test -s "${OUT_DIR}/agh.Ainv.id"
test -s "${OUT_DIR}/agh.G.bin"
test -s "${OUT_DIR}/agh.G.id"
test -s "${OUT_DIR}/agh.Hinv.bin"
test -s "${OUT_DIR}/agh.Hinv.id"
test -s "${OUT_DIR}/agh.agh.meta.txt"

echo "CosmicScanner release validation passed."
