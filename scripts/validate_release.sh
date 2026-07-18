#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${COSMICSCANNER_BUILD_DIR:-${ROOT_DIR}/build}"
OUT_DIR="${COSMICSCANNER_VALIDATION_OUT:-${ROOT_DIR}/test/results}"

if [[ "${COSMICSCANNER_SKIP_BUILD:-0}" != "1" ]]; then
  "${ROOT_DIR}/scripts/build.sh"
fi
EXE_SUFFIX=""
case "${OSTYPE:-}" in
  msys*|cygwin*) EXE_SUFFIX=".exe" ;;
esac
BIN="${BUILD_DIR}/cosmicscanner${EXE_SUFFIX}"
rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"
printf 'rs1\nrs2\nrs3\n' >"${OUT_DIR}/extract.txt"
printf 'rs2\n' >"${OUT_DIR}/exclude.txt"

run_case() {
  local name="$1"
  shift
  echo "[validate] START ${name}"
  "$@"
  echo "[validate] PASS  ${name}"
}

assert_header_has() {
  local file="$1"
  local col="$2"
  awk -v c="${col}" 'NR == 1 { for (i = 1; i <= NF; ++i) if ($i == c) found = 1; exit found ? 0 : 1 }' "${file}"
}

assert_meta() {
  local file="$1"
  local key="$2"
  local expected="$3"
  test "$(awk -F '\t' -v k="${key}" '$1 == k { print $2 }' "${file}")" = "${expected}"
}

"${BIN}" --help | grep -q -- "--mv-method <method>"
"${BIN}" --help | grep -q -- "Hybrid and GxE are unavailable"
if "${BIN}" --bfile "${ROOT_DIR}/examples/test_data" --pheno "${ROOT_DIR}/examples/test_data.pheno" \
  --mpheno-list 1,2 --multivariate --gwas --mv-method hybrid --out "${OUT_DIR}/bad_hybrid" \
  >"${OUT_DIR}/bad_hybrid.log" 2>&1; then
  echo "Unsupported --mv-method hybrid was not rejected" >&2
  exit 1
fi
grep -q "unavailable in this focused release" "${OUT_DIR}/bad_hybrid.log"
if "${BIN}" --bfile "${ROOT_DIR}/examples/test_data" --pheno "${ROOT_DIR}/examples/test_data.pheno" \
  --mpheno-list 1,2 --multivariate --gwas --mv-method gemma-gxe --out "${OUT_DIR}/bad_gxe" \
  >"${OUT_DIR}/bad_gxe.log" 2>&1; then
  echo "Unsupported --mv-method gemma-gxe was not rejected" >&2
  exit 1
fi
grep -q "unavailable in this focused release" "${OUT_DIR}/bad_gxe.log"

run_case "single-trait LMM" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/test_data.pheno" \
  --mpheno 1 --gwas --scan-tier semi --no-loco \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/single_trait"
test -s "${OUT_DIR}/single_trait.assoc.txt"
test "$(wc -l < "${OUT_DIR}/single_trait.assoc.txt")" -eq 12
assert_meta "${OUT_DIR}/single_trait.scan.meta.txt" "scan_tier" "semi"
assert_meta "${OUT_DIR}/single_trait.scan.meta.txt" "statistical_test" "legacy_p"

run_case "single-trait LMM approx" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/test_data.pheno" \
  --mpheno 1 --gwas --scan-tier approx --no-loco \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/single_approx"
test -s "${OUT_DIR}/single_approx.assoc.txt"
assert_meta "${OUT_DIR}/single_approx.scan.meta.txt" "scan_tier" "approx"

run_case "single-trait LM all tests" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/test_data.pheno" \
  --mpheno 1 --gwas --lm 4 \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/single_lm"
test -s "${OUT_DIR}/single_lm.lm.assoc.txt"
assert_header_has "${OUT_DIR}/single_lm.lm.assoc.txt" "P_WALD"
assert_header_has "${OUT_DIR}/single_lm.lm.assoc.txt" "P_LRT"
assert_header_has "${OUT_DIR}/single_lm.lm.assoc.txt" "P_SCORE"

run_case "README single-trait exact Wald/LRT/Score" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/test_data.pheno" \
  --mpheno 1 --gwas --scan-tier exact --lmm-gemma-mode 4 --no-loco \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/single_exact"
test -s "${OUT_DIR}/single_exact.assoc.txt"
assert_header_has "${OUT_DIR}/single_exact.assoc.txt" "P_WALD"
assert_header_has "${OUT_DIR}/single_exact.assoc.txt" "P_LRT"
assert_header_has "${OUT_DIR}/single_exact.assoc.txt" "P_SCORE"
assert_meta "${OUT_DIR}/single_exact.scan.meta.txt" "scan_tier" "exact"
assert_meta "${OUT_DIR}/single_exact.scan.meta.txt" "statistical_test" "wald,lrt,score"
assert_meta "${OUT_DIR}/single_exact.scan.meta.txt" "loco" "disabled_no_loco"

run_case "multivariate LMM" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 \
  --multivariate --gwas --scan-tier semi --no-loco \
  --out "${OUT_DIR}/multivariate"
test -s "${OUT_DIR}/multivariate.assoc.mv.txt"
assert_meta "${OUT_DIR}/multivariate.scan.meta.txt" "scan_tier" "semi"
assert_meta "${OUT_DIR}/multivariate.scan.meta.txt" "multivariate_method" "emmv"

run_case "multivariate fixed-null Wald" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 \
  --multivariate --gwas --mv-method wald --no-loco \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/multivariate_wald"
test -s "${OUT_DIR}/multivariate_wald.assoc.mv.txt"
assert_meta "${OUT_DIR}/multivariate_wald.scan.meta.txt" "multivariate_method" "wald"

run_case "multivariate GEMMA Wald" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 \
  --multivariate --gwas --mv-method gemma --no-loco \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/multivariate_gemma"
test -s "${OUT_DIR}/multivariate_gemma.assoc.mv.txt"
assert_meta "${OUT_DIR}/multivariate_gemma.scan.meta.txt" "multivariate_method" "gemma"
assert_meta "${OUT_DIR}/multivariate_gemma.scan.meta.txt" "release_status" "experimental"

run_case "multivariate GEMMA full" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 \
  --multivariate --gwas --mv-method gemma-full --no-loco \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/multivariate_gemma_full"
test -s "${OUT_DIR}/multivariate_gemma_full.assoc.mv.txt"
assert_header_has "${OUT_DIR}/multivariate_gemma_full.assoc.mv.txt" "P_WALD"
assert_header_has "${OUT_DIR}/multivariate_gemma_full.assoc.mv.txt" "P_LRT"
assert_header_has "${OUT_DIR}/multivariate_gemma_full.assoc.mv.txt" "P_SCORE"
assert_meta "${OUT_DIR}/multivariate_gemma_full.scan.meta.txt" "multivariate_method" "gemma-full"
assert_meta "${OUT_DIR}/multivariate_gemma_full.scan.meta.txt" "release_status" "experimental"

run_case "README multivariate GEMMA exact full" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 \
  --multivariate --gwas --mv-method gemma-exact-full --no-loco \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/multivariate_gemma_exact_full"
test -s "${OUT_DIR}/multivariate_gemma_exact_full.assoc.mv.txt"
assert_header_has "${OUT_DIR}/multivariate_gemma_exact_full.assoc.mv.txt" "P_WALD"
assert_header_has "${OUT_DIR}/multivariate_gemma_exact_full.assoc.mv.txt" "P_LRT"
assert_header_has "${OUT_DIR}/multivariate_gemma_exact_full.assoc.mv.txt" "P_SCORE"
assert_meta "${OUT_DIR}/multivariate_gemma_exact_full.scan.meta.txt" "multivariate_method" "gemma-exact-full"
assert_meta "${OUT_DIR}/multivariate_gemma_exact_full.scan.meta.txt" "release_status" "experimental"
assert_meta "${OUT_DIR}/multivariate_gemma_exact_full.scan.meta.txt" "statistical_test" "wald,lrt,score"

run_case "multivariate EB experimental" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 \
  --multivariate --gwas --mv-method eb --no-loco \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/multivariate_eb"
test -s "${OUT_DIR}/multivariate_eb.assoc.mv.txt"
assert_meta "${OUT_DIR}/multivariate_eb.scan.meta.txt" "multivariate_method" "eb"
assert_meta "${OUT_DIR}/multivariate_eb.scan.meta.txt" "release_status" "experimental"

run_case "multivariate EMmvGWAS compatibility experimental" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 \
  --emmvgwas --gwas \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/multivariate_emmvgwas"
test -s "${OUT_DIR}/multivariate_emmvgwas.assoc.mv.txt"
assert_meta "${OUT_DIR}/multivariate_emmvgwas.scan.meta.txt" "multivariate_method" "emmv"
assert_meta "${OUT_DIR}/multivariate_emmvgwas.scan.meta.txt" "compatibility_mode" "emmvgwas"
assert_meta "${OUT_DIR}/multivariate_emmvgwas.scan.meta.txt" "release_status" "experimental"

run_case "repeatability exact mode, scan-tier semi" \
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
assert_meta "${OUT_DIR}/repeatability_semi.scan.meta.txt" "release_status" "stable"
assert_meta "${OUT_DIR}/repeatability_semi.scan.meta.txt" "scan_tier" "semi"

run_case "repeatability fast comparison" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_repeatability.pheno" \
  --qcovar "${ROOT_DIR}/examples/demo_repeatability.qcovar" \
  --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode fast --no-loco --gwas \
  --out "${OUT_DIR}/repeatability_fast"
test -s "${OUT_DIR}/repeatability_fast.assoc.txt"

run_case "repeatability STCG filtered smoke" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_repeatability.pheno" \
  --qcovar "${ROOT_DIR}/examples/demo_repeatability.qcovar" \
  --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode stcg --no-loco --gwas \
  --chr 1 --from-bp 1000 --to-bp 3000 \
  --extract "${OUT_DIR}/extract.txt" --exclude "${OUT_DIR}/exclude.txt" \
  --out "${OUT_DIR}/repeatability_stcg"
test -s "${OUT_DIR}/repeatability_stcg.stcg.assoc.txt"

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
