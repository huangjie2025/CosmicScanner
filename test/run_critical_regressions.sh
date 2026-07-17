#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${1:-${ROOT_DIR}/build/cosmicscanner}"
OUT_DIR="${COSMICSCANNER_REGRESSION_OUT:-${ROOT_DIR}/build/critical-regressions}"

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"
printf 'rs1\nrs2\nrs3\n' >"${OUT_DIR}/extract.txt"
printf 'rs2\n' >"${OUT_DIR}/exclude.txt"

common_repeatability=(
  --bfile "${ROOT_DIR}/examples/test_data"
  --pheno "${ROOT_DIR}/examples/demo_repeatability.pheno"
  --qcovar "${ROOT_DIR}/examples/demo_repeatability.qcovar"
  --tpheno 1 --mpheno 2 --repeatability-model --gwas
  --chr 1 --from-bp 1000 --to-bp 3000 --threads 4
  --extract "${OUT_DIR}/extract.txt" --exclude "${OUT_DIR}/exclude.txt"
)

assert_region() {
  local file="$1"
  awk 'NR == 1 { next } $1 != 1 || $3 < 1000 || $3 > 3000 || ($2 != "rs1" && $2 != "rs3") { exit 1 } END { if (NR != 3) exit 1 }' "${file}"
}

assert_p_equivalent() {
  local expected="$1"
  local actual="$2"
  awk 'NR == FNR { if (FNR > 1) p[$2] = $10; next }
       FNR > 1 { d = $10 - p[$2]; if (d < 0) d = -d; if (d > 1e-9) exit 1; seen++ }
       END { if (seen != 2) exit 1 }' "${expected}" "${actual}"
}

test "$("${BIN}" --version)" = "CosmicScanner 1.1.3"
"${BIN}" --help | grep -q -- "--lmm-null-solver <solver>"
"${BIN}" --help | grep -q -- "--mv-method <method>"
if "${BIN}" --bfile x --pfile y --gwas >/dev/null 2>&1; then
  echo "Multiple genotype inputs were not rejected" >&2
  exit 1
fi
if "${BIN}" --bfile x --pheno y --mpheno-list 1,2 --multivariate --gwas --mv-method gemma-gxe >/dev/null 2>"${OUT_DIR}/bad_gxe.err"; then
  echo "Unsupported multivariate GxE was not rejected" >&2
  exit 1
fi
grep -q "unavailable in this focused release" "${OUT_DIR}/bad_gxe.err"
if "${BIN}" --bfile x --gwas --from-bp 2000 --to-bp 1000 >/dev/null 2>&1; then
  echo "Invalid base-pair range was not rejected" >&2
  exit 1
fi

"${BIN}" "${common_repeatability[@]}" --no-loco --repeatability-mode exact --scan-tier semi \
  --out "${OUT_DIR}/repeatability_semi" >"${OUT_DIR}/repeatability_semi.log" 2>&1
assert_region "${OUT_DIR}/repeatability_semi.assoc.txt"

"${BIN}" "${common_repeatability[@]}" --no-loco --repeatability-mode fast \
  --out "${OUT_DIR}/repeatability_fast" >"${OUT_DIR}/repeatability_fast.log" 2>&1
assert_region "${OUT_DIR}/repeatability_fast.assoc.txt"

"${BIN}" "${common_repeatability[@]}" --no-loco --repeatability-mode stcg \
  --out "${OUT_DIR}/repeatability_stcg" >"${OUT_DIR}/repeatability_stcg.log" 2>&1
assert_region "${OUT_DIR}/repeatability_stcg.stcg.assoc.txt"
if grep -q 'Computed stats for 240 SNPs' "${OUT_DIR}/repeatability_stcg.log"; then
  echo "STCG VCE ignored SNP list filtering" >&2
  exit 1
fi
grep -q 'R3 scan with filtered-GRM VCE fallback' "${OUT_DIR}/repeatability_stcg.log"

for format in pgen bgen; do
  if [[ "${format}" == "pgen" ]]; then
    genotype_args=(--pfile "${ROOT_DIR}/examples/test_data")
  else
    genotype_args=(--bgen "${ROOT_DIR}/examples/test_data.bgen" --sample "${ROOT_DIR}/examples/test_data.sample")
  fi
  common_alt=(
    "${genotype_args[@]}"
    --pheno "${ROOT_DIR}/examples/demo_repeatability.pheno"
    --qcovar "${ROOT_DIR}/examples/demo_repeatability.qcovar"
    --tpheno 1 --mpheno 2 --repeatability-model --gwas --no-loco
    --chr 1 --from-bp 1000 --to-bp 3000 --threads 4
    --extract "${OUT_DIR}/extract.txt" --exclude "${OUT_DIR}/exclude.txt"
  )
  "${BIN}" "${common_alt[@]}" --repeatability-mode exact --scan-tier semi \
    --out "${OUT_DIR}/${format}_semi" >"${OUT_DIR}/${format}_semi.log" 2>&1
  assert_region "${OUT_DIR}/${format}_semi.assoc.txt"
  assert_p_equivalent "${OUT_DIR}/repeatability_semi.assoc.txt" "${OUT_DIR}/${format}_semi.assoc.txt"

  "${BIN}" "${common_alt[@]}" --repeatability-mode fast \
    --out "${OUT_DIR}/${format}_fast" >"${OUT_DIR}/${format}_fast.log" 2>&1
  assert_region "${OUT_DIR}/${format}_fast.assoc.txt"
  assert_p_equivalent "${OUT_DIR}/repeatability_fast.assoc.txt" "${OUT_DIR}/${format}_fast.assoc.txt"

  "${BIN}" "${common_alt[@]}" --repeatability-mode stcg \
    --out "${OUT_DIR}/${format}_stcg" >"${OUT_DIR}/${format}_stcg.log" 2>&1
  assert_region "${OUT_DIR}/${format}_stcg.stcg.assoc.txt"
  assert_p_equivalent "${OUT_DIR}/repeatability_stcg.stcg.assoc.txt" "${OUT_DIR}/${format}_stcg.stcg.assoc.txt"
done

"${BIN}" "${common_repeatability[@]}" --no-loco --repeatability-mode exact --scan-tier exact \
  --out "${OUT_DIR}/repeatability_exact" >"${OUT_DIR}/repeatability_exact.log" 2>&1
test "$(awk -F '\t' '$1 == "variants_attempted" { print $2 }' "${OUT_DIR}/repeatability_exact.refit.meta.txt")" = "2"

"${BIN}" "${common_repeatability[@]}" --repeatability-mode exact --scan-tier semi \
  --out "${OUT_DIR}/repeatability_loco" >"${OUT_DIR}/repeatability_loco.log" 2>&1
assert_region "${OUT_DIR}/repeatability_loco.assoc.txt"

baseline=""
for threads in 1 2 4 8; do
  for repetition in 1 2; do
    prefix="${OUT_DIR}/agh_t${threads}_r${repetition}"
    "${BIN}" --bfile "${ROOT_DIR}/examples/test_data" --make-agh \
      --ped "${ROOT_DIR}/examples/demo.ped" \
      --genotyped-list "${ROOT_DIR}/examples/demo.genotyped.txt" \
      --threads "${threads}" --out "${prefix}" >"${prefix}.log" 2>&1
    if [[ -z "${baseline}" ]]; then
      baseline="${prefix}.G.bin"
    else
      case "${OSTYPE:-}" in
        linux*) cmp "${baseline}" "${prefix}.G.bin" ;;
        *) python3 "${ROOT_DIR}/test/compare_float_binary.py" "${baseline}" "${prefix}.G.bin" 1e-7 ;;
      esac
    fi
  done
done

echo "Critical regression tests passed."
