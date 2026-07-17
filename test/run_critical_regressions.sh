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
"${BIN}" --help | grep -q -- "gemma-compat"
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

"${BIN}" --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 --multivariate --gwas --mv-method gemma-compat \
  --chr 1 --from-bp 1000 --to-bp 3000 --threads 2 --no-loco \
  --out "${OUT_DIR}/mv_gemma_compat" >"${OUT_DIR}/mv_gemma_compat.log" 2>&1
test -s "${OUT_DIR}/mv_gemma_compat.assoc.mv.txt"
test -s "${OUT_DIR}/mv_gemma_compat.mv.trace.tsv"
grep -q $'compatibility_target\tgemma_default_finite_em' "${OUT_DIR}/mv_gemma_compat.scan.meta.txt"
grep -q $'per_snp_max_iterations\t1000' "${OUT_DIR}/mv_gemma_compat.scan.meta.txt"
grep -q $'per_snp_convergence_tolerance\t0.001' "${OUT_DIR}/mv_gemma_compat.scan.meta.txt"
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

"${BIN}" --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 --multivariate --gwas --mv-method gemma-exact-full \
  --chr 1 --from-bp 1000 --to-bp 1000 --threads 2 --no-loco \
  --out "${OUT_DIR}/mv_gemma_exact_alias" >"${OUT_DIR}/mv_gemma_exact_alias.log" 2>&1
grep -q "deprecated experimental GEMMA-style mode" "${OUT_DIR}/mv_gemma_exact_alias.log"
grep -q $'requested_method\tgemma-exact-full' "${OUT_DIR}/mv_gemma_exact_alias.scan.meta.txt"
grep -q $'resolved_method\tgemma-exact-full' "${OUT_DIR}/mv_gemma_exact_alias.scan.meta.txt"
grep -q $'deprecated_alias\ttrue' "${OUT_DIR}/mv_gemma_exact_alias.scan.meta.txt"

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
