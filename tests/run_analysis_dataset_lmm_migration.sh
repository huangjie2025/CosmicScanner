#!/usr/bin/env bash
set -euo pipefail

ROOT="${COSMICSCANNER_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
OUT="${COSMIC_ANALYSIS_DATASET_LMM_OUT:-${ROOT}/build/test-results/analysis_dataset_lmm}"
DATA="${ROOT}/tests/data/gemma_accuracy/gemma_accuracy"
EBFMT_DATA="${ROOT}/tests/data/eb_format_consistency/ebfmt"
PGEN_DATA="${ROOT}/tests/results/eb_format_consistency/pgen/ebfmt"
BGEN_DATA="${ROOT}/tests/results/eb_format_consistency/bgen/ebfmt"

mkdir -p "${OUT}"
rm -f "${OUT}"/*

PHENO_ORIG="${DATA}.pheno"
PHENO_SHUF="${OUT}/gemma_accuracy.shuffled.pheno"
PHENO_MISS="${OUT}/gemma_accuracy.missing.pheno"
PHENO_EXTRA="${OUT}/gemma_accuracy.extra.pheno"
PHENO_DROP="${OUT}/gemma_accuracy.drop.pheno"
QCOVAR_ORIG="${OUT}/gemma_accuracy.qcovar"
QCOVAR_SHUF="${OUT}/gemma_accuracy.shuffled.qcovar"
COVAR_ORIG="${OUT}/gemma_accuracy.covar"
COVAR_SHUF="${OUT}/gemma_accuracy.shuffled.covar"
PHENO_INF="${OUT}/gemma_accuracy.inf.pheno"
EB_PHENO_ORIG="${EBFMT_DATA}.pheno"
EB_PHENO_SHUF="${OUT}/ebfmt.shuffled.pheno"

tac "${PHENO_ORIG}" > "${PHENO_SHUF}"
awk 'NR<=3 {$3=-9} {print $1"\t"$2"\t"$3}' "${PHENO_ORIG}" > "${PHENO_MISS}"
awk '{print $1"\t"$2"\t"(NR == 1 ? "inf" : $3)}' "${PHENO_ORIG}" > "${PHENO_INF}"
{ cat "${PHENO_ORIG}"; printf "F9999\tI9999\t1.0\n"; } > "${PHENO_EXTRA}"
awk 'NR>1 {print}' "${PHENO_ORIG}" > "${PHENO_DROP}"
awk '{print $1"\t"$2"\t"(NR % 7)}' "${PHENO_ORIG}" > "${QCOVAR_ORIG}"
tac "${QCOVAR_ORIG}" > "${QCOVAR_SHUF}"
awk '{print $1"\t"$2"\t"((NR % 2) ? "A" : "B")}' "${PHENO_ORIG}" > "${COVAR_ORIG}"
tac "${COVAR_ORIG}" > "${COVAR_SHUF}"
tac "${EB_PHENO_ORIG}" > "${EB_PHENO_SHUF}"

run_lmm() {
    local prefix="$1"
    shift
    "${BIN}" --bfile "${DATA}" --pheno "$1" --gwas --no-loco --out "${prefix}" "${@:2}" \
        > "${prefix}.stdout" 2> "${prefix}.stderr"
}

run_bed_ebfmt_lmm() {
    local prefix="$1"
    local pheno="$2"
    "${BIN}" --bfile "${EBFMT_DATA}" --pheno "${pheno}" --mpheno 1 --gwas --no-loco --maf 0.0 --out "${prefix}" \
        > "${prefix}.stdout" 2> "${prefix}.stderr"
}

run_pgen_lmm() {
    local prefix="$1"
    local pheno="$2"
    "${BIN}" --pfile "${PGEN_DATA}" --pheno "${pheno}" --mpheno 1 --gwas --no-loco --maf 0.0 --out "${prefix}" \
        > "${prefix}.stdout" 2> "${prefix}.stderr"
}

run_bgen_lmm() {
    local prefix="$1"
    local pheno="$2"
    "${BIN}" --bgen "${BGEN_DATA}.bgen" --sample "${BGEN_DATA}.sample" --pheno "${pheno}" --mpheno 1 --gwas --no-loco --maf 0.0 --out "${prefix}" \
        > "${prefix}.stdout" 2> "${prefix}.stderr"
}

check_manifest() {
    local prefix="$1"
    local component="$2"
    local expected_format="$3"
    local expected_budget="$4"
    python3 - "${prefix}.run_manifest.json" "${component}" "${expected_format}" "${expected_budget}" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
component = sys.argv[2]
expected_format = sys.argv[3]
expected_budget = sys.argv[4]
if not path.exists() or path.stat().st_size == 0:
    raise SystemExit(f"missing manifest: {path}")
doc = json.loads(path.read_text(encoding="utf-8"))
required = ["schema_version", "software", "run", "inputs", "data", "model", "solver", "validation", "performance", "outputs"]
missing = [k for k in required if k not in doc]
if missing:
    raise SystemExit(f"manifest missing keys: {missing}")
if doc["schema_version"] != "1":
    raise SystemExit("schema_version mismatch")
if doc["software"]["component"] != component:
    raise SystemExit("component mismatch")
if doc["data"]["genotype"]["format"] != expected_format:
    raise SystemExit("genotype format mismatch")
if doc["solver"]["diagnostic_only"] is not True or doc["solver"]["planner_mode"] != "diagnostic-only":
    raise SystemExit("solver planner mode mismatch")
if expected_budget == "unknown":
    if doc["solver"]["memory_budget_bytes"] is not None:
        raise SystemExit("unknown memory budget serialized as non-null")
    if doc["solver"]["memory_budget_source"] != "UNKNOWN":
        raise SystemExit("unknown memory source mismatch")
else:
    if doc["solver"]["memory_budget_bytes"] != int(expected_budget):
        raise SystemExit("explicit memory budget did not propagate")
    if doc["solver"]["memory_budget_source"] != "USER_SPECIFIED":
        raise SystemExit("explicit memory source mismatch")
if doc["validation"]["status"] not in {"PASS", "WARNING"}:
    raise SystemExit("successful manifest has invalid validation status")
out_paths = {item["type"]: Path(item["path"]) for item in doc["outputs"]}
if "association" not in out_paths or not out_paths["association"].exists():
    raise SystemExit("association output inventory mismatch")
PY
}

run_lmm "${OUT}/baseline" "${PHENO_ORIG}"
check_manifest "${OUT}/baseline" "CosmicScanner" "PLINK_BED" "unknown"
run_lmm "${OUT}/baseline_mem" "${PHENO_ORIG}" --memory-gib 1
cmp -s "${OUT}/baseline.assoc.txt" "${OUT}/baseline_mem.assoc.txt"
check_manifest "${OUT}/baseline_mem" "CosmicScanner" "PLINK_BED" "1073741824"
run_lmm "${OUT}/pheno_shuffled" "${PHENO_SHUF}"
cmp -s "${OUT}/baseline.assoc.txt" "${OUT}/pheno_shuffled.assoc.txt"
grep -q "Validation : PASS" "${OUT}/baseline.stdout"
grep -q "AnalysisDataset production-contract mode" "${OUT}/baseline.stdout"

run_lmm "${OUT}/pheno_extra" "${PHENO_EXTRA}"
cmp -s "${OUT}/baseline.assoc.txt" "${OUT}/pheno_extra.assoc.txt"

run_lmm "${OUT}/covar_baseline" "${PHENO_ORIG}" --qcovar "${QCOVAR_ORIG}" --covar "${COVAR_ORIG}"
run_lmm "${OUT}/covar_shuffled" "${PHENO_ORIG}" --qcovar "${QCOVAR_SHUF}" --covar "${COVAR_SHUF}"
cmp -s "${OUT}/covar_baseline.assoc.txt" "${OUT}/covar_shuffled.assoc.txt"
grep -q "Using 3 fixed-effect columns (including intercept)." "${OUT}/covar_baseline.stdout"

run_lmm "${OUT}/missing" "${PHENO_MISS}"
grep -q "Validation : WARNING" "${OUT}/missing.stdout"
grep -q "\[IO_MISSING_PHENOTYPE\]" "${OUT}/missing.stdout"
grep -q "After filtering missing data, 157 individuals remain." "${OUT}/missing.stdout"
run_lmm "${OUT}/dropped_row" "${PHENO_DROP}"
grep -q "After filtering missing data, 159 individuals remain." "${OUT}/dropped_row.stdout"

set +e
"${BIN}" --bfile "${DATA}" --pheno "${PHENO_INF}" --gwas --no-loco --out "${OUT}/nonfinite_pheno" \
    > "${OUT}/nonfinite_pheno.stdout" 2> "${OUT}/nonfinite_pheno.stderr"
status=$?
set -e
if [[ ${status} -eq 0 ]]; then
    echo "non-finite phenotype validation unexpectedly succeeded" >&2
    exit 1
fi
cat "${OUT}/nonfinite_pheno.stdout" "${OUT}/nonfinite_pheno.stderr" > "${OUT}/nonfinite_pheno.combined"
grep -q "Validation : FAIL" "${OUT}/nonfinite_pheno.combined"
grep -q "\[PHENOTYPE_NONFINITE\]" "${OUT}/nonfinite_pheno.combined"

run_bed_ebfmt_lmm "${OUT}/ebfmt_bed" "${EB_PHENO_ORIG}"
check_manifest "${OUT}/ebfmt_bed" "CosmicScanner" "PLINK_BED" "unknown"
grep -q "Validation : PASS" "${OUT}/ebfmt_bed.stdout"
grep -q "AnalysisDataset production-contract mode" "${OUT}/ebfmt_bed.stdout"
grep -Eq "Genotype format[[:space:]]*: PLINK_BED" "${OUT}/ebfmt_bed.stdout"
grep -Eq "Samples[[:space:]]*: 100" "${OUT}/ebfmt_bed.stdout"
grep -Eq "Markers[[:space:]]*: 280" "${OUT}/ebfmt_bed.stdout"
grep -Eq "Genotype materialized[[:space:]]*: NO" "${OUT}/ebfmt_bed.stdout"

if [[ -f "${PGEN_DATA}.pgen" && -f "${PGEN_DATA}.psam" && -f "${PGEN_DATA}.pvar" ]]; then
    run_pgen_lmm "${OUT}/ebfmt_pgen" "${EB_PHENO_ORIG}"
    check_manifest "${OUT}/ebfmt_pgen" "CosmicScanner" "PLINK_PGEN" "unknown"
    run_pgen_lmm "${OUT}/ebfmt_pgen_shuffled" "${EB_PHENO_SHUF}"
    cmp -s "${OUT}/ebfmt_pgen.assoc.txt" "${OUT}/ebfmt_pgen_shuffled.assoc.txt"
    grep -q "Validation : PASS" "${OUT}/ebfmt_pgen.stdout"
    grep -q "AnalysisDataset production-contract mode" "${OUT}/ebfmt_pgen.stdout"
    grep -Eq "Genotype format[[:space:]]*: PLINK_PGEN" "${OUT}/ebfmt_pgen.stdout"
    grep -Eq "Samples[[:space:]]*: 100" "${OUT}/ebfmt_pgen.stdout"
    grep -Eq "Markers[[:space:]]*: 280" "${OUT}/ebfmt_pgen.stdout"
    grep -Eq "Genotype materialized[[:space:]]*: NO" "${OUT}/ebfmt_pgen.stdout"
fi

if [[ -f "${BGEN_DATA}.bgen" && -f "${BGEN_DATA}.sample" ]]; then
    run_bgen_lmm "${OUT}/ebfmt_bgen" "${EB_PHENO_ORIG}"
    check_manifest "${OUT}/ebfmt_bgen" "CosmicScanner" "BGEN" "unknown"
    run_bgen_lmm "${OUT}/ebfmt_bgen_shuffled" "${EB_PHENO_SHUF}"
    cmp -s "${OUT}/ebfmt_bgen.assoc.txt" "${OUT}/ebfmt_bgen_shuffled.assoc.txt"
    grep -q "Validation : PASS" "${OUT}/ebfmt_bgen.stdout"
    grep -q "AnalysisDataset production-contract mode" "${OUT}/ebfmt_bgen.stdout"
    grep -Eq "Genotype format[[:space:]]*: BGEN" "${OUT}/ebfmt_bgen.stdout"
    grep -Eq "Samples[[:space:]]*: 100" "${OUT}/ebfmt_bgen.stdout"
    grep -Eq "Markers[[:space:]]*: 280" "${OUT}/ebfmt_bgen.stdout"
    grep -Eq "Genotype materialized[[:space:]]*: NO" "${OUT}/ebfmt_bgen.stdout"
fi

echo "analysis_dataset_lmm_migration: passed"
