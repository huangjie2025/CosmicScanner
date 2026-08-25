#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPT_DIR="${ROOT}/scripts/cfw_mouse_benchmark"
DATA_DIR="${DATA_DIR:-${ROOT}/tests/CFW_mouse}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/results/cfw_mouse_benchmark}"
OUT_DIR="${RESULT_DIR}/scaled_harmonized"
PAPER_DIR="${DATA_DIR}/processed/paper_like"
PLINK2_BIN="${PLINK2_BIN:-/disk228/huangjie/software/miniconda3/envs/benchmark/bin/plink2}"
GEMMA_BIN="${GEMMA_BIN:-/disk228/huangjie/software/miniconda3/envs/benchmark/bin/gemma}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript}"
COSMICSCANNER_BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
THREADS="${THREADS:-1}"
RUN_LOWMISS_5K="${RUN_LOWMISS_5K:-0}"
TRAIT_LIST="${TRAIT_LIST:-2 4 10}"
METHOD_TIMEOUT="${METHOD_TIMEOUT:-900}"

mkdir -p "${OUT_DIR}" "${RESULT_DIR}/plots" "${OUT_DIR}/logs"

run_logged() {
  local name="$1"; shift
  local log="${OUT_DIR}/logs/${name}.log"
  printf '[scaled-harmonized] %s\n' "$*" | tee "${log}"
  "$@" >>"${log}" 2>&1
}

run_timed() {
  local log_dir="$1"; local name="$2"; shift 2
  mkdir -p "${log_dir}"
  /usr/bin/time -f 'seconds\t%e\nmax_rss_kb\t%M' -o "${log_dir}/${name}.time" timeout "${METHOD_TIMEOUT}" "$@" >"${log_dir}/${name}.log" 2>&1
}

RAW_MISSING_PREFIX="${OUT_DIR}/raw_missing"
if [[ ! -s "${RAW_MISSING_PREFIX}.vmiss" ]]; then
  run_logged raw_missing "${PLINK2_BIN}" --bfile "${DATA_DIR}/cfw" --missing --out "${RAW_MISSING_PREFIX}"
fi

build_set() {
  local set_name="$1"; local n_snps="$2"; local max_miss="$3"
  if [[ ! -s "${OUT_DIR}/${set_name}/manifest.json" ]]; then
    run_logged "build_${set_name}" "${PYTHON_BIN}" "${SCRIPT_DIR}/scaled_harmonized.py" build-set \
      --raw-prefix "${DATA_DIR}/cfw" --vmiss "${RAW_MISSING_PREFIX}.vmiss" \
      --samples "${PAPER_DIR}/mouse_10traits_complete_samples.txt" \
      --pheno "${PAPER_DIR}/cfw_mouse_traits.tsv" --out-dir "${OUT_DIR}" \
      --set-name "${set_name}" --n-snps "${n_snps}" --max-missing-rate "${max_miss}" --min-maf 0.01
  fi
}

run_combo() {
  local set_name="$1"; local traits="$2"
  local set_dir="${OUT_DIR}/${set_name}"
  local combo_dir="${set_dir}/traits_${traits}"
  local res_dir="${combo_dir}/results"
  local log_dir="${combo_dir}/logs"
  mkdir -p "${res_dir}" "${log_dir}" "${res_dir}/gemma" "${res_dir}/emmvgwas" "${res_dir}/cosmicscanner_gemma_full" "${res_dir}/cosmicscanner_emmv"
  local prefix="${set_dir}/plink"
  local kin_dir="${res_dir}/gemma/kinship"
  mkdir -p "${kin_dir}"
  run_timed "${log_dir}" gemma_kinship bash -c 'cd "$1" && "$2" -bfile "$3" -p "$4" -n 1 -gk 1 -o kinship' _ \
    "${kin_dir}" "${GEMMA_BIN}" "${prefix}" "${combo_dir}/gemma.no_header.tsv"
  local kinship="${kin_dir}/output/kinship.cXX.txt"
  run_timed "${log_dir}" gemma bash -c 'cd "$1" && "$2" -bfile "$3" -p "$4" -n $5 -k "$6" -lmm -o gemma' _ \
    "${res_dir}/gemma" "${GEMMA_BIN}" "${prefix}" "${combo_dir}/gemma.no_header.tsv" "$(seq -s' ' 1 "${traits}")" "${kinship}"
  run_timed "${log_dir}" emmvgwas "${RSCRIPT_BIN}" "${SCRIPT_DIR}/run_emmvgwas_cfw.R" \
    "${combo_dir}/pheno.tsv" "${set_dir}/emmv.Z.csv" "${res_dir}/emmvgwas/emmv"
  run_timed "${log_dir}" cosmic_gemma_full "${COSMICSCANNER_BIN}" --bfile "${prefix}" \
    --pheno "${combo_dir}/pheno.tsv" --mpheno-list "$(seq -s, 1 "${traits}")" \
    --multivariate --mv-method gemma-full --mv-test wald,lrt,score --gwas --no-loco --threads "${THREADS}" \
    --out "${res_dir}/cosmicscanner_gemma_full/cosmic_gemma_full"
  run_timed "${log_dir}" cosmic_emmv "${COSMICSCANNER_BIN}" --bfile "${prefix}" \
    --pheno "${combo_dir}/pheno.tsv" --mpheno-list "$(seq -s, 1 "${traits}")" \
    --multivariate --mv-method emmv --mv-compatible emmvgwas --gwas --no-loco --threads "${THREADS}" \
    --out "${res_dir}/cosmicscanner_emmv/cosmic_emmv"
  local input_snps
  input_snps="$(${PYTHON_BIN} - "${set_dir}/manifest.json" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))['snps'])
PY
)"
  run_logged "summarize_${set_name}_t${traits}" "${PYTHON_BIN}" "${SCRIPT_DIR}/scaled_harmonized.py" summarize-combo \
    --out-dir "${res_dir}" --set-name "${set_name}" --traits "${traits}" --input-snps "${input_snps}" \
    --gemma "${res_dir}/gemma/output/gemma.assoc.txt" \
    --emmv "${res_dir}/emmvgwas/emmv_association_with_ids.csv" \
    --cosmic-gemma-full "${res_dir}/cosmicscanner_gemma_full/cosmic_gemma_full.assoc.mv.txt" \
    --cosmic-emmv "${res_dir}/cosmicscanner_emmv/cosmic_emmv.assoc.mv.txt"
}

build_set missing0_1k 1000 0
for t in ${TRAIT_LIST}; do
  run_combo missing0_1k "${t}"
done

"${PYTHON_BIN}" "${SCRIPT_DIR}/scaled_harmonized.py" collect --out-dir "${OUT_DIR}"

if [[ "${RUN_LOWMISS_5K}" == "1" ]]; then
  build_set lowmiss_5k 5000 0.01
  for t in ${TRAIT_LIST}; do
    run_combo lowmiss_5k "${t}"
  done
  "${PYTHON_BIN}" "${SCRIPT_DIR}/scaled_harmonized.py" collect --out-dir "${OUT_DIR}"
fi

echo "Scaled harmonized benchmark complete: ${OUT_DIR}"
