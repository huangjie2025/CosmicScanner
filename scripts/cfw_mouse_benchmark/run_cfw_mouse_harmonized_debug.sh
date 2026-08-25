#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPT_DIR="${ROOT}/scripts/cfw_mouse_benchmark"
DATA_DIR="${DATA_DIR:-${ROOT}/tests/CFW_mouse}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/results/cfw_mouse_benchmark}"
OUT_DIR="${RESULT_DIR}/debug/harmonized_2traits"
LOG_DIR="${RESULT_DIR}/logs/harmonized_2traits"
PLOT_DIR="${RESULT_DIR}/plots"
PAPER_DIR="${DATA_DIR}/processed/paper_like"
PLINK2_BIN="${PLINK2_BIN:-/disk228/huangjie/software/miniconda3/envs/benchmark/bin/plink2}"
GEMMA_BIN="${GEMMA_BIN:-/disk228/huangjie/software/miniconda3/envs/benchmark/bin/gemma}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript}"
COSMICSCANNER_BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
THREADS="${THREADS:-1}"

mkdir -p "${OUT_DIR}" "${LOG_DIR}" "${PLOT_DIR}" "${OUT_DIR}/gemma" "${OUT_DIR}/emmvgwas" "${OUT_DIR}/cosmicscanner_gemma_full" "${OUT_DIR}/cosmicscanner_emmv" "${OUT_DIR}/mv_p1"

run_logged() {
  local name="$1"; shift
  local log="${LOG_DIR}/${name}.log"
  printf '[harmonized-debug] %s\n' "$*" | tee "${log}"
  "$@" >>"${log}" 2>&1
}

run_timed() {
  local name="$1"; shift
  /usr/bin/time -f 'seconds\t%e\nmax_rss_kb\t%M' -o "${LOG_DIR}/${name}.time" "$@" >"${LOG_DIR}/${name}.log" 2>&1
}

RAW_MISSING_PREFIX="${OUT_DIR}/raw_missing"
if [[ ! -s "${RAW_MISSING_PREFIX}.vmiss" ]]; then
  run_logged raw_missing "${PLINK2_BIN}" --bfile "${DATA_DIR}/cfw" --missing --out "${RAW_MISSING_PREFIX}"
fi

if [[ ! -s "${PAPER_DIR}/mouse_10traits_complete_samples.txt" ]]; then
  echo "Missing paper_like samples; run scripts/cfw_mouse_benchmark/run_cfw_mouse_debug.sh first." >&2
  exit 2
fi

run_logged build_harmonized "${PYTHON_BIN}" "${SCRIPT_DIR}/harmonized_debug.py" select \
  --raw-prefix "${DATA_DIR}/cfw" --vmiss "${RAW_MISSING_PREFIX}.vmiss" \
  --samples "${PAPER_DIR}/mouse_10traits_complete_samples.txt" \
  --pheno "${PAPER_DIR}/cfw_mouse_traits.tsv" --out-dir "${OUT_DIR}" --n-snps 200

H_PREFIX="${OUT_DIR}/harmonized"
KIN_DIR="${OUT_DIR}/gemma/kinship"
mkdir -p "${KIN_DIR}"
run_logged gemma_kinship bash -c 'cd "$1" && "$2" -bfile "$3" -p "$4" -n 1 -gk 1 -o harmonized_kinship' _ \
  "${KIN_DIR}" "${GEMMA_BIN}" "${H_PREFIX}" "${OUT_DIR}/harmonized_pheno_TA.gemma.no_header.tsv"
KINSHIP="${KIN_DIR}/output/harmonized_kinship.cXX.txt"

run_timed gemma_2traits bash -c 'cd "$1" && "$2" -bfile "$3" -p "$4" -n 1 2 -k "$5" -lmm -o harmonized_2traits' _ \
  "${OUT_DIR}/gemma" "${GEMMA_BIN}" "${H_PREFIX}" "${OUT_DIR}/harmonized_pheno_TA_EDL.gemma.no_header.tsv" "${KINSHIP}"

run_timed emmvgwas_2traits "${RSCRIPT_BIN}" "${SCRIPT_DIR}/run_emmvgwas_cfw.R" \
  "${OUT_DIR}/harmonized_pheno_TA_EDL.tsv" "${OUT_DIR}/harmonized_emmv.Z.csv" "${OUT_DIR}/emmvgwas/harmonized_2traits"

run_timed cosmic_gemma_full "${COSMICSCANNER_BIN}" --bfile "${H_PREFIX}" \
  --pheno "${OUT_DIR}/harmonized_pheno_TA_EDL.tsv" --mpheno-list 1,2 \
  --multivariate --mv-method gemma-full --mv-test wald,lrt,score --gwas --no-loco --threads "${THREADS}" \
  --out "${OUT_DIR}/cosmicscanner_gemma_full/harmonized_2traits"

run_timed cosmic_emmv "${COSMICSCANNER_BIN}" --bfile "${H_PREFIX}" \
  --pheno "${OUT_DIR}/harmonized_pheno_TA_EDL.tsv" --mpheno-list 1,2 \
  --multivariate --mv-method emmv --mv-compatible emmvgwas --gwas --no-loco --threads "${THREADS}" \
  --out "${OUT_DIR}/cosmicscanner_emmv/harmonized_2traits"

run_logged compare_2traits "${PYTHON_BIN}" "${SCRIPT_DIR}/harmonized_debug.py" compare --out-dir "${OUT_DIR}" \
  --gemma "${OUT_DIR}/gemma/output/harmonized_2traits.assoc.txt" \
  --emmv "${OUT_DIR}/emmvgwas/harmonized_2traits_association_with_ids.csv" \
  --cosmic-gemma-full "${OUT_DIR}/cosmicscanner_gemma_full/harmonized_2traits.assoc.mv.txt" \
  --cosmic-emmv "${OUT_DIR}/cosmicscanner_emmv/harmonized_2traits.assoc.mv.txt" \
  --plot "${PLOT_DIR}/harmonized_2traits_comparison.png"

run_timed gemma_univariate_TA bash -c 'cd "$1" && "$2" -bfile "$3" -p "$4" -n 1 -k "$5" -lmm 4 -o harmonized_TA' _ \
  "${OUT_DIR}/mv_p1" "${GEMMA_BIN}" "${H_PREFIX}" "${OUT_DIR}/harmonized_pheno_TA.gemma.no_header.tsv" "${KINSHIP}"

run_timed cosmic_univariate_TA "${COSMICSCANNER_BIN}" --bfile "${H_PREFIX}" \
  --pheno "${OUT_DIR}/harmonized_pheno_TA.tsv" --mpheno 1 --gwas --no-loco --threads "${THREADS}" \
  --out "${OUT_DIR}/mv_p1/cosmic_univariate_TA"

run_timed cosmic_mv_p1_TA "${COSMICSCANNER_BIN}" --bfile "${H_PREFIX}" \
  --pheno "${OUT_DIR}/harmonized_pheno_TA.tsv" --mpheno-list 1 \
  --multivariate --mv-method gemma-full --mv-test wald,lrt,score --gwas --no-loco --threads "${THREADS}" \
  --out "${OUT_DIR}/mv_p1/cosmic_mv_p1_TA"

run_logged compare_p1 "${PYTHON_BIN}" "${SCRIPT_DIR}/harmonized_debug.py" p1 --out-dir "${OUT_DIR}" \
  --gemma "${OUT_DIR}/mv_p1/output/harmonized_TA.assoc.txt" \
  --cosmic-uni "${OUT_DIR}/mv_p1/cosmic_univariate_TA.assoc.txt" \
  --cosmic-mv-p1 "${OUT_DIR}/mv_p1/cosmic_mv_p1_TA.assoc.mv.txt" \
  --plot "${PLOT_DIR}/mv_p1_parity.png"

run_logged final_report "${PYTHON_BIN}" "${SCRIPT_DIR}/harmonized_debug.py" report --out-dir "${OUT_DIR}" --results-dir "${RESULT_DIR}"

echo "Harmonized 2-trait debug complete: ${OUT_DIR}"
