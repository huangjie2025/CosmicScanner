#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPT_DIR="${ROOT}/scripts/cfw_mouse_benchmark"
DATA_DIR="${DATA_DIR:-${ROOT}/tests/CFW_mouse}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/results/cfw_mouse_benchmark}"
DEBUG_DIR="${RESULT_DIR}/debug"
LOG_DIR="${RESULT_DIR}/logs/debug"
TABLE_DIR="${RESULT_DIR}/tables"
PLOT_DIR="${RESULT_DIR}/plots"
PROCESSED="${DATA_DIR}/processed"
PAPER_DIR="${PROCESSED}/paper_like"
STRICT_DIR="${PROCESSED}/strict_qc"

COSMICSCANNER_BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
PLINK2_BIN="${PLINK2_BIN:-/disk228/huangjie/software/miniconda3/envs/benchmark/bin/plink2}"
GEMMA_BIN="${GEMMA_BIN:-/disk228/huangjie/software/miniconda3/envs/benchmark/bin/gemma}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
THREADS="${THREADS:-1}"

mkdir -p "${DEBUG_DIR}" "${LOG_DIR}" "${TABLE_DIR}" "${PLOT_DIR}" "${PAPER_DIR}" "${STRICT_DIR}"

run_logged() {
  local name="$1"
  shift
  local log="${LOG_DIR}/${name}.log"
  printf '[cfw-debug] %s\n' "$*" | tee "${log}"
  "$@" >>"${log}" 2>&1
}

run_timed() {
  local name="$1"
  shift
  local log="${LOG_DIR}/${name}.log"
  local metric="${LOG_DIR}/${name}.time"
  /usr/bin/time -f 'seconds\t%e\nmax_rss_kb\t%M' -o "${metric}" "$@" >"${log}" 2>&1
}

RAW_PREFIX="${DATA_DIR}/cfw"
RAW_MISSING_PREFIX="${DEBUG_DIR}/raw_missing"
run_logged raw_missing "${PLINK2_BIN}" --bfile "${RAW_PREFIX}" --missing --out "${RAW_MISSING_PREFIX}"
run_logged summarize_missing "${PYTHON_BIN}" "${SCRIPT_DIR}/debug_cfw_mouse.py" missing --raw-prefix "${RAW_PREFIX}" --missing-prefix "${RAW_MISSING_PREFIX}" --results-dir "${RESULT_DIR}"

run_logged build_profiles "${PYTHON_BIN}" "${SCRIPT_DIR}/debug_cfw_mouse.py" profiles --data-dir "${DATA_DIR}" --results-dir "${RESULT_DIR}"

PAPER_1000_PREFIX="${PAPER_DIR}/cfw.paper_like.complete10.first1000"
if [[ ! -s "${PAPER_1000_PREFIX}.bed" ]]; then
  run_logged paper_1000_bfile "${PLINK2_BIN}" --bfile "${RAW_PREFIX}" \
    --keep "${PAPER_DIR}/mouse_10traits_complete_samples.txt" \
    --extract "${PAPER_DIR}/first_1000_raw_snps.txt" \
    --make-bed --out "${PAPER_1000_PREFIX}"
fi

STRICT_PREFIX="${PROCESSED}/cfw.qc"
if [[ -s "${STRICT_PREFIX}.bed" ]]; then
  ln -sfn "../cfw.qc.bed" "${STRICT_DIR}/cfw.strict_qc.bed"
  ln -sfn "../cfw.qc.bim" "${STRICT_DIR}/cfw.strict_qc.bim"
  ln -sfn "../cfw.qc.fam" "${STRICT_DIR}/cfw.strict_qc.fam"
fi

run_logged decode_emmv_z "${PYTHON_BIN}" "${SCRIPT_DIR}/debug_cfw_mouse.py" decode-emmv-z \
  --bfile "${PAPER_1000_PREFIX}" --keep "${PAPER_DIR}/mouse_10traits_complete_samples.txt" \
  --out "${PAPER_DIR}/cfw.paper_like.first1000.emmv.Z.csv" --max-snps 1000

KIN_DIR="${DEBUG_DIR}/kinship_univariate"
mkdir -p "${KIN_DIR}"
if [[ ! -s "${KIN_DIR}/output/paper_1000.cXX.txt" ]]; then
  run_logged gemma_kinship_univariate bash -c 'cd "$1" && "$2" -bfile "$3" -p "$4" -n 1 -gk 1 -o paper_1000' _ \
    "${KIN_DIR}" "${GEMMA_BIN}" "${PAPER_1000_PREFIX}" "${PAPER_DIR}/TA.complete.gemma.no_header.tsv"
fi
KINSHIP="${KIN_DIR}/output/paper_1000.cXX.txt"

mkdir -p "${DEBUG_DIR}/gemma_univariate" "${DEBUG_DIR}/gemma_2traits"

run_timed gemma_univariate_TA bash -c 'cd "$1" && "$2" -bfile "$3" -p "$4" -n 1 -k "$5" -lmm 4 -o debug_univariate_TA' _ \
  "${DEBUG_DIR}/gemma_univariate" "${GEMMA_BIN}" "${PAPER_1000_PREFIX}" "${PAPER_DIR}/TA.complete.gemma.no_header.tsv" "${KINSHIP}"

run_timed cosmic_univariate_TA "${COSMICSCANNER_BIN}" --bfile "${PAPER_1000_PREFIX}" \
  --pheno "${PAPER_DIR}/TA.complete.pheno.tsv" --mpheno 1 --gwas --no-loco --threads "${THREADS}" \
  --out "${DEBUG_DIR}/cosmic_univariate_TA"

run_logged join_univariate "${PYTHON_BIN}" "${SCRIPT_DIR}/debug_cfw_mouse.py" join-univariate \
  --gemma "${DEBUG_DIR}/gemma_univariate/output/debug_univariate_TA.assoc.txt" \
  --cosmic "${DEBUG_DIR}/cosmic_univariate_TA.assoc.txt" \
  --out "${TABLE_DIR}/debug_univariate_joined.tsv" \
  --summary "${TABLE_DIR}/debug_univariate_correlation.tsv" \
  --plot "${PLOT_DIR}/debug_univariate_cosmic_vs_gemma.png"

UNIV_R="$(${PYTHON_BIN} - "${TABLE_DIR}/debug_univariate_correlation.tsv" <<'PY'
import pandas as pd, sys
df = pd.read_csv(sys.argv[1], sep='\t')
row = df[df['metric'] == 'minus_log10_p']
print(float(row['pearson'].iloc[0]) if len(row) else float('nan'))
PY
)"

if "${PYTHON_BIN}" - "${UNIV_R}" <<'PY'
import math, sys
r = float(sys.argv[1])
raise SystemExit(0 if math.isfinite(r) and r >= 0.95 else 1)
PY
then
  run_timed gemma_2traits bash -c 'cd "$1" && "$2" -bfile "$3" -p "$4" -n 1 2 -k "$5" -lmm -o debug_2traits' _ \
    "${DEBUG_DIR}/gemma_2traits" "${GEMMA_BIN}" "${PAPER_1000_PREFIX}" "${PAPER_DIR}/traits_2.complete.gemma.no_header.tsv" "${KINSHIP}"
  run_timed emmv_2traits "${RSCRIPT_BIN}" "${SCRIPT_DIR}/run_emmvgwas_cfw.R" \
    "${PAPER_DIR}/traits_2.complete.pheno.tsv" "${PAPER_DIR}/cfw.paper_like.first1000.emmv.Z.csv" "${DEBUG_DIR}/emmv_2traits"
  run_timed cosmic_2traits "${COSMICSCANNER_BIN}" --bfile "${PAPER_1000_PREFIX}" \
    --pheno "${PAPER_DIR}/traits_2.complete.pheno.tsv" --mpheno-list 1,2 \
    --multivariate --mv-method emmv --mv-compatible emmvgwas --gwas --no-loco --threads "${THREADS}" \
    --out "${DEBUG_DIR}/cosmic_2traits"
  run_logged join_2traits "${PYTHON_BIN}" "${SCRIPT_DIR}/debug_cfw_mouse.py" join-multi \
    --gemma "${DEBUG_DIR}/gemma_2traits/output/debug_2traits.assoc.txt" \
    --emmv "${DEBUG_DIR}/emmv_2traits_association_with_ids.csv" \
    --cosmic "${DEBUG_DIR}/cosmic_2traits.assoc.mv.txt" \
    --out "${TABLE_DIR}/debug_2traits_joined.tsv" \
    --summary "${TABLE_DIR}/debug_2traits_correlation.tsv" \
    --plot-emmv "${PLOT_DIR}/debug_2traits_emmvgwas_vs_gemma.png" \
    --plot-cosmic "${PLOT_DIR}/debug_2traits_cosmic_vs_gemma.png"
else
  printf 'Univariate CosmicScanner vs GEMMA -log10(p) Pearson r=%s; stopping before multivariate debug.\n' "${UNIV_R}" | tee "${DEBUG_DIR}/stop_after_univariate.txt"
fi

if [[ -s "${DEBUG_DIR}/gemma_2traits/output/debug_2traits.assoc.txt" && -s "${DEBUG_DIR}/cosmic_2traits.assoc.mv.txt" ]]; then
  run_logged schema_dump "${PYTHON_BIN}" "${SCRIPT_DIR}/debug_cfw_mouse.py" schema --results-dir "${RESULT_DIR}" \
    --gemma "${DEBUG_DIR}/gemma_2traits/output/debug_2traits.assoc.txt" \
    --emmv "${DEBUG_DIR}/emmv_2traits_association_with_ids.csv" \
    --cosmic "${DEBUG_DIR}/cosmic_2traits.assoc.mv.txt"
  run_logged alignment "${PYTHON_BIN}" "${SCRIPT_DIR}/debug_cfw_mouse.py" alignment --results-dir "${RESULT_DIR}" \
    --raw-prefix "${RAW_PREFIX}" --debug-prefix "${PAPER_1000_PREFIX}" \
    --gemma "${DEBUG_DIR}/gemma_2traits/output/debug_2traits.assoc.txt" \
    --emmv "${DEBUG_DIR}/emmv_2traits_association_with_ids.csv" \
    --cosmic "${DEBUG_DIR}/cosmic_2traits.assoc.mv.txt"
fi

"${PYTHON_BIN}" - "${RESULT_DIR}" <<'PY'
from pathlib import Path
import pandas as pd
res = Path(__import__('sys').argv[1])
def txt(p): return p.read_text(encoding='utf-8') if p.exists() else ''
univ = txt(res / 'tables/debug_univariate_correlation.tsv')
multi = txt(res / 'tables/debug_2traits_correlation.tsv')
missing = txt(res / 'tables/missingness_summary.tsv')
lines = ['# CFW Mouse Debug Summary', '', '## Missingness', '', '```tsv', missing.strip(), '```', '', '## Univariate TA Debug', '', '```tsv', univ.strip(), '```', '', '## 2-Trait Debug', '', '```tsv', multi.strip() if multi else 'not run; see debug/stop_after_univariate.txt', '```', '', '## Current Decision', '']
if multi:
    lines.append('- Single-trait debug passed the configured threshold and 2-trait debug was run. Inspect correlations before any full benchmark.')
else:
    lines.append('- Single-trait debug did not pass the configured threshold or failed; do not enter full benchmark.')
lines += ['', '## Statistical Columns', '', '- GEMMA multivariate debug uses `p_wald` from GEMMA `-lmm` output.', '- EMmvGWAS debug uses `p` from `gwas_all()` output, with SNP IDs carried in `SNP_ID`.', '- CosmicScanner multivariate debug uses `P` from `.assoc.mv.txt`, produced by `--mv-method emmv --mv-compatible emmvgwas`.', '- If correlations remain low, these outputs should be treated as not yet calibrated even though they are all Wald-like columns.']
(res / 'debug_summary.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
report = res / 'CFW_mouse_benchmark_report.md'
prev = txt(report)
append = ['','## Debug Calibration Update','', '- Smoke test success is not equivalent to paper reproduction.', f"- Debug summary: `{res / 'debug_summary.md'}`", f"- Preprocessing profiles: `{res / 'preprocessing_profiles.md'}`", f"- Trait matching debug: `{res / 'trait_matching_debug.md'}`", f"- Alignment debug: `{res / 'alignment_debug.md'}`", '- Full benchmark remains blocked until single-trait and 2-trait calibration are both satisfactory.']
report.write_text(prev.rstrip() + '\n' + '\n'.join(append) + '\n', encoding='utf-8')
PY

echo "CFW mouse debug pipeline complete: ${RESULT_DIR}"
