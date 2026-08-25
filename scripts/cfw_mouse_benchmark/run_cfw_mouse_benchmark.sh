#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPT_DIR="${ROOT}/scripts/cfw_mouse_benchmark"
DATA_DIR="${DATA_DIR:-${ROOT}/tests/CFW_mouse}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/results/cfw_mouse_benchmark}"
MODE="${MODE:-smoke}"
RUN_METHODS="${RUN_METHODS:-1}"
RUN_PLOTS="${RUN_PLOTS:-1}"
RUN_NULL="${RUN_NULL:-0}"
THREADS="${THREADS:-1}"
ALLOW_QC_FALLBACK="${ALLOW_QC_FALLBACK:-1}"

COSMICSCANNER_BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
PLINK2_BIN="${PLINK2_BIN:-/disk228/huangjie/software/miniconda3/envs/benchmark/bin/plink2}"
GEMMA_BIN="${GEMMA_BIN:-/disk228/huangjie/software/miniconda3/envs/benchmark/bin/gemma}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

LOG_DIR="${RESULT_DIR}/logs"
mkdir -p "${LOG_DIR}" "${RESULT_DIR}/gemma" "${RESULT_DIR}/emmvgwas" "${RESULT_DIR}/cosmicscanner" "${RESULT_DIR}/plots" "${RESULT_DIR}/tables"

RAW_SPACE_DIR="${ROOT}/tests/CFW mouse"
if [[ ! -e "${DATA_DIR}" && -d "${RAW_SPACE_DIR}" ]]; then
  ln -s "CFW mouse" "${ROOT}/tests/CFW_mouse"
fi

run_logged() {
  local name="$1"
  shift
  local log="${LOG_DIR}/${name}.log"
  printf '[cfw-benchmark] %s\n' "$*" | tee "${log}"
  "$@" >>"${log}" 2>&1
}

run_timed() {
  local method="$1"
  local traits="$2"
  local prefix="$3"
  shift 3
  local log="${LOG_DIR}/${method}_t${traits}.log"
  local metric="${LOG_DIR}/${method}_t${traits}.time"
  set +e
  /usr/bin/time -f 'seconds\t%e\nmax_rss_kb\t%M' -o "${metric}" "$@" >"${log}" 2>&1
  local rc=$?
  set -e
  "${PYTHON_BIN}" - "${RESULT_DIR}/tables/runtime.tsv" "${method}" "${traits}" "${prefix}" "${log}" "${metric}" "${rc}" <<'PY'
import sys
from pathlib import Path
table, method, traits, prefix, log, metric, rc = sys.argv[1:]
vals = {"seconds": "nan", "max_rss_kb": "nan"}
p = Path(metric)
if p.exists():
    for line in p.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) == 2:
            vals[parts[0]] = parts[1]
out = Path(table)
if not out.exists():
    out.write_text("method\ttraits\tseconds\tmax_rss_kb\tstatus\tprefix\tlog\n", encoding="utf-8")
with out.open("a", encoding="utf-8") as f:
    f.write(f"{method}\t{traits}\t{vals['seconds']}\t{vals['max_rss_kb']}\t{rc}\t{prefix}\t{log}\n")
PY
  if [[ "${rc}" != "0" ]]; then
    echo "Method failed: ${method} traits=${traits}. See ${log}" >&2
    return "${rc}"
  fi
}

tool_status_json="${RESULT_DIR}/tables/tool_status.json"
"${PYTHON_BIN}" - "${tool_status_json}" "${COSMICSCANNER_BIN}" "${PLINK2_BIN}" "${GEMMA_BIN}" "${RSCRIPT_BIN}" "${PYTHON_BIN}" <<'PY'
import json, shutil, subprocess, sys
from pathlib import Path
out, cosmic, plink2, gemma, rscript, python = sys.argv[1:]
def status(path, probe=None):
    exe = path if Path(path).exists() else shutil.which(path)
    if not exe:
        return "missing"
    if probe:
        try:
            txt = subprocess.run([exe] + probe, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=20).stdout.splitlines()[0]
            return f"{exe} ({txt})"
        except Exception as e:
            return f"{exe} (probe_failed: {e})"
    return exe
d = {
  "CosmicScanner": status(cosmic, ["--help"]),
  "plink": status("plink", ["--version"]),
  "plink2": status(plink2, ["--version"]),
  "gemma": status(gemma, ["-h"]),
  "Rscript": status(rscript, ["--version"]),
  "python3": status(python, ["--version"]),
}
exe = rscript if Path(rscript).exists() else shutil.which(rscript)
if exe:
    try:
        d["R_EMmvGWAS"] = subprocess.run([exe, "-e", "cat(requireNamespace('EMmvGWAS', quietly=TRUE),'\\n')"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30).stdout.strip()
    except Exception as e:
        d["R_EMmvGWAS"] = f"probe_failed: {e}"
Path(out).parent.mkdir(parents=True, exist_ok=True)
Path(out).write_text(json.dumps(d, indent=2), encoding="utf-8")
PY

run_logged audit "${PYTHON_BIN}" "${SCRIPT_DIR}/prepare_cfw_mouse.py" \
  --data-dir "${DATA_DIR}" --qc-prefix "${DATA_DIR}/cfw" --results-dir "${RESULT_DIR}" \
  --tool-status-json "${tool_status_json}" --audit-only

PROCESSED="${DATA_DIR}/processed"
QC_PREFIX="${PROCESSED}/cfw.qc"
if [[ ! -s "${QC_PREFIX}.bed" ]]; then
  set +e
  run_logged plink_qc "${PLINK2_BIN}" --bfile "${DATA_DIR}/cfw" --chr 1-19 --maf 0.01 --geno 0.05 --mind 0.05 --make-bed --out "${QC_PREFIX}"
  qc_rc=$?
  set -e
  if [[ "${qc_rc}" != "0" ]]; then
    if [[ "${ALLOW_QC_FALLBACK}" != "1" ]]; then
      exit "${qc_rc}"
    fi
    printf 'Requested PLINK QC with --mind 0.05 failed; running fallback without --mind. See %s\n' "${LOG_DIR}/plink_qc.log" | tee "${LOG_DIR}/plink_qc_fallback.log"
    rm -f "${QC_PREFIX}.bed" "${QC_PREFIX}.bim" "${QC_PREFIX}.fam"
    run_logged plink_qc_fallback "${PLINK2_BIN}" --bfile "${DATA_DIR}/cfw" --chr 1-19 --maf 0.01 --geno 0.05 --make-bed --out "${QC_PREFIX}"
  fi
fi

run_logged preprocess "${PYTHON_BIN}" "${SCRIPT_DIR}/prepare_cfw_mouse.py" \
  --data-dir "${DATA_DIR}" --qc-prefix "${QC_PREFIX}" --results-dir "${RESULT_DIR}" \
  --tool-status-json "${tool_status_json}" --export-emmv-z --max-emmv-snps "$([[ "${MODE}" == "smoke" ]] && echo 1000 || echo 0)"

RUN_PREFIX="${QC_PREFIX}"
if [[ "${MODE}" == "smoke" ]]; then
  RUN_PREFIX="${PROCESSED}/cfw.qc.smoke1000"
  if [[ ! -s "${RUN_PREFIX}.bed" ]]; then
    run_logged plink_smoke "${PLINK2_BIN}" --bfile "${QC_PREFIX}" --extract "${PROCESSED}/first_1000_snps.txt" --make-bed --out "${RUN_PREFIX}"
  fi
fi

KIN_DIR="${PROCESSED}/kinship"
mkdir -p "${KIN_DIR}"
if [[ -x "${GEMMA_BIN}" && ! -s "${KIN_DIR}/output/cfw_qc.cXX.txt" ]]; then
  run_logged gemma_kinship bash -c 'cd "$1" && "$2" -bfile "$3" -p "$4" -n 1 -gk 1 -o cfw_qc' _ "${KIN_DIR}" "${GEMMA_BIN}" "${QC_PREFIX}" "${PROCESSED}/trait_sets/traits_2.gemma.no_header.tsv"
fi

if [[ "${RUN_METHODS}" == "1" ]]; then
  rm -f "${RESULT_DIR}/tables/runtime.tsv"
  for t in 2 4 6 8 10; do
    pheno="${PROCESSED}/trait_sets/traits_${t}.pheno.tsv"
    gemma_pheno="${PROCESSED}/trait_sets/traits_${t}.gemma.no_header.tsv"
    trait_list="$(seq -s, 1 "${t}")"
    gemma_ns=( $(seq 1 "${t}") )

    if [[ -x "${COSMICSCANNER_BIN}" ]]; then
      run_timed "CosmicScanner" "${t}" "${RESULT_DIR}/cosmicscanner/cfw_t${t}" \
        "${COSMICSCANNER_BIN}" --bfile "${RUN_PREFIX}" --pheno "${pheno}" --mpheno-list "${trait_list}" \
        --multivariate --mv-method emmv --mv-compatible emmvgwas --gwas --no-loco --threads "${THREADS}" --out "${RESULT_DIR}/cosmicscanner/cfw_t${t}"
    fi

    if [[ -x "${RSCRIPT_BIN}" && -s "${PROCESSED}/cfw.qc.emmv.Z.csv" ]]; then
      run_timed "EMmvGWAS" "${t}" "${RESULT_DIR}/emmvgwas/cfw_t${t}" \
        "${RSCRIPT_BIN}" "${SCRIPT_DIR}/run_emmvgwas_cfw.R" "${pheno}" "${PROCESSED}/cfw.qc.emmv.Z.csv" "${RESULT_DIR}/emmvgwas/cfw_t${t}"
    fi

    if [[ -x "${GEMMA_BIN}" && -s "${KIN_DIR}/output/cfw_qc.cXX.txt" ]]; then
      run_timed "GEMMA" "${t}" "${RESULT_DIR}/gemma/output/cfw_t${t}" bash -c \
        'work="$1"; exe="$2"; bfile="$3"; pheno="$4"; kinship="$5"; shift 5; cd "${work}" && "${exe}" -bfile "${bfile}" -p "${pheno}" -n "$@" -k "${kinship}" -lmm -o "cfw_t'$t'"' _ \
        "${RESULT_DIR}/gemma" "${GEMMA_BIN}" "${RUN_PREFIX}" "${gemma_pheno}" "${KIN_DIR}/output/cfw_qc.cXX.txt" "${gemma_ns[@]}"
    fi
  done
fi

if [[ "${RUN_NULL}" == "1" ]]; then
  run_logged simulate_null "${PYTHON_BIN}" "${SCRIPT_DIR}/simulate_null.py" --template-pheno "${PROCESSED}/cfw_mouse_traits.tsv" --out "${RESULT_DIR}/null_sim" --replicates 10
  run_logged simulate_power "${PYTHON_BIN}" "${SCRIPT_DIR}/simulate_power.py" --template-pheno "${PROCESSED}/cfw_mouse_traits.tsv" --out "${RESULT_DIR}/power_sim"
fi

if [[ "${RUN_PLOTS}" == "1" ]]; then
  "${PYTHON_BIN}" "${SCRIPT_DIR}/plot_runtime.py" --runtime-tsv "${RESULT_DIR}/tables/runtime.tsv" --out "${RESULT_DIR}/plots/runtime_mouse_3methods.png"
  "${PYTHON_BIN}" "${SCRIPT_DIR}/plot_consistency.py" \
    --gemma "${RESULT_DIR}/gemma/output/cfw_t10.assoc.txt" \
    --emmv "${RESULT_DIR}/emmvgwas/cfw_t10_association_with_ids.csv" \
    --cosmic "${RESULT_DIR}/cosmicscanner/cfw_t10.assoc.mv.txt" \
    --out "${RESULT_DIR}/plots/consistency_mouse_10traits.png" \
    --summary "${RESULT_DIR}/tables/consistency_r.tsv"
  "${PYTHON_BIN}" "${SCRIPT_DIR}/plot_manhattan.py" \
    --bim "${RUN_PREFIX}.bim" \
    --gemma "${RESULT_DIR}/gemma/output/cfw_t10.assoc.txt" \
    --emmv "${RESULT_DIR}/emmvgwas/cfw_t10_association_with_ids.csv" \
    --cosmic "${RESULT_DIR}/cosmicscanner/cfw_t10.assoc.mv.txt" \
    --out "${RESULT_DIR}/plots/manhattan_mouse_10traits_3methods.png"
  "${PYTHON_BIN}" "${SCRIPT_DIR}/plot_qq.py" \
    --inputs "${RESULT_DIR}/gemma/output/cfw_t10.assoc.txt" "${RESULT_DIR}/emmvgwas/cfw_t10_association_with_ids.csv" "${RESULT_DIR}/cosmicscanner/cfw_t10.assoc.mv.txt" \
    --labels GEMMA EMmvGWAS CosmicScanner \
    --out "${RESULT_DIR}/plots/qq_mouse_null_3methods.png" \
    --summary "${RESULT_DIR}/tables/qq_lambda_gc.tsv"
fi

"${PYTHON_BIN}" - "${RESULT_DIR}" "${DATA_DIR}" "${MODE}" <<'PY'
import json, sys
from pathlib import Path
res, data, mode = map(Path, sys.argv[1:3]) + [sys.argv[3]] if False else (Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3])
def read(p):
    return Path(p).read_text(encoding="utf-8") if Path(p).exists() else ""
runtime = read(res / "tables/runtime.tsv")
cons = read(res / "tables/consistency_r.tsv")
tools = json.loads(read(res / "tables/tool_status.json") or "{}")
lines = [
    "# CFW Mouse Benchmark Report", "",
    f"Mode: `{mode}`", "",
    "## Data", "",
    f"Data directory: `{data}`",
    f"Audit report: `{res / 'data_audit.md'}`",
    f"Preprocess report: `{data / 'processed/preprocess_report.md'}`", "",
    "## Tool Status", "",
]
for k, v in tools.items():
    lines.append(f"- {k}: {v}")
lines += ["", "## Runtime Table", "", "```tsv", runtime.strip(), "```", "", "## Consistency", "", "```tsv", cons.strip(), "```", "", "## Figures", "",
          f"- Runtime: `{res / 'plots/runtime_mouse_3methods.png'}`",
          f"- Consistency: `{res / 'plots/consistency_mouse_10traits.png'}`",
          f"- Manhattan: `{res / 'plots/manhattan_mouse_10traits_3methods.png'}`",
          f"- QQ/null smoke: `{res / 'plots/qq_mouse_null_3methods.png'}`", "",
          "## Reproduction Differences", "",
          "- This report is a smoke benchmark by default, not a full reproduction of the EMmvGWAS mouse panels.",
          "- The requested PLINK QC with `--mind 0.05` removed all samples on the hard-call BED files; the default pipeline records that failure and falls back to `--maf 0.01 --geno 0.05` unless `ALLOW_QC_FALLBACK=0` is set.",
          "- Smoke mode scans the first 1000 QC SNPs. Runtime and consistency values should not be compared directly to the paper's full 92734-SNP benchmark.",
          "- The QQ plot produced by the default run uses smoke association p-values. Use `RUN_NULL=1` plus method runs on generated null phenotypes for a true null-simulation QQ benchmark.",
          "- Current smoke results show low CosmicScanner-vs-GEMMA and CosmicScanner-vs-EMmvGWAS correlation under `--mv-method emmv --mv-compatible emmvgwas`; this is retained as a diagnostic finding for follow-up input/coding compatibility checks.", "",
          "## Notes", "",
          "- Default smoke mode uses the first 1000 QC SNPs for method runs.",
          "- Full benchmark command: `MODE=full RUN_METHODS=1 RUN_PLOTS=1 bash scripts/cfw_mouse_benchmark/run_cfw_mouse_benchmark.sh`.",
          "- If any method is missing or fails, inspect `results/cfw_mouse_benchmark/logs/` and rerun after fixing the dependency."]
(res / "CFW_mouse_benchmark_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
PY

echo "CFW mouse benchmark pipeline complete: ${RESULT_DIR}"
