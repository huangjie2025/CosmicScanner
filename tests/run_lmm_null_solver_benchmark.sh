#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
N="${N:-1000}"
M="${M:-10000}"
RUN_LOCO="${RUN_LOCO:-0}"
RUN_GEMMA="${RUN_GEMMA:-0}"
SOLVERS="${SOLVERS:-brent idul idul_plus aidul}"
OUT_DIR="${ROOT}/tests/results/lmm_null_solver_n${N}_m${M}_loco${RUN_LOCO}"
mkdir -p "${OUT_DIR}"

SUMMARY="${OUT_DIR}/solver_summary.tsv"
printf 'solver\tn\tm\tloco\tseconds\tmax_rss_kb\tassoc_rows\tgrm_s\tnull_s\tscan_s\teta_iter\tlog\n' > "${SUMMARY}"

for solver in ${SOLVERS}; do
  echo "[LMM-null-benchmark] solver=${solver}, N=${N}, M=${M}, RUN_LOCO=${RUN_LOCO}"
  LMM_NULL_SOLVER="${solver}" N="${N}" M="${M}" RUN_GEMMA="${RUN_GEMMA}" RUN_LOCO="${RUN_LOCO}" \
    bash "${ROOT}/tests/run_gemma_speed_benchmark.sh" >/"${OUT_DIR}/${solver}.stdout" 2>&1

  SPEED_DIR="${ROOT}/tests/results/gwas_speed_n${N}_m${M}"
  cp "${SPEED_DIR}/cosmic.log" "${OUT_DIR}/cosmic_${solver}.log"
  cp "${SPEED_DIR}/speed_summary.tsv" "${OUT_DIR}/speed_summary_${solver}.tsv"

  python3 - "${solver}" "${RUN_LOCO}" "${OUT_DIR}" "${SUMMARY}" <<'PY'
import re
import sys
from pathlib import Path

solver, run_loco, out_dir, summary = sys.argv[1:]
out_dir = Path(out_dir)
summary = Path(summary)
speed = (out_dir / f"speed_summary_{solver}.tsv").read_text().strip().splitlines()[1].split("\t")
log_path = out_dir / f"cosmic_{solver}.log"
log = log_path.read_text().splitlines()
times = {}
eta = ""
null_sum = 0.0
scan_sum = 0.0
for line in log:
    if "[time]" in line:
        part = line.split("[time]", 1)[1].strip()
        label, val = part.rsplit(":", 1)
        label = label.strip()
        sec = val.strip().rstrip("s")
        times[label] = sec
        try:
            if label.startswith("null model"):
                null_sum += float(sec)
            if label.startswith("SNP scan"):
                scan_sum += float(sec)
        except ValueError:
            pass
    if " eta=" in line and "iterations=" in line:
        eta = line.strip()
grm = times.get("PLINK global GRM build", times.get("PLINK split GRM build", ""))
null = times.get("null model", f"{null_sum:.2f}" if null_sum else "")
scan = times.get("SNP scan", f"{scan_sum:.2f}" if scan_sum else "")
with summary.open("a", encoding="utf-8") as f:
    f.write("\t".join([solver, speed[1], speed[2], run_loco, speed[4], speed[5], speed[6], grm, null, scan, eta, str(log_path)]) + "\n")
PY
done

cat "${SUMMARY}"
