#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_PREFIX="${ROOT}/tests/data/psr_gwas_smoke/psr_smoke"
RESULT_DIR="${ROOT}/tests/results/psr_gwas_smoke"
ASSOC="${RESULT_DIR}/cosmic_poisson.assoc.txt"

resolve_rscript() {
  if [[ -n "${RSCRIPT_BIN:-}" ]]; then printf '%s\n' "${RSCRIPT_BIN}"; return 0; fi
  if command -v Rscript >/dev/null 2>&1; then command -v Rscript; return 0; fi
  local conda_r="/disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript"
  if [[ -x "${conda_r}" ]]; then printf '%s\n' "${conda_r}"; return 0; fi
  return 1
}

if ! RSCRIPT_EXE="$(resolve_rscript)"; then
  echo "Rscript not found; skipping PSR candidate refinement suite." >&2
  exit 0
fi

mkdir -p "${RESULT_DIR}"

check_columns() {
  local file="$1"
  python3 - "$file" <<'PY'
import csv, sys
path = sys.argv[1]
required = ["trait_type", "method", "snp", "beta", "se", "wald_or_stat", "p_value", "runtime", "converged", "error_message"]
with open(path, encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    missing = [x for x in required if x not in (reader.fieldnames or [])]
    rows = list(reader)
if missing:
    raise SystemExit(f"Missing columns in {path}: {missing}")
if not rows:
    raise SystemExit(f"No refinement rows in {path}")
PY
}

for trait in binary poisson binomial; do
  echo "[PSR-GWAS] Candidate refinement suite: ${trait} glm..."
  assoc="${RESULT_DIR}/cosmic_${trait}.assoc.txt"
  "${RSCRIPT_EXE}" "${ROOT}/tests/psr_candidate_refine.R" "${DATA_PREFIX}" "${assoc}" "${trait}" "${RESULT_DIR}/candidate_suite.${trait}.glm.tsv" 5 glm
  check_columns "${RESULT_DIR}/candidate_suite.${trait}.glm.tsv"
done

for method in firth glmer glmmtmb; do
  trait="binary"
  assoc="${RESULT_DIR}/cosmic_binary.assoc.txt"
  if [[ "${method}" != "firth" ]]; then
    trait="poisson"
    assoc="${ASSOC}"
  fi
  out="${RESULT_DIR}/candidate_suite.${trait}.${method}.tsv"
  echo "[PSR-GWAS] Candidate refinement suite: ${method} optional hook..."
  "${RSCRIPT_EXE}" "${ROOT}/tests/psr_candidate_refine.R" "${DATA_PREFIX}" "${assoc}" "${trait}" "${out}" 3 "${method}"
  check_columns "${out}"
done

echo "[PSR-GWAS] Candidate refinement suite passed."
