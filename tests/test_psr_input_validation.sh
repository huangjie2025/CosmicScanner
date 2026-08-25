#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/results/psr_input_validation"
PFX="${DATA_DIR}/invalid_psr"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

resolve_rscript() {
  if [[ -n "${RSCRIPT_BIN:-}" ]]; then printf '%s\n' "${RSCRIPT_BIN}"; return 0; fi
  if command -v Rscript >/dev/null 2>&1; then command -v Rscript; return 0; fi
  local conda_r="/disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript"
  if [[ -x "${conda_r}" ]]; then printf '%s\n' "${conda_r}"; return 0; fi
  return 1
}

mkdir -p "${DATA_DIR}"

python3 - "${PFX}" <<'PY'
import sys
from pathlib import Path
prefix = Path(sys.argv[1])
n, m = 12, 16
G = [[(i + j) % 3 for j in range(m)] for i in range(n)]
prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for i in range(n): f.write(f"F{i+1:04d}\tI{i+1:04d}\t0\t0\t1\t-9\n")
with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for j in range(m): f.write(f"1\tiv{j+1:06d}\t0\t{1000+j}\tA\tC\n")
def code(g): return {2:0b00,1:0b10,0:0b11}.get(g,0b01)
bed = bytearray(b"\x6c\x1b\x01")
bytes_per_snp = (n + 3) // 4
for j in range(m):
    block = bytearray(bytes_per_snp)
    for i in range(n): block[i//4] |= code(G[i][j]) << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix(".bed").write_bytes(bed)
with prefix.with_suffix(".geno.tsv").open("w", encoding="utf-8") as f:
    f.write("FID\tIID" + "".join(f"\tiv{j+1:06d}" for j in range(m)) + "\n")
    for i in range(n): f.write(f"F{i+1:04d}\tI{i+1:04d}" + "".join(f"\t{G[i][j]}" for j in range(m)) + "\n")
with prefix.with_suffix(".assoc.txt").open("w", encoding="utf-8") as f:
    f.write("SNP\tP\n")
    for j in range(4): f.write(f"iv{j+1:06d}\t0.0{j+1}\n")
def pheno(name, vals):
    with prefix.with_suffix(name).open("w", encoding="utf-8") as f:
        for i, v in enumerate(vals): f.write(f"F{i+1:04d}\tI{i+1:04d}\t{v}\n")
pheno(".binary_oneclass.pheno", [1]*n)
pheno(".poisson_negative.pheno", [0,1,2,-1] + [1]*(n-4))
pheno(".binary.pheno", [0,1]*(n//2))
pheno(".poisson.pheno", [1]*n)
with prefix.with_suffix(".binomial_event_gt_trial.pheno").open("w", encoding="utf-8") as f:
    for i in range(n): f.write(f"F{i+1:04d}\tI{i+1:04d}\t2\t1\n")
with prefix.with_suffix(".binomial_trial_zero.pheno").open("w", encoding="utf-8") as f:
    for i in range(n): f.write(f"F{i+1:04d}\tI{i+1:04d}\t0\t0\n")
with prefix.with_suffix(".binomial.pheno").open("w", encoding="utf-8") as f:
    for i in range(n): f.write(f"F{i+1:04d}\tI{i+1:04d}\t{ i % 2 }\t1\n")
with prefix.with_suffix(".offset_short").open("w", encoding="utf-8") as f:
    for i in range(n-2): f.write("0\n")
with prefix.with_suffix(".offset_nan").open("w", encoding="utf-8") as f:
    for i in range(n): f.write("nan\n" if i == 3 else "0\n")
PY

expect_fail() {
  local label="$1"; shift
  if "$@" >"${DATA_DIR}/${label}.stdout" 2>"${DATA_DIR}/${label}.stderr"; then
    echo "Expected failure but command passed: ${label}" >&2
    exit 1
  fi
  echo "[PSR-GWAS] Input validation expected failure passed: ${label}"
}

expect_fail binary_oneclass "${BIN}" --bfile "${PFX}" --pheno "${PFX}.binary_oneclass.pheno" --mpheno 1 --gwas --threshold-model --no-loco --maf 0 --out "${DATA_DIR}/binary_oneclass"
expect_fail binomial_event_gt_trial "${BIN}" --bfile "${PFX}" --gwas --count-model --family binomial --binomial-pheno "${PFX}.binomial_event_gt_trial.pheno" --no-loco --maf 0 --out "${DATA_DIR}/binomial_event_gt_trial"
expect_fail binomial_trial_zero "${BIN}" --bfile "${PFX}" --gwas --count-model --family binomial --binomial-pheno "${PFX}.binomial_trial_zero.pheno" --no-loco --maf 0 --out "${DATA_DIR}/binomial_trial_zero"
expect_fail poisson_negative "${BIN}" --bfile "${PFX}" --pheno "${PFX}.poisson_negative.pheno" --mpheno 1 --gwas --count-model --family poisson --no-loco --maf 0 --out "${DATA_DIR}/poisson_negative"
expect_fail offset_short "${BIN}" --bfile "${PFX}" --pheno "${PFX}.poisson.pheno" --mpheno 1 --gwas --count-model --family poisson --offset "${PFX}.offset_short" --no-loco --maf 0 --out "${DATA_DIR}/offset_short"
expect_fail offset_nan "${BIN}" --bfile "${PFX}" --pheno "${PFX}.poisson.pheno" --mpheno 1 --gwas --count-model --family poisson --offset "${PFX}.offset_nan" --no-loco --maf 0 --out "${DATA_DIR}/offset_nan"

if RSCRIPT_EXE="$(resolve_rscript)"; then
  "${RSCRIPT_EXE}" "${ROOT}/tests/psr_candidate_refine.R" "${PFX}" "${PFX}.assoc.txt" binary "${DATA_DIR}/firth_missing_package.tsv" 2 firth
  python3 - "${DATA_DIR}/firth_missing_package.tsv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t"))
if not rows:
    raise SystemExit("No firth refinement rows")
if "error_message" not in rows[0]:
    raise SystemExit("Missing error_message column")
PY
fi

echo "[PSR-GWAS] Input validation regression passed."
