#!/usr/bin/env bash
set -euo pipefail

ROOT="${COSMICSCANNER_ROOT:-/disk228/huangjie/PanOmics/CosmicScanner}"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
RESULT_DIR="${COSMIC_BSLMM_REF_OUT:-${ROOT}/tests/results/gemma_bslmm_reference_validation}"
DATA_PREFIX="${RESULT_DIR}/compat_data/data/smoke"
GEMMA_WORK="${RESULT_DIR}/gemma"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

resolve_gemma() {
  if [[ -n "${GEMMA_BIN:-}" ]]; then
    printf '%s\n' "${GEMMA_BIN}"
    return 0
  fi
  if command -v gemma >/dev/null 2>&1; then
    command -v gemma
    return 0
  fi
  local conda_bin="/disk228/huangjie/software/miniconda3/envs/benchmark/bin/gemma"
  if [[ -x "${conda_bin}" ]]; then
    printf '%s\n' "${conda_bin}"
    return 0
  fi
  return 1
}

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi
if ! GEMMA_EXE="$(resolve_gemma)"; then
  echo "Missing GEMMA executable; skipping external BSLMM reference validation." >&2
  exit 2
fi

mkdir -p "${RESULT_DIR}" "${GEMMA_WORK}"
COSMIC_SMOKE_OUT="${RESULT_DIR}/compat_data" "${ROOT}/tests/run_gemma_emmvgwas_compat_smoke.sh" >/dev/null

rm -rf "${GEMMA_WORK}/output"
mkdir -p "${GEMMA_WORK}/output"

echo "[GEMMA-BSLMM-ref] Running GEMMA -bslmm 2..."
(
  cd "${GEMMA_WORK}"
  "${GEMMA_EXE}" -bfile "${DATA_PREFIX}" -bslmm 2 -o gemma_bslmm2
)

echo "[GEMMA-BSLMM-ref] Running CosmicScanner -bslmm 2..."
"${BIN}" \
  --bfile "${DATA_PREFIX}" \
  -p "${DATA_PREFIX}.pheno" \
  -n 1 \
  -bslmm 2 \
  --bslmm-h2 0.01 \
  --maf 0.0 \
  -o "${RESULT_DIR}/cosmic_bslmm2" >/dev/null

python3 - "${GEMMA_WORK}/output/gemma_bslmm2.param.txt" \
          "${GEMMA_WORK}/output/gemma_bslmm2.bv.txt" \
          "${RESULT_DIR}/cosmic_bslmm2.param.txt" \
          "${RESULT_DIR}/cosmic_bslmm2.bv.txt" \
          "${RESULT_DIR}/cosmic_bslmm2.bslmm.meta.txt" \
          "${RESULT_DIR}/bslmm_reference.summary.tsv" <<'PY'
import csv
import math
import sys
from pathlib import Path

gemma_param, gemma_bv, cosmic_param, cosmic_bv, cosmic_meta, summary = map(Path, sys.argv[1:])

def require(condition, message):
    if not condition:
        raise SystemExit(message)

def read_param(path):
    rows = list(csv.DictReader(path.open("r", encoding="utf-8"), delimiter="\t"))
    require(rows and {"chr", "rs", "ps", "n_miss", "alpha", "beta", "gamma"}.issubset(rows[0]), f"bad param header: {path}")
    return {row["rs"]: float(row["alpha"]) for row in rows}

def read_vec(path):
    vals = []
    for line in path.open("r", encoding="utf-8"):
        line = line.strip()
        if not line or line == "NA":
            continue
        vals.append(float(line))
    return vals

def corr(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    require(len(pairs) >= 3, "not enough finite pairs")
    mx = sum(x for x, _ in pairs) / len(pairs)
    my = sum(y for _, y in pairs) / len(pairs)
    vx = sum((x - mx) ** 2 for x, _ in pairs)
    vy = sum((y - my) ** 2 for _, y in pairs)
    require(vx > 0 and vy > 0, "zero variance for correlation")
    return sum((x - mx) * (y - my) for x, y in pairs) / math.sqrt(vx * vy)

gp = read_param(gemma_param)
cp = read_param(cosmic_param)
ids = [snp for snp in gp if snp in cp]
require(len(ids) == len(gp) == len(cp), "GEMMA/CosmicScanner param SNP sets differ")
alpha_corr = corr([gp[i] for i in ids], [cp[i] for i in ids])
alpha_abs_corr = corr([abs(gp[i]) for i in ids], [abs(cp[i]) for i in ids])

gbv = read_vec(gemma_bv)
cbv = read_vec(cosmic_bv)
require(len(gbv) == len(cbv) and len(gbv) > 0, "GEMMA/CosmicScanner BV lengths differ")
bv_corr = corr(gbv, cbv)

meta = {}
for line in cosmic_meta.read_text(encoding="utf-8").splitlines():
    if "\t" in line:
        k, v = line.split("\t", 1)
        meta[k] = v
require(meta.get("mode") == "2", "CosmicScanner BSLMM metadata mode mismatch")
require(meta.get("implementation") == "deterministic-ridge-gblup", "CosmicScanner BSLMM metadata implementation mismatch")

summary.parent.mkdir(parents=True, exist_ok=True)
with summary.open("w", encoding="utf-8") as out:
    out.write("metric\tvalue\tthreshold\tpass\n")
    checks = [
        ("n_snps", len(ids), len(ids), True),
        ("alpha_corr", alpha_corr, 0.85, alpha_corr >= 0.85),
        ("alpha_abs_corr", alpha_abs_corr, 0.70, alpha_abs_corr >= 0.70),
        ("bv_corr", bv_corr, 0.95, bv_corr >= 0.95),
    ]
    for name, value, threshold, ok in checks:
        out.write(f"{name}\t{value}\t{threshold}\t{int(ok)}\n")
        require(ok, f"{name} failed: {value} < {threshold}")
PY

cat "${RESULT_DIR}/bslmm_reference.summary.tsv"
echo "GEMMA BSLMM2 reference validation passed. Outputs: ${RESULT_DIR}"
