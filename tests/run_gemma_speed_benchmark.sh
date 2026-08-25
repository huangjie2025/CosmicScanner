#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
N="${N:-1000}"
M="${M:-10000}"
SEED="${SEED:-20260601}"
RUN_GEMMA="${RUN_GEMMA:-1}"
RUN_LOCO="${RUN_LOCO:-0}"
LMM_NULL_SOLVER="${LMM_NULL_SOLVER:-brent}"
if [[ -z "${PYTHON_BIN:-}" ]]; then
  if [[ -x "/disk228/huangjie/software/miniconda3/envs/benchmark/bin/python" ]]; then
    PYTHON_BIN="/disk228/huangjie/software/miniconda3/envs/benchmark/bin/python"
  else
    PYTHON_BIN="python3"
  fi
fi
DATA_DIR="${ROOT}/tests/data/gwas_speed_n${N}_m${M}"
RESULT_DIR="${ROOT}/tests/results/gwas_speed_n${N}_m${M}"
DATA_PREFIX="${DATA_DIR}/speed"
COSMIC_PREFIX="${RESULT_DIR}/cosmic/speed"
GEMMA_WORK="${RESULT_DIR}/gemma"
SUMMARY="${RESULT_DIR}/speed_summary.tsv"

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
  echo "Build it first or set COSMICSCANNER_BIN=/path/to/cosmicscanner." >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}/cosmic" "${GEMMA_WORK}"
rm -f "${SUMMARY}"

"${PYTHON_BIN}" - "${DATA_PREFIX}" "${N}" "${M}" "${SEED}" <<'PY'
import sys
from pathlib import Path

try:
    import numpy as np
except Exception as exc:
    raise SystemExit(f"This benchmark generator requires numpy. Set PYTHON_BIN to a Python with numpy. Error: {exc}")

prefix = Path(sys.argv[1])
n = int(sys.argv[2])
m = int(sys.argv[3])
seed = int(sys.argv[4])
rng = np.random.default_rng(seed)

prefix.parent.mkdir(parents=True, exist_ok=True)
chrom_count = 22 if m >= 22 else max(1, m)

samples = []
pop_arr = np.empty(n, dtype=np.int16)
for i in range(n):
    pop = i % 4
    pop_arr[i] = pop
    fid = f"F{i+1:07d}"
    iid = f"I{i+1:07d}"
    samples.append((fid, iid, pop))

causal = {}
for idx, eff in [(m // 13, 0.35), (m // 3, -0.28), ((2 * m) // 3, 0.22)]:
    if 0 <= idx < m:
        causal[idx] = eff

phenotypes = 0.2 * pop_arr.astype(np.float64) + rng.normal(0.0, 0.8, size=n)
bytes_per_snp = (n + 3) // 4
pad_n = bytes_per_snp * 4

with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for j in range(m):
        chrom = (j % chrom_count) + 1
        f.write(f"{chrom}\trs{j+1:09d}\t0\t{1000 + j * 7}\tA\tC\n")

chunk = 2048
pop_shift_base = (pop_arr.astype(np.float32) - 1.5)[None, :]
with prefix.with_suffix(".bed").open("wb") as bed:
    bed.write(b"\x6c\x1b\x01")
    for start in range(0, m, chunk):
        end = min(m, start + chunk)
        idx = np.arange(start, end, dtype=np.float64)
        base_p = 0.05 + 0.45 * ((idx.astype(np.int64) * 7919) % 1000) / 999.0
        shift = 0.04 * np.sin(idx * 0.031)
        p = np.clip(base_p[:, None].astype(np.float32) + shift[:, None].astype(np.float32) * pop_shift_base,
                    0.01, 0.95)
        g = (rng.random((end - start, n), dtype=np.float32) < p).astype(np.uint8)
        g += (rng.random((end - start, n), dtype=np.float32) < p).astype(np.uint8)

        for snp_idx, eff in causal.items():
            if start <= snp_idx < end:
                local = snp_idx - start
                phenotypes += eff * (g[local].astype(np.float64) - 2.0 * p[local].astype(np.float64))

        codes = np.empty((end - start, pad_n), dtype=np.uint8)
        codes[:, :n] = np.where(g == 2, 0, np.where(g == 1, 2, 3)).astype(np.uint8)
        if pad_n > n:
            codes[:, n:] = 0
        packed = (codes[:, 0::4]
                  | (codes[:, 1::4] << 2)
                  | (codes[:, 2::4] << 4)
                  | (codes[:, 3::4] << 6)).astype(np.uint8)
        bed.write(packed.tobytes(order="C"))

with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for (fid, iid, _), y in zip(samples, phenotypes):
        f.write(f"{fid}\t{iid}\t0\t0\t1\t{float(y):.8f}\n")

with prefix.with_suffix(".pheno").open("w", encoding="utf-8") as f:
    for (fid, iid, _), y in zip(samples, phenotypes):
        f.write(f"{fid}\t{iid}\t{float(y):.8f}\n")

with prefix.with_suffix(".truth.txt").open("w", encoding="utf-8") as f:
    f.write("SNP\tTRUE_BETA\n")
    for idx, eff in sorted(causal.items()):
        f.write(f"rs{idx+1:09d}\t{eff}\n")
PY

write_summary_header() {
  if [[ ! -f "${SUMMARY}" ]]; then
    printf 'tool\tn\tm\tmode\tseconds\tmax_rss_kb\tassoc_rows\tlog\n' > "${SUMMARY}"
  fi
}

count_assoc_rows() {
  local file="$1"
  if [[ ! -f "${file}" ]]; then
    printf '0'
    return
  fi
  python3 - "$file" <<'PY'
import sys
from pathlib import Path
p = Path(sys.argv[1])
print(max(0, sum(1 for _ in p.open()) - 1))
PY
}

run_timed() {
  local tool="$1"
  local mode="$2"
  local log="$3"
  shift 3
  local metric="${log}.time"
  /usr/bin/time -f 'seconds\t%e\nmax_rss_kb\t%M' -o "${metric}" "$@" >"${log}" 2>&1
  local seconds rss
  seconds="$(awk '$1=="seconds"{print $2}' "${metric}")"
  rss="$(awk '$1=="max_rss_kb"{print $2}' "${metric}")"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t' "${tool}" "${N}" "${M}" "${mode}" "${seconds}" "${rss}" >> "${SUMMARY}"
}

write_summary_header

COSMIC_MODE="no_loco"
COSMIC_ARGS=("${BIN}" --bfile "${DATA_PREFIX}" --pheno "${DATA_PREFIX}.pheno" --mpheno 1 --gwas --maf 0.0 --lmm-null-solver "${LMM_NULL_SOLVER}" --out "${COSMIC_PREFIX}")
if [[ "${RUN_LOCO}" != "1" ]]; then
  COSMIC_ARGS+=(--no-loco)
else
  COSMIC_MODE="loco"
fi

echo "[GWAS-speed] Running CosmicScanner (${COSMIC_MODE}) for N=${N}, M=${M}..."
run_timed "CosmicScanner" "${COSMIC_MODE}" "${RESULT_DIR}/cosmic.log" "${COSMIC_ARGS[@]}"
printf '%s\t%s\n' "$(count_assoc_rows "${COSMIC_PREFIX}.assoc.txt")" "${RESULT_DIR}/cosmic.log" >> "${SUMMARY}"

if [[ "${RUN_GEMMA}" == "1" ]]; then
  if ! GEMMA_EXE="$(resolve_gemma)"; then
    echo "Missing GEMMA executable; rerun with RUN_GEMMA=0 or set GEMMA_BIN=/path/to/gemma." >&2
    exit 2
  fi
  rm -rf "${GEMMA_WORK}/output"
  echo "[GWAS-speed] Running GEMMA kinship..."
  run_timed "GEMMA" "kinship" "${RESULT_DIR}/gemma_kinship.log" \
    bash -c 'cd "$1" && "$2" -bfile "$3" -gk 1 -o speed_kinship' _ "${GEMMA_WORK}" "${GEMMA_EXE}" "${DATA_PREFIX}"
  printf '0\t%s\n' "${RESULT_DIR}/gemma_kinship.log" >> "${SUMMARY}"

  KINSHIP="${GEMMA_WORK}/output/speed_kinship.cXX.txt"
  echo "[GWAS-speed] Running GEMMA LMM..."
  run_timed "GEMMA" "lmm" "${RESULT_DIR}/gemma_lmm.log" \
    bash -c 'cd "$1" && "$2" -bfile "$3" -k "$4" -lmm 4 -o speed_lmm' _ "${GEMMA_WORK}" "${GEMMA_EXE}" "${DATA_PREFIX}" "${KINSHIP}"
  printf '%s\t%s\n' "$(count_assoc_rows "${GEMMA_WORK}/output/speed_lmm.assoc.txt")" "${RESULT_DIR}/gemma_lmm.log" >> "${SUMMARY}"
fi

echo "[GWAS-speed] Summary: ${SUMMARY}"
cat "${SUMMARY}"
