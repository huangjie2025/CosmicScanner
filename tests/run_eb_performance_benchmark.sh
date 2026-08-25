#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/eb_performance_benchmark"
RESULT_DIR="${ROOT}/tests/results/eb_performance_benchmark"
PREFIX="${DATA_DIR}/ebperf"
N="${N:-180}"
M="${M:-1200}"
THREADS_LIST="${THREADS_LIST:-1 2 4}"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}"

python3 - "${PREFIX}" "${N}" "${M}" <<'PY'
import math, random, sys
from pathlib import Path

prefix = Path(sys.argv[1])
n = int(sys.argv[2]); m = int(sys.argv[3])
rng = random.Random(20260628)
samples = [(f'F{i+1:05d}', f'I{i+1:05d}', 0 if i < n//2 else 1) for i in range(n)]
genos, variants = [], []
causal = min(m - 1, max(0, m // 3))
for j in range(m):
    chrom = str(1 + j // max(1, m // 8))
    if int(chrom) > 8:
        chrom = '8'
    base = 0.08 + 0.42 * ((j * 37) % 101) / 100.0
    shift = 0.05 * math.sin(j * 0.11)
    row = []
    for _, _, pop in samples:
        p = min(0.95, max(0.02, base + (shift if pop else -shift)))
        row.append((1 if rng.random() < p else 0) + (1 if rng.random() < p else 0))
    genos.append(row)
    variants.append((chrom, f'rs{j+1:07d}', 1000 + 10*j, 'A', 'C'))

prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix('.fam').open('w') as f:
    for fid, iid, _ in samples:
        f.write(f'{fid}\t{iid}\t0\t0\t1\t-9\n')
with prefix.with_suffix('.bim').open('w') as f:
    for chrom, snp, bp, a1, a2 in variants:
        f.write(f'{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n')
with prefix.with_suffix('.pheno').open('w') as f:
    for i, (fid, iid, pop) in enumerate(samples):
        g = genos[causal][i] - 1.0
        y = 0.20 * pop + 0.60 * g + rng.gauss(0, 0.80)
        f.write(f'{fid}\t{iid}\t{y:.8f}\n')

def code(g):
    return 0b00 if g == 2 else (0b10 if g == 1 else 0b11)
bps = (n + 3) // 4
bed = bytearray(b'\x6c\x1b\x01')
for row in genos:
    block = bytearray(bps)
    for i, g in enumerate(row):
        block[i // 4] |= code(g) << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix('.bed').write_bytes(bed)
PY

run_case() {
  local label="$1" threads="$2" model="$3" out="$4"
  local log="${out}.log"
  local start end elapsed rows rss
  start=$(date +%s%N)
  /usr/bin/time -f 'max_rss_kb\t%M' -o "${out}.time" \
    "${BIN}" --bfile "${PREFIX}" --pheno "${PREFIX}.pheno" --mpheno 1 --gwas --model "${model}" --no-loco --maf 0.0 --threads "${threads}" --out "${out}" >"${log}" 2>&1
  end=$(date +%s%N)
  elapsed=$(python3 - <<PY
print((${end} - ${start}) / 1e9)
PY
)
  rows=$(python3 - "${out}.assoc.txt" <<'PY'
import sys
print(max(0, sum(1 for _ in open(sys.argv[1])) - 1))
PY
)
  rss=$(python3 - "${out}.time" <<'PY'
import sys
vals = {}
for line in open(sys.argv[1]):
    k, _, v = line.strip().partition('\t')
    vals[k] = v
print(vals.get('max_rss_kb', 'NA'))
PY
)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${label}" "${model}" "${threads}" "${rows}" "${elapsed}" "${rss}" "$(python3 - <<PY
elapsed=${elapsed}; rows=${rows}
print(rows / elapsed if elapsed > 0 else 0.0)
PY
)" >>"${RESULT_DIR}/eb_performance.summary.tsv"
}

printf 'case\tmodel\tthreads\trows\telapsed_sec\tmax_rss_kb\trows_per_sec\n' >"${RESULT_DIR}/eb_performance.summary.tsv"
run_case standard 1 lmm "${RESULT_DIR}/standard_t1"
for t in ${THREADS_LIST}; do
  run_case "eb_t${t}" "${t}" lmm-eb "${RESULT_DIR}/eb_t${t}"
done

python3 - "${RESULT_DIR}/eb_performance.summary.tsv" "${RESULT_DIR}/eb_performance.summary.md" <<'PY'
import csv, sys
tsv, md = sys.argv[1:]
rows = list(csv.DictReader(open(tsv), delimiter='\t'))
if not rows:
    raise SystemExit('empty performance summary')
with open(md, 'w') as f:
    f.write('# EB-LMM Performance Benchmark\n\n')
    f.write('| Case | Model | Threads | Rows | Elapsed sec | Max RSS KB | Rows/sec |\n')
    f.write('|---|---|---:|---:|---:|---:|---:|\n')
    for r in rows:
        f.write(f"| {r['case']} | {r['model']} | {r['threads']} | {r['rows']} | {float(r['elapsed_sec']):.4g} | {r['max_rss_kb']} | {float(r['rows_per_sec']):.4g} |\n")
print(f'EB performance benchmark passed: {tsv}')
PY
