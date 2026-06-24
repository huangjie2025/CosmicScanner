#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
N="${N:-300}"
M="${M:-1500}"
TIMES="${TIMES:-7}"
CHROMS="${CHROMS:-6}"
SEED="${SEED:-20260615}"
RRM_ORDER="${RRM_ORDER:-3}"
RRM_GRID="${RRM_GRID:-11}"
RRM_TOP="${RRM_TOP:-50}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/tests/results/rrm_bspline_benchmark_n${N}_m${M}_t${TIMES}}"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

mkdir -p "${RESULT_DIR}"

run_one() {
  local basis="$1"
  local outdir="${RESULT_DIR}/${basis}"
  RESULT_DIR="${outdir}" \
  N="${N}" M="${M}" TIMES="${TIMES}" CHROMS="${CHROMS}" SEED="${SEED}" \
  RRM_ORDER="${RRM_ORDER}" RRM_BASIS="${basis}" RRM_GRID="${RRM_GRID}" RRM_TOP="${RRM_TOP}" \
  bash "${ROOT}/tests/run_rrm_dynamic_benchmark.sh" > "${RESULT_DIR}/${basis}.log" 2>&1
}

run_one legendre
run_one bspline

"${PYTHON_BIN}" - "${RESULT_DIR}" <<'PY'
from pathlib import Path
import sys
import pandas as pd

root = Path(sys.argv[1])
rows = []
for basis in ['legendre', 'bspline']:
    summary = root / basis / 'summary.tsv'
    if not summary.exists():
        raise SystemExit(f'Missing summary: {summary}')
    kv = {}
    for line in summary.read_text(encoding='utf-8').strip().splitlines()[1:]:
        key, value = line.split('\t', 1)
        kv[key] = value
    rows.append({
        'Basis': basis,
        'N': kv.get('N'),
        'M': kv.get('M'),
        'TimePoints': kv.get('TimePoints'),
        'RrmOrder': kv.get('RrmOrder'),
        'ElapsedSeconds': kv.get('ElapsedSeconds'),
        'MaxRSSKB': kv.get('MaxRSSKB'),
        'DynamicRank': kv.get('rrmb_dyn_causal.Rank'),
        'DynamicPattern': kv.get('rrmb_dyn_causal.Pattern'),
        'StableRank': kv.get('rrmb_stable_causal.Rank'),
        'EarlyRank': kv.get('rrmb_early_causal.Rank'),
        'LateRank': kv.get('rrmb_late_causal.Rank'),
        'TransientRank': kv.get('rrmb_transient_causal.Rank'),
    })

df = pd.DataFrame(rows)
df.to_csv(root / 'bspline_benchmark_summary.tsv', sep='\t', index=False)
with (root / 'bspline_benchmark_summary.md').open('w', encoding='utf-8') as handle:
    handle.write('# RRM B-spline Benchmark\n\n')
    cols = list(df.columns)
    handle.write('| ' + ' | '.join(cols) + ' |\n')
    handle.write('| ' + ' | '.join(['---'] * len(cols)) + ' |\n')
    for row in df.astype(str).itertuples(index=False):
        handle.write('| ' + ' | '.join(row) + ' |\n')
print(root / 'bspline_benchmark_summary.tsv')
PY

echo "Done. RRM B-spline benchmark results: ${RESULT_DIR}"
