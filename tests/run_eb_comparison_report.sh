#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/eb_lmm_smoke"
RESULT_DIR="${ROOT}/tests/results/eb_comparison_report"
PREFIX="${DATA_DIR}/eb_lmm_smoke"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${RESULT_DIR}"
bash "${ROOT}/tests/run_eb_lmm_smoke.sh" >/dev/null

"${BIN}" --bfile "${PREFIX}" --pheno "${PREFIX}.pheno" --mpheno 1 --gwas --no-loco --model lmm --maf 0.0 --out "${RESULT_DIR}/standard" >/dev/null
"${BIN}" --bfile "${PREFIX}" --pheno "${PREFIX}.pheno" --mpheno 1 --gwas --no-loco --model lmm-eb --adjust eb-mbonf --maf 0.0 --out "${RESULT_DIR}/eb" >/dev/null

python3 - "${RESULT_DIR}" <<'PY'
import csv, math, sys
from pathlib import Path

out = Path(sys.argv[1])
std = list(csv.DictReader(open(out/'standard.assoc.txt'), delimiter='\t'))
eb = list(csv.DictReader(open(out/'eb.assoc.txt'), delimiter='\t'))
if len(std) != len(eb):
    raise SystemExit(f'row mismatch: standard={len(std)} eb={len(eb)}')
std_by = {r['SNP']: r for r in std}
eb_by = {r['SNP']: r for r in eb}
overlap = sorted(set(std_by) & set(eb_by))
if len(overlap) != len(std):
    raise SystemExit(f'overlap mismatch: overlap={len(overlap)} standard={len(std)}')

def top(rows, k=10):
    return [r['SNP'] for r in sorted(rows, key=lambda r: float(r['P']))[:k]]

std_top = top(std)
eb_top = top(eb)
top10_overlap = len(set(std_top) & set(eb_top))
rows = []
for snp in overlap:
    s = std_by[snp]
    e = eb_by[snp]
    rows.append({
        'SNP': snp,
        'CHR': s['CHR'],
        'BP': s['BP'],
        'standard_P': float(s['P']),
        'standard_BETA': float(s['BETA']),
        'standard_SE': float(s['SE']),
        'EB_P': float(e['P']),
        'EB_BETA': float(e['EB_BETA']),
        'EB_SE': float(e['EB_SE']),
        'MARKER_LAMBDA': float(e['MARKER_LAMBDA']),
        'EB_CONFIDENCE': float(e['EB_CONFIDENCE']),
    })
for r in rows:
    for k, v in r.items():
        if isinstance(v, float) and not math.isfinite(v):
            raise SystemExit(f'non-finite {k} for {r["SNP"]}')

with (out/'standard_vs_eb.tsv').open('w', newline='') as f:
    fields = ['SNP','CHR','BP','standard_P','standard_BETA','standard_SE','EB_P','EB_BETA','EB_SE','MARKER_LAMBDA','EB_CONFIDENCE']
    w = csv.DictWriter(f, fields, delimiter='\t')
    w.writeheader(); w.writerows(sorted(rows, key=lambda r: r['EB_P']))

summary = {
    'markers': len(overlap),
    'standard_top1': top(std, 1)[0],
    'eb_top1': top(eb, 1)[0],
    'top10_overlap': top10_overlap,
    'mean_eb_confidence': sum(r['EB_CONFIDENCE'] for r in rows) / len(rows),
    'median_marker_lambda': sorted(r['MARKER_LAMBDA'] for r in rows)[len(rows)//2],
}
with (out/'standard_vs_eb.summary.tsv').open('w') as f:
    for k, v in summary.items():
        f.write(f'{k}\t{v}\n')
with (out/'standard_vs_eb.summary.md').open('w') as f:
    f.write('# Standard LMM vs EB-LMM Summary\n\n')
    f.write('| Metric | Value |\n|---|---:|\n')
    for k, v in summary.items():
        f.write(f'| {k} | {v} |\n')
print(f'EB comparison report passed: {out / "standard_vs_eb.summary.tsv"}')
PY
