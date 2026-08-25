#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/eb_loco_validation"
RESULT_DIR="${ROOT}/tests/results/eb_loco_validation"
PREFIX="${DATA_DIR}/eb_loco"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}"

python3 - "${PREFIX}" <<'PY'
import math, random, sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260625)
n, m = 90, 300
causal = 73
samples = [(f'F{i+1:04d}', f'I{i+1:04d}', 0 if i < n//2 else 1) for i in range(n)]
genos, variants = [], []
for j in range(m):
    chrom = str((j % 5) + 1)
    base = 0.10 + 0.38 * ((j * 31) % 101) / 100.0
    shift = 0.05 * math.sin(j * 0.17)
    row = []
    for _, _, pop in samples:
        p = min(0.95, max(0.02, base + (shift if pop else -shift)))
        row.append((1 if rng.random() < p else 0) + (1 if rng.random() < p else 0))
    genos.append(row)
    variants.append((chrom, f'rs{j+1:06d}', 1000 + 10*j, 'A', 'C'))

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
        y = 0.18 * pop + 0.85 * g + rng.gauss(0, 0.65)
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
with prefix.with_suffix('.truth.txt').open('w') as f:
    f.write(f'causal_snp\trs{causal+1:06d}\n')
PY

"${BIN}" --bfile "${PREFIX}" --pheno "${PREFIX}.pheno" --mpheno 1 --gwas --model lmm-eb --adjust eb-mbonf --maf 0.0 --out "${RESULT_DIR}/loco" >/dev/null

python3 - "${RESULT_DIR}" <<'PY'
import csv, math, sys
from pathlib import Path

result_dir = Path(sys.argv[1])
assoc = result_dir/'loco.assoc.txt'
summary = result_dir/'loco.eb.summary.txt'
rows = list(csv.DictReader(open(assoc), delimiter='\t'))
if len(rows) != 300:
    raise SystemExit(f'expected 300 rows, got {len(rows)}')
chroms = {r['CHR'] for r in rows}
if chroms != {'1','2','3','4','5'}:
    raise SystemExit(f'unexpected chromosome set: {sorted(chroms)}')
conf_sum = 0.0
for r in rows:
    for col in ['P','LOG10P','MARKER_LAMBDA','EB_CONFIDENCE']:
        v = float(r[col])
        if not math.isfinite(v):
            raise SystemExit(f'non-finite {col}')
    c = float(r['EB_CONFIDENCE'])
    if c < -1e-12 or c > 1.0 + 1e-8:
        raise SystemExit(f'EB_CONFIDENCE out of range: {c}')
    conf_sum += max(0.0, min(1.0, c))
vals = {}
for line in open(summary):
    k, _, v = line.rstrip('\n').partition('\t')
    vals[k] = v
meff = float(vals['effective_tests'])
if int(float(vals['tested_markers'])) != 300:
    raise SystemExit('tested_markers mismatch')
if abs(conf_sum - meff) > max(1e-6, 1e-6 * max(1.0, meff)):
    raise SystemExit(f'effective_tests mismatch summary={meff} assoc_sum={conf_sum}')
top10 = [r['SNP'] for r in sorted(rows, key=lambda r: float(r['P']))[:10]]
if 'rs000074' not in top10:
    raise SystemExit(f'causal rs000074 missing from top10: {top10}')
print(f'EB LOCO validation passed: rows={len(rows)} effective_tests={meff:.6g} top1={top10[0]}')
PY
