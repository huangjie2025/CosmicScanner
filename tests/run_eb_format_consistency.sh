#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
PLINK2="${PLINK2_BIN:-/disk228/huangjie/software/miniconda3/envs/benchmark/bin/plink2}"
DATA_DIR="${ROOT}/tests/data/eb_format_consistency"
RESULT_DIR="${ROOT}/tests/results/eb_format_consistency"
PREFIX="${DATA_DIR}/ebfmt"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi
if [[ ! -x "${PLINK2}" ]]; then
  echo "Missing PLINK2 executable: ${PLINK2}" >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}/plink" "${RESULT_DIR}/pgen" "${RESULT_DIR}/bgen"

python3 - "${PREFIX}" <<'PY'
import math, random, sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260627)
n, m = 100, 280
causal = 64
samples = [(f'F{i+1:04d}', f'I{i+1:04d}', 0 if i < n//2 else 1) for i in range(n)]
genos, variants = [], []
for j in range(m):
    chrom = str(1 + j // max(1, m // 4))
    if int(chrom) > 4:
        chrom = '4'
    base = 0.10 + 0.36 * ((j * 41) % 101) / 100.0
    shift = 0.06 * math.sin(j * 0.21)
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
        y = 0.20 * pop + 0.75 * g + rng.gauss(0, 0.70)
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

"${PLINK2}" --bfile "${PREFIX}" --make-pgen --out "${RESULT_DIR}/pgen/ebfmt" >/dev/null 2>&1
"${PLINK2}" --bfile "${PREFIX}" --export bgen-1.2 --out "${RESULT_DIR}/bgen/ebfmt" >/dev/null 2>&1

"${BIN}" --bfile "${PREFIX}" --pheno "${PREFIX}.pheno" --mpheno 1 --gwas --model lmm-eb --adjust eb-mbonf --no-loco --maf 0.0 --out "${RESULT_DIR}/plink/ebfmt" >/dev/null
"${BIN}" --pfile "${RESULT_DIR}/pgen/ebfmt" --pheno "${PREFIX}.pheno" --mpheno 1 --gwas --model lmm-eb --adjust eb-mbonf --no-loco --maf 0.0 --out "${RESULT_DIR}/pgen/ebfmt" >/dev/null
"${BIN}" --bgen "${RESULT_DIR}/bgen/ebfmt.bgen" --sample "${RESULT_DIR}/bgen/ebfmt.sample" --pheno "${PREFIX}.pheno" --mpheno 1 --gwas --model lmm-eb --adjust eb-mbonf --no-loco --maf 0.0 --out "${RESULT_DIR}/bgen/ebfmt" >/dev/null

python3 - "${RESULT_DIR}" <<'PY'
import csv, math, sys
from pathlib import Path

result = Path(sys.argv[1])
cols = ['EB_BETA','EB_SE','EB_WALD','P','LOG10P','MARKER_LAMBDA','PRIOR_VAR','POST_VAR','EB_CONFIDENCE']

def read_assoc(path):
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    return {r['SNP']: r for r in rows}

def read_summary(path):
    vals = {}
    for line in open(path):
        k, _, v = line.rstrip('\n').partition('\t')
        vals[k] = v
    return vals

plink = read_assoc(result/'plink/ebfmt.assoc.txt')
pgen = read_assoc(result/'pgen/ebfmt.assoc.txt')
bgen = read_assoc(result/'bgen/ebfmt.assoc.txt')
if len(plink) != 280 or len(pgen) != 280 or len(bgen) != 280:
    raise SystemExit(f'row count mismatch: plink={len(plink)} pgen={len(pgen)} bgen={len(bgen)}')

def compare(ref, test, label, tol):
    snps = sorted(set(ref) & set(test))
    if len(snps) != len(ref):
        raise SystemExit(f'{label}: overlap mismatch {len(snps)} vs {len(ref)}')
    max_diff = {c: 0.0 for c in cols}
    for snp in snps:
        for c in cols:
            rv = float(ref[snp][c]); tv = float(test[snp][c])
            if not math.isfinite(rv) or not math.isfinite(tv):
                raise SystemExit(f'{label}: non-finite {c} at {snp}')
            d = abs(rv - tv)
            max_diff[c] = max(max_diff[c], d)
            if d > tol:
                raise SystemExit(f'{label}: {c} diff too large at {snp}: ref={rv} test={tv} diff={d} tol={tol}')
    return max_diff

pgen_diff = compare(plink, pgen, 'PGEN', 2e-6)
bgen_diff = compare(plink, bgen, 'BGEN', 2e-5)

plink_sum = read_summary(result/'plink/ebfmt.eb.summary.txt')
pgen_sum = read_summary(result/'pgen/ebfmt.eb.summary.txt')
bgen_sum = read_summary(result/'bgen/ebfmt.eb.summary.txt')
for label, vals in [('PGEN', pgen_sum), ('BGEN', bgen_sum)]:
    if int(float(vals['tested_markers'])) != 280:
        raise SystemExit(f'{label}: tested_markers mismatch')
    d = abs(float(vals['effective_tests']) - float(plink_sum['effective_tests']))
    if d > (2e-6 if label == 'PGEN' else 2e-5):
        raise SystemExit(f'{label}: effective_tests diff too large: {d}')

out = result/'eb_format_consistency.summary.tsv'
with out.open('w') as f:
    f.write('format\tmetric\tmax_abs_diff\n')
    for fmt, diffs in [('PGEN', pgen_diff), ('BGEN', bgen_diff)]:
        for k, v in diffs.items():
            f.write(f'{fmt}\t{k}\t{v}\n')
print(f'EB format consistency passed: {out}')
PY
