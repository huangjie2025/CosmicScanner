#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/repeatability_eb_smoke"
RESULT_DIR="${ROOT}/tests/results/repeatability_eb_smoke"
PREFIX="${DATA_DIR}/repeatability_eb"
OUT="${RESULT_DIR}/repeatability_eb"

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
prefix=Path(sys.argv[1]); rng=random.Random(20260622)
n=72; m=180; causal=22
samples=[(f"F{i+1:04d}", f"I{i+1:04d}") for i in range(n)]
genos=[]; variants=[]
for j in range(m):
    chrom=str((j%3)+1); p=0.12+0.38*((j*17)%101)/100.0
    row=[(1 if rng.random()<p else 0)+(1 if rng.random()<p else 0) for _ in range(n)]
    genos.append(row); variants.append((chrom,f"rs{j+1:06d}",1000+10*j,"A","C"))
u=[rng.gauss(0,0.25) for _ in range(n)]
with prefix.with_suffix('.fam').open('w') as f:
    for fid,iid in samples: f.write(f"{fid}\t{iid}\t0\t0\t1\t-9\n")
with prefix.with_suffix('.bim').open('w') as f:
    for v in variants: f.write('\t'.join(map(str,v))+'\n')
def code(g): return 0b00 if g==2 else (0b10 if g==1 else 0b11)
bps=(n+3)//4; bed=bytearray(b"\x6c\x1b\x01")
for row in genos:
    block=bytearray(bps)
    for i,g in enumerate(row): block[i//4]|=code(g)<<((i%4)*2)
    bed.extend(block)
prefix.with_suffix('.bed').write_bytes(bed)
with prefix.with_suffix('.long.pheno').open('w') as f:
    for i,(fid,iid) in enumerate(samples):
        for r in range(3):
            y=0.85*(genos[causal][i]-1.0)+u[i]+0.08*r+rng.gauss(0,0.35)
            f.write(f"{fid}\t{iid}\t{y:.8f}\n")
PY

"${BIN}" \
  --bfile "${PREFIX}" \
  --pheno "${PREFIX}.long.pheno" \
  --mpheno 1 \
  --repeatability-model \
  --gwas \
  --model lmm-eb \
  --adjust eb-mbonf \
  --no-loco \
  --maf 0.0 \
  --out "${OUT}"

python3 - "${OUT}.assoc.txt" "${OUT}.repeatability.eb.summary.txt" <<'PY'
import csv, math, sys
assoc, summary = sys.argv[1:]
rows=list(csv.DictReader(open(assoc), delimiter='\t'))
if len(rows)!=180: raise SystemExit(f"Expected 180 rows, got {len(rows)}")
for col in ['EB_BETA','EB_SE','EB_WALD','P','MARKER_LAMBDA','PRIOR_VAR','POST_VAR','EB_CONFIDENCE']:
    bad=sum(1 for r in rows if not math.isfinite(float(r[col])))
    if bad: raise SystemExit(f"{col} non-finite count={bad}")
top10=[r['SNP'] for r in sorted(rows, key=lambda r: float(r['P']))[:10]]
if 'rs000023' not in top10: raise SystemExit(f"causal rs000023 missing top10: {top10}")
vals={}
for line in open(summary):
    k,_,v=line.rstrip('\n').partition('\t'); vals[k]=v
meff=float(vals.get('effective_tests','nan'))
if not (0 <= meff <= 180): raise SystemExit(f"invalid meff {meff}")
print(f"Repeatability EB smoke passed: rows={len(rows)} effective_tests={meff:.6g} top1={top10[0]}")
PY
