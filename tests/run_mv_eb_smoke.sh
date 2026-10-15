#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/mv_eb_smoke"
RESULT_DIR="${ROOT}/tests/results/mv_eb_smoke"
PREFIX="${DATA_DIR}/mv_eb"
OUT="${RESULT_DIR}/mv_eb"

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
prefix=Path(sys.argv[1]); rng=random.Random(20260623)
n=90; m=220; causal=44
samples=[(f"F{i+1:04d}", f"I{i+1:04d}") for i in range(n)]
genos=[]; vars=[]
for j in range(m):
    chrom=str((j%4)+1); p=0.10+0.40*((j*31)%101)/100.0
    row=[(1 if rng.random()<p else 0)+(1 if rng.random()<p else 0) for _ in range(n)]
    genos.append(row); vars.append((chrom,f"rs{j+1:06d}",1000+10*j,"A","C"))
with prefix.with_suffix('.fam').open('w') as f:
    for fid,iid in samples: f.write(f"{fid}\t{iid}\t0\t0\t1\t-9\n")
with prefix.with_suffix('.bim').open('w') as f:
    for v in vars: f.write('\t'.join(map(str,v))+'\n')
def code(g): return 0b00 if g==2 else (0b10 if g==1 else 0b11)
bps=(n+3)//4; bed=bytearray(b"\x6c\x1b\x01")
for row in genos:
    block=bytearray(bps)
    for i,g in enumerate(row): block[i//4]|=code(g)<<((i%4)*2)
    bed.extend(block)
prefix.with_suffix('.bed').write_bytes(bed)
with prefix.with_suffix('.pheno').open('w') as f:
    for i,(fid,iid) in enumerate(samples):
        g=genos[causal][i]-1.0
        e1=rng.gauss(0,0.55); e2=0.35*e1+rng.gauss(0,0.50)
        y1=0.75*g+e1; y2=-0.55*g+e2
        f.write(f"{fid}\t{iid}\t{y1:.8f}\t{y2:.8f}\n")
PY

"${BIN}" \
  --bfile "${PREFIX}" \
  --pheno "${PREFIX}.pheno" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method eb \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${OUT}"

python3 - "${OUT}.assoc.mv.txt" "${OUT}.eb.mv.summary.txt" <<'PY'
import csv, math, sys
assoc, summary = sys.argv[1:]
rows=list(csv.DictReader(open(assoc), delimiter='\t'))
if len(rows)!=220: raise SystemExit(f"Expected 220 rows, got {len(rows)}")
for col in ['EB_WALD','P','EB_BETA_1','EB_SE_1','EB_BETA_2','EB_SE_2','MARKER_LAMBDA','EB_CONFIDENCE']:
    bad=sum(1 for r in rows if not math.isfinite(float(r[col])))
    if bad: raise SystemExit(f"{col} non-finite count={bad}")
top10=[r['SNP'] for r in sorted(rows, key=lambda r: float(r['P']))[:10]]
if 'rs000045' not in top10: raise SystemExit(f"causal rs000045 missing top10: {top10}")
vals={}
for line in open(summary):
    k,_,v=line.rstrip('\n').partition('\t'); vals[k]=v
meff=float(vals.get('effective_tests','nan'))
if not (0 <= meff <= 2*220): raise SystemExit(f"invalid meff {meff}")
print(f"Multivariate EB smoke passed: rows={len(rows)} effective_tests={meff:.6g} top1={top10[0]}")
PY
