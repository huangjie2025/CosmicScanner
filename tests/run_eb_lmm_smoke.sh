#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/eb_lmm_smoke"
RESULT_DIR="${ROOT}/tests/results/eb_lmm_smoke"
PREFIX="${DATA_DIR}/eb_lmm_smoke"
OUT="${RESULT_DIR}/eb_lmm_smoke"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}"

python3 - "${PREFIX}" <<'PY'
import math
import random
import sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260621)
n = 96
m = 240
samples = [(f"F{i+1:04d}", f"I{i+1:04d}", 0 if i < n // 2 else 1) for i in range(n)]
genotypes = []
variants = []
for j in range(m):
    chrom = str((j % 3) + 1)
    p_base = 0.12 + 0.38 * ((j * 29) % 101) / 100.0
    shift = 0.07 * math.sin(j * 0.23)
    row = []
    for _, _, pop in samples:
        p = min(0.95, max(0.02, p_base + (shift if pop else -shift)))
        row.append((1 if rng.random() < p else 0) + (1 if rng.random() < p else 0))
    variants.append((chrom, f"rs{j+1:06d}", 1000 + 10 * j, "A", "C"))
    genotypes.append(row)

causal_idx = 36
phenotypes = []
for i, (_, _, pop) in enumerate(samples):
    y = 0.25 * pop + 0.9 * (genotypes[causal_idx][i] - 1.0) + rng.gauss(0, 0.55)
    phenotypes.append(y)

prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix(".fam").open("w") as f:
    for (fid, iid, _), y in zip(samples, phenotypes):
        f.write(f"{fid}\t{iid}\t0\t0\t1\t{y:.8f}\n")
with prefix.with_suffix(".pheno").open("w") as f:
    for (fid, iid, _), y in zip(samples, phenotypes):
        f.write(f"{fid}\t{iid}\t{y:.8f}\n")
with prefix.with_suffix(".bim").open("w") as f:
    for chrom, snp, bp, a1, a2 in variants:
        f.write(f"{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n")

def plink_code(g):
    return 0b00 if g == 2 else (0b10 if g == 1 else 0b11)

bytes_per_snp = (n + 3) // 4
bed = bytearray(b"\x6c\x1b\x01")
for row in genotypes:
    block = bytearray(bytes_per_snp)
    for i, g in enumerate(row):
        block[i // 4] |= plink_code(g) << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix(".bed").write_bytes(bed)
with prefix.with_suffix(".truth.txt").open("w") as f:
    f.write("causal_snp\trs000037\n")
PY

"${BIN}" \
  --bfile "${PREFIX}" \
  --pheno "${PREFIX}.pheno" \
  --mpheno 1 \
  --gwas \
  --no-loco \
  --model lmm-eb \
  --maf 0.0 \
  --out "${OUT}"

python3 - "${OUT}.assoc.txt" "${OUT}.eb.summary.txt" <<'PY'
import csv
import math
import sys

assoc, summary = sys.argv[1:]
rows = list(csv.DictReader(open(assoc), delimiter="\t"))
if len(rows) != 240:
    raise SystemExit(f"Expected 240 EB rows, got {len(rows)}")
required = ["EB_BETA", "EB_SE", "EB_WALD", "P", "MARKER_LAMBDA", "PRIOR_VAR", "POST_VAR", "EB_CONFIDENCE"]
for col in required:
    bad = [r["SNP"] for r in rows if not math.isfinite(float(r[col]))]
    if bad:
        raise SystemExit(f"Column {col} has non-finite values, first={bad[0]}")
meff = sum(float(r["EB_CONFIDENCE"]) for r in rows)
if not (0.0 <= meff <= len(rows) + 1e-8):
    raise SystemExit(f"Invalid effective test sum: {meff}")
top10 = [r["SNP"] for r in sorted(rows, key=lambda r: float(r["P"]))[:10]]
if "rs000037" not in top10:
    raise SystemExit(f"Causal SNP rs000037 missing from EB top10: {top10}")
summary_vals = {}
for line in open(summary):
    k, _, v = line.rstrip("\n").partition("\t")
    summary_vals[k] = v
if "effective_tests" not in summary_vals:
    raise SystemExit("EB summary missing effective_tests")
print(f"EB-LMM smoke passed: rows={len(rows)} effective_tests={meff:.6g} top1={top10[0]}")
PY
