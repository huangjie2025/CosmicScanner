#!/usr/bin/env bash
set -euo pipefail

ROOT="/disk228/huangjie/PanOmics/CosmicScanner"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
RESULT_DIR="${ROOT}/tests/results/multivariate_validation_matrix"
DATA_PREFIX="${RESULT_DIR}/data/mv_matrix"
SUMMARY="${RESULT_DIR}/validation_matrix.summary.txt"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${RESULT_DIR}/data"
rm -f "${SUMMARY}"

python3 - "${DATA_PREFIX}" <<'PY'
import math
import random
import sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260604)
n = 128
m = 720
n_traits = 3
n_pop = 4
missing_rate = 0.025

samples = []
for i in range(n):
    pop = i % n_pop
    fid = f"F{i+1:04d}"
    iid = f"I{i+1:04d}"
    age = 20.0 + (i % 31) + rng.gauss(0.0, 1.5)
    pc1 = (pop - 1.5) / 2.0 + rng.gauss(0.0, 0.05)
    samples.append((fid, iid, pop, age, pc1))

variants = []
genotypes = []
for j in range(m):
    chrom = str((j % 8) + 1)
    snp = f"mvm{j+1:06d}"
    bp = 2000 + j * 13
    base = 0.06 + 0.44 * ((j * 37) % 103) / 102.0
    row = []
    for _, _, pop, _, _ in samples:
        shift = (pop - 1.5) * 0.035 * math.sin(j * 0.07)
        p = min(0.96, max(0.015, base + shift))
        if rng.random() < missing_rate:
            row.append(None)
        else:
            row.append((1 if rng.random() < p else 0) + (1 if rng.random() < p else 0))
    variants.append((chrom, snp, bp, "A", "C"))
    genotypes.append(row)

causal = {
    23: (0.35, 0.12, -0.08),
    118: (-0.25, -0.28, 0.16),
    263: (0.08, 0.34, 0.25),
    511: (-0.20, 0.18, 0.38),
}

phenos = []
for i, (_, _, pop, age, pc1) in enumerate(samples):
    shared1 = rng.gauss(0.0, 0.45)
    shared2 = rng.gauss(0.0, 0.25)
    y = [0.0] * n_traits
    y[0] = 0.015 * age + 0.30 * pc1 + 0.18 * pop + shared1 + rng.gauss(0.0, 0.50)
    y[1] = -0.010 * age - 0.20 * pc1 + 0.12 * pop + 0.70 * shared1 + rng.gauss(0.0, 0.52)
    y[2] = 0.006 * age + 0.45 * pc1 - 0.05 * pop + 0.35 * shared1 + shared2 + rng.gauss(0.0, 0.48)
    for snp_idx, effects in causal.items():
        g = genotypes[snp_idx][i]
        centered = 0.0 if g is None else g - 1.0
        for t, b in enumerate(effects):
            y[t] += b * centered
    phenos.append(y)

prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for (fid, iid, _, _, _), y in zip(samples, phenos):
        f.write(f"{fid}\t{iid}\t0\t0\t1\t{y[0]:.8f}\n")
with prefix.with_suffix(".pheno").open("w", encoding="utf-8") as f:
    f.write("FID\tIID\tTrait1\tTrait2\tTrait3\n")
    for (fid, iid, _, _, _), y in zip(samples, phenos):
        f.write(f"{fid}\t{iid}\t" + "\t".join(f"{v:.8f}" for v in y) + "\n")
with prefix.with_suffix(".qcovar").open("w", encoding="utf-8") as f:
    f.write("FID\tIID\tAge\tPC1\n")
    for fid, iid, _, age, pc1 in samples:
        f.write(f"{fid}\t{iid}\t{age:.8f}\t{pc1:.8f}\n")
with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for chrom, snp, bp, a1, a2 in variants:
        f.write(f"{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n")

def plink_code(dosage):
    if dosage is None:
        return 0b01
    if dosage == 2:
        return 0b00
    if dosage == 1:
        return 0b10
    return 0b11

bytes_per_snp = (n + 3) // 4
bed = bytearray(b"\x6c\x1b\x01")
for row in genotypes:
    block = bytearray(bytes_per_snp)
    for i, g in enumerate(row):
        block[i // 4] |= plink_code(g) << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix(".bed").write_bytes(bed)
PY

run_case() {
  local label="$1"
  shift
  python3 - "$label" "$SUMMARY" "$@" <<'PY'
import subprocess
import sys
import time
label, summary = sys.argv[1], sys.argv[2]
cmd = sys.argv[3:]
t0 = time.time()
proc = subprocess.run(cmd)
elapsed = time.time() - t0
with open(summary, "a", encoding="utf-8") as f:
    f.write(f"{label}_status\t{'ok' if proc.returncode == 0 else 'failed'}\n")
    f.write(f"time_{label}\t{elapsed:.6f}\n")
if proc.returncode != 0:
    raise SystemExit(proc.returncode)
PY
}

run_case emmv_3trait_cov_missing_no_loco "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${DATA_PREFIX}.pheno" \
  --qcovar "${DATA_PREFIX}.qcovar" \
  --mpheno-list 1,2,3 \
  --multivariate \
  --mv-method emmv \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${RESULT_DIR}/emmv_3trait_cov_missing_no_loco"

run_case emmv_3trait_cov_missing_loco "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${DATA_PREFIX}.pheno" \
  --qcovar "${DATA_PREFIX}.qcovar" \
  --mpheno-list 1,2,3 \
  --multivariate \
  --mv-method emmv \
  --gwas \
  --maf 0.0 \
  --out "${RESULT_DIR}/emmv_3trait_cov_missing_loco"

run_case gemma_full_fast_3trait_cov_missing "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${DATA_PREFIX}.pheno" \
  --qcovar "${DATA_PREFIX}.qcovar" \
  --mpheno-list 1,2,3 \
  --multivariate \
  --mv-method gemma-full \
  --mv-test wald,lrt,score \
  --mv-refit-p 0.02 \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${RESULT_DIR}/gemma_full_fast_3trait_cov_missing"

python3 - "${RESULT_DIR}" "${SUMMARY}" <<'PY'
import csv
import math
import sys
from pathlib import Path
result_dir = Path(sys.argv[1])
summary = Path(sys.argv[2])
with summary.open("a", encoding="utf-8") as f:
    for stem in [
        "emmv_3trait_cov_missing_no_loco.assoc.mv.txt",
        "emmv_3trait_cov_missing_loco.assoc.mv.txt",
        "gemma_full_fast_3trait_cov_missing.assoc.mv.txt",
    ]:
        path = result_dir / stem
        rows = list(csv.DictReader(path.open("r", encoding="utf-8"), delimiter="\t")) if path.exists() else []
        finite_p = sum(1 for row in rows if row.get("P") not in (None, "", "nan", "NA"))
        f.write(f"{stem}_rows\t{len(rows)}\n")
        f.write(f"{stem}_finite_p\t{finite_p}\n")
        if rows and "BETA_3" not in rows[0]:
            raise SystemExit(f"Missing BETA_3 in {stem}")
print(summary.read_text(encoding="utf-8"))
PY

echo "Done. Multivariate validation matrix complete: ${SUMMARY}"
