#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/psr_normal_special_case"
RESULT_DIR="${ROOT}/tests/results/psr_normal_special_case"
PREFIX="${DATA_DIR}/psr_normal"

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
rng = random.Random(20260623)
n = 80
m = 120
samples = [(f"F{i+1:04d}", f"I{i+1:04d}", 0 if i < n // 2 else 1) for i in range(n)]
variants = []
genotypes = []
for j in range(m):
    chrom = str((j % 3) + 1)
    snp = f"norm{j+1:06d}"
    bp = 2000 + j * 10
    p = 0.12 + 0.36 * ((j * 13) % 97) / 96.0
    row = []
    for _, _, pop in samples:
        pp = min(0.95, max(0.02, p + (0.04 if pop else -0.02)))
        row.append((1 if rng.random() < pp else 0) + (1 if rng.random() < pp else 0))
    variants.append((chrom, snp, bp, "A", "C"))
    genotypes.append(row)

causal = {17: 0.55, 88: -0.45}
y = []
for i, (_, _, pop) in enumerate(samples):
    val = 0.25 * pop + rng.gauss(0.0, 0.7)
    for idx, eff in causal.items():
        val += eff * (genotypes[idx][i] - 1.0)
    y.append(val)

prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for fid, iid, _ in samples:
        f.write(f"{fid}\t{iid}\t0\t0\t1\t-9\n")
with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for chrom, snp, bp, a1, a2 in variants:
        f.write(f"{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n")
with prefix.with_suffix(".pheno").open("w", encoding="utf-8") as f:
    for (fid, iid, _), val in zip(samples, y):
        f.write(f"{fid}\t{iid}\t{val:.12g}\n")

def plink_code(dosage):
    if dosage == 2:
        return 0b00
    if dosage == 1:
        return 0b10
    if dosage == 0:
        return 0b11
    return 0b01

bytes_per_snp = (n + 3) // 4
bed = bytearray(b"\x6c\x1b\x01")
for row in genotypes:
    block = bytearray(bytes_per_snp)
    for i, g in enumerate(row):
        block[i // 4] |= plink_code(g) << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix(".bed").write_bytes(bed)
PY

standard_out="${RESULT_DIR}/normal_standard"
special_out="${RESULT_DIR}/normal_special_case"

"${BIN}" --bfile "${PREFIX}" --gwas --no-loco --maf 0.0 --pheno "${PREFIX}.pheno" --out "${standard_out}" >/dev/null
"${BIN}" --bfile "${PREFIX}" --gwas --no-loco --maf 0.0 --pheno "${PREFIX}.pheno" --out "${special_out}" >/dev/null

python3 - "${standard_out}.assoc.txt" "${special_out}.assoc.txt" "${RESULT_DIR}" <<'PY'
import csv
import math
import sys
from pathlib import Path

std_path = Path(sys.argv[1])
sp_path = Path(sys.argv[2])
out_dir = Path(sys.argv[3])

def read_assoc(path):
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return {r["SNP"]: r for r in rows}

std = read_assoc(std_path)
sp = read_assoc(sp_path)
common = sorted(set(std) & set(sp))
if not common:
    raise SystemExit("No overlapping SNPs in normal special-case comparison")

fields = ["BETA", "SE", "P", "LOG10P"]
max_diff = {k: 0.0 for k in fields}
detail = out_dir / "normal_special_case_detail.tsv"
with detail.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["trait_type", "engine", "snp", "effect", "se", "wald_or_stat", "p_value", "log10p", "runtime"])
    for snp in common:
        a = std[snp]
        b = sp[snp]
        for field in fields:
            max_diff[field] = max(max_diff[field], abs(float(a[field]) - float(b[field])))
        beta = float(b["BETA"])
        se = float(b["SE"])
        wald = (beta / se) ** 2 if se > 0 else float("nan")
        w.writerow(["normal", "cosmicscanner_lmm_special_case", snp, beta, se, wald, b["P"], b["LOG10P"], "NA"])

summary = out_dir / "normal_special_case_summary.tsv"
with summary.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["trait_type", "engine", "n_snps", "max_abs_beta_diff", "max_abs_se_diff", "max_abs_p_diff", "max_abs_log10p_diff", "validation_status"])
    status = "PASS" if max(max_diff.values()) < 1e-8 else "FAIL"
    w.writerow(["normal", "cosmicscanner_lmm_special_case", len(common), max_diff["BETA"], max_diff["SE"], max_diff["P"], max_diff["LOG10P"], status])
if max(max_diff.values()) >= 1e-8:
    raise SystemExit(f"Normal special-case mismatch: {max_diff}")
print(f"[PSR-GWAS] Normal special-case summary written to {summary}")
PY

echo "[PSR-GWAS] Normal special-case validation passed."
