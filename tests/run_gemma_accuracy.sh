#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/gemma_accuracy"
RESULT_DIR="${ROOT}/tests/results/gemma_accuracy"
DATA_PREFIX="${DATA_DIR}/gemma_accuracy"
COSMIC_PREFIX="${RESULT_DIR}/cosmic/gemma_accuracy"
GEMMA_WORK="${RESULT_DIR}/gemma"
LMM_NULL_SOLVER="${LMM_NULL_SOLVER:-brent}"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

resolve_gemma() {
  if [[ -n "${GEMMA_BIN:-}" ]]; then
    printf '%s\n' "${GEMMA_BIN}"
    return 0
  fi
  if command -v gemma >/dev/null 2>&1; then
    command -v gemma
    return 0
  fi
  local conda_bin="/disk228/huangjie/software/miniconda3/envs/benchmark/bin/gemma"
  if [[ -x "${conda_bin}" ]]; then
    printf '%s\n' "${conda_bin}"
    return 0
  fi
  local repo_bin="/disk228/huangjie/PanOmics/References_and_Tools/GEMMA-master/bin/gemma"
  if [[ -x "${repo_bin}" ]]; then
    printf '%s\n' "${repo_bin}"
    return 0
  fi
  return 1
}

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  echo "Build it first or set COSMICSCANNER_BIN=/path/to/cosmicscanner." >&2
  exit 2
fi

if ! GEMMA_EXE="$(resolve_gemma)"; then
  echo "Missing GEMMA executable." >&2
  echo "Install GEMMA or set GEMMA_BIN=/path/to/gemma." >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}/cosmic" "${GEMMA_WORK}"
rm -rf "${GEMMA_WORK}/output"

python3 - "${DATA_PREFIX}" <<'PY'
import math
import random
import struct
import sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260531)
n = 160
m = 800
n_pop = 2

samples = []
for i in range(n):
    pop = 0 if i < n // 2 else 1
    fid = f"F{i+1:04d}"
    iid = f"I{i+1:04d}"
    samples.append((fid, iid, pop))

variants = []
genotypes = []
for j in range(m):
    chrom = str((j % 4) + 1)
    snp = f"rs{j+1:06d}"
    bp = 1000 + j * 10
    base_p = 0.08 + 0.42 * ((j * 37) % 100) / 99.0
    pop_shift = 0.10 * math.sin(j * 0.17)
    p0 = min(0.95, max(0.02, base_p - pop_shift / 2.0))
    p1 = min(0.95, max(0.02, base_p + pop_shift / 2.0))
    row = []
    for _, _, pop in samples:
        p = p0 if pop == 0 else p1
        g = (1 if rng.random() < p else 0) + (1 if rng.random() < p else 0)
        row.append(g)
    variants.append((chrom, snp, bp, "A", "C"))
    genotypes.append(row)

causal = {12: 0.45, 211: -0.35, 509: 0.30}
phenotypes = []
for i, (_, _, pop) in enumerate(samples):
    y = 0.65 * pop
    y += 0.25 * math.sin(i * 0.13)
    for snp_idx, eff in causal.items():
        y += eff * (genotypes[snp_idx][i] - 1.0)
    y += rng.gauss(0.0, 0.75)
    phenotypes.append(y)

prefix.parent.mkdir(parents=True, exist_ok=True)

with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for (fid, iid, _), y in zip(samples, phenotypes):
        f.write(f"{fid}\t{iid}\t0\t0\t1\t{y:.8f}\n")

with prefix.with_suffix(".pheno").open("w", encoding="utf-8") as f:
    for (fid, iid, _), y in zip(samples, phenotypes):
        f.write(f"{fid}\t{iid}\t{y:.8f}\n")

with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for chrom, snp, bp, a1, a2 in variants:
        f.write(f"{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n")

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

with prefix.with_suffix(".truth.txt").open("w", encoding="utf-8") as f:
    f.write("SNP\tTRUE_BETA\n")
    for idx, eff in causal.items():
        f.write(f"rs{idx+1:06d}\t{eff}\n")
PY

echo "[GEMMA-accuracy] Running CosmicScanner exact no-LOCO GWAS..."
"${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${DATA_PREFIX}.pheno" \
  --mpheno 1 \
  --gwas \
  --no-loco \
  --lmm-null-solver "${LMM_NULL_SOLVER}" \
  --maf 0.0 \
  --out "${COSMIC_PREFIX}"

echo "[GEMMA-accuracy] Running GEMMA kinship..."
(
  cd "${GEMMA_WORK}"
  "${GEMMA_EXE}" -bfile "${DATA_PREFIX}" -gk 1 -o gemma_accuracy_kinship
)

KINSHIP="${GEMMA_WORK}/output/gemma_accuracy_kinship.cXX.txt"
if [[ ! -f "${KINSHIP}" ]]; then
  echo "Missing GEMMA kinship output: ${KINSHIP}" >&2
  exit 3
fi

echo "[GEMMA-accuracy] Running GEMMA LMM..."
(
  cd "${GEMMA_WORK}"
  "${GEMMA_EXE}" -bfile "${DATA_PREFIX}" -k "${KINSHIP}" -lmm 4 -o gemma_accuracy_lmm
)

COSMIC_ASSOC="${COSMIC_PREFIX}.assoc.txt"
GEMMA_ASSOC="${GEMMA_WORK}/output/gemma_accuracy_lmm.assoc.txt"
REPORT="${RESULT_DIR}/gemma_accuracy.compare.tsv"
SUMMARY="${RESULT_DIR}/gemma_accuracy.summary.txt"

python3 - "${COSMIC_ASSOC}" "${GEMMA_ASSOC}" "${REPORT}" "${SUMMARY}" <<'PY'
import csv
import math
import os
import sys

cosmic_path, gemma_path, report_path, summary_path = sys.argv[1:]

def read_table(path):
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise SystemExit(f"No rows in {path}")
    return rows

def fval(row, *names):
    for name in names:
        if name in row and row[name] not in ("", "NA", "nan", "-nan"):
            return float(row[name])
    lower = {k.lower(): k for k in row}
    for name in names:
        key = lower.get(name.lower())
        if key and row[key] not in ("", "NA", "nan", "-nan"):
            return float(row[key])
    return math.nan

def sval(row, *names):
    for name in names:
        if name in row:
            return row[name]
    lower = {k.lower(): k for k in row}
    for name in names:
        key = lower.get(name.lower())
        if key:
            return row[key]
    return ""

def corr(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 3:
        return math.nan
    mx = sum(x for x, _ in pairs) / len(pairs)
    my = sum(y for _, y in pairs) / len(pairs)
    vx = sum((x - mx) ** 2 for x, _ in pairs)
    vy = sum((y - my) ** 2 for _, y in pairs)
    if vx <= 0 or vy <= 0:
        return math.nan
    cov = sum((x - mx) * (y - my) for x, y in pairs)
    return cov / math.sqrt(vx * vy)

cosmic = {sval(r, "SNP", "rs"): r for r in read_table(cosmic_path)}
gemma = {sval(r, "rs", "SNP"): r for r in read_table(gemma_path)}
snps = sorted(set(cosmic) & set(gemma))
if not snps:
    raise SystemExit("No overlapping SNPs between CosmicScanner and GEMMA outputs")

rows = []
unaligned = 0
for snp in snps:
    c = cosmic[snp]
    g = gemma[snp]
    c_a1 = sval(c, "A1", "allele1")
    c_a2 = sval(c, "A2", "allele0")
    g_a1 = sval(g, "allele1", "A1")
    g_a0 = sval(g, "allele0", "A0", "A2")
    g_beta = fval(g, "beta")
    sign = 1.0
    if g_a1 == c_a1:
        sign = 1.0
    elif g_a0 == c_a1 or g_a1 == c_a2:
        sign = -1.0
    else:
        unaligned += 1
    c_beta = fval(c, "BETA", "beta")
    c_se = fval(c, "SE", "se")
    c_p = fval(c, "P", "p", "p_value")
    g_se = fval(g, "se")
    g_p = fval(g, "p_wald", "p_lrt", "p_score", "p")
    g_beta_a1 = sign * g_beta
    rows.append({
        "SNP": snp,
        "CHR": sval(c, "CHR", "chr"),
        "BP": sval(c, "BP", "ps"),
        "cosmic_A1": c_a1,
        "cosmic_A2": c_a2,
        "gemma_allele1": g_a1,
        "gemma_allele0": g_a0,
        "cosmic_beta": c_beta,
        "gemma_beta_a1": g_beta_a1,
        "beta_diff": c_beta - g_beta_a1,
        "abs_beta_diff": abs(c_beta - g_beta_a1),
        "cosmic_se": c_se,
        "gemma_se": g_se,
        "se_diff": c_se - g_se,
        "abs_se_diff": abs(c_se - g_se),
        "cosmic_p": c_p,
        "gemma_p": g_p,
        "log10p_cosmic": -math.log10(max(c_p, 1e-300)) if c_p > 0 else 300.0,
        "log10p_gemma": -math.log10(max(g_p, 1e-300)) if g_p > 0 else 300.0,
        "p_abs_diff": abs(c_p - g_p),
    })

fieldnames = list(rows[0].keys())
with open(report_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

beta_corr = corr([r["cosmic_beta"] for r in rows], [r["gemma_beta_a1"] for r in rows])
se_corr = corr([r["cosmic_se"] for r in rows], [r["gemma_se"] for r in rows])
logp_corr = corr([r["log10p_cosmic"] for r in rows], [r["log10p_gemma"] for r in rows])
mean_abs_beta = sum(r["abs_beta_diff"] for r in rows) / len(rows)
max_abs_beta = max(r["abs_beta_diff"] for r in rows)
mean_abs_se = sum(r["abs_se_diff"] for r in rows) / len(rows)
max_abs_se = max(r["abs_se_diff"] for r in rows)
mean_abs_p = sum(r["p_abs_diff"] for r in rows) / len(rows)
max_abs_p = max(r["p_abs_diff"] for r in rows)

thresholds = {
    "BETA_CORR_MIN": float(os.environ.get("BETA_CORR_MIN", "0.99")),
    "SE_CORR_MIN": float(os.environ.get("SE_CORR_MIN", "0.99")),
    "LOGP_CORR_MIN": float(os.environ.get("LOGP_CORR_MIN", "0.99")),
    "MEAN_ABS_BETA_MAX": float(os.environ.get("MEAN_ABS_BETA_MAX", "0.05")),
}
passed = (
    beta_corr >= thresholds["BETA_CORR_MIN"] and
    se_corr >= thresholds["SE_CORR_MIN"] and
    logp_corr >= thresholds["LOGP_CORR_MIN"] and
    mean_abs_beta <= thresholds["MEAN_ABS_BETA_MAX"]
)

lines = [
    f"overlap_snps\t{len(rows)}",
    f"unaligned_alleles\t{unaligned}",
    f"beta_corr\t{beta_corr:.12g}",
    f"se_corr\t{se_corr:.12g}",
    f"log10p_corr\t{logp_corr:.12g}",
    f"mean_abs_beta_diff\t{mean_abs_beta:.12g}",
    f"max_abs_beta_diff\t{max_abs_beta:.12g}",
    f"mean_abs_se_diff\t{mean_abs_se:.12g}",
    f"max_abs_se_diff\t{max_abs_se:.12g}",
    f"mean_abs_p_diff\t{mean_abs_p:.12g}",
    f"max_abs_p_diff\t{max_abs_p:.12g}",
    f"threshold_beta_corr_min\t{thresholds['BETA_CORR_MIN']}",
    f"threshold_se_corr_min\t{thresholds['SE_CORR_MIN']}",
    f"threshold_log10p_corr_min\t{thresholds['LOGP_CORR_MIN']}",
    f"threshold_mean_abs_beta_max\t{thresholds['MEAN_ABS_BETA_MAX']}",
    f"pass\t{str(passed).lower()}",
]
with open(summary_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")

print("\n".join(lines))
if not passed:
    raise SystemExit(1)
PY

echo "[GEMMA-accuracy] Comparison report: ${REPORT}"
echo "[GEMMA-accuracy] Summary: ${SUMMARY}"
