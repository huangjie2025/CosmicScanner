#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
PLINK2="${PLINK2_BIN:-/disk228/huangjie/software/miniconda3/envs/benchmark/bin/plink2}"
DATA_DIR="${ROOT}/tests/data/repeatability_format_consistency"
RESULT_DIR="${ROOT}/tests/results/repeatability_format_consistency"
DATA_PREFIX="${DATA_DIR}/repfmt"
PGEN_PREFIX="${RESULT_DIR}/pgen/repfmt"
BGEN_PREFIX="${RESULT_DIR}/bgen/repfmt"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner: ${BIN}" >&2
  exit 2
fi
if [[ ! -x "${PLINK2}" ]]; then
  echo "Missing PLINK2: ${PLINK2}" >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}/plink" "${RESULT_DIR}/pgen" "${RESULT_DIR}/bgen"

echo "[repeatability-format] Generating PLINK repeated-measure data..."
python3 - "${DATA_PREFIX}" <<'PY'
import math, random, sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260608)
n = 36
m = 144
samples = [(f"F{i+1:04d}", f"I{i+1:04d}") for i in range(n)]

genotypes = []
variants = []
for j in range(m):
    chrom = "1" if j < m // 2 else "2"
    snp = f"rrs{j+1:05d}"
    bp = 1000 + j * 17
    p = 0.18 + 0.55 * ((j * 37) % 101) / 100.0
    row = []
    for i in range(n):
        pop_shift = 0.10 if i >= n // 2 else -0.05
        pi = min(0.92, max(0.02, p + pop_shift))
        g = (1 if rng.random() < pi else 0) + (1 if rng.random() < pi else 0)
        row.append(g)
    genotypes.append(row)
    variants.append((chrom, snp, bp, "A", "G"))

causal = {9: 0.55, 83: -0.45}
subject_effect = [rng.gauss(0.0, 0.45) for _ in range(n)]

prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for fid, iid in samples:
        f.write(f"{fid}\t{iid}\t0\t0\t1\t-9\n")

with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for chrom, snp, bp, a1, a2 in variants:
        f.write(f"{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n")

with prefix.with_suffix(".pheno").open("w", encoding="utf-8") as f, \
     Path(str(prefix) + ".weights").open("w", encoding="utf-8") as wf:
    for i, (fid, iid) in enumerate(samples):
        reps = 2 + (i % 3)
        base = subject_effect[i] + (0.35 if i >= n // 2 else -0.15)
        for j, eff in causal.items():
            base += eff * (genotypes[j][i] - 1.0)
        for t in range(1, reps + 1):
            year = float(2018 + t)
            y = base + 0.20 * t + 0.10 * math.sin(i + t) + rng.gauss(0.0, 0.35 + 0.03 * t)
            w = 1.0 / (1.0 + 0.15 * t)
            f.write(f"{fid}\t{iid}\t{year:.1f}\t{y:.8f}\n")
            wf.write(f"{fid}\t{iid}\t{year:.1f}\t{w:.8f}\n")

def plink_code(dosage):
    if dosage == 2: return 0b00
    if dosage == 1: return 0b10
    if dosage == 0: return 0b11
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

echo "[repeatability-format] Converting to PGEN and BGEN..."
"${PLINK2}" --bfile "${DATA_PREFIX}" --make-pgen --out "${PGEN_PREFIX}" >/dev/null 2>&1
"${PLINK2}" --bfile "${DATA_PREFIX}" --export bgen-1.2 --out "${BGEN_PREFIX}" >/dev/null 2>&1

run_scan() {
  local label="$1"; shift
  local out_prefix="$1"; shift
  echo "[repeatability-format] Running ${label}..."
  "${BIN}" "$@" \
    --pheno "${DATA_PREFIX}.pheno" \
    --tpheno 1 \
    --mpheno 2 \
    --repeatability-model \
    --repeatability-mode exact \
    --reml-vg 0.35 \
    --reml-vpe 0.25 \
    --reml-ve 0.40 \
    --maf 0.0 \
    --gwas \
    --out "${out_prefix}" >/dev/null
}

run_group() {
  local suffix="$1"
  local extra="$2"
  local loco_flag="$3"
  local plink_out="${RESULT_DIR}/plink/repfmt_${suffix}"
  local pgen_out="${RESULT_DIR}/pgen/repfmt_${suffix}"
  local bgen_out="${RESULT_DIR}/bgen/repfmt_${suffix}"

  # shellcheck disable=SC2086
  run_scan "PLINK ${suffix}" "${plink_out}" --bfile "${DATA_PREFIX}" ${extra} ${loco_flag}
  # shellcheck disable=SC2086
  run_scan "PGEN ${suffix}" "${pgen_out}" --pfile "${PGEN_PREFIX}" ${extra} ${loco_flag}
  # shellcheck disable=SC2086
  run_scan "BGEN ${suffix}" "${bgen_out}" --bgen "${BGEN_PREFIX}.bgen" --sample "${BGEN_PREFIX}.sample" ${extra} ${loco_flag}

  python3 - "${plink_out}.assoc.txt" "${pgen_out}.assoc.txt" "${bgen_out}.assoc.txt" "${suffix}" <<'PY'
import csv, math, sys

ref_path, pgen_path, bgen_path, label = sys.argv[1:]

def read(path):
    with open(path, "r", encoding="utf-8") as f:
        return {r["SNP"]: r for r in csv.DictReader(f, delimiter="\t")}

def f(row, key):
    try:
        v = float(row[key])
        return v if math.isfinite(v) else math.nan
    except Exception:
        return math.nan

def corr(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 3:
        return math.nan
    xs, ys = zip(*pairs)
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx <= 0 or vy <= 0:
        return math.nan
    return sum((x - mx) * (y - my) for x, y in pairs) / math.sqrt(vx * vy)

def compare(test_path, name):
    ref = read(ref_path)
    test = read(test_path)
    snps = sorted(set(ref) & set(test))
    if len(snps) < 100:
        raise SystemExit(f"{label} {name}: too few overlapping SNPs: {len(snps)}")
    max_beta = max(abs(f(ref[s], "BETA") - f(test[s], "BETA")) for s in snps)
    max_se = max(abs(f(ref[s], "SE") - f(test[s], "SE")) for s in snps)
    max_logp = max(abs(f(ref[s], "LOG10P") - f(test[s], "LOG10P")) for s in snps)
    beta_corr = corr([f(ref[s], "BETA") for s in snps], [f(test[s], "BETA") for s in snps])
    se_corr = corr([f(ref[s], "SE") for s in snps], [f(test[s], "SE") for s in snps])
    logp_corr = corr([f(ref[s], "LOG10P") for s in snps], [f(test[s], "LOG10P") for s in snps])
    print(f"  {label} {name}: overlap={len(snps)} max_beta={max_beta:.3e} max_se={max_se:.3e} max_log10p={max_logp:.3e} beta_corr={beta_corr:.12f} se_corr={se_corr:.12f} logp_corr={logp_corr:.12f}")
    if max_beta > 2e-6 or max_se > 2e-6 or max_logp > 1e-5:
        raise SystemExit(f"{label} {name}: format consistency failed")

compare(pgen_path, "PGEN")
compare(bgen_path, "BGEN")
PY
}

run_group "global" "" "--no-loco"
run_group "weighted_global" "--repeatability-weight ${DATA_PREFIX}.weights" "--no-loco"
run_group "loco" "" ""

echo "[repeatability-format] PLINK/PGEN/BGEN repeatability consistency passed."
