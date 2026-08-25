#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/psr_gwas_smoke"
RESULT_DIR="${ROOT}/tests/results/psr_gwas_smoke"
DATA_PREFIX="${DATA_DIR}/psr_smoke"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  echo "Build it first or set COSMICSCANNER_BIN=/path/to/cosmicscanner." >&2
  exit 2
fi

resolve_rscript() {
  if [[ -n "${RSCRIPT_BIN:-}" ]]; then
    printf '%s\n' "${RSCRIPT_BIN}"
    return 0
  fi
  if command -v Rscript >/dev/null 2>&1; then
    command -v Rscript
    return 0
  fi
  local conda_r="/disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript"
  if [[ -x "${conda_r}" ]]; then
    printf '%s\n' "${conda_r}"
    return 0
  fi
  return 1
}

mkdir -p "${DATA_DIR}" "${RESULT_DIR}"

python3 - "${DATA_PREFIX}" <<'PY'
import math
import random
import sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260622)
n = 96
m = 240

samples = []
for i in range(n):
    fid = f"F{i+1:04d}"
    iid = f"I{i+1:04d}"
    pop = 0 if i < n // 2 else 1
    samples.append((fid, iid, pop))

variants = []
genotypes = []
for j in range(m):
    chrom = str((j % 3) + 1)
    snp = f"psr{j+1:06d}"
    bp = 1000 + j * 10
    p = 0.08 + 0.42 * ((j * 17) % 100) / 99.0
    row = []
    for _, _, pop in samples:
        pp = min(0.95, max(0.02, p + (0.06 if pop else -0.03) * math.sin(j * 0.11)))
        g = (1 if rng.random() < pp else 0) + (1 if rng.random() < pp else 0)
        row.append(g)
    variants.append((chrom, snp, bp, "A", "C"))
    genotypes.append(row)

binary_causal = {18: 0.95, 119: -0.75}
poisson_causal = {42: 0.35, 177: -0.30}
binary = []
counts = []
offsets = []
for i, (_, _, pop) in enumerate(samples):
    eta = -0.35 + 0.45 * pop
    for snp_idx, eff in binary_causal.items():
        eta += eff * (genotypes[snp_idx][i] - 1.0)
    pr = 1.0 / (1.0 + math.exp(-eta))
    binary.append(1 if rng.random() < pr else 0)

    exposure = 0.8 + 0.4 * ((i * 7) % 11) / 10.0
    offsets.append(math.log(exposure))
    log_mu = 0.7 + math.log(exposure) + 0.15 * pop
    for snp_idx, eff in poisson_causal.items():
        log_mu += eff * (genotypes[snp_idx][i] - 1.0)
    mu = math.exp(log_mu)
    # Knuth Poisson sampler; mu is small in this smoke data.
    limit = math.exp(-mu)
    k = 0
    prod = 1.0
    while prod > limit:
        k += 1
        prod *= rng.random()
    counts.append(k - 1)

prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for fid, iid, _ in samples:
        f.write(f"{fid}\t{iid}\t0\t0\t1\t-9\n")
with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for chrom, snp, bp, a1, a2 in variants:
        f.write(f"{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n")
with prefix.with_suffix(".binary.pheno").open("w", encoding="utf-8") as f:
    for (fid, iid, _), y in zip(samples, binary):
        f.write(f"{fid}\t{iid}\t{y}\n")
with prefix.with_suffix(".poisson.pheno").open("w", encoding="utf-8") as f:
    for (fid, iid, _), y in zip(samples, counts):
        f.write(f"{fid}\t{iid}\t{y}\n")
with prefix.with_suffix(".binomial.pheno").open("w", encoding="utf-8") as f:
    for (fid, iid, _), y in zip(samples, binary):
        # Deterministic two-column event/trial representation of the binary data.
        f.write(f"{fid}\t{iid}\t{y}\t1\n")
with prefix.with_suffix(".offset").open("w", encoding="utf-8") as f:
    for x in offsets:
        f.write(f"{x:.12g}\n")
with prefix.with_suffix(".geno.tsv").open("w", encoding="utf-8") as f:
    f.write("FID\tIID" + "".join(f"\t{snp}" for _, snp, _, _, _ in variants) + "\n")
    for i, (fid, iid, _) in enumerate(samples):
        f.write(f"{fid}\t{iid}" + "".join(f"\t{genotypes[j][i]}" for j in range(m)) + "\n")

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

echo "[PSR-GWAS] Running CosmicScanner binary/probit PSR smoke..."
/usr/bin/time -f "cosmic_binary_seconds\t%e\ncosmic_binary_maxrss_kb\t%M" \
  -o "${RESULT_DIR}/timing.binary.tsv" \
  "${BIN}" \
    --bfile "${DATA_PREFIX}" \
    --pheno "${DATA_PREFIX}.binary.pheno" \
    --mpheno 1 \
    --gwas \
    --threshold-model \
    --link probit \
    --no-loco \
    --maf 0.0 \
    --out "${RESULT_DIR}/cosmic_binary"

echo "[PSR-GWAS] Running CosmicScanner Poisson PSR smoke..."
/usr/bin/time -f "cosmic_poisson_seconds\t%e\ncosmic_poisson_maxrss_kb\t%M" \
  -o "${RESULT_DIR}/timing.poisson.tsv" \
  "${BIN}" \
    --bfile "${DATA_PREFIX}" \
    --pheno "${DATA_PREFIX}.poisson.pheno" \
    --mpheno 1 \
    --gwas \
    --count-model \
    --family poisson \
    --offset "${DATA_PREFIX}.offset" \
    --no-loco \
    --maf 0.0 \
    --out "${RESULT_DIR}/cosmic_poisson"

echo "[PSR-GWAS] Running CosmicScanner binomial event/trial PSR smoke..."
/usr/bin/time -f "cosmic_binomial_seconds\t%e\ncosmic_binomial_maxrss_kb\t%M" \
  -o "${RESULT_DIR}/timing.binomial.tsv" \
  "${BIN}" \
    --bfile "${DATA_PREFIX}" \
    --gwas \
    --count-model \
    --family binomial \
    --binomial-pheno "${DATA_PREFIX}.binomial.pheno" \
    --no-loco \
    --maf 0.0 \
    --out "${RESULT_DIR}/cosmic_binomial"

python3 - "${RESULT_DIR}" <<'PY'
import csv
import math
import sys
from pathlib import Path

result_dir = Path(sys.argv[1])
summary = result_dir / "psr_gwas_smoke.summary.tsv"
rows = []
for label in ("binary", "poisson", "binomial"):
    assoc = result_dir / f"cosmic_{label}.assoc.txt"
    with assoc.open("r", encoding="utf-8") as f:
        data = list(csv.DictReader(f, delimiter="\t"))
    pvals = []
    for row in data:
        try:
            p = float(row.get("P", row.get("P_VALUE", "nan")))
        except ValueError:
            p = math.nan
        if math.isfinite(p):
            pvals.append(p)
    if not pvals:
        raise SystemExit(f"No finite P values in {assoc}")
    rows.append((label, len(data), min(pvals), sum(1 for p in pvals if p < 0.05)))
with summary.open("w", encoding="utf-8") as f:
    f.write("trait\tassoc_rows\tmin_p\tn_p_lt_0.05\n")
    for row in rows:
        f.write("%s\t%d\t%.12g\t%d\n" % row)
print(f"[PSR-GWAS] Summary written to {summary}")
PY

if [[ -n "${PSR_GLMM_R_SCRIPT:-}" ]]; then
  if [[ ! -f "${PSR_GLMM_R_SCRIPT}" ]]; then
    echo "PSR_GLMM_R_SCRIPT is set but not found: ${PSR_GLMM_R_SCRIPT}" >&2
    exit 3
  fi
  echo "[PSR-GWAS] Running external PSR-GLMM/R benchmark hook..."
  if ! RSCRIPT_EXE="$(resolve_rscript)"; then
    echo "Rscript not found; cannot run PSR_GLMM_R_SCRIPT." >&2
    exit 4
  fi
  "${RSCRIPT_EXE}" "${PSR_GLMM_R_SCRIPT}" "${DATA_PREFIX}" "${RESULT_DIR}"
else
  echo "[PSR-GWAS] PSR_GLMM_R_SCRIPT not set; skipped R reference comparison hook."
fi

if [[ -z "${PSR_GLMM_R_SCRIPT:-}" ]] && RSCRIPT_EXE="$(resolve_rscript)"; then
  echo "[PSR-GWAS] Running local R GLM baseline adapter..."
  "${RSCRIPT_EXE}" "${ROOT}/tests/psr_glmm_r_adapter.R" "${DATA_PREFIX}" "${RESULT_DIR}"
fi

if RSCRIPT_EXE="$(resolve_rscript)"; then
  echo "[PSR-GWAS] Running candidate SNP GLM refinement smoke..."
  "${RSCRIPT_EXE}" "${ROOT}/tests/psr_candidate_refine.R" \
    "${DATA_PREFIX}" \
    "${RESULT_DIR}/cosmic_poisson.assoc.txt" \
    "poisson" \
    "${RESULT_DIR}/cosmic_poisson.refined.glm.tsv" \
    20 \
    glm
else
  echo "[PSR-GWAS] Rscript not found; skipped candidate refinement smoke."
fi

echo "[PSR-GWAS] Running candidate refinement suite..."
bash "${ROOT}/tests/run_psr_candidate_refinement_suite.sh"

python3 "${ROOT}/tests/test_weighted_lmm_equivalence.py"
python3 "${ROOT}/tests/run_psr_simulation_summary.py"
python3 "${ROOT}/tests/summarize_psr_gwas_benchmark.py" "${RESULT_DIR}"

echo "[PSR-GWAS] Running input validation regression..."
bash "${ROOT}/tests/test_psr_input_validation.sh"
