#!/usr/bin/env bash
# Biobank full-flow calibration against REGENIE
# Tests Step1+Step2 end-to-end for PLINK, PGEN, BGEN paths
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
PLINK2="${PLINK2_BIN:-plink2}"
REGENIE="${REGENIE_BIN:-regenie}"
PLINK2="$(command -v "${PLINK2}" || printf '%s' "${PLINK2}")"
REGENIE="$(command -v "${REGENIE}" || printf '%s' "${REGENIE}")"
DATA_DIR="${ROOT}/tests/data/biobank_calibration"
RESULT_DIR="${ROOT}/tests/results/biobank_calibration"
DATA_PREFIX="${DATA_DIR}/biobank_calibration"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

for cmd_var in BIN PLINK2 REGENIE; do
    cmd_path="${!cmd_var}"
    if [[ ! -x "${cmd_path}" ]]; then
        echo "Missing ${cmd_var}: ${cmd_path}" >&2; exit 2
    fi
done

mkdir -p "${DATA_DIR}" "${RESULT_DIR}/cosmic_plink" "${RESULT_DIR}/cosmic_pgen" \
         "${RESULT_DIR}/cosmic_bgen" "${RESULT_DIR}/regenie"

echo "[biobank-calibration] Generating test data..."
python3 - "${DATA_PREFIX}" <<'PY'
import math, random, struct, sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260601)
n = 500
m = 2000

samples = []
for i in range(n):
    pop = 0 if i < n // 2 else 1
    samples.append((f"F{i+1:04d}", f"I{i+1:04d}", pop))

# Generate with realistic LD blocks
variants = []
genotypes = []
for j in range(m):
    chrom = str(1 + j // (m // 5))
    snp = f"rs{j+1:06d}"
    bp = 1000 + j * 10
    # Layer 1: population structure
    p = 0.08 + 0.42 * ((j * 37) % 100) / 99.0
    # Layer 2: haplotype blocks of 10 for LD
    block_id = j // 10
    haplo = [1 if rng.random() < 0.3 else 0 for _ in range(n)]
    row = []
    for i in range(n):
        _, _, pop = samples[i]
        p_adj = p + 0.10 * (pop - 0.5) * math.sin(j * 0.17)
        p_adj = min(0.95, max(0.02, p_adj))
        # Mix population-based with haplotype-based
        use_haplo = 1 if rng.random() < 0.2 else 0
        if use_haplo:
            g = haplo[i] + (1 if rng.random() < 0.3 else 0)
        else:
            g = (1 if rng.random() < p_adj else 0) + (1 if rng.random() < p_adj else 0)
        row.append(min(2, g))
    variants.append((chrom, snp, bp, "A", "C"))
    genotypes.append(row)

# Causal variants: 20 across the genome
causal = {}
for idx in [50, 150, 320, 480, 610, 790, 920, 1050, 1180, 1310,
            1420, 1530, 1640, 1720, 1810, 1890, 5, 205, 405, 705]:
    causal[idx] = rng.gauss(0, 0.3)

phenotypes = []
for i in range(n):
    _, _, pop = samples[i]
    y = 0.5 * pop + rng.gauss(0, 0.5)
    for idx, eff in causal.items():
        y += eff * (genotypes[idx][i] - 1.0)
    y += rng.gauss(0.0, 1.0)
    phenotypes.append(y)

prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix(".fam").open("w") as f:
    for (fid, iid, _), y in zip(samples, phenotypes):
        f.write(f"{fid}\t{iid}\t0\t0\t0\t{y:.8f}\n")
with prefix.with_suffix(".pheno").open("w") as f:
    f.write("FID\tIID\tPHENO\n")
    for (fid, iid, _), y in zip(samples, phenotypes):
        f.write(f"{fid}\t{iid}\t{y:.8f}\n")
with prefix.with_suffix(".bim").open("w") as f:
    for chrom, snp, bp, a1, a2 in variants:
        f.write(f"{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n")

def plink_code(d):
    return 0b00 if d == 2 else 0b10 if d == 1 else 0b11 if d == 0 else 0b01
bed = bytearray(b"\x6c\x1b\x01")
for row in genotypes:
    block = bytearray((n + 3) // 4)
    for i, g in enumerate(row):
        block[i // 4] |= plink_code(g) << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix(".bed").write_bytes(bed)

# Covariate file for REGENIE (sex from FAM col 5)
with (prefix.parent / "covar.txt").open("w") as f:
    f.write("FID\tIID\tCOVAR1\n")
    for (fid, iid, _) in samples:
        f.write(f"{fid}\t{iid}\t0\n")

# Model SNPs (first 200) and test SNPs (remaining 1800)
with (prefix.parent / "model_snps.txt").open("w") as f:
    f.write("SNP\n")
    for j in range(200):
        f.write(f"rs{j+1:06d}\n")

# Truth
with (prefix.parent / "causal_truth.txt").open("w") as f:
    f.write("SNP\tTRUE_BETA\n")
    for idx, eff in sorted(causal.items()):
        f.write(f"rs{idx+1:06d}\t{eff}\n")
PY

MODEL_SNPS="${DATA_DIR}/model_snps.txt"

echo "[biobank-calibration] Converting to PGEN + BGEN..."
"${PLINK2}" --bfile "${DATA_PREFIX}" --make-pgen --out "${RESULT_DIR}/cosmic_pgen/biobank" 2>&1 | tail -2
"${PLINK2}" --bfile "${DATA_PREFIX}" --export bgen-1.2 --out "${RESULT_DIR}/cosmic_bgen/biobank" 2>&1 | tail -2

echo ""
echo "[biobank-calibration] ===== REGENIE Step 1 ====="
# REGENIE: no covariates, just phenotype
"${REGENIE}" \
    --step 1 \
    --bed "${DATA_PREFIX}" \
    --phenoFile "${DATA_PREFIX}.pheno" \
    --covarFile "${DATA_DIR}/covar.txt" \
    --bsize 200 \
    --extract "${MODEL_SNPS}" \
    --pred "${RESULT_DIR}/regenie/step1_pred" \
    --out "${RESULT_DIR}/regenie/step1" 2>&1 | tail -10

echo ""
echo "[biobank-calibration] ===== REGENIE Step 2 ====="
"${REGENIE}" \
    --step 2 \
    --bed "${DATA_PREFIX}" \
    --phenoFile "${DATA_PREFIX}.pheno" \
    --covarFile "${DATA_DIR}/covar.txt" \
    --bsize 200 \
    --pred "${RESULT_DIR}/regenie/step1_pred.list" \
    --out "${RESULT_DIR}/regenie/step2" 2>&1 | tail -10

echo ""
echo "[biobank-calibration] ===== CosmicScanner PLINK biobank mode ====="
"${BIN}" \
    --bfile "${DATA_PREFIX}" \
    --pheno "${DATA_PREFIX}.pheno" \
    --mpheno 1 \
    --gwas \
    --biobank-mode \
    --block-size 200 \
    --regenie-model-snps "${MODEL_SNPS}" \
    --maf 0.0 \
    --out "${RESULT_DIR}/cosmic_plink/biobank" 2>&1 | tail -10

echo ""
echo "[biobank-calibration] ===== CosmicScanner PGEN biobank mode ====="
"${BIN}" \
    --pfile "${RESULT_DIR}/cosmic_pgen/biobank" \
    --pheno "${DATA_PREFIX}.pheno" \
    --mpheno 1 \
    --gwas \
    --biobank-mode \
    --block-size 200 \
    --regenie-model-snps "${MODEL_SNPS}" \
    --maf 0.0 \
    --out "${RESULT_DIR}/cosmic_pgen/biobank" 2>&1 | tail -10

echo ""
echo "[biobank-calibration] ===== CosmicScanner BGEN biobank mode ====="
"${BIN}" \
    --bgen "${RESULT_DIR}/cosmic_bgen/biobank.bgen" \
    --sample "${RESULT_DIR}/cosmic_bgen/biobank.sample" \
    --pheno "${DATA_PREFIX}.pheno" \
    --mpheno 1 \
    --gwas \
    --biobank-mode \
    --block-size 200 \
    --regenie-model-snps "${MODEL_SNPS}" \
    --maf 0.0 \
    --out "${RESULT_DIR}/cosmic_bgen/biobank" 2>&1 | tail -10

echo ""
echo "[biobank-calibration] ===== Cross-Comparison ====="

python3 - "${RESULT_DIR}" <<'PY'
import csv, glob, math, os, sys

res_dir = sys.argv[1]

# Find REGENIE output (appends trait name like _PHENO)
regenie_files = sorted(glob.glob(f"{res_dir}/regenie/step2_*.regenie"))
regenie_path = regenie_files[0] if regenie_files else f"{res_dir}/regenie/step2_PHENO.regenie"

results = {
    "REGENIE": regenie_path,
    "Cosmic_PLINK": f"{res_dir}/cosmic_plink/biobank.biobank.assoc.txt",
    "Cosmic_PGEN": f"{res_dir}/cosmic_pgen/biobank.biobank.assoc.txt",
    "Cosmic_BGEN": f"{res_dir}/cosmic_bgen/biobank.biobank.assoc.txt",
}

def read_table(path):
    rows = {}
    # Auto-detect delimiter: check if header has tabs
    with open(path) as f:
        first_line = f.readline().rstrip('\n')
    delim = '\t' if '\t' in first_line else None  # None = csv.Sniffer or split
    if delim is None:
        # space-separated (REGENIE)
        delim = None  # let split handle
    with open(path) as f:
        if delim:
            reader = csv.DictReader(f, delimiter=delim)
        else:
            # For space-separated: read header manually
            header = f.readline().rstrip('\n').split()
            for line in f:
                parts = line.rstrip('\n').split()
                if len(parts) != len(header): continue
                row = dict(zip(header, parts))
                snp = row.get("ID", "")
                beta = float(row.get("BETA", "nan") or "nan")
                se = float(row.get("SE", "nan") or "nan")
                log10p = float(row.get("LOG10P", "nan") or "nan")
                p = 10.0 ** -log10p if math.isfinite(log10p) else float("nan")
                rows[snp] = {"BETA": beta, "SE": se, "P": p}
    if not rows:
        # Fallback: try tab-separated
        with open(path) as f:
            if delim:
                reader = csv.DictReader(f, delimiter='\t')
            else:
                reader = csv.DictReader(f, delimiter=' ')
            for r in reader:
                snp = r.get("SNP", r.get("ID", ""))
                rows[snp] = {
                    "BETA": float(r.get("BETA", "nan") or "nan"),
                    "SE": float(r.get("SE", "nan") or "nan"),
                    "P": float(r.get("P", "nan") or "nan"),
                }
    return rows

def corr(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 3: return math.nan
    mx = sum(x for x,_ in pairs)/len(pairs)
    my = sum(y for _,y in pairs)/len(pairs)
    vx = sum((x-mx)**2 for x,_ in pairs)
    vy = sum((y-my)**2 for _,y in pairs)
    if vx <= 0 or vy <= 0: return math.nan
    cov = sum((x-mx)*(y-my) for x,y in pairs)
    return cov / math.sqrt(vx*vy)

tables = {}
for name, path in results.items():
    if os.path.exists(path):
        tables[name] = read_table(path)
        print(f"  {name}: {len(tables[name])} SNPs loaded from {path}")
    else:
        print(f"  {name}: NOT FOUND at {path}")

if "REGENIE" not in tables:
    print("ERROR: REGENIE results missing")
    exit(1)

ref = tables["REGENIE"]
for name in ["Cosmic_PLINK", "Cosmic_PGEN", "Cosmic_BGEN"]:
    if name not in tables: continue
    test = tables[name]
    snps = sorted(set(ref) & set(test))
    if not snps: continue

    betas_r = [ref[s]["BETA"] for s in snps]
    betas_t = [test[s]["BETA"] for s in snps]
    ses_r = [ref[s]["SE"] for s in snps]
    ses_t = [test[s]["SE"] for s in snps]
    logps_r = [-math.log10(max(ref[s]["P"], 1e-300)) for s in snps]
    logps_t = [-math.log10(max(test[s]["P"], 1e-300)) for s in snps]

    bc = corr(betas_r, betas_t)
    sc = corr(ses_r, ses_t)
    lc = corr(logps_r, logps_t)
    mab = sum(abs(betas_r[i] - betas_t[i]) for i in range(len(snps))) / len(snps)
    mas = sum(abs(ses_r[i] - ses_t[i]) for i in range(len(snps))) / len(snps)
    map_ = sum(abs(ref[s]["P"] - test[s]["P"]) for s in snps) / len(snps)
    maxab = max(abs(betas_r[i] - betas_t[i]) for i in range(len(snps)))
    maxas = max(abs(ses_r[i] - ses_t[i]) for i in range(len(snps)))
    maxap = max(abs(ref[s]["P"] - test[s]["P"]) for s in snps)

    print(f"\n  {name} vs REGENIE ({len(snps)} SNPs):")
    print(f"    beta_corr:      {bc:.6f}")
    print(f"    se_corr:        {sc:.6f}")
    print(f"    log10p_corr:    {lc:.6f}")
    print(f"    mean_abs_beta:  {mab:.6f}")
    print(f"    max_abs_beta:   {maxab:.6f}")
    print(f"    mean_abs_se:    {mas:.6f}")
    print(f"    max_abs_se:     {maxas:.6f}")
    print(f"    mean_abs_p:     {map_:.6e}")
    print(f"    max_abs_p:      {maxap:.6e}")
PY

echo ""
echo "[biobank-calibration] Done. Results in ${RESULT_DIR}/"
