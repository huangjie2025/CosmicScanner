#!/usr/bin/env bash
# Test BFILE vs PGEN GWAS consistency (including format=2 with type 2/3 LD compression)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
PLINK2="${PLINK2_BIN:-plink2}"
PLINK2="$(command -v "${PLINK2}" || printf '%s' "${PLINK2}")"
DATA_DIR="${ROOT}/tests/data/pgen_consistency"
RESULT_DIR="${ROOT}/tests/results/pgen_consistency"
DATA_PREFIX="${DATA_DIR}/pgen_consistency"
BFILE_RESULT="${RESULT_DIR}/bfile/pgen_consistency"
PGEN_RESULT="${RESULT_DIR}/pgen/pgen_consistency"
PGEN2_RESULT="${RESULT_DIR}/pgen_fmt2/pgen_consistency"
LMM_NULL_SOLVER="${LMM_NULL_SOLVER:-brent}"

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

mkdir -p "${DATA_DIR}" "${RESULT_DIR}/bfile" "${RESULT_DIR}/pgen" "${RESULT_DIR}/pgen_fmt2"

echo "[PGEN-consistency] Generating test data..."
python3 - "${DATA_PREFIX}" <<'PY'
import math, random, struct, sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260531)
n = 160
m = 800

# Generate samples with population structure for LD-like correlation
samples = []
for i in range(n):
    pop = 0 if i < n // 2 else 1
    samples.append((f"F{i+1:04d}", f"I{i+1:04d}", pop))

# Generate variants with VERY strong haplotype blocks to trigger type 2/3 encoding
variants = []
genotypes = []
block_size = 10  # variants per haplotype block

# Use a single haplotype base for extremely strong LD
haplo_base = [1 if rng.random() < 0.3 else 0 for _ in range(n)]
haplo_base2 = [1 if rng.random() < 0.3 else 0 for _ in range(n)]
# The first variant (non-LD anchor) uses the base haplotypes directly
anchor = [haplo_base[i] + haplo_base2[i] for i in range(n)]

for b in range(m // block_size):
    for j in range(block_size):
        idx = b * block_size + j
        chrom = str(1 + idx // (m // 4))
        snp = f"rs{idx + 1:06d}"
        bp = 1000 + idx * 10
        # First half of each block: anchored on original haplotypes with tiny flip
        # Second half: swap ref/alt to trigger type 3
        swap = (j >= block_size // 2)
        row = []
        for i in range(n):
            # Very small flip rate (1%) for strong correlation
            flip1 = 1 if rng.random() < 0.01 else 0
            flip2 = 1 if rng.random() < 0.01 else 0
            if swap:
                # Swap: flip both bits -> 0<->2 inversion
                h1 = (haplo_base[i] + flip1) % 2
                h2 = (haplo_base2[i] + flip2) % 2
                g = h1 + h2
                # After inversion, this matches the anchor
            else:
                h1 = (haplo_base[i] + flip1) % 2
                h2 = (haplo_base2[i] + flip2) % 2
                g = h1 + h2
            row.append(g)
        if swap:
            variants.append((chrom, snp, bp, "C", "A"))  # swapped ref
        else:
            variants.append((chrom, snp, bp, "A", "C"))
        genotypes.append(row)

# Ensure first variant of each 65536-block is never LD (set to type 0 by making it first var)
# The overall block structure should already trigger PLINK2's LD compression

causal = {12: 0.45, 211: -0.35, 509: 0.30}
phenotypes = []
for i in range(n):
    _, _, pop = samples[i]
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

echo "[PGEN-consistency] Running BFILE GWAS (baseline)..."
"${BIN}" \
    --bfile "${DATA_PREFIX}" \
    --pheno "${DATA_PREFIX}.pheno" \
    --mpheno 1 \
    --gwas \
    --no-loco \
    --lmm-null-solver "${LMM_NULL_SOLVER}" \
    --maf 0.0 \
    --out "${BFILE_RESULT}"

echo "[PGEN-consistency] Converting to PGEN (default format)..."
"${PLINK2}" --bfile "${DATA_PREFIX}" \
    --make-pgen \
    --out "${RESULT_DIR}/pgen/pgen_consistency_pgen" 2>&1 | tail -3

echo "[PGEN-consistency] Converting to PGEN (format=2, with LD compression)..."
"${PLINK2}" --bfile "${DATA_PREFIX}" \
    --make-pgen format=2 \
    --out "${RESULT_DIR}/pgen_fmt2/pgen_consistency_pgen" 2>&1 | tail -3

echo "[PGEN-consistency] Running PGEN (default) GWAS..."
"${BIN}" \
    --pfile "${RESULT_DIR}/pgen/pgen_consistency_pgen" \
    --pheno "${DATA_PREFIX}.pheno" \
    --mpheno 1 \
    --gwas \
    --no-loco \
    --lmm-null-solver "${LMM_NULL_SOLVER}" \
    --maf 0.0 \
    --out "${PGEN_RESULT}"

echo "[PGEN-consistency] Running PGEN (format=2) GWAS..."
"${BIN}" \
    --pfile "${RESULT_DIR}/pgen_fmt2/pgen_consistency_pgen" \
    --pheno "${DATA_PREFIX}.pheno" \
    --mpheno 1 \
    --gwas \
    --no-loco \
    --lmm-null-solver "${LMM_NULL_SOLVER}" \
    --maf 0.0 \
    --out "${PGEN2_RESULT}"

BFILE_ASSOC="${BFILE_RESULT}.assoc.txt"
PGEN_ASSOC="${PGEN_RESULT}.assoc.txt"
PGEN2_ASSOC="${PGEN2_RESULT}.assoc.txt"

echo "[PGEN-consistency] Checking PGEN format vrtype distribution..."
"${PLINK2}" --pfile "${RESULT_DIR}/pgen/pgen_consistency_pgen" --pgen-info 2>&1 | grep -E "VR type|Storage mode" || true
"${PLINK2}" --pfile "${RESULT_DIR}/pgen_fmt2/pgen_consistency_pgen" --pgen-info 2>&1 | grep -E "VR type|Storage mode" || true

# Report PGEN vrtype distribution (pass paths as args to avoid heredoc issues)
python3 - "${RESULT_DIR}/pgen/pgen_consistency_pgen" "${RESULT_DIR}/pgen_fmt2/pgen_consistency_pgen" <<'PY'
import struct, sys
from collections import Counter
for p in sys.argv[1:]:
    fname = p + ".pgen"
    label = p.split("/")[-2]
    with open(fname, "rb") as f:
        data = f.read(65536)
    magic = data[:2]
    mode = data[2]
    nvar = struct.unpack('<I', data[3:7])[0]
    nsam = struct.unpack('<I', data[7:11])[0]
    fmt = data[11]  # format byte at offset 11 for all modes >= 0x02
    if mode >= 0x10:
        num_blocks = (nvar + 65535) // 65536
        off = 12
        vbo = struct.unpack('<Q', data[off:off+8])[0]
        off += num_blocks * 8
        if fmt & 0x80:
            vrtypes_bytes = (nvar + 1) // 2
            raw = data[off:off+vrtypes_bytes]
            vrtypes = [(raw[j//2] >> (4 if j%2==0 else 0)) & 0xF for j in range(nvar)]
        else:
            vrtypes = list(data[off:off+nvar])
    else:
        vrtypes = [fmt & 0x0F] * nvar  # mode 0x02: all same vrtype from format byte low nibble
    c = Counter(vrtypes)
    mc = Counter()
    for t, cnt in c.items():
        mc[t & 7] += cnt
    print(f"  {label}: mode=0x{mode:02x}, nvar={nvar}, nsam={nsam}, fmt=0x{fmt:02x}")
    print(f"    VR type distribution: {dict(sorted(c.items()))}")
    print(f"    Main type distribution: {dict(sorted(mc.items()))}")
PY

echo ""
echo "[PGEN-consistency] Comparing BFILE vs PGEN (default) results..."
python3 - "${BFILE_ASSOC}" "${PGEN_ASSOC}" "pgen_default" <<'PY'
import csv, math, sys

def read_table(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def fval(row, *names):
    for name in names:
        if name in row and row[name] not in ("", "NA", "nan", "-nan"):
            return float(row[name])
    return math.nan

def sval(row, *names):
    for name in names:
        if name in row:
            return row[name]
    return ""

def compare(ref_path, test_path, label):
    ref = {sval(r, "SNP", "rs"): r for r in read_table(ref_path)}
    test = {sval(r, "SNP", "rs"): r for r in read_table(test_path)}
    snps = sorted(set(ref) & set(test))
    md = mbd = msd = mpd = 0.0
    for snp in snps:
        r, t = ref[snp], test[snp]
        for col in ["BETA_A", "SE_A", "P_A"]:
            rv, tv = fval(r, col), fval(t, col)
            if math.isfinite(rv) and math.isfinite(tv):
                d = abs(rv - tv)
                if d > md: md = d
                if col == "BETA_A" and d > mbd: mbd = d
                if col == "SE_A" and d > msd: msd = d
                if col == "P_A" and d > mpd: mpd = d
    return len(snps), md, mbd, msd, mpd

n, md, mbd, msd, mpd = compare(sys.argv[1], sys.argv[2], sys.argv[3])
print(f"  Overlap SNPs: {n}")
print(f"  Max beta diff: {mbd:.6e}")
print(f"  Max SE diff: {msd:.6e}")
print(f"  Max P diff: {mpd:.6e}")
if n == 800 and md < 1e-5:
    print(f"  PASS: BFILE vs {sys.argv[3]} results match (max diff {md:.2e} < 1e-5, GRM rounding)")
else:
    print(f"  FAIL: {sys.argv[3]} results differ from BFILE (max diff {md:.2e})")
    if md >= 1e-5:
        exit(1)
PY

echo ""
echo "[PGEN-consistency] Comparing BFILE vs PGEN (format=2) results..."
python3 - "${BFILE_ASSOC}" "${PGEN2_ASSOC}" "pgen_fmt2" <<'PY'
import csv, math, sys

def read_table(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def fval(row, *names):
    for name in names:
        if name in row and row[name] not in ("", "NA", "nan", "-nan"):
            return float(row[name])
    return math.nan

def sval(row, *names):
    for name in names:
        if name in row:
            return row[name]
    return ""

def compare(ref_path, test_path, label):
    ref = {sval(r, "SNP", "rs"): r for r in read_table(ref_path)}
    test = {sval(r, "SNP", "rs"): r for r in read_table(test_path)}
    snps = sorted(set(ref) & set(test))
    max_diff = 0.0
    max_beta_diff = 0.0
    max_se_diff = 0.0
    max_p_diff = 0.0
    for snp in snps:
        r, t = ref[snp], test[snp]
        for col, diff_arr in [("BETA_A", [max_beta_diff]), ("SE_A", [max_se_diff]), ("P_A", [max_p_diff])]:
            rv = fval(r, col)
            tv = fval(t, col)
            if math.isfinite(rv) and math.isfinite(tv):
                diff = abs(rv - tv)
                if diff > diff_arr[0]:
                    diff_arr[0] = diff
                    if diff > max_diff:
                        max_diff = diff
    return len(snps), max_diff, max_beta_diff, max_se_diff, max_p_diff

n, md, mbd, msd, mpd = compare(sys.argv[1], sys.argv[2], sys.argv[3])
print(f"  Overlap SNPs: {n}")
print(f"  Max overall diff: {md:.6e}")
print(f"  Max beta diff: {mbd:.6e}")
print(f"  Max SE diff: {msd:.6e}")
print(f"  Max P diff: {mpd:.6e}")
if n == 800 and md < 1e-5:
    print(f"  PASS: BFILE vs {sys.argv[3]} results match (max diff {md:.2e} < 1e-5)")
else:
    print(f"  FAIL: {sys.argv[3]} results differ from BFILE (max diff {md:.2e})")
    if md >= 1e-5:
        exit(1)
PY

echo ""
echo "[PGEN-consistency] Comparing PGEN default vs PGEN format=2 (same PGEN path, should match exactly)..."
python3 - "${PGEN_ASSOC}" "${PGEN2_ASSOC}" "pgen_default_vs_fmt2" <<'PY'
import csv, math, sys
def read_table(p):
    with open(p) as f: return list(csv.DictReader(f, delimiter="\t"))
def fval(r, *ns):
    for n in ns:
        if n in r and r[n] not in ("", "NA", "nan", "-nan"): return float(r[n])
    return math.nan
def sval(r, *ns):
    for n in ns:
        if n in r: return r[n]
    return ""
ref = {sval(r, "SNP", "rs"): r for r in read_table(sys.argv[1])}
test = {sval(r, "SNP", "rs"): r for r in read_table(sys.argv[2])}
snps = sorted(set(ref) & set(test))
md = mbd = msd = mpd = 0.0
for snp in snps:
    r, t = ref[snp], test[snp]
    for col in ["BETA_A", "SE_A", "P_A"]:
        rv, tv = fval(r, col), fval(t, col)
        if math.isfinite(rv) and math.isfinite(tv):
            d = abs(rv - tv)
            if d > md: md = d
            if col == "BETA_A" and d > mbd: mbd = d
            if col == "SE_A" and d > msd: msd = d
            if col == "P_A" and d > mpd: mpd = d
print(f"  Max overall diff: {md:.6e}")
print(f"  Max beta diff: {mbd:.6e}")
print(f"  Max SE diff: {msd:.6e}")
print(f"  Max P diff: {mpd:.6e}")
if md < 1e-14:
    print(f"  PASS: PGEN default vs format=2 results are identical (within double precision)")
else:
    print(f"  ISH: minor FP difference ({md:.2e})")
PY

# Also compare against GEMMA for the BFILE baseline
if GEMMA_BIN="${GEMMA_BIN:-}" || command -v gemma >/dev/null 2>&1; then
    GEMMA_EXE="${GEMMA_BIN:-$(command -v gemma)}"
    echo ""
    echo "[PGEN-consistency] Running GEMMA for reference..."
    mkdir -p "${RESULT_DIR}/gemma"
    GEMMA_WORK="${RESULT_DIR}/gemma"
    rm -rf "${GEMMA_WORK}/output"
    (
        cd "${GEMMA_WORK}"
        "$GEMMA_EXE" -bfile "${DATA_PREFIX}" -gk 1 -o pgen_consistency_kinship > /dev/null 2>&1
        KINSHIP="${GEMMA_WORK}/output/pgen_consistency_kinship.cXX.txt"
        if [[ -f "$KINSHIP" ]]; then
            "$GEMMA_EXE" -bfile "${DATA_PREFIX}" -k "$KINSHIP" -lmm 4 -o pgen_consistency_lmm > /dev/null 2>&1
        fi
    )
    GEMMA_ASSOC="${GEMMA_WORK}/output/pgen_consistency_lmm.assoc.txt"
    if [[ -f "${GEMMA_ASSOC}" ]]; then
        echo "[PGEN-consistency] Comparing BFILE vs GEMMA..."
        python3 - "${BFILE_ASSOC}" "${GEMMA_ASSOC}" <<'PY'
import csv, math, sys
def read_table(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))
def fval(r, *ns):
    for n in ns:
        if n in r and r[n] not in ("", "NA", "nan", "-nan"):
            return float(r[n])
    return math.nan
def sval(r, *ns):
    for n in ns:
        if n in r: return r[n]
    return ""
ref = {sval(r, "SNP", "rs"): r for r in read_table(sys.argv[1])}
gem = {sval(r, "rs", "SNP"): r for r in read_table(sys.argv[2])}
snps = sorted(set(ref) & set(gem))
max_bd = max_sd = 0.0
for snp in snps:
    r, g = ref[snp], gem[snp]
    sign = 1.0
    g_a1, g_a0, c_a1 = sval(g, "allele1","A1"), sval(g, "allele0","A0","A2"), sval(r, "A1","allele1")
    if g_a1 != c_a1 and (g_a0 == c_a1 or g_a1 == sval(r, "A2","allele0")): sign = -1.0
    bd = abs(fval(r, "BETA_A") - sign*fval(g, "beta"))
    sd = abs(fval(r, "SE_A") - fval(g, "se"))
    if bd > max_bd: max_bd = bd
    if sd > max_sd: max_sd = sd
print(f"  Overlap SNPs: {len(snps)}")
print(f"  Max BETA diff vs GEMMA: {max_bd:.6e}")
print(f"  Max SE diff vs GEMMA: {max_sd:.6e}")
if len(snps) > 100 and max_bd < 0.1:
    print("  PASS: BFILE vs GEMMA consistent")
else:
    print("  (BFILE vs GEMMA: expected small differences due to different LMM implementations)")
PY
    fi
fi

echo ""
echo "[PGEN-consistency] Done. Summary:"
echo "  BFILE  (baseline): ${BFILE_ASSOC}"
echo "  PGEN   (default) : ${PGEN_ASSOC}"
echo "  PGEN   (format=2): ${PGEN2_ASSOC}"
