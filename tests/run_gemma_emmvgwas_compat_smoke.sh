#!/usr/bin/env bash
set -euo pipefail

ROOT="${COSMICSCANNER_ROOT:-/disk228/huangjie/PanOmics/CosmicScanner}"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
OUT_DIR="${COSMIC_SMOKE_OUT:-${ROOT}/tests/results/gemma_emmvgwas_compat_smoke}"
DATA_PREFIX="${OUT_DIR}/data/smoke"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${OUT_DIR}/data"

python3 - "${DATA_PREFIX}" <<'PY'
import math
import random
import sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260704)
n = 96
m = 360
samples = [(f"F{i+1:04d}", f"I{i+1:04d}", i % 3) for i in range(n)]
genotypes = []
variants = []
for j in range(m):
    chrom = str((j % 6) + 1)
    snp = f"gems{j+1:06d}"
    bp = 1000 + 7 * j
    base = 0.08 + 0.42 * ((j * 37) % 101) / 100.0
    row = []
    for _, _, pop in samples:
        p = min(0.95, max(0.02, base + (pop - 1) * 0.035 * math.sin(j * 0.09)))
        row.append((1 if rng.random() < p else 0) + (1 if rng.random() < p else 0))
    variants.append((chrom, snp, bp, "A", "C"))
    genotypes.append(row)

causal = {17: (0.45, 0.20), 91: (-0.35, -0.25), 177: (0.10, 0.42), 251: (-0.25, 0.30)}
phenos = []
for i, (_, _, pop) in enumerate(samples):
    shared = rng.gauss(0.0, 0.45)
    y1 = 0.25 * pop + shared + rng.gauss(0.0, 0.55)
    y2 = -0.15 * pop + 0.65 * shared + rng.gauss(0.0, 0.55)
    for snp_idx, (b1, b2) in causal.items():
        centered = genotypes[snp_idx][i] - 1.0
        y1 += b1 * centered
        y2 += b2 * centered
    phenos.append((y1, y2))

prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for (fid, iid, _), (y1, _) in zip(samples, phenos):
        f.write(f"{fid}\t{iid}\t0\t0\t1\t{y1:.8f}\n")
with prefix.with_suffix(".pheno").open("w", encoding="utf-8") as f:
    f.write("FID\tIID\tTrait1\tTrait2\n")
    for (fid, iid, _), (y1, y2) in zip(samples, phenos):
        f.write(f"{fid}\t{iid}\t{y1:.8f}\t{y2:.8f}\n")
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

bed = bytearray(b"\x6c\x1b\x01")
bytes_per_snp = (n + 3) // 4
for row in genotypes:
    block = bytearray(bytes_per_snp)
    for i, g in enumerate(row):
        block[i // 4] |= plink_code(g) << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix(".bed").write_bytes(bed)

with (prefix.parent / "snps40.txt").open("w", encoding="utf-8") as f:
    for _, snp, _, _, _ in variants[:40]:
        f.write(snp + "\n")
with (prefix.parent / "cat.tsv").open("w", encoding="utf-8") as f:
    f.write("SNP C1 C2\n")
    for idx, (_, snp, _, _, _) in enumerate(variants):
        f.write(f"{snp} {1 if idx < m // 2 else 0} {0 if idx < m // 2 else 1}\n")
with (prefix.parent / "sumstats.tsv").open("w", encoding="utf-8") as f:
    f.write("SNP N Z\n")
    for idx, (_, snp, _, _, _) in enumerate(variants):
        z = ((idx % 17) - 8) / 3.0
        f.write(f"{snp} {n} {z:.8f}\n")
PY

PHENO="${DATA_PREFIX}.pheno"
SNP40="${OUT_DIR}/data/snps40.txt"
CAT="${OUT_DIR}/data/cat.tsv"
SUMSTATS="${OUT_DIR}/data/sumstats.tsv"

run() {
  echo "[smoke] $*"
  "$@"
}

assert_file() {
  [[ -s "$1" ]] || { echo "Missing or empty output: $1" >&2; exit 1; }
}

assert_no_bad_number() {
  local file="$1"
  if grep -Ei '(^|[^A-Za-z])(nan|inf)([^A-Za-z]|$)' "$file" >/dev/null; then
    echo "Found NaN/Inf in ${file}" >&2
    exit 1
  fi
}

run "${BIN}" --bfile "${DATA_PREFIX}" -p "${PHENO}" -n 1 -lm 4 --maf 0.0 --extract "${SNP40}" -o "${OUT_DIR}/lm4"
assert_file "${OUT_DIR}/lm4.lm.assoc.txt"

run "${BIN}" --bfile "${DATA_PREFIX}" -p "${PHENO}" -n 1 -lmm 4 --maf 0.0 --extract "${SNP40}" --no-loco -o "${OUT_DIR}/lmm4"
assert_file "${OUT_DIR}/lmm4.assoc.txt"

run "${BIN}" --bfile "${DATA_PREFIX}" -p "${PHENO}" -n 1 -lmm 5 --maf 0.0 -o "${OUT_DIR}/lmm5_vc"
assert_file "${OUT_DIR}/lmm5_vc.vc.txt"

run "${BIN}" --bfile "${DATA_PREFIX}" -p "${PHENO}" -n 1 -lmm 9 --maf 0.0 --extract "${SNP40}" --no-loco -o "${OUT_DIR}/lmm9"
assert_file "${OUT_DIR}/lmm9.assoc.txt"

run "${BIN}" --bfile "${DATA_PREFIX}" --pheno "${PHENO}" --mpheno-list 1,2 --emmvgwas-vc-only --maf 0.0 -o "${OUT_DIR}/emmv_vc"
assert_file "${OUT_DIR}/emmv_vcG.csv"
assert_file "${OUT_DIR}/emmv_vcR.csv"
assert_file "${OUT_DIR}/emmv_vcbeta.csv"

run "${BIN}" --bfile "${DATA_PREFIX}" --pheno "${PHENO}" --mpheno-list 1,2 --emmvgwas \
  --emmvgwas-G "${OUT_DIR}/emmv_vcG.csv" --emmvgwas-R "${OUT_DIR}/emmv_vcR.csv" \
  --emmvgwas-output --maf 0.0 --extract "${SNP40}" -o "${OUT_DIR}/emmv_semi"
assert_file "${OUT_DIR}/emmv_semi.assoc.mv.txt"
assert_file "${OUT_DIR}/emmv_semi_association_with_ids.csv"

run "${BIN}" --bfile "${DATA_PREFIX}" --gs 1 --cat "${CAT}" --maf 0.0 -o "${OUT_DIR}/gs_cat"
assert_file "${OUT_DIR}/gs_cat.S.txt"
assert_file "${OUT_DIR}/gs_cat.size.txt"
assert_no_bad_number "${OUT_DIR}/gs_cat.S.txt"

run "${BIN}" --bfile "${DATA_PREFIX}" --gq 1 --beta "${SUMSTATS}" --cat "${CAT}" -o "${OUT_DIR}/gq_cat"
assert_file "${OUT_DIR}/gq_cat.Vq.txt"
assert_file "${OUT_DIR}/gq_cat.q.txt"
assert_no_bad_number "${OUT_DIR}/gq_cat.Vq.txt"

run "${BIN}" --bfile "${DATA_PREFIX}" --study "${OUT_DIR}/gq_cat" --ref "${OUT_DIR}/gs_cat" --vc 1 -o "${OUT_DIR}/summary_vc"
assert_file "${OUT_DIR}/summary_vc.summary.vc.txt"
assert_file "${OUT_DIR}/summary_vc.summary.vc.meta.txt"
assert_no_bad_number "${OUT_DIR}/summary_vc.summary.vc.txt"

run "${BIN}" --bfile "${DATA_PREFIX}" --study "${OUT_DIR}/gq_cat" --ref "${OUT_DIR}/gs_cat" --ci 1 -o "${OUT_DIR}/summary_ci"
assert_file "${OUT_DIR}/summary_ci.summary.ci.txt"

run "${BIN}" --bfile "${DATA_PREFIX}" -p "${PHENO}" -n 1 -bslmm 4 --extract "${SNP40}" --bslmm-hgrid 2 --bslmm-rho-grid 2 \
  -hmin 0.2 -hmax 0.8 -rmin 0.3 -rmax 0.7 -pmin -2 -pmax -0.7 -gmean 4 \
  -hscale 0.5 -rscale 0.5 -pscale 0.5 -mh 2 -smin 2 -smax 8 -o "${OUT_DIR}/bslmm4"
assert_file "${OUT_DIR}/bslmm4.bf.txt"
assert_file "${OUT_DIR}/bslmm4.hyp.txt"
assert_file "${OUT_DIR}/bslmm4.bslmm.meta.txt"

run "${BIN}" --bfile "${DATA_PREFIX}" -bslmm 5 -bf "${OUT_DIR}/bslmm4.bf.txt" -hyp "${OUT_DIR}/bslmm4.hyp.txt" \
  -pmin -2 -pmax -0.7 -gmean 4 -smin 2 -smax 8 -o "${OUT_DIR}/bslmm5"
assert_file "${OUT_DIR}/bslmm5.gamma.txt"
assert_file "${OUT_DIR}/bslmm5.coef.txt"
assert_file "${OUT_DIR}/bslmm5.hyp.txt"
assert_file "${OUT_DIR}/bslmm5.bslmm.meta.txt"

python3 - "${OUT_DIR}" <<'PY'
import csv
import math
import sys
from pathlib import Path

out = Path(sys.argv[1])

def read_tsv(path):
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def require(condition, message):
    if not condition:
        raise SystemExit(message)

def count_data_lines(path):
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip()) - 1

def assert_finite_table(path):
    text = path.read_text(encoding="utf-8")
    for token in text.replace(",", "\t").split():
        try:
            value = float(token)
        except ValueError:
            continue
        require(math.isfinite(value), f"{path} contains non-finite value: {token}")

lm_header = (out / "lm4.lm.assoc.txt").read_text(encoding="utf-8").splitlines()[0].split("\t")
require({"P_WALD", "P_LRT", "P_SCORE"}.issubset(set(lm_header)), "LM mode 4 missing Wald/LRT/Score columns")
require(count_data_lines(out / "lm4.lm.assoc.txt") == 40, "LM mode 4 did not test 40 variants")

lmm_header = (out / "lmm4.assoc.txt").read_text(encoding="utf-8").splitlines()[0].split("\t")
require({"P_WALD", "P_LRT", "P_SCORE"}.issubset(set(lmm_header)), "LMM mode 4 missing Wald/LRT/Score columns")
require(count_data_lines(out / "lmm4.assoc.txt") == 40, "LMM mode 4 did not test 40 variants")

vc_rows = read_tsv(out / "lmm5_vc.vc.txt")
require(vc_rows and {"component", "variance"}.issubset(vc_rows[0]), "LMM mode 5 VC output header mismatch")

for matrix_name in ["emmv_vcG.csv", "emmv_vcR.csv"]:
    rows = list(csv.reader((out / matrix_name).open("r", encoding="utf-8")))
    require(len(rows) == 3 and all(len(r) == 2 for r in rows), f"{matrix_name} is not a 2x2 CSV matrix with header")

require(count_data_lines(out / "emmv_semi.assoc.mv.txt") == 40, "EMmvGWAS semi route did not test 40 variants")
assoc_header = next(csv.reader((out / "emmv_semi_association_with_ids.csv").open("r", encoding="utf-8")))
require({"snp", "SNP_ID", "wald", "p", "log10p"}.issubset(set(assoc_header)), "EMmvGWAS association IDs CSV header mismatch")

s_rows = [line.split() for line in (out / "gs_cat.S.txt").read_text(encoding="utf-8").splitlines() if line.strip()]
require(len(s_rows) == 4 and all(len(row) == 2 for row in s_rows), "multi-category S file should be 4x2 for two components")
size_vals = [float(x) for x in (out / "gs_cat.size.txt").read_text(encoding="utf-8").split()]
require(size_vals == [180.0, 180.0, 96.0], f"unexpected gs_cat size vector: {size_vals}")

vq_rows = [line.split() for line in (out / "gq_cat.Vq.txt").read_text(encoding="utf-8").splitlines() if line.strip()]
q_rows = [line.split() for line in (out / "gq_cat.q.txt").read_text(encoding="utf-8").splitlines() if line.strip()]
require(len(vq_rows) == 2 and all(len(row) == 2 for row in vq_rows), "multi-category Vq file should be 2x2")
require(len(q_rows) == 2, "multi-category q file should have two rows")

meta = {}
for line in (out / "summary_vc.summary.vc.meta.txt").read_text(encoding="utf-8").splitlines():
    if "\t" in line:
        k, v = line.split("\t", 1)
        meta[k] = v
require(meta.get("n_components") == "2", "summary VC meta should report two components")
require(meta.get("S_solver") in {"full-piv-lu", "svd-pseudoinverse"}, f"unexpected S_solver: {meta.get('S_solver')}")

ci_rows = read_tsv(out / "summary_ci.summary.ci.txt")
require(len(ci_rows) >= 3, "summary CI output should include component and total rows")

bf_header = (out / "bslmm4.bf.txt").read_text(encoding="utf-8").splitlines()[0].split("\t")
require(bf_header[:4] == ["chr", "rs", "ps", "n_miss"] and len(bf_header) == 8, "BSLMM mode 4 BF header mismatch")
hyp_rows = (out / "bslmm4.hyp.txt").read_text(encoding="utf-8").splitlines()
require(len(hyp_rows) == 5, "BSLMM mode 4 HYP should contain 4 grid rows plus header")
gamma_header = (out / "bslmm5.gamma.txt").read_text(encoding="utf-8").splitlines()[0].split("\t")
require(gamma_header == ["rs", "gamma"], "BSLMM mode 5 gamma header mismatch")
for stem, mode in [("bslmm4", "4"), ("bslmm5", "5")]:
    meta = {}
    for line in (out / f"{stem}.bslmm.meta.txt").read_text(encoding="utf-8").splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            meta[k] = v
    require(meta.get("mode") == mode, f"{stem} metadata mode mismatch")
    require(abs(float(meta.get("pmin_log10_pi", "nan")) - (-2.0)) < 1e-12, f"{stem} pmin metadata mismatch")
    require(abs(float(meta.get("pmax_log10_pi", "nan")) - (-0.7)) < 1e-12, f"{stem} pmax metadata mismatch")
    require(0.0 < float(meta.get("pi_min", "nan")) <= float(meta.get("pi_max", "nan")) <= 1.0, f"{stem} pi range metadata invalid")

for path in [
    out / "lm4.lm.assoc.txt",
    out / "lmm4.assoc.txt",
    out / "emmv_semi.assoc.mv.txt",
    out / "gs_cat.S.txt",
    out / "gq_cat.Vq.txt",
    out / "summary_vc.summary.vc.txt",
    out / "summary_ci.summary.ci.txt",
    out / "bslmm4.bf.txt",
    out / "bslmm5.gamma.txt",
    out / "bslmm4.bslmm.meta.txt",
    out / "bslmm5.bslmm.meta.txt",
]:
    assert_finite_table(path)
PY

cat > "${OUT_DIR}/compat_smoke_manifest.tsv" <<EOF
route	output
gemma_lm4	${OUT_DIR}/lm4.lm.assoc.txt
gemma_lmm4	${OUT_DIR}/lmm4.assoc.txt
gemma_lmm5_vc	${OUT_DIR}/lmm5_vc.vc.txt
gemma_lmm9	${OUT_DIR}/lmm9.assoc.txt
emmvgwas_vc_estimation	${OUT_DIR}/emmv_vcG.csv;${OUT_DIR}/emmv_vcR.csv;${OUT_DIR}/emmv_vcbeta.csv
emmvgwas_semi_gwas	${OUT_DIR}/emmv_semi.assoc.mv.txt;${OUT_DIR}/emmv_semi_association_with_ids.csv
summary_stat_gs_category	${OUT_DIR}/gs_cat.S.txt;${OUT_DIR}/gs_cat.size.txt
summary_stat_gq_category	${OUT_DIR}/gq_cat.Vq.txt;${OUT_DIR}/gq_cat.q.txt
summary_stat_vc_ci	${OUT_DIR}/summary_vc.summary.vc.txt;${OUT_DIR}/summary_ci.summary.ci.txt
bslmm_dap_mode4	${OUT_DIR}/bslmm4.bf.txt;${OUT_DIR}/bslmm4.hyp.txt
bslmm_dap_mode5	${OUT_DIR}/bslmm5.gamma.txt;${OUT_DIR}/bslmm5.coef.txt;${OUT_DIR}/bslmm5.hyp.txt
bslmm_metadata	${OUT_DIR}/bslmm4.bslmm.meta.txt;${OUT_DIR}/bslmm5.bslmm.meta.txt
EOF
assert_file "${OUT_DIR}/compat_smoke_manifest.tsv"

echo "GEMMA/EMmvGWAS compatibility smoke passed. Outputs: ${OUT_DIR}"
