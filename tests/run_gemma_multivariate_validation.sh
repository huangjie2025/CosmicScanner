#!/usr/bin/env bash
set -euo pipefail

ROOT="/disk228/huangjie/PanOmics/CosmicScanner"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
RESULT_DIR="${ROOT}/tests/results/multivariate_gemma_validation"
DATA_PREFIX="${RESULT_DIR}/data/sim_mv"
PHENO_SRC="${DATA_PREFIX}.pheno"
PHENO_GXE="${DATA_PREFIX}.gxe.pheno"
GEMMA_WORK="${RESULT_DIR}/gemma_work"
PHENO_NO_HEADER="${RESULT_DIR}/mouse_multivariate.no_header.pheno"
SUMMARY="${RESULT_DIR}/multivariate_validation.summary.txt"
EMMV_PREFIX="${RESULT_DIR}/emmvgwas/emmv"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

resolve_gemma() {
  if [[ -n "${GEMMA_BIN:-}" ]]; then
    printf '%s\n' "${GEMMA_BIN}"
    return 0
  fi
  local conda_bin="/disk228/huangjie/software/miniconda3/envs/benchmark/bin/gemma"
  if [[ -x "${conda_bin}" ]]; then
    printf '%s\n' "${conda_bin}"
    return 0
  fi
  if command -v gemma >/dev/null 2>&1; then
    command -v gemma
    return 0
  fi
  return 1
}

resolve_rscript() {
  if [[ -n "${RSCRIPT_BIN:-}" ]]; then
    printf '%s\n' "${RSCRIPT_BIN}"
    return 0
  fi
  local rgs="/disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript"
  if [[ -x "${rgs}" ]]; then
    printf '%s\n' "${rgs}"
    return 0
  fi
  if command -v Rscript >/dev/null 2>&1; then
    command -v Rscript
    return 0
  fi
  return 1
}

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi
if ! GEMMA_EXE="$(resolve_gemma)"; then
  echo "Missing GEMMA executable." >&2
  exit 2
fi
if RSCRIPT_EXE="$(resolve_rscript)" && "${RSCRIPT_EXE}" -e 'ok <- requireNamespace("EMmvGWAS", quietly = TRUE); if (ok) { cl <- try(parallel::makeCluster(1L), silent = TRUE); ok <- !inherits(cl, "try-error"); if (ok) parallel::stopCluster(cl) }; quit(status = ifelse(ok, 0L, 1L))' >/dev/null 2>&1; then
  HAVE_EMMVGWAS_R=1
else
  HAVE_EMMVGWAS_R=0
fi

mkdir -p "${RESULT_DIR}" "${RESULT_DIR}/data" "${GEMMA_WORK}/output"
rm -rf "${GEMMA_WORK}/output"
rm -rf "${RESULT_DIR}/emmvgwas"
mkdir -p "${GEMMA_WORK}/output"
mkdir -p "${RESULT_DIR}/emmvgwas"

python3 - "${DATA_PREFIX}" <<'PY'
import math
import random
import struct
import sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260603)
n = 96
m = 360
n_pop = 3

samples = []
for i in range(n):
    pop = i % n_pop
    fid = f"F{i+1:04d}"
    iid = f"I{i+1:04d}"
    samples.append((fid, iid, pop))

variants = []
genotypes = []
for j in range(m):
    chrom = str((j % 6) + 1)
    snp = f"mvs{j+1:06d}"
    bp = 1000 + j * 11
    base = 0.08 + 0.42 * ((j * 29) % 101) / 100.0
    row = []
    for _, _, pop in samples:
        shift = (pop - 1) * 0.045 * math.sin(j * 0.11)
        p = min(0.95, max(0.02, base + shift))
        g = (1 if rng.random() < p else 0) + (1 if rng.random() < p else 0)
        row.append(g)
    variants.append((chrom, snp, bp, "A", "C"))
    genotypes.append(row)

causal = {17: (0.45, 0.20), 91: (-0.35, -0.25), 177: (0.10, 0.42), 251: (-0.25, 0.30)}
phenos = []
for i, (_, _, pop) in enumerate(samples):
    shared = rng.gauss(0.0, 0.45)
    e1 = rng.gauss(0.0, 0.55)
    e2 = rng.gauss(0.0, 0.55)
    y1 = 0.25 * pop + shared + e1
    y2 = -0.15 * pop + 0.65 * shared + e2
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
with prefix.with_suffix(".gxe.pheno").open("w", encoding="utf-8") as f:
    f.write("FID\tIID\tEnv\tTrait1\tTrait2\n")
    for (idx, ((fid, iid, pop), (y1, y2))) in enumerate(zip(samples, phenos)):
        env = -1.0 + pop + 0.15 * math.sin(idx * 0.17)
        f.write(f"{fid}\t{iid}\t{env:.8f}\t{y1:.8f}\t{y2:.8f}\n")
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

with prefix.with_suffix(".emmv.P.csv").open("w", encoding="utf-8") as f:
    f.write("Trait1,Trait2\n")
    for y1, y2 in phenos:
        f.write(f"{y1:.8f},{y2:.8f}\n")
with prefix.with_suffix(".emmv.Z.csv").open("w", encoding="utf-8") as f:
    f.write("SNP," + ",".join(iid for _, iid, _ in samples) + "\n")
    for (chrom, snp, bp, a1, a2), row in zip(variants, genotypes):
        centered = [str(g - 1) for g in row]
        f.write(snp + "," + ",".join(centered) + "\n")
PY

python3 - "${PHENO_SRC}" "${PHENO_NO_HEADER}" <<'PY'
import sys
src, dst = sys.argv[1:]
with open(src, "r", encoding="utf-8") as fin, open(dst, "w", encoding="utf-8") as fout:
    first = True
    for line in fin:
        if first:
            first = False
            if line.lower().startswith("fid"):
                continue
        parts = line.strip().split()
        if len(parts) >= 4:
            fout.write(f"{parts[2]}\t{parts[3]}\n")
PY

run_timed() {
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
    f.write(f"time_{label}\t{elapsed:.6f}\n")
if proc.returncode != 0:
    raise SystemExit(proc.returncode)
PY
}

rm -f "${SUMMARY}"

run_timed cosmic_emmv "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${PHENO_SRC}" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method emmv \
  --mv-compatible emmvgwas \
  --mv-emmvgwas-output \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${RESULT_DIR}/cosmic_emmv"

run_timed cosmic_emmv_semi "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${PHENO_SRC}" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method emmv \
  --mv-compatible emmvgwas \
  --mv-emmvgwas-G "${RESULT_DIR}/cosmic_emmvG.csv" \
  --mv-emmvgwas-R "${RESULT_DIR}/cosmic_emmvR.csv" \
  --mv-emmvgwas-output \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${RESULT_DIR}/cosmic_emmv_semi"

run_timed cosmic_gemma_full "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${PHENO_SRC}" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method gemma-full \
  --mv-test wald,lrt,score \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${RESULT_DIR}/cosmic_gemma_full"

run_timed cosmic_gemma_exact "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${PHENO_SRC}" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method gemma-exact \
  --mv-gemma-mode 4 \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${RESULT_DIR}/cosmic_gemma_exact"

run_timed cosmic_gemma_exact_crt "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${PHENO_SRC}" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method gemma-exact \
  --mv-gemma-mode 4 \
  --mv-crt \
  --mv-pnr 1 \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${RESULT_DIR}/cosmic_gemma_exact_crt"

run_timed cosmic_gemma_gxe "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${PHENO_GXE}" \
  --mpheno-list 2,3 \
  --env-covar 3 \
  --multivariate \
  --mv-method gemma-gxe \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${RESULT_DIR}/cosmic_gemma_gxe"

run_timed cosmic_gemma_gxe_full "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${PHENO_GXE}" \
  --mpheno-list 2,3 \
  --env-covar 3 \
  --multivariate \
  --mv-method gemma-gxe-full \
  --mv-gemma-mode 4 \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${RESULT_DIR}/cosmic_gemma_gxe_full"

run_timed cosmic_gemma_full_fast "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${PHENO_SRC}" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method gemma-full \
  --mv-test wald,lrt,score \
  --mv-refit-p 0.05 \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${RESULT_DIR}/cosmic_gemma_full_fast"

for mode in 1 2 3; do
  run_timed "cosmic_gemma_mode${mode}" "${BIN}" \
    --bfile "${DATA_PREFIX}" \
    --pheno "${PHENO_SRC}" \
    --mpheno-list 1,2 \
    --multivariate \
    --mv-method gemma-full \
    --mv-gemma-mode "${mode}" \
    --gwas \
    --no-loco \
    --maf 0.0 \
    --out "${RESULT_DIR}/cosmic_gemma_mode${mode}"
done

run_timed cosmic_hybrid "${BIN}" \
  --bfile "${DATA_PREFIX}" \
  --pheno "${PHENO_SRC}" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method hybrid \
  --mv-test wald,lrt,score \
  --mv-refine-p 1 \
  --mv-refine-top 5 \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${RESULT_DIR}/cosmic_hybrid"

if [[ "${HAVE_EMMVGWAS_R}" -eq 1 ]]; then
  run_timed emmvgwas "${RSCRIPT_EXE}" "${ROOT}/tests/run_emmvgwas_validation.R" "${DATA_PREFIX}.emmv.P.csv" "${DATA_PREFIX}.emmv.Z.csv" "${EMMV_PREFIX}"
else
  printf 'emmvgwas_preflight\tblocked_no_Rscript_EMmvGWAS_package_or_R_socket\n' >> "${SUMMARY}"
fi

(
  cd "${GEMMA_WORK}"
  python3 - "gemma_kinship" "${SUMMARY}" "${GEMMA_EXE}" -bfile "${DATA_PREFIX}" -gk 1 -o gemma_mv_kinship <<'PY'
import subprocess
import sys
import time
label, summary = sys.argv[1], sys.argv[2]
cmd = sys.argv[3:]
t0 = time.time()
proc = subprocess.run(cmd)
elapsed = time.time() - t0
with open(summary, "a", encoding="utf-8") as f:
    f.write(f"time_{label}\t{elapsed:.6f}\n")
if proc.returncode != 0:
    raise SystemExit(proc.returncode)
PY
)

KINSHIP="${GEMMA_WORK}/output/gemma_mv_kinship.cXX.txt"
if [[ ! -f "${KINSHIP}" ]]; then
  echo "Missing GEMMA kinship output: ${KINSHIP}" >&2
  exit 3
fi

(
  cd "${GEMMA_WORK}"
  python3 - "gemma_mv" "${SUMMARY}" "${GEMMA_EXE}" -bfile "${DATA_PREFIX}" -p "${PHENO_NO_HEADER}" -n 1 2 -k "${KINSHIP}" -lmm -o gemma_mv_lmm <<'PY'
import subprocess
import sys
import time
label, summary = sys.argv[1], sys.argv[2]
cmd = sys.argv[3:]
t0 = time.time()
proc = subprocess.run(cmd)
elapsed = time.time() - t0
with open(summary, "a", encoding="utf-8") as f:
    f.write(f"time_{label}\t{elapsed:.6f}\n")
if proc.returncode != 0:
    raise SystemExit(proc.returncode)
PY
)

python3 - "${RESULT_DIR}" "${GEMMA_WORK}/output/gemma_mv_lmm.assoc.txt" "${SUMMARY}" <<'PY'
import csv
import math
import sys
from pathlib import Path

result_dir = Path(sys.argv[1])
gemma_path = Path(sys.argv[2])
summary_path = Path(sys.argv[3])

def read_tsv(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def read_csv(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def fval(row, *names):
    lower = {k.lower(): k for k in row}
    for name in names:
        key = lower.get(name.lower())
        if key and row[key] not in ("", "NA", "nan", "-nan"):
            return float(row[key])
    return math.nan

def corr(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 3:
        return math.nan
    xs, ys = zip(*pairs)
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx <= 0 or vy <= 0:
        return math.nan
    return sum((x - mx) * (y - my) for x, y in pairs) / math.sqrt(vx * vy)

gemma = {row.get("rs", row.get("SNP")): row for row in read_tsv(gemma_path)}
emmvgwas_path = result_dir / "emmvgwas" / "emmv_association_with_ids.csv"

comparisons = []
for name, rel, pcol in [
    ("emmv", "cosmic_emmv.assoc.mv.txt", "P"),
    ("gemma_full", "cosmic_gemma_full.assoc.mv.txt", "P_WALD"),
    ("gemma_exact", "cosmic_gemma_exact.assoc.mv.txt", "P_WALD"),
]:
    cosmic_rows = read_tsv(result_dir / rel)
    overlap = []
    p_cosmic, p_gemma = [], []
    b1_cosmic, b1_gemma = [], []
    b2_cosmic, b2_gemma = [], []
    for row in cosmic_rows:
        snp = row.get("SNP")
        if snp not in gemma:
            continue
        grow = gemma[snp]
        overlap.append(snp)
        pc = fval(row, pcol)
        pg = fval(grow, "p_wald")
        p_cosmic.append(-math.log10(max(pc, 1e-300)) if math.isfinite(pc) else math.nan)
        p_gemma.append(-math.log10(max(pg, 1e-300)) if math.isfinite(pg) else math.nan)
        b1_cosmic.append(fval(row, "BETA_1"))
        b1_gemma.append(fval(grow, "beta_1"))
        b2_cosmic.append(fval(row, "BETA_2"))
        b2_gemma.append(fval(grow, "beta_2"))
    comparisons.append((name, len(overlap), corr(p_cosmic, p_gemma), corr(b1_cosmic, b1_gemma), corr(b2_cosmic, b2_gemma)))

reference_thresholds = {
    # These are intentionally just below the empirical validation values to catch
    # route regressions while allowing tiny BLAS/platform drift.
    "emmv": (360, 0.98, 0.98, 0.98),
    "gemma_full": (360, 0.995, 0.995, 0.995),
    "gemma_exact": (360, 0.995, 0.995, 0.995),
}

with open(summary_path, "a", encoding="utf-8") as f:
    for name, overlap, pcorr, b1corr, b2corr in comparisons:
        f.write(f"{name}_overlap\t{overlap}\n")
        f.write(f"{name}_log10p_corr_vs_gemma\t{pcorr:.12g}\n")
        f.write(f"{name}_beta1_corr_vs_gemma\t{b1corr:.12g}\n")
        f.write(f"{name}_beta2_corr_vs_gemma\t{b2corr:.12g}\n")
        min_overlap, min_pcorr, min_b1corr, min_b2corr = reference_thresholds[name]
        if overlap != min_overlap or pcorr < min_pcorr or b1corr < min_b1corr or b2corr < min_b2corr:
            raise SystemExit(
                f"{name} reference validation failed: overlap={overlap}, "
                f"log10p_corr={pcorr}, beta1_corr={b1corr}, beta2_corr={b2corr}"
            )

    fast_path = result_dir / "cosmic_gemma_full_fast.assoc.mv.txt"
    if fast_path.exists():
        fast_rows = read_tsv(fast_path)
        full_refit = sum(1 for row in fast_rows if math.isfinite(fval(row, "P_WALD")))
        score_only = sum(1 for row in fast_rows if not math.isfinite(fval(row, "P_WALD")) and math.isfinite(fval(row, "P_SCORE")))
        f.write(f"gemma_full_fast_rows\t{len(fast_rows)}\n")
        f.write(f"gemma_full_fast_full_refit\t{full_refit}\n")
        f.write(f"gemma_full_fast_score_only\t{score_only}\n")

    for mode, finite_col, nan_cols in [
        (1, "P_WALD", ["P_LRT", "P_SCORE"]),
        (2, "P_LRT", ["P_WALD", "P_SCORE"]),
        (3, "P_SCORE", ["P_WALD", "P_LRT"]),
    ]:
        mode_path = result_dir / f"cosmic_gemma_mode{mode}.assoc.mv.txt"
        rows = read_tsv(mode_path)
        finite_count = sum(1 for row in rows if math.isfinite(fval(row, finite_col)))
        nan_ok = all(not math.isfinite(fval(row, col)) for row in rows for col in nan_cols)
        f.write(f"gemma_mode{mode}_rows\t{len(rows)}\n")
        f.write(f"gemma_mode{mode}_{finite_col}_finite\t{finite_count}\n")
        f.write(f"gemma_mode{mode}_other_p_nan\t{str(nan_ok).lower()}\n")
        if finite_count != len(rows) or not nan_ok:
            raise SystemExit(f"GEMMA mode {mode} output columns failed validation")

    exact_rows = read_tsv(result_dir / "cosmic_gemma_exact.assoc.mv.txt")
    exact_finite = {
        col: sum(1 for row in exact_rows if math.isfinite(fval(row, col)))
        for col in ["P_WALD", "P_LRT", "P_SCORE"]
    }
    f.write(f"gemma_exact_rows\t{len(exact_rows)}\n")
    for col, count in exact_finite.items():
        f.write(f"gemma_exact_{col}_finite\t{count}\n")
        if count != len(exact_rows):
            raise SystemExit(f"GEMMA exact {col} finite-count validation failed")

    exact_crt_rows = read_tsv(result_dir / "cosmic_gemma_exact_crt.assoc.mv.txt")
    exact_crt_finite = {
        col: sum(1 for row in exact_crt_rows if math.isfinite(fval(row, col)))
        for col in ["P_WALD_CRT", "P_LRT_CRT", "P_SCORE_CRT"]
    }
    f.write(f"gemma_exact_crt_rows\t{len(exact_crt_rows)}\n")
    for col, count in exact_crt_finite.items():
        f.write(f"gemma_exact_crt_{col}_finite\t{count}\n")
        if count != len(exact_crt_rows):
            raise SystemExit(f"GEMMA exact CRT {col} finite-count validation failed")

    gxe_rows = read_tsv(result_dir / "cosmic_gemma_gxe.assoc.mv.gxe.txt")
    gxe_p = sum(1 for row in gxe_rows if math.isfinite(fval(row, "P_GXE")))
    gxe_wald = sum(1 for row in gxe_rows if math.isfinite(fval(row, "WALD_GXE")))
    gxe_b1 = sum(1 for row in gxe_rows if math.isfinite(fval(row, "BETA_GXE_1")))
    gxe_b2 = sum(1 for row in gxe_rows if math.isfinite(fval(row, "BETA_GXE_2")))
    f.write(f"gemma_gxe_rows\t{len(gxe_rows)}\n")
    f.write(f"gemma_gxe_P_GXE_finite\t{gxe_p}\n")
    f.write(f"gemma_gxe_WALD_GXE_finite\t{gxe_wald}\n")
    f.write(f"gemma_gxe_BETA_GXE_1_finite\t{gxe_b1}\n")
    f.write(f"gemma_gxe_BETA_GXE_2_finite\t{gxe_b2}\n")
    if min(gxe_p, gxe_wald, gxe_b1, gxe_b2) != len(gxe_rows) or len(gxe_rows) != 360:
        raise SystemExit("CosmicScanner GEMMA-style multivariate GxE validation failed")

    gxe_full_rows = read_tsv(result_dir / "cosmic_gemma_gxe_full.assoc.mv.gxe.txt")
    gxe_full_counts = {
        col: sum(1 for row in gxe_full_rows if math.isfinite(fval(row, col)))
        for col in ["P_WALD_GXE", "P_LRT_GXE", "P_SCORE_GXE", "BETA_GXE_1", "BETA_GXE_2"]
    }
    f.write(f"gemma_gxe_full_rows\t{len(gxe_full_rows)}\n")
    for col, count in gxe_full_counts.items():
        f.write(f"gemma_gxe_full_{col}_finite\t{count}\n")
        if count != len(gxe_full_rows):
            raise SystemExit(f"CosmicScanner GEMMA-style multivariate GxE full {col} validation failed")
    if len(gxe_full_rows) != 360:
        raise SystemExit("CosmicScanner GEMMA-style multivariate GxE full row-count validation failed")

    full_emmv_rows = {row.get("SNP"): row for row in read_tsv(result_dir / "cosmic_emmv.assoc.mv.txt")}
    semi_emmv_rows = read_tsv(result_dir / "cosmic_emmv_semi.assoc.mv.txt")
    semi_p_full, semi_p_fixed = [], []
    semi_b1_full, semi_b1_fixed = [], []
    semi_b2_full, semi_b2_fixed = [], []
    semi_overlap = 0
    for row in semi_emmv_rows:
        snp = row.get("SNP")
        if snp not in full_emmv_rows:
            continue
        full = full_emmv_rows[snp]
        semi_overlap += 1
        pf = fval(full, "P")
        ps = fval(row, "P")
        semi_p_full.append(-math.log10(max(pf, 1e-300)) if math.isfinite(pf) else math.nan)
        semi_p_fixed.append(-math.log10(max(ps, 1e-300)) if math.isfinite(ps) else math.nan)
        semi_b1_full.append(fval(full, "BETA_1"))
        semi_b1_fixed.append(fval(row, "BETA_1"))
        semi_b2_full.append(fval(full, "BETA_2"))
        semi_b2_fixed.append(fval(row, "BETA_2"))
    semi_pcorr = corr(semi_p_full, semi_p_fixed)
    semi_b1corr = corr(semi_b1_full, semi_b1_fixed)
    semi_b2corr = corr(semi_b2_full, semi_b2_fixed)
    f.write(f"emmv_semi_overlap_vs_full\t{semi_overlap}\n")
    f.write(f"emmv_semi_log10p_corr_vs_full\t{semi_pcorr:.12g}\n")
    f.write(f"emmv_semi_beta1_corr_vs_full\t{semi_b1corr:.12g}\n")
    f.write(f"emmv_semi_beta2_corr_vs_full\t{semi_b2corr:.12g}\n")
    if semi_overlap != len(full_emmv_rows) or min(semi_pcorr, semi_b1corr, semi_b2corr) < 0.999999:
        raise SystemExit("CosmicScanner EMmvGWAS semi_gwas fixed G/R route failed validation")

    if emmvgwas_path.exists():
        emmvgwas = {row.get("SNP_ID"): row for row in read_csv(emmvgwas_path)}
        cosmic_rows = read_tsv(result_dir / "cosmic_emmv.assoc.mv.txt")
        cosmic_emmv_csv = read_csv(result_dir / "cosmic_emmv_association.csv")
        cosmic_emmv_ids_csv = read_csv(result_dir / "cosmic_emmv_association_with_ids.csv")
        cosmic_g_rows = read_csv(result_dir / "cosmic_emmvG.csv")
        cosmic_r_rows = read_csv(result_dir / "cosmic_emmvR.csv")
        cosmic_beta_rows = read_csv(result_dir / "cosmic_emmvbeta.csv")
        p_cosmic, p_emmv = [], []
        b1_cosmic, b1_emmv = [], []
        b2_cosmic, b2_emmv = [], []
        overlap = 0
        for row in cosmic_rows:
            snp = row.get("SNP")
            if snp not in emmvgwas:
                continue
            erow = emmvgwas[snp]
            overlap += 1
            pc = fval(row, "P")
            pe = fval(erow, "p")
            p_cosmic.append(-math.log10(max(pc, 1e-300)) if math.isfinite(pc) else math.nan)
            p_emmv.append(-math.log10(max(pe, 1e-300)) if math.isfinite(pe) else math.nan)
            b1_cosmic.append(fval(row, "BETA_1"))
            b1_emmv.append(fval(erow, "Trait1", "trait1"))
            b2_cosmic.append(fval(row, "BETA_2"))
            b2_emmv.append(fval(erow, "Trait2", "trait2"))
        f.write("emmvgwas_status\tok\n")
        f.write(f"cosmic_emmvgwas_assoc_rows\t{len(cosmic_emmv_csv)}\n")
        f.write(f"cosmic_emmvgwas_assoc_with_ids_rows\t{len(cosmic_emmv_ids_csv)}\n")
        f.write(f"cosmic_emmvgwas_G_rows\t{len(cosmic_g_rows)}\n")
        f.write(f"cosmic_emmvgwas_R_rows\t{len(cosmic_r_rows)}\n")
        f.write(f"cosmic_emmvgwas_beta_rows\t{len(cosmic_beta_rows)}\n")
        f.write(f"emmv_overlap_vs_emmvgwas\t{overlap}\n")
        emmv_pcorr = corr(p_cosmic, p_emmv)
        emmv_b1corr = corr(b1_cosmic, b1_emmv)
        emmv_b2corr = corr(b2_cosmic, b2_emmv)
        f.write(f"emmv_log10p_corr_vs_emmvgwas\t{emmv_pcorr:.12g}\n")
        f.write(f"emmv_beta1_corr_vs_emmvgwas\t{emmv_b1corr:.12g}\n")
        f.write(f"emmv_beta2_corr_vs_emmvgwas\t{emmv_b2corr:.12g}\n")
        if len(cosmic_emmv_csv) != len(cosmic_rows) or len(cosmic_emmv_ids_csv) != len(cosmic_rows):
            raise SystemExit("CosmicScanner EMmvGWAS-style association CSV row count mismatch")
        if len(cosmic_g_rows) != 2 or len(cosmic_r_rows) != 2 or len(cosmic_beta_rows) != 2:
            raise SystemExit("CosmicScanner EMmvGWAS-style G/R/beta CSV shape mismatch")
        if overlap != len(cosmic_rows) or min(emmv_pcorr, emmv_b1corr, emmv_b2corr) < 0.999999:
            raise SystemExit("CosmicScanner EMmvGWAS R-package reference validation failed")
    else:
        f.write("emmvgwas_status\tblocked_no_Rscript_or_output\n")

print(summary_path.read_text(encoding="utf-8"))
PY

echo "Done. Multivariate GEMMA validation complete: ${SUMMARY}"
