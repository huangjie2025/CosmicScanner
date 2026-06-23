#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/psr_ordinal_smoke"
RESULT_DIR="${ROOT}/tests/results/psr_ordinal_smoke"
PREFIX="${DATA_DIR}/psr_ordinal"
OUT="${RESULT_DIR}/cosmic_ordinal"

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
rng = random.Random(20260624)
n = 90
m = 150
samples = [(f"F{i+1:04d}", f"I{i+1:04d}", 0 if i < n // 2 else 1) for i in range(n)]
variants = []
genotypes = []
for j in range(m):
    p = 0.10 + 0.40 * ((j * 19) % 101) / 100.0
    row = []
    for _, _, pop in samples:
        pp = min(0.95, max(0.02, p + (0.05 if pop else -0.025)))
        row.append((1 if rng.random() < pp else 0) + (1 if rng.random() < pp else 0))
    variants.append((str((j % 3) + 1), f"ord{j+1:06d}", 3000 + j * 10, "A", "C"))
    genotypes.append(row)

causal = {21: 0.85, 99: -0.65}
ordinal = []
for i, (_, _, pop) in enumerate(samples):
    eta = -0.15 + 0.35 * pop + rng.gauss(0.0, 0.8)
    for idx, eff in causal.items():
        eta += eff * (genotypes[idx][i] - 1.0)
    if eta < -0.8:
        y = 1
    elif eta < 0.0:
        y = 2
    elif eta < 0.8:
        y = 3
    else:
        y = 4
    ordinal.append(y)

prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for fid, iid, _ in samples:
        f.write(f"{fid}\t{iid}\t0\t0\t1\t-9\n")
with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for chrom, snp, bp, a1, a2 in variants:
        f.write(f"{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n")
with prefix.with_suffix(".ordinal.pheno").open("w", encoding="utf-8") as f:
    for (fid, iid, _), y in zip(samples, ordinal):
        f.write(f"{fid}\t{iid}\t{y}\n")

def plink_code(dosage):
    return 0b00 if dosage == 2 else (0b10 if dosage == 1 else 0b11)

bytes_per_snp = (n + 3) // 4
bed = bytearray(b"\x6c\x1b\x01")
for row in genotypes:
    block = bytearray(bytes_per_snp)
    for i, g in enumerate(row):
        block[i // 4] |= plink_code(g) << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix(".bed").write_bytes(bed)
PY

"${BIN}" --bfile "${PREFIX}" --gwas --no-loco --maf 0.0 \
  --ordinal-model --ordinal-pheno "${PREFIX}.ordinal.pheno" --link logit \
  --psr-maxit 25 --out "${OUT}"

python3 - "${OUT}" "${RESULT_DIR}" <<'PY'
import csv
import sys
from pathlib import Path

prefix = Path(sys.argv[1])
result_dir = Path(sys.argv[2])
gwas = prefix.with_suffix(".gwas.tsv")
conv = Path(str(prefix) + ".psr.convergence.tsv")
thr = Path(str(prefix) + ".ordinal.thresholds.tsv")
for path in (gwas, conv, thr):
    if not path.exists() or path.stat().st_size == 0:
        raise SystemExit(f"Missing ordinal output: {path}")

with gwas.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
required = {"trait_type", "ordinal_mode", "snp", "effect", "se", "wald_or_stat", "p_value", "log10p", "n_levels", "n_thresholds", "converged", "status_class", "runtime"}
missing = required - set(rows[0])
if missing:
    raise SystemExit(f"Ordinal GWAS schema missing: {sorted(missing)}")
if not rows or any(r["ordinal_mode"] != "experimental_cumulative_binary" for r in rows[:5]):
    raise SystemExit("Ordinal output did not mark experimental cumulative-binary mode")
if int(rows[0]["n_levels"]) < 3 or int(rows[0]["n_thresholds"]) < 2:
    raise SystemExit("Ordinal output has insufficient level/threshold metadata")

summary = result_dir / "ordinal_validation_summary.tsv"
with summary.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["trait_type", "ordinal_mode", "n_rows", "n_levels", "n_thresholds", "validation_status"])
    w.writerow(["ordinal", "experimental_cumulative_binary", len(rows), rows[0]["n_levels"], rows[0]["n_thresholds"], "PASS"])
print(f"[PSR-GWAS] Ordinal validation summary written to {summary}")
PY

RSCRIPT_EXE="${RSCRIPT_BIN:-}"
if [[ -z "${RSCRIPT_EXE}" ]] && command -v Rscript >/dev/null 2>&1; then
  RSCRIPT_EXE="$(command -v Rscript)"
fi
if [[ -z "${RSCRIPT_EXE}" && -x /disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript ]]; then
  RSCRIPT_EXE="/disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript"
fi

if [[ -n "${RSCRIPT_EXE}" ]]; then
  "${RSCRIPT_EXE}" - "${PREFIX}.ordinal.pheno" "${RESULT_DIR}/ordinal_r_baseline_availability.tsv" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
out <- args[2]
available <- requireNamespace("MASS", quietly=TRUE)
err <- ""
if (available) {
  dat <- tryCatch(read.table(args[1], col.names=c("FID","IID","ORDINAL")), error=function(e) e)
  if (inherits(dat, "error")) {
    available <- FALSE
    err <- conditionMessage(dat)
  } else {
    fit <- tryCatch(MASS::polr(as.ordered(ORDINAL) ~ 1, data=dat, Hess=TRUE), error=function(e) e)
    if (inherits(fit, "error")) {
      available <- FALSE
      err <- conditionMessage(fit)
    }
  }
} else {
  err <- "R package MASS is not available"
}
write.table(data.frame(trait_type="ordinal", method="MASS::polr", method_available=available, error_message=err), out, sep="\t", quote=FALSE, row.names=FALSE)
RS
else
  printf 'trait_type\tmethod\tmethod_available\terror_message\nordinal\tMASS::polr\tFALSE\tRscript not found\n' > "${RESULT_DIR}/ordinal_r_baseline_availability.tsv"
fi

echo "[PSR-GWAS] Ordinal validation passed."
