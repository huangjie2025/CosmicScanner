#!/bin/bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export COSMICSCANNER_ROOT="${ROOT}"
RESULT_DIR="${ROOT}/tests/results/multivariate_mouse_clean"
DATA_DIR="${ROOT}/tests/data/multivariate_mouse_clean"
OUT_PREFIX="${RESULT_DIR}/multivariate_mouse_clean"
OUT_PREFIX_GEMMA="${RESULT_DIR}/multivariate_mouse_clean.gemma"
OUT_PREFIX_GEMMA_FULL="${RESULT_DIR}/multivariate_mouse_clean.gemma_full"
OUT_PREFIX_HYBRID="${RESULT_DIR}/multivariate_mouse_clean.hybrid"
BIN="${ROOT}/build/cosmicscanner"
SUBSET_PREFIX="${DATA_DIR}/mouse_multivariate_subset"
PHENO="${RESULT_DIR}/mouse_multivariate.pheno"

mkdir -p "${RESULT_DIR}" "${DATA_DIR}"

python3 - <<'PY'
from pathlib import Path

import os
root = Path(os.environ["COSMICSCANNER_ROOT"])
src_prefix = root / "examples" / "mouse_clean" / "mouse_clean"
if not src_prefix.with_suffix(".bed").exists():
    src_prefix = root / "tests" / "data" / "mv_eb_smoke" / "mv_eb"
dst_prefix = root / "tests" / "data" / "multivariate_mouse_clean" / "mouse_multivariate_subset"
pheno_out = root / "tests" / "results" / "multivariate_mouse_clean" / "mouse_multivariate.pheno"

pheno_map = {}
with src_prefix.with_suffix(".pheno").open("r", encoding="utf-8") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 3 or parts[2] == "NA":
            continue
        try:
            y1 = float(parts[2])
            y2 = float(parts[3]) if len(parts) > 3 and parts[3] != "NA" else None
            pheno_map[(parts[0], parts[1])] = (y1, y2)
        except ValueError:
            continue

fam_rows = []
keep_indices = []
with src_prefix.with_suffix(".fam").open("r", encoding="utf-8") as f:
    for idx, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) < 6:
            continue
        key = (parts[0], parts[1])
        if key in pheno_map:
            fam_rows.append(parts)
            keep_indices.append(idx)
        if len(fam_rows) >= 24:
            break

if len(fam_rows) < 10:
    raise SystemExit("Not enough samples with phenotype for multivariate subset.")

bim_rows = []
with src_prefix.with_suffix(".bim").open("r", encoding="utf-8") as f:
    chr1 = []
    chr2 = []
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 5:
            if parts[0] == "1" and len(chr1) < 48:
                chr1.append(parts)
            elif parts[0] == "2" and len(chr2) < 48:
                chr2.append(parts)
        if len(chr1) >= 48 and len(chr2) >= 48:
            break
    bim_rows = chr1 + chr2

if len(bim_rows) < 32:
    raise SystemExit("Not enough SNPs for multivariate subset.")

src_bed = src_prefix.with_suffix(".bed")
raw = src_bed.read_bytes()
if raw[:3] != b"\x6c\x1b\x01":
    raise SystemExit("Unexpected PLINK .bed header.")

n_src = sum(1 for _ in src_prefix.with_suffix(".fam").open("r", encoding="utf-8"))
n_keep = len(keep_indices)
src_bps = (n_src + 3) // 4
dst_bps = (n_keep + 3) // 4

out = bytearray(b"\x6c\x1b\x01")
for snp_idx in range(len(bim_rows)):
    block = raw[3 + snp_idx * src_bps: 3 + (snp_idx + 1) * src_bps]
    new_block = bytearray(dst_bps)
    for new_i, old_i in enumerate(keep_indices):
        code = (block[old_i // 4] >> ((old_i % 4) * 2)) & 0b11
        new_block[new_i // 4] |= code << ((new_i % 4) * 2)
    out.extend(new_block)

dst_prefix.with_suffix(".bed").write_bytes(out)
with dst_prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for row in bim_rows:
        f.write("\t".join(row) + "\n")
with dst_prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for row in fam_rows:
        f.write("\t".join(row) + "\n")

with pheno_out.open("w", encoding="utf-8") as f:
    f.write("FID\tIID\tTrait1\tTrait2\n")
    for idx, row in enumerate(fam_rows):
        fid, iid = row[0], row[1]
        y1, y2_src = pheno_map[(fid, iid)]
        y2 = y2_src if y2_src is not None else 0.65 * y1 + (0.20 if idx % 2 == 0 else -0.15) + 0.01 * idx
        f.write(f"{fid}\t{iid}\t{y1:.6f}\t{y2:.6f}\n")
PY

"${BIN}" \
  --bfile "${SUBSET_PREFIX}" \
  --pheno "${PHENO}" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method emmv \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${OUT_PREFIX}"

python3 - <<'PY'
import csv
import math
from pathlib import Path

import os
assoc = Path(os.environ["COSMICSCANNER_ROOT"]) / "tests/results/multivariate_mouse_clean/multivariate_mouse_clean.assoc.mv.txt"
rows = list(csv.DictReader(assoc.open("r", encoding="utf-8"), delimiter="\t"))
if not rows:
    raise SystemExit("Missing multivariate GWAS result rows")

required = ["CHR", "SNP", "BP", "WALD_STAT", "P", "BETA_1", "SE_1", "BETA_2", "SE_2"]
missing = [name for name in required if name not in rows[0]]
if missing:
    raise SystemExit(f"Missing columns: {missing}")

for row in rows[:10]:
    for key in ["WALD_STAT", "P", "BETA_1", "SE_1", "BETA_2", "SE_2"]:
        value = float(row[key])
        if not math.isfinite(value):
            raise SystemExit(f"Non-finite {key}: {row[key]}")
    p = float(row["P"])
    if not 0.0 <= p <= 1.0:
        raise SystemExit(f"Invalid p-value: {p}")

print("Multivariate GWAS smoke test passed.")
PY

"${BIN}" \
  --bfile "${SUBSET_PREFIX}" \
  --pheno "${PHENO}" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method gemma \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${OUT_PREFIX_GEMMA}"

python3 - <<'PY'
import csv
import math
from pathlib import Path

import os
assoc = Path(os.environ["COSMICSCANNER_ROOT"]) / "tests/results/multivariate_mouse_clean/multivariate_mouse_clean.gemma.assoc.mv.txt"
rows = list(csv.DictReader(assoc.open("r", encoding="utf-8"), delimiter="\t"))
if not rows:
    raise SystemExit("Missing GEMMA-style multivariate GWAS result rows")
for row in rows[:10]:
    for key in ["WALD_STAT", "P", "BETA_1", "SE_1", "BETA_2", "SE_2"]:
        value = float(row[key])
        if not math.isfinite(value):
            raise SystemExit(f"Non-finite GEMMA-style {key}: {row[key]}")
    p = float(row["P"])
    if not 0.0 <= p <= 1.0:
        raise SystemExit(f"Invalid GEMMA-style p-value: {p}")
print("GEMMA-style multivariate refit smoke test passed.")
PY

"${BIN}" \
  --bfile "${SUBSET_PREFIX}" \
  --pheno "${PHENO}" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method gemma-full \
  --mv-test wald,lrt,score \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${OUT_PREFIX_GEMMA_FULL}"

python3 - <<'PY'
import csv
import math
from pathlib import Path

import os
assoc = Path(os.environ["COSMICSCANNER_ROOT"]) / "tests/results/multivariate_mouse_clean/multivariate_mouse_clean.gemma_full.assoc.mv.txt"
rows = list(csv.DictReader(assoc.open("r", encoding="utf-8"), delimiter="\t"))
if not rows:
    raise SystemExit("Missing GEMMA-full multivariate GWAS result rows")
required = ["P_WALD", "P_LRT", "P_SCORE", "P_WALD_GC", "P_LRT_GC", "P_SCORE_GC"]
missing = [name for name in required if name not in rows[0]]
if missing:
    raise SystemExit(f"Missing GEMMA-full columns: {missing}")
for row in rows[:10]:
    for key in ["P_WALD", "P_LRT", "P_SCORE", "P_WALD_GC", "P_LRT_GC", "P_SCORE_GC"]:
        value = float(row[key])
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise SystemExit(f"Invalid GEMMA-full {key}: {row[key]}")
print("GEMMA-full multivariate Wald/LRT/Score smoke test passed.")
PY

"${BIN}" \
  --bfile "${SUBSET_PREFIX}" \
  --pheno "${PHENO}" \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method hybrid \
  --mv-test wald,lrt,score \
  --mv-refine-p 1 \
  --mv-refine-top 5 \
  --gwas \
  --no-loco \
  --maf 0.0 \
  --out "${OUT_PREFIX_HYBRID}"

python3 - <<'PY'
import csv
import math
from pathlib import Path

import os
base = Path(os.environ["COSMICSCANNER_ROOT"]) / "tests/results/multivariate_mouse_clean/multivariate_mouse_clean.hybrid"
scan = Path(str(base) + ".assoc.mv.txt")
snplist = Path(str(base) + ".refine.snplist")
refined = Path(str(base) + ".refined.assoc.mv.txt")
for path in (scan, snplist, refined):
    if not path.exists() or path.stat().st_size == 0:
        raise SystemExit(f"Missing hybrid output: {path}")
rows = list(csv.DictReader(refined.open("r", encoding="utf-8"), delimiter="\t"))
if not rows:
    raise SystemExit("Missing hybrid refined rows")
if len(rows) > 5:
    raise SystemExit(f"Hybrid refined too many SNPs: {len(rows)}")
for row in rows:
    for key in ["P_WALD", "P_LRT", "P_SCORE", "P_WALD_GC", "P_LRT_GC", "P_SCORE_GC"]:
        value = float(row[key])
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise SystemExit(f"Invalid hybrid {key}: {row[key]}")
print("Hybrid EMmvGWAS + GEMMA-full refinement smoke test passed.")
PY

echo "Done. 多性状 GWAS smoke test 完成。"
