#!/bin/bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="/disk228/huangjie/PanOmics/CosmicScanner"
RESULT_DIR="${ROOT}/tests/results/repeatability_mouse_clean"
DATA_DIR="${ROOT}/tests/data/repeatability_mouse_clean"
OUT_PREFIX_EXACT="${RESULT_DIR}/repeatability_mouse_clean.exact"
OUT_PREFIX_EXACT_CACHED="${RESULT_DIR}/repeatability_mouse_clean.exact_cached"
OUT_PREFIX_LOCO="${RESULT_DIR}/repeatability_mouse_clean.loco"
OUT_PREFIX_LOCO_CACHED="${RESULT_DIR}/repeatability_mouse_clean.loco_cached"
OUT_PREFIX_FAST="${RESULT_DIR}/repeatability_mouse_clean.fast"
BIN="${ROOT}/build/cosmicscanner"
LONG_PHENO="${RESULT_DIR}/mouse_clean.repeatability.pheno"
SUBSET_PREFIX="${DATA_DIR}/mouse_repeatability_subset"
GRM_CACHE_PREFIX="${DATA_DIR}/repeatability_cache/grm"

mkdir -p "${RESULT_DIR}"
mkdir -p "${DATA_DIR}"
rm -rf "${DATA_DIR}/repeatability_cache"

python3 - <<'PY'
from pathlib import Path

root = Path("/disk228/huangjie/PanOmics/CosmicScanner")
src_prefix = root / "examples" / "mouse_clean" / "mouse_clean"
dst_prefix = root / "tests" / "data" / "repeatability_mouse_clean" / "mouse_repeatability_subset"
dst = root / "tests" / "results" / "repeatability_mouse_clean" / "mouse_clean.repeatability.pheno"

pheno_map = {}
with (src_prefix.with_suffix(".pheno")).open("r", encoding="utf-8") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 3 or parts[2] == "NA":
            continue
        try:
            pheno_map[(parts[0], parts[1])] = float(parts[2])
        except ValueError:
            continue

fam_rows = []
keep_indices = []
with (src_prefix.with_suffix(".fam")).open("r", encoding="utf-8") as f:
    for idx, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) < 6:
            continue
        key = (parts[0], parts[1])
        if key in pheno_map:
            fam_rows.append(parts)
            keep_indices.append(idx)
        if len(fam_rows) >= 20:
            break

if len(fam_rows) < 10:
    raise SystemExit("Not enough samples with non-missing phenotype for repeatability subset.")

bim_rows = []
with (src_prefix.with_suffix(".bim")).open("r", encoding="utf-8") as f:
    chr1 = []
    chr2 = []
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 6:
            if parts[0] == "1" and len(chr1) < 64:
                chr1.append(parts)
            elif parts[0] == "2" and len(chr2) < 64:
                chr2.append(parts)
        if len(chr1) >= 64 and len(chr2) >= 64:
            break
    bim_rows = chr1 + chr2

if len(bim_rows) < 32:
    raise SystemExit("Not enough SNPs for repeatability subset.")
chroms = {row[0] for row in bim_rows}
if len(chroms) < 2:
    raise SystemExit("Repeatability subset must contain at least two chromosomes for LOCO testing.")

src_bed = src_prefix.with_suffix(".bed")
dst_bed = dst_prefix.with_suffix(".bed")
dst_bim = dst_prefix.with_suffix(".bim")
dst_fam = dst_prefix.with_suffix(".fam")

raw = src_bed.read_bytes()
if raw[:3] != b"\x6c\x1b\x01":
    raise SystemExit("Unexpected PLINK .bed header.")

n_src = 1940
n_keep = len(keep_indices)
m_keep = len(bim_rows)
src_bps = (n_src + 3) // 4
dst_bps = (n_keep + 3) // 4

out = bytearray(b"\x6c\x1b\x01")
for snp_idx in range(m_keep):
    block = raw[3 + snp_idx * src_bps: 3 + (snp_idx + 1) * src_bps]
    new_block = bytearray(dst_bps)
    for new_i, old_i in enumerate(keep_indices):
        code = (block[old_i // 4] >> ((old_i % 4) * 2)) & 0b11
        new_block[new_i // 4] |= (code << ((new_i % 4) * 2))
    out.extend(new_block)

dst_bed.write_bytes(out)
with dst_bim.open("w", encoding="utf-8") as f:
    for row in bim_rows:
        f.write("\t".join(row) + "\n")
with dst_fam.open("w", encoding="utf-8") as f:
    for row in fam_rows:
        f.write("\t".join(row) + "\n")

rows = []
for row in fam_rows:
    fid, iid = row[0], row[1]
    rows.append((fid, iid, pheno_map[(fid, iid)]))

with dst.open("w", encoding="utf-8") as f:
    for idx, (fid, iid, y) in enumerate(rows):
        f.write(f"{fid}\t{iid}\t{y:.6f}\n")
        shifted = y + (0.15 if idx % 2 == 0 else -0.10)
        f.write(f"{fid}\t{iid}\t{shifted:.6f}\n")
PY

"${BIN}" \
  --bfile "${SUBSET_PREFIX}" \
  --pheno "${LONG_PHENO}" \
  --mpheno 1 \
  --repeatability-model \
  --repeatability-mixed-precision \
  --repeatability-grm-cache "${GRM_CACHE_PREFIX}" \
  --gwas \
  --no-loco \
  --out "${OUT_PREFIX_EXACT}"

"${BIN}" \
  --bfile "${SUBSET_PREFIX}" \
  --pheno "${LONG_PHENO}" \
  --mpheno 1 \
  --repeatability-model \
  --repeatability-mixed-precision \
  --repeatability-grm-cache "${GRM_CACHE_PREFIX}" \
  --gwas \
  --no-loco \
  --out "${OUT_PREFIX_EXACT_CACHED}"

"${BIN}" \
  --bfile "${SUBSET_PREFIX}" \
  --pheno "${LONG_PHENO}" \
  --mpheno 1 \
  --repeatability-model \
  --repeatability-mixed-precision \
  --repeatability-grm-cache "${GRM_CACHE_PREFIX}" \
  --gwas \
  --out "${OUT_PREFIX_LOCO}"

"${BIN}" \
  --bfile "${SUBSET_PREFIX}" \
  --pheno "${LONG_PHENO}" \
  --mpheno 1 \
  --repeatability-model \
  --repeatability-mixed-precision \
  --repeatability-grm-cache "${GRM_CACHE_PREFIX}" \
  --gwas \
  --out "${OUT_PREFIX_LOCO_CACHED}"

"${BIN}" \
  --bfile "${SUBSET_PREFIX}" \
  --pheno "${LONG_PHENO}" \
  --mpheno 1 \
  --repeatability-model \
  --gwas \
  --fast-lmm \
  --no-loco \
  --out "${OUT_PREFIX_FAST}"

python3 - <<'PY'
import csv
from pathlib import Path

def check_run(base_str, expected_mode, extra_required=None):
    base = Path(base_str)
    assoc = base.with_suffix(".assoc.txt")
    nullf = Path(str(base) + ".repeatability.null.txt")
    if not assoc.exists():
        raise SystemExit(f"Missing assoc output file: {assoc}")
    if not nullf.exists():
        raise SystemExit(f"Missing repeatability null summary: {nullf}")
    with assoc.open("r", encoding="utf-8") as f:
        assoc_rows = [line.strip() for line in f if line.strip()]
    if len(assoc_rows) < 2:
        raise SystemExit(f"Association output does not contain any result rows: {assoc}")

    summary = {}
    with nullf.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            summary[row["Item"]] = row["Value"]

    required = {"n_records", "n_subjects", "Vg", "Vpe", "Ve", "h2", "r2", "scan_mode"}
    if extra_required:
        required.update(extra_required)
    missing = required.difference(summary)
    if missing:
        raise SystemExit(f"Missing keys in null summary {nullf}: {sorted(missing)}")
    if int(float(summary["n_records"])) != 40:
        raise SystemExit(f"Unexpected record count: {summary['n_records']}")
    if int(float(summary["n_subjects"])) != 20:
        raise SystemExit(f"Unexpected subject count: {summary['n_subjects']}")
    if summary["scan_mode"] != expected_mode:
        raise SystemExit(f"Unexpected scan mode: {summary['scan_mode']}")
    for key in ("Vg", "Vpe", "Ve", "h2", "r2"):
        float(summary[key])
    return summary

def read_cache_summary(base_str):
    cachef = Path(str(Path(base_str)) + ".repeatability.cache.txt")
    if not cachef.exists():
        raise SystemExit(f"Missing cache summary file: {cachef}")
    rows = {}
    with cachef.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows[row["Label"]] = row["Status"]
    return rows

check_run(
    "/disk228/huangjie/PanOmics/CosmicScanner/tests/results/repeatability_mouse_clean/repeatability_mouse_clean.exact",
    "exact_repeatability_global",
    {"kernel_precision", "grm_cache_prefix", "grm_cache_status"},
)
cached_summary = check_run(
    "/disk228/huangjie/PanOmics/CosmicScanner/tests/results/repeatability_mouse_clean/repeatability_mouse_clean.exact_cached",
    "exact_repeatability_global",
    {"kernel_precision", "grm_cache_prefix", "grm_cache_status"},
)
fast_summary = check_run(
    "/disk228/huangjie/PanOmics/CosmicScanner/tests/results/repeatability_mouse_clean/repeatability_mouse_clean.fast",
    "fast_repeatability_mean_model",
    {"mean_records", "min_records", "max_records"},
)
loco_summary = check_run(
    "/disk228/huangjie/PanOmics/CosmicScanner/tests/results/repeatability_mouse_clean/repeatability_mouse_clean.loco",
    "exact_repeatability_loco",
    {"kernel_precision", "grm_cache_prefix", "grm_cache_status"},
)
loco_cached_summary = check_run(
    "/disk228/huangjie/PanOmics/CosmicScanner/tests/results/repeatability_mouse_clean/repeatability_mouse_clean.loco_cached",
    "exact_repeatability_loco",
    {"kernel_precision", "grm_cache_prefix", "grm_cache_status"},
)
agg = Path("/disk228/huangjie/PanOmics/CosmicScanner/tests/results/repeatability_mouse_clean/repeatability_mouse_clean.fast.agg.pheno")
if not agg.exists():
    raise SystemExit("Missing aggregated phenotype output for fast repeatability mode.")
if float(fast_summary["mean_records"]) <= 1.0:
    raise SystemExit("Fast repeatability summary reports invalid mean_records.")
if cached_summary["kernel_precision"] != "float32_mixed":
    raise SystemExit(f"Unexpected kernel precision: {cached_summary['kernel_precision']}")
if cached_summary["grm_cache_prefix"] == "NA":
    raise SystemExit("Expected grm_cache_prefix in cached exact mode.")
if check_run(
    "/disk228/huangjie/PanOmics/CosmicScanner/tests/results/repeatability_mouse_clean/repeatability_mouse_clean.exact",
    "exact_repeatability_global",
    {"kernel_precision", "grm_cache_prefix", "grm_cache_status"},
)["grm_cache_status"] != "computed_save":
    raise SystemExit("Expected first exact run to compute and save GRM cache.")
if cached_summary["grm_cache_status"] != "mmap_load":
    raise SystemExit(f"Expected cached exact run to use mmap_load, got {cached_summary['grm_cache_status']}")
if loco_summary["grm_cache_status"] not in {"mmap_load", "computed_save"}:
    raise SystemExit(f"Unexpected LOCO global cache status: {loco_summary['grm_cache_status']}")
if loco_cached_summary["grm_cache_status"] != "mmap_load":
    raise SystemExit(f"Expected cached LOCO run to use mmap_load, got {loco_cached_summary['grm_cache_status']}")

loco_cache_rows = read_cache_summary("/disk228/huangjie/PanOmics/CosmicScanner/tests/results/repeatability_mouse_clean/repeatability_mouse_clean.loco")
loco_cached_rows = read_cache_summary("/disk228/huangjie/PanOmics/CosmicScanner/tests/results/repeatability_mouse_clean/repeatability_mouse_clean.loco_cached")
for label in ("global", "chr_1", "chr_2"):
    if label not in loco_cache_rows:
        raise SystemExit(f"Missing LOCO cache entry: {label}")
    if label not in loco_cached_rows:
        raise SystemExit(f"Missing cached LOCO cache entry: {label}")
if loco_cached_rows["chr_1"] != "mmap_load" or loco_cached_rows["chr_2"] != "mmap_load":
    raise SystemExit(f"Expected chr caches to load via mmap, got: {loco_cached_rows}")

cache_bin = Path("/disk228/huangjie/PanOmics/CosmicScanner/tests/data/repeatability_mouse_clean/repeatability_cache/grm.global.f32grm.bin")
cache_meta = Path("/disk228/huangjie/PanOmics/CosmicScanner/tests/data/repeatability_mouse_clean/repeatability_cache/grm.global.meta.txt")
if not cache_bin.exists() or not cache_meta.exists():
    raise SystemExit("Missing GRM cache files for repeatability exact mode.")
for chrom in ("1", "2"):
    chr_bin = Path(f"/disk228/huangjie/PanOmics/CosmicScanner/tests/data/repeatability_mouse_clean/repeatability_cache/grm.chr_{chrom}.f32grm.bin")
    chr_meta = Path(f"/disk228/huangjie/PanOmics/CosmicScanner/tests/data/repeatability_mouse_clean/repeatability_cache/grm.chr_{chrom}.meta.txt")
    if not chr_bin.exists() or not chr_meta.exists():
        raise SystemExit(f"Missing chromosome cache files for chr{chrom}.")

print("Repeatability GWAS exact/cached/loco/loco_cached/fast smoke tests passed.")
PY

echo "Done. 重复力 GWAS exact/cached/loco/loco_cached/fast smoke tests 完成。"
