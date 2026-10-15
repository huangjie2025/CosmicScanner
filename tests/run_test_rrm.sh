#!/bin/bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export COSMICSCANNER_ROOT="${ROOT}"
RESULT_DIR="${ROOT}/tests/results/rrm_mouse_clean"
DATA_DIR="${ROOT}/tests/data/rrm_mouse_clean"
OUT_PREFIX="${RESULT_DIR}/rrm_mouse_clean"
BIN="${ROOT}/build/cosmicscanner"
LONG_PHENO="${RESULT_DIR}/mouse_clean.rrm.long.pheno"
SUBSET_PREFIX="${DATA_DIR}/mouse_rrm_subset"

mkdir -p "${RESULT_DIR}"
mkdir -p "${DATA_DIR}"

python3 - <<'PY'
from pathlib import Path

import os
root = Path(os.environ["COSMICSCANNER_ROOT"])
src_prefix = root / "examples" / "mouse_clean" / "mouse_clean"
dst_prefix = root / "tests" / "data" / "rrm_mouse_clean" / "mouse_rrm_subset"
dst = root / "tests" / "results" / "rrm_mouse_clean" / "mouse_clean.rrm.long.pheno"

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
    raise SystemExit("Not enough samples with phenotype for RRM subset.")

bim_rows = []
with (src_prefix.with_suffix(".bim")).open("r", encoding="utf-8") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 6:
            bim_rows.append(parts)
        if len(bim_rows) >= 64:
            break

if len(bim_rows) < 32:
    raise SystemExit("Not enough SNPs for RRM subset.")

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

with dst.open("w", encoding="utf-8") as f:
    f.write("FID IID TIME Y\n")
    for idx, row in enumerate(fam_rows):
        fid, iid = row[0], row[1]
        y = pheno_map[(fid, iid)]
        slope = 0.10 if idx % 2 == 0 else -0.06
        f.write(f"{fid}\t{iid}\t0.0\t{y:.6f}\n")
        f.write(f"{fid}\t{iid}\t1.0\t{(y + slope):.6f}\n")
PY

"${BIN}" \
  --bfile "${SUBSET_PREFIX}" \
  --long-pheno "${LONG_PHENO}" \
  --tpheno 1 \
  --mpheno 2 \
  --rrm-order 1 \
  --rrm-curve-grid 5 \
  --rrm-top 7 \
  --rrm-sig-threshold 0.1 \
  --gwas \
  --no-loco \
  --out "${OUT_PREFIX}"

python3 - <<'PY'
import csv
from pathlib import Path

import os
base = Path(os.environ["COSMICSCANNER_ROOT"]) / "tests/results/rrm_mouse_clean/rrm_mouse_clean"

meta = {}
with Path(str(base) + ".rrm.meta.txt").open("r", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        meta[row["Key"]] = row["Value"]

joint = list(csv.DictReader(Path(str(base) + ".rrm_joint.assoc.txt").open("r", encoding="utf-8"), delimiter="\t"))
null_rows = list(csv.DictReader(Path(str(base) + ".rrm_joint.null.txt").open("r", encoding="utf-8"), delimiter="\t"))
coef_rows = list(csv.DictReader(Path(str(base) + ".rrmcoef.txt").open("r", encoding="utf-8"), delimiter="\t"))
summary_rows = list(csv.DictReader(Path(str(base) + ".rrm.assoc.summary.txt").open("r", encoding="utf-8"), delimiter="\t"))
curve_rows = list(csv.DictReader(Path(str(base) + ".rrm.top_curve.txt").open("r", encoding="utf-8"), delimiter="\t"))
dynamic_curve_rows = list(csv.DictReader(Path(str(base) + ".rrm.dynamic_curve.txt").open("r", encoding="utf-8"), delimiter="\t"))
dynamic_summary_rows = list(csv.DictReader(Path(str(base) + ".rrm.dynamic_summary.txt").open("r", encoding="utf-8"), delimiter="\t"))

expected = {
    "Mode": "GWAS",
    "Order": "1",
    "BasisCount": "2",
    "SubjectCount": "20",
    "ObservationCount": "40",
}
for key, val in expected.items():
    if meta.get(key) != val:
        raise SystemExit(f"Unexpected metadata {key}: {meta.get(key)!r} != {val!r}")

if abs(float(meta["TimeMin"]) - 0.0) > 1e-8 or abs(float(meta["TimeMax"]) - 1.0) > 1e-8:
    raise SystemExit(f"Unexpected time range: {meta['TimeMin']}..{meta['TimeMax']}")

if len(coef_rows) != 20:
    raise SystemExit(f"Unexpected coefficient row count: {len(coef_rows)}")
if not joint:
    raise SystemExit("Missing joint RRM GWAS results")
if len(null_rows) != 1 or null_rows[0]["Scope"] != "GLOBAL":
    raise SystemExit(f"Unexpected null-model rows: {null_rows}")
for key in ("Sigma2E", "LogL", "K_0_0", "K_0_1", "K_1_0", "K_1_1"):
    float(null_rows[0][key])
if not summary_rows:
    raise SystemExit("Missing summary rows")
if not curve_rows:
    raise SystemExit("Missing top-curve rows")
if not dynamic_curve_rows:
    raise SystemExit("Missing dynamic curve rows")
if not dynamic_summary_rows:
    raise SystemExit("Missing dynamic summary rows")
if len(dynamic_summary_rows) != 7:
    raise SystemExit(f"Unexpected dynamic summary row count: {len(dynamic_summary_rows)}")
times = sorted({float(row["Time"]) for row in curve_rows})
if times != [0.0, 1.0]:
    raise SystemExit(f"Unexpected curve times: {times}")
dynamic_times = sorted({float(row["Time"]) for row in dynamic_curve_rows})
if dynamic_times != [0.0, 1.0]:
    raise SystemExit(f"Unexpected dynamic curve times: {dynamic_times}")
if len(dynamic_curve_rows) != 14:
    raise SystemExit(f"Unexpected dynamic curve row count: {len(dynamic_curve_rows)}")
for key in ("BestTime", "BestTimeP", "BetaAtBestTime", "SEAtBestTime", "BetaStart", "BetaEnd", "DeltaBeta", "BetaRange", "MaxAbsBeta", "EffectRatio", "PeakTime"):
    float(dynamic_summary_rows[0][key])
if abs(float(dynamic_summary_rows[0]["SigThreshold"]) - 0.1) > 1e-12:
    raise SystemExit(f"Unexpected SigThreshold: {dynamic_summary_rows[0]['SigThreshold']}")
for key in ("SigWindowStart", "SigWindowEnd"):
    val = dynamic_summary_rows[0][key]
    if val.lower() != "nan":
        float(val)
int(dynamic_summary_rows[0]["SigGridPoints"])
int(dynamic_summary_rows[0]["GridPoints"])
if dynamic_summary_rows[0]["Pattern"] not in {"stable", "increasing", "decreasing", "stage_specific", "nonlinear_dynamic", "early_specific", "late_specific", "transient", "sign_switch"}:
    raise SystemExit(f"Unexpected dynamic pattern: {dynamic_summary_rows[0]['Pattern']}")

print("RRM GWAS demo checks passed.")
PY

echo "Done. GWAS 随机回归 demo 验证完成。"
