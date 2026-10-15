#!/bin/bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="${COSMICSCANNER_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BIN="${ROOT}/build/cosmicscanner"
DATA_DIR="${ROOT}/tests/data/rrm_dynamic_simulation"
RESULT_DIR="${ROOT}/tests/results/rrm_dynamic_simulation"
PREFIX="${DATA_DIR}/rrm_dynamic_sim"
LONG_PHENO="${DATA_DIR}/rrm_dynamic_sim.long.pheno"
OUT_PREFIX="${RESULT_DIR}/rrm_dynamic_sim"
MEAN_PHENO="${DATA_DIR}/rrm_dynamic_sim.mean.pheno"
T0_PHENO="${DATA_DIR}/rrm_dynamic_sim.t0.pheno"
T1_PHENO="${DATA_DIR}/rrm_dynamic_sim.t1.pheno"
WIDE_PHENO="${DATA_DIR}/rrm_dynamic_sim.wide.pheno"
RUN_LOCO_SMOKE=0
RUN_ORDER_SWEEP=0
RUN_PGEN_SMOKE=0
RUN_BGEN_SMOKE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --loco-smoke) RUN_LOCO_SMOKE=1 ;;
    --order-sweep) RUN_ORDER_SWEEP=1 ;;
    --pgen-smoke) RUN_PGEN_SMOKE=1 ;;
    --bgen-smoke) RUN_BGEN_SMOKE=1 ;;
    -h|--help)
      echo "Usage: tests/run_rrm_dynamic_simulation.sh [--loco-smoke] [--order-sweep] [--pgen-smoke] [--bgen-smoke]"
      exit 0
      ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

mkdir -p "${DATA_DIR}" "${RESULT_DIR}"

COSMICSCANNER_ROOT="${ROOT}" python3 - <<'PY'
import os
import math
import random
from pathlib import Path

random.seed(20260610)

root = Path(os.environ["COSMICSCANNER_ROOT"])
data_dir = root / "tests" / "data" / "rrm_dynamic_simulation"
prefix = data_dir / "rrm_dynamic_sim"
long_pheno = data_dir / "rrm_dynamic_sim.long.pheno"
mean_pheno = data_dir / "rrm_dynamic_sim.mean.pheno"
t0_pheno = data_dir / "rrm_dynamic_sim.t0.pheno"
t1_pheno = data_dir / "rrm_dynamic_sim.t1.pheno"
wide_pheno = data_dir / "rrm_dynamic_sim.wide.pheno"

n = 140
m = 360
times = [0.0, 0.25, 0.5, 0.75, 1.0]
causal_dynamic = 12
causal_stable = 37
causal_early = 82
causal_late = 126
causal_transient = 171

genotypes = []
for j in range(m):
    maf = 0.18 + 0.32 * ((j * 37) % 101) / 100.0
    row = []
    for i in range(n):
        u = random.random()
        p0 = (1.0 - maf) ** 2
        p1 = p0 + 2.0 * maf * (1.0 - maf)
        if u < p0:
            g = 0
        elif u < p1:
            g = 1
        else:
            g = 2
        row.append(g)
    genotypes.append(row)

def mean_center(values):
    mu = sum(values) / len(values)
    return [v - mu for v in values]

g_dyn = mean_center(genotypes[causal_dynamic])
g_stable = mean_center(genotypes[causal_stable])
g_early = mean_center(genotypes[causal_early])
g_late = mean_center(genotypes[causal_late])
g_transient = mean_center(genotypes[causal_transient])
subject_effect = [random.gauss(0.0, 0.12) for _ in range(n)]

with (prefix.with_suffix(".fam")).open("w", encoding="utf-8") as f:
    for i in range(n):
        fid = f"F{i+1:04d}"
        iid = f"I{i+1:04d}"
        f.write(f"{fid}\t{iid}\t0\t0\t0\t-9\n")

with (prefix.with_suffix(".bim")).open("w", encoding="utf-8") as f:
    for j in range(m):
        chrom = min(3, 1 + j // (m // 3))
        if j == causal_dynamic:
            snp = "rs_dyn_causal"
        elif j == causal_stable:
            snp = "rs_stable_causal"
        elif j == causal_early:
            snp = "rs_early_causal"
        elif j == causal_late:
            snp = "rs_late_causal"
        elif j == causal_transient:
            snp = "rs_transient_causal"
        else:
            snp = f"rs_bg_{j+1:04d}"
        f.write(f"{chrom}\t{snp}\t0\t{1000 + j * 10}\tA\tG\n")

bed = bytearray(b"\x6c\x1b\x01")
code = {0: 0b00, 1: 0b10, 2: 0b11}
bytes_per_snp = (n + 3) // 4
for j in range(m):
    block = bytearray(bytes_per_snp)
    for i, g in enumerate(genotypes[j]):
        block[i // 4] |= code[g] << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix(".bed").write_bytes(bed)

phenotypes = []
with long_pheno.open("w", encoding="utf-8") as f:
    f.write("FID IID TIME Y\n")
    for i in range(n):
        fid = f"F{i+1:04d}"
        iid = f"I{i+1:04d}"
        values = []
        for t in times:
            centered_t = 2.0 * t - 1.0
            transient_shape = max(0.0, 1.0 - 4.0 * (t - 0.5) * (t - 0.5))
            y = (
                1.0
                + subject_effect[i]
                + 0.35 * centered_t
                + 0.75 * g_stable[i]
                + 1.35 * g_dyn[i] * centered_t
                + 1.10 * g_early[i] * (1.0 - t) * (1.0 - t)
                + 1.10 * g_late[i] * t * t
                + 1.20 * g_transient[i] * transient_shape
                + random.gauss(0.0, 0.08)
            )
            values.append(y)
            f.write(f"{fid}\t{iid}\t{t:.2f}\t{y:.8f}\n")
        phenotypes.append((fid, iid, values))

with mean_pheno.open("w", encoding="utf-8") as f_mean, \
     t0_pheno.open("w", encoding="utf-8") as f_t0, \
     t1_pheno.open("w", encoding="utf-8") as f_t1, \
     wide_pheno.open("w", encoding="utf-8") as f_wide:
    f_mean.write("FID IID Y\n")
    f_t0.write("FID IID Y\n")
    f_t1.write("FID IID Y\n")
    f_wide.write("FID IID Y_T0 Y_T025 Y_T05 Y_T075 Y_T1\n")
    for fid, iid, values in phenotypes:
        f_mean.write(f"{fid}\t{iid}\t{sum(values) / len(values):.8f}\n")
        f_t0.write(f"{fid}\t{iid}\t{values[0]:.8f}\n")
        f_t1.write(f"{fid}\t{iid}\t{values[-1]:.8f}\n")
        f_wide.write(f"{fid}\t{iid}\t" + "\t".join(f"{v:.8f}" for v in values) + "\n")

with (data_dir / "causal_snps.txt").open("w", encoding="utf-8") as f:
    f.write("rs_dyn_causal\tsign_switch_slope\n")
    f.write("rs_stable_causal\tstable\n")
    f.write("rs_early_causal\tearly_specific\n")
    f.write("rs_late_causal\tlate_specific\n")
    f.write("rs_transient_causal\ttransient\n")

print(f"Generated RRM dynamic simulation: n={n}, m={m}, observations={n * len(times)}")
PY

"${BIN}" \
  --bfile "${PREFIX}" \
  --long-pheno "${LONG_PHENO}" \
  --tpheno 1 \
  --mpheno 2 \
  --rrm-order 2 \
  --rrm-curve-grid 5 \
  --rrm-top 30 \
  --rrm-sig-threshold 0.01 \
  --maf 0 \
  --gwas \
  --no-loco \
  --out "${OUT_PREFIX}"

"${BIN}" \
  --bfile "${PREFIX}" \
  --pheno "${MEAN_PHENO}" \
  --mpheno 1 \
  --maf 0 \
  --gwas \
  --no-loco \
  --out "${RESULT_DIR}/mean_gwas"

"${BIN}" \
  --bfile "${PREFIX}" \
  --pheno "${T0_PHENO}" \
  --mpheno 1 \
  --maf 0 \
  --gwas \
  --no-loco \
  --out "${RESULT_DIR}/time0_gwas"

"${BIN}" \
  --bfile "${PREFIX}" \
  --pheno "${T1_PHENO}" \
  --mpheno 1 \
  --maf 0 \
  --gwas \
  --no-loco \
  --out "${RESULT_DIR}/time1_gwas"

"${BIN}" \
  --bfile "${PREFIX}" \
  --pheno "${WIDE_PHENO}" \
  --mpheno-list 1,2,3,4,5 \
  --multivariate \
  --mv-method emmv \
  --maf 0 \
  --gwas \
  --no-loco \
  --out "${RESULT_DIR}/multivariate_gwas"

COSMICSCANNER_ROOT="${ROOT}" python3 - <<'PY'
import os
import csv
from pathlib import Path

base = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation/rrm_dynamic_sim")
summary_path = Path(str(base) + ".rrm.dynamic_summary.txt")
curve_path = Path(str(base) + ".rrm.dynamic_curve.txt")
joint_path = Path(str(base) + ".rrm_joint.assoc.txt")
mean_path = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation/mean_gwas.assoc.txt")
t0_path = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation/time0_gwas.assoc.txt")
t1_path = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation/time1_gwas.assoc.txt")
mv_path = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation/multivariate_gwas.assoc.mv.txt")
comparison_path = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation/comparison.tsv")
plot_dir = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation/figures")

summary = list(csv.DictReader(summary_path.open("r", encoding="utf-8"), delimiter="\t"))
curve = list(csv.DictReader(curve_path.open("r", encoding="utf-8"), delimiter="\t"))
joint = list(csv.DictReader(joint_path.open("r", encoding="utf-8"), delimiter="\t"))

def assoc_rank(path, snp):
    rows = list(csv.DictReader(path.open("r", encoding="utf-8"), delimiter="\t"))
    rows = [row for row in rows if row.get("SNP") and row.get("P") not in (None, "", "nan")]
    rows.sort(key=lambda row: float(row["P"]))
    for rank, row in enumerate(rows, 1):
        if row["SNP"] == snp:
            return rank, float(row["P"])
    return None, float("nan")

def mv_rank(path, snp):
    rows = list(csv.DictReader(path.open("r", encoding="utf-8"), delimiter="\t"))
    rows = [row for row in rows if row.get("SNP") and row.get("P") not in (None, "", "nan")]
    rows.sort(key=lambda row: float(row["P"]))
    for rank, row in enumerate(rows, 1):
        if row["SNP"] == snp:
            return rank, float(row["P"])
    return None, float("nan")

def summary_rank(snp):
    for rank, row in enumerate(summary, 1):
        if row["SNP"] == snp:
            return rank, float(row["JointP"])
    return None, float("nan")

if len(summary) != 30:
    raise SystemExit(f"Expected 30 dynamic summary rows, got {len(summary)}")
if len(curve) != 150:
    raise SystemExit(f"Expected 150 dynamic curve rows, got {len(curve)}")
if not joint:
    raise SystemExit("Missing joint association results")

top10 = {row["SNP"] for row in summary[:10]}
if "rs_dyn_causal" not in top10:
    raise SystemExit(f"Dynamic causal SNP not recovered in top10: {sorted(top10)}")
top30 = {row["SNP"] for row in summary}
expected_causal = {"rs_dyn_causal", "rs_stable_causal", "rs_early_causal", "rs_late_causal", "rs_transient_causal"}
missing_causal = expected_causal - top30
if missing_causal:
    raise SystemExit(f"Causal SNPs missing from dynamic top30: {sorted(missing_causal)}")

dyn_row = next(row for row in summary if row["SNP"] == "rs_dyn_causal")
if dyn_row["Pattern"] == "stable":
    raise SystemExit(f"Dynamic causal SNP classified as stable: {dyn_row}")
if abs(float(dyn_row["SigThreshold"]) - 0.01) > 1e-12:
    raise SystemExit(f"Unexpected dynamic SigThreshold: {dyn_row['SigThreshold']}")
if int(dyn_row["GridPoints"]) != 5:
    raise SystemExit(f"Unexpected grid count: {dyn_row['GridPoints']}")

patterns = {row["SNP"]: row["Pattern"] for row in summary if row["SNP"] in expected_causal}
if patterns.get("rs_dyn_causal") != "sign_switch":
    raise SystemExit(f"Expected sign-switch dynamic causal pattern, got: {patterns.get('rs_dyn_causal')}")
if patterns.get("rs_stable_causal") != "stable":
    raise SystemExit(f"Expected stable causal pattern, got: {patterns.get('rs_stable_causal')}")
if patterns.get("rs_early_causal") not in {"early_specific", "decreasing", "stage_specific"}:
    raise SystemExit(f"Unexpected early causal pattern: {patterns.get('rs_early_causal')}")
if patterns.get("rs_late_causal") not in {"late_specific", "increasing", "stage_specific"}:
    raise SystemExit(f"Unexpected late causal pattern: {patterns.get('rs_late_causal')}")
if patterns.get("rs_transient_causal") not in {"transient", "nonlinear_dynamic", "stage_specific"}:
    raise SystemExit(f"Unexpected transient causal pattern: {patterns.get('rs_transient_causal')}")

times = sorted({float(row["Time"]) for row in curve if row["SNP"] == "rs_dyn_causal"})
if times != [0.0, 0.25, 0.5, 0.75, 1.0]:
    raise SystemExit(f"Unexpected dynamic causal curve times: {times}")

rrm_dyn_rank, rrm_dyn_p = summary_rank("rs_dyn_causal")
rrm_stable_rank, rrm_stable_p = summary_rank("rs_stable_causal")
rrm_early_rank, rrm_early_p = summary_rank("rs_early_causal")
rrm_late_rank, rrm_late_p = summary_rank("rs_late_causal")
rrm_transient_rank, rrm_transient_p = summary_rank("rs_transient_causal")
mean_dyn_rank, mean_dyn_p = assoc_rank(mean_path, "rs_dyn_causal")
mean_stable_rank, mean_stable_p = assoc_rank(mean_path, "rs_stable_causal")
t0_dyn_rank, t0_dyn_p = assoc_rank(t0_path, "rs_dyn_causal")
t1_dyn_rank, t1_dyn_p = assoc_rank(t1_path, "rs_dyn_causal")
mv_dyn_rank, mv_dyn_p = mv_rank(mv_path, "rs_dyn_causal")
mv_stable_rank, mv_stable_p = mv_rank(mv_path, "rs_stable_causal")
mv_early_rank, mv_early_p = mv_rank(mv_path, "rs_early_causal")
mv_late_rank, mv_late_p = mv_rank(mv_path, "rs_late_causal")
mv_transient_rank, mv_transient_p = mv_rank(mv_path, "rs_transient_causal")

if mean_dyn_rank is None or mean_stable_rank is None:
    raise SystemExit("Missing causal SNPs in mean GWAS results")
if rrm_dyn_rank > mean_dyn_rank:
    raise SystemExit(f"RRM did not improve dynamic causal rank: RRM={rrm_dyn_rank}, mean={mean_dyn_rank}")
if mean_stable_rank > 10:
    raise SystemExit(f"Stable causal SNP not recovered by mean GWAS top10: rank={mean_stable_rank}")
if min(t0_dyn_rank or 999999, t1_dyn_rank or 999999) > 10:
    raise SystemExit(f"Dynamic causal SNP not recovered by endpoint GWAS top10: t0={t0_dyn_rank}, t1={t1_dyn_rank}")
if mv_stable_rank is None or mv_stable_rank > 10:
    raise SystemExit(f"Stable causal SNP not recovered by multivariate GWAS top10: rank={mv_stable_rank}")

def fmt_rank(value):
    return "NA" if value is None else str(value)

def fmt_p(value):
    return "NA" if value != value else str(value)

with comparison_path.open("w", encoding="utf-8") as f:
    f.write("Method\tSNP\tRank\tP\tNote\n")
    f.write(f"RRM_dynamic\trs_dyn_causal\t{rrm_dyn_rank}\t{rrm_dyn_p}\tjoint_dynamic_effect\n")
    f.write(f"RRM_dynamic\trs_stable_causal\t{rrm_stable_rank}\t{rrm_stable_p}\tjoint_dynamic_effect\n")
    f.write(f"RRM_dynamic\trs_early_causal\t{rrm_early_rank}\t{rrm_early_p}\tearly_specific_effect\n")
    f.write(f"RRM_dynamic\trs_late_causal\t{rrm_late_rank}\t{rrm_late_p}\tlate_specific_effect\n")
    f.write(f"RRM_dynamic\trs_transient_causal\t{rrm_transient_rank}\t{rrm_transient_p}\ttransient_effect\n")
    f.write(f"Mean_GWAS\trs_dyn_causal\t{mean_dyn_rank}\t{mean_dyn_p}\tmean_effect_can_cancel_dynamic_signal\n")
    f.write(f"Mean_GWAS\trs_stable_causal\t{mean_stable_rank}\t{mean_stable_p}\tstable_effect_expected\n")
    f.write(f"Time0_GWAS\trs_dyn_causal\t{t0_dyn_rank}\t{t0_dyn_p}\tendpoint_specific_effect\n")
    f.write(f"Time1_GWAS\trs_dyn_causal\t{t1_dyn_rank}\t{t1_dyn_p}\tendpoint_specific_effect\n")
    f.write(f"MvGWAS\trs_dyn_causal\t{fmt_rank(mv_dyn_rank)}\t{fmt_p(mv_dyn_p)}\tdiscrete_multi_time_joint_test_no_curve\n")
    f.write(f"MvGWAS\trs_stable_causal\t{fmt_rank(mv_stable_rank)}\t{fmt_p(mv_stable_p)}\tdiscrete_multi_time_joint_test_no_curve\n")
    f.write(f"MvGWAS\trs_early_causal\t{fmt_rank(mv_early_rank)}\t{fmt_p(mv_early_p)}\tdiscrete_multi_time_joint_test_no_curve\n")
    f.write(f"MvGWAS\trs_late_causal\t{fmt_rank(mv_late_rank)}\t{fmt_p(mv_late_p)}\tdiscrete_multi_time_joint_test_no_curve\n")
    f.write(f"MvGWAS\trs_transient_causal\t{fmt_rank(mv_transient_rank)}\t{fmt_p(mv_transient_p)}\tdiscrete_multi_time_joint_test_no_curve\n")

print("RRM dynamic simulation checks passed.")
print(f"Dynamic causal summary: rank={dyn_row['Rank']} pattern={dyn_row['Pattern']} best_time={dyn_row['BestTime']} joint_p={dyn_row['JointP']}")
print(f"Comparison ranks: RRM dynamic={rrm_dyn_rank}, mean={mean_dyn_rank}, time0={t0_dyn_rank}, time1={t1_dyn_rank}")
print(f"Multivariate ranks: dynamic={mv_dyn_rank}, stable={mv_stable_rank}, early={mv_early_rank}, late={mv_late_rank}, transient={mv_transient_rank}")
print(f"Causal trajectory patterns: {patterns}")
PY

python3 "${ROOT}/tools/plot_rrm_dynamic.py" \
  --curve "${OUT_PREFIX}.rrm.dynamic_curve.txt" \
  --summary "${OUT_PREFIX}.rrm.dynamic_summary.txt" \
  --comparison "${RESULT_DIR}/comparison.tsv" \
  --top 5 \
  --outdir "${RESULT_DIR}/figures"

COSMICSCANNER_ROOT="${ROOT}" python3 - <<'PY'
import os
from pathlib import Path

fig_dir = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation/figures")
expected = [
    "rrm_dynamic_effect_curves.png",
    "rrm_dynamic_effect_heatmap.png",
    "rrm_dynamic_method_comparison.png",
    "rrm_dynamic_plot_summary.tsv",
]
missing = [name for name in expected if not (fig_dir / name).exists()]
if missing:
    raise SystemExit(f"Missing RRM dynamic figures: {missing}")
print("RRM dynamic plotting checks passed.")
PY

if [[ ${RUN_LOCO_SMOKE} -eq 1 ]]; then
  "${BIN}" \
    --bfile "${PREFIX}" \
    --long-pheno "${LONG_PHENO}" \
    --tpheno 1 \
    --mpheno 2 \
    --rrm-order 1 \
    --rrm-curve-grid 5 \
    --rrm-top 10 \
    --rrm-sig-threshold 0.01 \
    --maf 0 \
    --gwas \
    --out "${RESULT_DIR}/rrm_dynamic_loco"

  COSMICSCANNER_ROOT="${ROOT}" python3 - <<'PY'
import os
import csv
from pathlib import Path

base = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation/rrm_dynamic_loco")
null_rows = list(csv.DictReader(Path(str(base) + ".rrm_joint.null.txt").open("r", encoding="utf-8"), delimiter="\t"))
summary_rows = list(csv.DictReader(Path(str(base) + ".rrm.dynamic_summary.txt").open("r", encoding="utf-8"), delimiter="\t"))
curve_rows = list(csv.DictReader(Path(str(base) + ".rrm.dynamic_curve.txt").open("r", encoding="utf-8"), delimiter="\t"))
scopes = {row["Scope"] for row in null_rows}
expected = {"LOCO_1", "LOCO_2", "LOCO_3"}
if not expected.issubset(scopes):
    raise SystemExit(f"Missing LOCO null scopes: expected {expected}, got {scopes}")
if not summary_rows:
    raise SystemExit("Missing LOCO dynamic summary rows")
if not curve_rows:
    raise SystemExit("Missing LOCO dynamic curve rows")
print("RRM dynamic LOCO smoke checks passed.")
PY
fi

if [[ ${RUN_PGEN_SMOKE} -eq 1 ]]; then
  if ! command -v plink2 >/dev/null 2>&1; then
    echo "plink2 is required for --pgen-smoke" >&2
    exit 1
  fi

  PGEN_PREFIX="${DATA_DIR}/rrm_dynamic_sim_pgen"
  plink2 --bfile "${PREFIX}" --make-pgen --out "${PGEN_PREFIX}" >/dev/null

  "${BIN}" \
    --pfile "${PGEN_PREFIX}" \
    --long-pheno "${LONG_PHENO}" \
    --tpheno 1 \
    --mpheno 2 \
    --rrm-order 1 \
    --rrm-curve-grid 5 \
    --rrm-top 10 \
    --rrm-sig-threshold 0.01 \
    --maf 0 \
    --gwas \
    --no-loco \
    --out "${RESULT_DIR}/rrm_dynamic_pgen"

  COSMICSCANNER_ROOT="${ROOT}" python3 - <<'PY'
import os
import csv
from pathlib import Path

base = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation/rrm_dynamic_pgen")
summary = list(csv.DictReader(Path(str(base) + ".rrm.dynamic_summary.txt").open("r", encoding="utf-8"), delimiter="\t"))
curve = list(csv.DictReader(Path(str(base) + ".rrm.dynamic_curve.txt").open("r", encoding="utf-8"), delimiter="\t"))
joint = list(csv.DictReader(Path(str(base) + ".rrm_joint.assoc.txt").open("r", encoding="utf-8"), delimiter="\t"))
if len(summary) != 10:
    raise SystemExit(f"PGEN smoke: expected 10 summary rows, got {len(summary)}")
if len(curve) != 50:
    raise SystemExit(f"PGEN smoke: expected 50 curve rows, got {len(curve)}")
if not joint:
    raise SystemExit("PGEN smoke: missing joint association rows")
dyn = next((row for row in summary if row["SNP"] == "rs_dyn_causal"), None)
if dyn is None:
    raise SystemExit("PGEN smoke: rs_dyn_causal missing from top10 summary")
if dyn["Pattern"] == "stable":
    raise SystemExit("PGEN smoke: rs_dyn_causal classified as stable")
print("RRM dynamic PGEN smoke checks passed.")
PY
fi

if [[ ${RUN_BGEN_SMOKE} -eq 1 ]]; then
  if ! command -v plink2 >/dev/null 2>&1; then
    echo "plink2 is required for --bgen-smoke" >&2
    exit 1
  fi

  BGEN_PREFIX="${DATA_DIR}/rrm_dynamic_sim_bgen"
  plink2 --bfile "${PREFIX}" --export bgen-1.2 --out "${BGEN_PREFIX}" >/dev/null

  "${BIN}" \
    --bgen "${BGEN_PREFIX}.bgen" \
    --sample "${BGEN_PREFIX}.sample" \
    --long-pheno "${LONG_PHENO}" \
    --tpheno 1 \
    --mpheno 2 \
    --rrm-order 1 \
    --rrm-curve-grid 5 \
    --rrm-top 10 \
    --rrm-sig-threshold 0.01 \
    --maf 0 \
    --gwas \
    --no-loco \
    --out "${RESULT_DIR}/rrm_dynamic_bgen"

  COSMICSCANNER_ROOT="${ROOT}" python3 - <<'PY'
import os
import csv
from pathlib import Path

base = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation/rrm_dynamic_bgen")
summary = list(csv.DictReader(Path(str(base) + ".rrm.dynamic_summary.txt").open("r", encoding="utf-8"), delimiter="\t"))
curve = list(csv.DictReader(Path(str(base) + ".rrm.dynamic_curve.txt").open("r", encoding="utf-8"), delimiter="\t"))
joint = list(csv.DictReader(Path(str(base) + ".rrm_joint.assoc.txt").open("r", encoding="utf-8"), delimiter="\t"))
if len(summary) != 10:
    raise SystemExit(f"BGEN smoke: expected 10 summary rows, got {len(summary)}")
if len(curve) != 50:
    raise SystemExit(f"BGEN smoke: expected 50 curve rows, got {len(curve)}")
if not joint:
    raise SystemExit("BGEN smoke: missing joint association rows")
dyn = next((row for row in summary if row["SNP"] == "rs_dyn_causal"), None)
if dyn is None:
    raise SystemExit("BGEN smoke: rs_dyn_causal missing from top10 summary")
if dyn["Pattern"] == "stable":
    raise SystemExit("BGEN smoke: rs_dyn_causal classified as stable")
print("RRM dynamic BGEN smoke checks passed.")
PY
fi

if [[ ${RUN_ORDER_SWEEP} -eq 1 ]]; then
  for order in 1 2 3; do
    "${BIN}" \
      --bfile "${PREFIX}" \
      --long-pheno "${LONG_PHENO}" \
      --tpheno 1 \
      --mpheno 2 \
      --rrm-order "${order}" \
      --rrm-curve-grid 5 \
      --rrm-top 30 \
      --rrm-sig-threshold 0.01 \
      --maf 0 \
      --gwas \
      --no-loco \
      --out "${RESULT_DIR}/rrm_dynamic_order${order}"
  done

  COSMICSCANNER_ROOT="${ROOT}" python3 - <<'PY'
import os
import csv
from pathlib import Path

result_dir = Path(os.environ["COSMICSCANNER_ROOT"] + "/tests/results/rrm_dynamic_simulation")
causal_snps = [
    "rs_dyn_causal",
    "rs_stable_causal",
    "rs_early_causal",
    "rs_late_causal",
    "rs_transient_causal",
]
out_path = result_dir / "rrm_order_sweep.tsv"

with out_path.open("w", encoding="utf-8") as out:
    out.write("Order\tSNP\tRank\tJointP\tPattern\tGridPoints\tSummaryRows\tCurveRows\n")
    for order in (1, 2, 3):
        base = result_dir / f"rrm_dynamic_order{order}"
        summary_path = Path(str(base) + ".rrm.dynamic_summary.txt")
        curve_path = Path(str(base) + ".rrm.dynamic_curve.txt")
        summary = list(csv.DictReader(summary_path.open("r", encoding="utf-8"), delimiter="\t"))
        curve = list(csv.DictReader(curve_path.open("r", encoding="utf-8"), delimiter="\t"))
        if len(summary) != 30:
            raise SystemExit(f"order {order}: expected 30 summary rows, got {len(summary)}")
        if len(curve) != 150:
            raise SystemExit(f"order {order}: expected 150 curve rows, got {len(curve)}")
        rows_by_snp = {row["SNP"]: row for row in summary}
        dyn = rows_by_snp.get("rs_dyn_causal")
        if dyn is None:
            raise SystemExit(f"order {order}: rs_dyn_causal missing from top30 summary")
        if int(dyn["Rank"]) > 10:
            raise SystemExit(f"order {order}: rs_dyn_causal rank > 10: {dyn['Rank']}")
        if dyn["Pattern"] == "stable":
            raise SystemExit(f"order {order}: rs_dyn_causal classified as stable")
        for snp in causal_snps:
            row = rows_by_snp.get(snp)
            if row is None:
                out.write(f"{order}\t{snp}\tNA\tNA\tNA\tNA\t{len(summary)}\t{len(curve)}\n")
            else:
                out.write(
                    f"{order}\t{snp}\t{row['Rank']}\t{row['JointP']}\t{row['Pattern']}\t"
                    f"{row['GridPoints']}\t{len(summary)}\t{len(curve)}\n"
                )

print(f"RRM order sweep checks passed: {out_path}")
PY
fi

echo "Done. RRM dynamic simulation validation completed."
