#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
N="${N:-120}"
M="${M:-300}"
TIMES="${TIMES:-5}"
REPS="${REPS:-5}"
SEED="${SEED:-20260612}"
THRESHOLD="${THRESHOLD:-0.01}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/tests/results/rrm_method_comparison_n${N}_m${M}_r${REPS}}"
METRICS_TSV="${RESULT_DIR}/method_metrics.tsv"
SUMMARY_TSV="${RESULT_DIR}/summary.tsv"
SUMMARY_MD="${RESULT_DIR}/summary.md"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${RESULT_DIR}/data" "${RESULT_DIR}/logs"
rm -f "${METRICS_TSV}" "${SUMMARY_TSV}" "${SUMMARY_MD}"
printf 'Rep\tMethod\tTarget\tRank\tP\tRuntimeSec\tMaxRSSKB\tCurveRMSE\tPattern\n' > "${METRICS_TSV}"

for rep in $(seq 1 "${REPS}"); do
  REP_DIR="${RESULT_DIR}/data/rep${rep}"
  PREFIX="${REP_DIR}/sim"
  LONG_PHENO="${PREFIX}.long.pheno"
  MEAN_PHENO="${PREFIX}.mean.pheno"
  T0_PHENO="${PREFIX}.t0.pheno"
  T1_PHENO="${PREFIX}.t1.pheno"
  WIDE_PHENO="${PREFIX}.wide.pheno"
  TRUE_CURVE="${PREFIX}.true_curve.tsv"
  mkdir -p "${REP_DIR}" "${RESULT_DIR}/rep${rep}"

  "${PYTHON_BIN}" - "${PREFIX}" "${N}" "${M}" "${TIMES}" "$((SEED + rep))" <<'PY'
import sys
from pathlib import Path
import numpy as np

prefix = Path(sys.argv[1])
n = int(sys.argv[2])
m = int(sys.argv[3])
time_count = int(sys.argv[4])
seed = int(sys.argv[5])
rng = np.random.default_rng(seed)
prefix.parent.mkdir(parents=True, exist_ok=True)
times = np.linspace(0.0, 1.0, time_count)
samples = [(f"F{i+1:05d}", f"I{i+1:05d}") for i in range(n)]
bytes_per_snp = (n + 3) // 4
pad_n = bytes_per_snp * 4

idx_dyn = max(0, min(m - 1, m // 11))
idx_stable = max(0, min(m - 1, m // 4))
idx_null = max(0, min(m - 1, (3 * m) // 4))
idx_early = max(0, min(m - 1, m // 2))
idx_late = max(0, min(m - 1, (2 * m) // 3))
causal_idx = {idx_dyn, idx_stable, idx_early, idx_late}

maf = 0.10 + 0.38 * ((np.arange(m) * 37) % 997) / 996.0
G = np.empty((m, n), dtype=np.uint8)
for j in range(m):
    p = maf[j]
    G[j] = (rng.random(n) < p).astype(np.uint8) + (rng.random(n) < p).astype(np.uint8)

centered = {j: G[j].astype(float) - G[j].mean() for j in causal_idx}
subject = rng.normal(0.0, 0.18, size=n)
Y = np.zeros((n, time_count), dtype=float)
true_dyn = []
for k, t in enumerate(times):
    centered_t = 2.0 * t - 1.0
    beta_dyn = 0.90 * centered_t
    beta_stable = 0.55
    beta_early = 0.75 * (1.0 - t) * (1.0 - t)
    beta_late = 0.75 * t * t
    true_dyn.append(beta_dyn)
    Y[:, k] = (
        1.0 + subject + 0.25 * centered_t
        + beta_dyn * centered[idx_dyn]
        + beta_stable * centered[idx_stable]
        + beta_early * centered[idx_early]
        + beta_late * centered[idx_late]
        + rng.normal(0.0, 0.10, size=n)
    )

with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for fid, iid in samples:
        f.write(f"{fid}\t{iid}\t0\t0\t0\t-9\n")

with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for j in range(m):
        chrom = min(3, 1 + j // max(1, m // 3))
        if j == idx_dyn:
            snp = "rrmc_dyn_causal"
        elif j == idx_stable:
            snp = "rrmc_stable_causal"
        elif j == idx_null:
            snp = "rrmc_null_snp"
        elif j == idx_early:
            snp = "rrmc_early_causal"
        elif j == idx_late:
            snp = "rrmc_late_causal"
        else:
            snp = f"rrmc_bg_{j+1:06d}"
        f.write(f"{chrom}\t{snp}\t0\t{1000 + j * 17}\tA\tG\n")

with prefix.with_suffix(".bed").open("wb") as bed:
    bed.write(b"\x6c\x1b\x01")
    codes = np.empty((m, pad_n), dtype=np.uint8)
    codes[:, :n] = np.where(G == 2, 0, np.where(G == 1, 2, 3)).astype(np.uint8)
    if pad_n > n:
        codes[:, n:] = 0
    packed = (codes[:, 0::4] | (codes[:, 1::4] << 2) | (codes[:, 2::4] << 4) | (codes[:, 3::4] << 6)).astype(np.uint8)
    bed.write(packed.tobytes(order="C"))

phenos = [(fid, iid, row) for (fid, iid), row in zip(samples, Y)]
with Path(str(prefix) + ".long.pheno").open("w", encoding="utf-8") as f:
    f.write("FID IID TIME Y\n")
    for fid, iid, row in phenos:
        for t, y in zip(times, row):
            f.write(f"{fid}\t{iid}\t{t:.8f}\t{y:.8f}\n")

with Path(str(prefix) + ".mean.pheno").open("w", encoding="utf-8") as fm, \
     Path(str(prefix) + ".t0.pheno").open("w", encoding="utf-8") as f0, \
     Path(str(prefix) + ".t1.pheno").open("w", encoding="utf-8") as f1, \
     Path(str(prefix) + ".wide.pheno").open("w", encoding="utf-8") as fw:
    fm.write("FID IID Y\n")
    f0.write("FID IID Y\n")
    f1.write("FID IID Y\n")
    fw.write("FID IID " + " ".join(f"Y{k+1}" for k in range(time_count)) + "\n")
    for fid, iid, row in phenos:
        fm.write(f"{fid}\t{iid}\t{row.mean():.8f}\n")
        f0.write(f"{fid}\t{iid}\t{row[0]:.8f}\n")
        f1.write(f"{fid}\t{iid}\t{row[-1]:.8f}\n")
        fw.write(f"{fid}\t{iid}\t" + "\t".join(f"{x:.8f}" for x in row) + "\n")

with Path(str(prefix) + ".true_curve.tsv").open("w", encoding="utf-8") as f:
    f.write("SNP\tTime\tBeta\n")
    for t, b in zip(times, true_dyn):
        f.write(f"rrmc_dyn_causal\t{t:.8f}\t{b:.8f}\n")
PY

  run_timed() {
    local method="$1"; shift
    local out_prefix="${RESULT_DIR}/rep${rep}/${method}"
    local time_file="${out_prefix}.time"
    /usr/bin/time -f '%e\t%M' -o "${time_file}" "$@" --out "${out_prefix}" > "${out_prefix}.log" 2>&1
  }

  run_timed rrm_fast "${BIN}" --bfile "${PREFIX}" --long-pheno "${LONG_PHENO}" --tpheno 1 --mpheno 2 --rrm-order 2 --rrm-ai-iter 0 --rrm-curve-grid "${TIMES}" --rrm-top 20 --rrm-sig-threshold "${THRESHOLD}" --maf 0 --gwas --no-loco
  run_timed rrm_exact "${BIN}" --bfile "${PREFIX}" --long-pheno "${LONG_PHENO}" --tpheno 1 --mpheno 2 --rrm-order 2 --rrm-ai-iter 8 --rrm-curve-grid "${TIMES}" --rrm-top 20 --rrm-sig-threshold "${THRESHOLD}" --maf 0 --gwas --no-loco
  run_timed mean_gwas "${BIN}" --bfile "${PREFIX}" --pheno "${MEAN_PHENO}" --mpheno 1 --maf 0 --gwas --no-loco
  run_timed time0_gwas "${BIN}" --bfile "${PREFIX}" --pheno "${T0_PHENO}" --mpheno 1 --maf 0 --gwas --no-loco
  run_timed time1_gwas "${BIN}" --bfile "${PREFIX}" --pheno "${T1_PHENO}" --mpheno 1 --maf 0 --gwas --no-loco
  run_timed mv_gwas "${BIN}" --bfile "${PREFIX}" --pheno "${WIDE_PHENO}" --mpheno-list "$(seq -s, 1 "${TIMES}")" --multivariate --mv-method emmv --maf 0 --gwas --no-loco

  "${PYTHON_BIN}" - "${RESULT_DIR}" "${rep}" "${TRUE_CURVE}" >> "${METRICS_TSV}" <<'PY'
import csv
import math
import sys
from pathlib import Path

result_dir = Path(sys.argv[1])
rep = sys.argv[2]
true_curve_path = Path(sys.argv[3])
targets = ["rrmc_dyn_causal", "rrmc_stable_causal", "rrmc_null_snp"]

def read_time(prefix):
    vals = Path(str(prefix) + ".time").read_text(encoding="utf-8").strip().split("\t")
    return vals[0], vals[1]

def rank_assoc(path, snp, pcol="P"):
    rows = list(csv.DictReader(path.open("r", encoding="utf-8"), delimiter="\t"))
    rows = [r for r in rows if r.get("SNP") and r.get(pcol) not in (None, "", "nan")]
    rows.sort(key=lambda r: float(r[pcol]))
    for i, r in enumerate(rows, 1):
        if r["SNP"] == snp:
            return str(i), r[pcol], r
    return "NA", "NA", None

true_curve = {}
for row in csv.DictReader(true_curve_path.open("r", encoding="utf-8"), delimiter="\t"):
    true_curve[float(row["Time"])] = float(row["Beta"])

def curve_rmse(prefix, snp):
    path = Path(str(prefix) + ".rrm.dynamic_curve.txt")
    if not path.exists():
        return "NA"
    sq = []
    for row in csv.DictReader(path.open("r", encoding="utf-8"), delimiter="\t"):
        if row["SNP"] != snp:
            continue
        t = float(row["Time"])
        if t in true_curve:
            sq.append((float(row["Beta"]) - true_curve[t]) ** 2)
    return "NA" if not sq else str(math.sqrt(sum(sq) / len(sq)))

def pattern(prefix, snp):
    path = Path(str(prefix) + ".rrm.dynamic_summary.txt")
    if not path.exists():
        return "NA"
    for row in csv.DictReader(path.open("r", encoding="utf-8"), delimiter="\t"):
        if row["SNP"] == snp:
            return row.get("Pattern", "NA")
    return "NA"

methods = {
    "rrm_fast": ("rrm_joint.assoc.txt", "JOINT_P"),
    "rrm_exact": ("rrm_joint.assoc.txt", "JOINT_P"),
    "mean_gwas": ("assoc.txt", "P"),
    "time0_gwas": ("assoc.txt", "P"),
    "time1_gwas": ("assoc.txt", "P"),
    "mv_gwas": ("assoc.mv.txt", "P"),
}

for method, (suffix, pcol) in methods.items():
    prefix = result_dir / f"rep{rep}" / method
    runtime, rss = read_time(prefix)
    assoc = Path(str(prefix) + "." + suffix)
    for snp in targets:
        rank, p, _ = rank_assoc(assoc, snp, pcol)
        rmse = curve_rmse(prefix, snp) if method.startswith("rrm_") and snp == "rrmc_dyn_causal" else "NA"
        pat = pattern(prefix, snp) if method.startswith("rrm_") else "NA"
        print(f"{rep}\t{method}\t{snp}\t{rank}\t{p}\t{runtime}\t{rss}\t{rmse}\t{pat}")
PY
done

"${PYTHON_BIN}" - "${METRICS_TSV}" "${SUMMARY_TSV}" "${SUMMARY_MD}" "${THRESHOLD}" <<'PY'
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

metrics_path = Path(sys.argv[1])
summary_path = Path(sys.argv[2])
md_path = Path(sys.argv[3])
threshold = float(sys.argv[4])
rows = list(csv.DictReader(metrics_path.open("r", encoding="utf-8"), delimiter="\t"))
by_method = defaultdict(list)
for row in rows:
    by_method[row["Method"]].append(row)

summary_rows = []
for method, mrows in sorted(by_method.items()):
    reps = sorted({r["Rep"] for r in mrows})
    dyn = [r for r in mrows if r["Target"] == "rrmc_dyn_causal"]
    stable = [r for r in mrows if r["Target"] == "rrmc_stable_causal"]
    null = [r for r in mrows if r["Target"] == "rrmc_null_snp"]
    def valid_p(row):
        try:
            return float(row["P"])
        except Exception:
            return math.nan
    def power(target_rows):
        vals = [valid_p(r) for r in target_rows]
        vals = [v for v in vals if math.isfinite(v)]
        return sum(v < threshold for v in vals) / len(vals) if vals else math.nan
    def mean_float(col, target_rows=mrows):
        vals = []
        for r in target_rows:
            try:
                v = float(r[col])
                if math.isfinite(v): vals.append(v)
            except Exception:
                pass
        return sum(vals) / len(vals) if vals else math.nan
    def median_rank(target_rows):
        vals = []
        for r in target_rows:
            try:
                vals.append(float(r["Rank"]))
            except Exception:
                pass
        if not vals: return math.nan
        vals.sort()
        return vals[len(vals)//2]
    summary_rows.append({
        "Method": method,
        "Replicates": len(reps),
        "DynamicPower": power(dyn),
        "StablePower": power(stable),
        "NullFPR": power(null),
        "MedianDynamicRank": median_rank(dyn),
        "MeanRuntimeSec": mean_float("RuntimeSec"),
        "MeanMaxRSSKB": mean_float("MaxRSSKB"),
        "MeanDynamicCurveRMSE": mean_float("CurveRMSE", dyn),
    })

with summary_path.open("w", encoding="utf-8") as out:
    fields = ["Method", "Replicates", "DynamicPower", "StablePower", "NullFPR", "MedianDynamicRank", "MeanRuntimeSec", "MeanMaxRSSKB", "MeanDynamicCurveRMSE"]
    out.write("\t".join(fields) + "\n")
    for r in summary_rows:
        vals = []
        for f in fields:
            v = r[f]
            vals.append("NA" if isinstance(v, float) and not math.isfinite(v) else str(v))
        out.write("\t".join(vals) + "\n")

with md_path.open("w", encoding="utf-8") as out:
    out.write("# RRM Method Comparison\n\n")
    out.write(f"Threshold: P < {threshold}\n\n")
    out.write("| Method | Dynamic power | Stable power | Null FPR | Median dynamic rank | Runtime(s) | Curve RMSE |\n")
    out.write("|---|---:|---:|---:|---:|---:|---:|\n")
    for r in summary_rows:
        def fmt(v, digits=3):
            return "NA" if not math.isfinite(v) else f"{v:.{digits}g}"
        out.write(f"| {r['Method']} | {fmt(r['DynamicPower'])} | {fmt(r['StablePower'])} | {fmt(r['NullFPR'])} | {fmt(r['MedianDynamicRank'])} | {fmt(r['MeanRuntimeSec'])} | {fmt(r['MeanDynamicCurveRMSE'], 4)} |\n")

print(f"Method comparison summary: {summary_path}")
PY

echo "Done. RRM method comparison results: ${RESULT_DIR}"
