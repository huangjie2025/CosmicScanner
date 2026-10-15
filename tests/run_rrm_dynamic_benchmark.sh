#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
N="${N:-500}"
M="${M:-3000}"
TIMES="${TIMES:-5}"
CHROMS="${CHROMS:-8}"
SEED="${SEED:-20260611}"
RRM_ORDER="${RRM_ORDER:-2}"
RRM_BASIS="${RRM_BASIS:-legendre}"
RRM_GRID="${RRM_GRID:-11}"
RRM_TOP="${RRM_TOP:-50}"
RRM_SIG_THRESHOLD="${RRM_SIG_THRESHOLD:-0.01}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/tests/results/rrm_dynamic_benchmark_n${N}_m${M}_t${TIMES}}"
DATA_PREFIX="${RESULT_DIR}/data/rrm_dynamic_bench"
OUT_PREFIX="${RESULT_DIR}/rrm_dynamic_bench"
LONG_PHENO="${DATA_PREFIX}.long.pheno"
SUMMARY_TSV="${RESULT_DIR}/summary.tsv"
SUMMARY_MD="${RESULT_DIR}/summary.md"
TIME_LOG="${RESULT_DIR}/runtime.tsv"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${RESULT_DIR}/data" "${RESULT_DIR}/logs"
rm -f "${SUMMARY_TSV}" "${SUMMARY_MD}" "${TIME_LOG}"

"${PYTHON_BIN}" - "${DATA_PREFIX}" "${N}" "${M}" "${TIMES}" "${CHROMS}" "${SEED}" <<'PY'
import sys
from pathlib import Path

try:
    import numpy as np
except Exception as exc:
    raise SystemExit(f"This benchmark requires numpy. Set PYTHON_BIN to a Python with numpy. Error: {exc}")

prefix = Path(sys.argv[1])
n = int(sys.argv[2])
m = int(sys.argv[3])
time_count = int(sys.argv[4])
chroms = max(1, int(sys.argv[5]))
seed = int(sys.argv[6])
rng = np.random.default_rng(seed)
prefix.parent.mkdir(parents=True, exist_ok=True)

times = np.linspace(0.0, 1.0, time_count)
samples = [(f"F{i+1:07d}", f"I{i+1:07d}") for i in range(n)]
pop = (np.arange(n, dtype=np.float64) % 4) - 1.5
subject_effect = rng.normal(0.0, 0.25, size=n)
bytes_per_snp = (n + 3) // 4
pad_n = bytes_per_snp * 4

causal_dynamic = max(0, min(m - 1, m // 13))
causal_stable = max(0, min(m - 1, m // 5))
causal_early = max(0, min(m - 1, m // 3))
causal_late = max(0, min(m - 1, (2 * m) // 3))
causal_transient = max(0, min(m - 1, (4 * m) // 5))
causal = {causal_dynamic, causal_stable, causal_early, causal_late, causal_transient}

Y = np.zeros((n, time_count), dtype=np.float64)
for k, t in enumerate(times):
    Y[:, k] = 1.0 + subject_effect + 0.18 * pop + 0.25 * (2.0 * t - 1.0) + rng.normal(0.0, 0.12, size=n)

with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for fid, iid in samples:
        f.write(f"{fid}\t{iid}\t0\t0\t0\t-9\n")

with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    block = max(1, (m + chroms - 1) // chroms)
    for j in range(m):
        chrom = min(chroms, 1 + j // block)
        if j == causal_dynamic:
            snp = "rrmb_dyn_causal"
        elif j == causal_stable:
            snp = "rrmb_stable_causal"
        elif j == causal_early:
            snp = "rrmb_early_causal"
        elif j == causal_late:
            snp = "rrmb_late_causal"
        elif j == causal_transient:
            snp = "rrmb_transient_causal"
        else:
            snp = f"rrmb{j+1:09d}"
        f.write(f"{chrom}\t{snp}\t0\t{1000 + j * 11}\tA\tG\n")

with prefix.with_suffix(".bed").open("wb") as bed:
    bed.write(b"\x6c\x1b\x01")
    chunk = 1024
    pop_shift = pop[None, :].astype(np.float32)
    for start in range(0, m, chunk):
        end = min(m, start + chunk)
        idx = np.arange(start, end, dtype=np.float64)
        base_p = 0.08 + 0.42 * ((idx.astype(np.int64) * 3571) % 1000) / 999.0
        p = np.clip(base_p[:, None].astype(np.float32) + 0.025 * np.sin(idx[:, None] * 0.019).astype(np.float32) * pop_shift, 0.01, 0.95)
        g = (rng.random((end - start, n), dtype=np.float32) < p).astype(np.uint8)
        g += (rng.random((end - start, n), dtype=np.float32) < p).astype(np.uint8)

        for snp_idx in causal:
            if start <= snp_idx < end:
                local = snp_idx - start
                centered = g[local].astype(np.float64) - np.mean(g[local].astype(np.float64))
                for k, t in enumerate(times):
                    centered_t = 2.0 * t - 1.0
                    transient = max(0.0, 1.0 - 4.0 * (t - 0.5) * (t - 0.5))
                    if snp_idx == causal_dynamic:
                        beta = 0.65 * centered_t
                    elif snp_idx == causal_stable:
                        beta = 0.45
                    elif snp_idx == causal_early:
                        beta = 0.55 * (1.0 - t) * (1.0 - t)
                    elif snp_idx == causal_late:
                        beta = 0.55 * t * t
                    else:
                        beta = 0.55 * transient
                    Y[:, k] += centered * beta

        codes = np.empty((end - start, pad_n), dtype=np.uint8)
        codes[:, :n] = np.where(g == 2, 0, np.where(g == 1, 2, 3)).astype(np.uint8)
        if pad_n > n:
            codes[:, n:] = 0
        packed = (codes[:, 0::4] | (codes[:, 1::4] << 2) | (codes[:, 2::4] << 4) | (codes[:, 3::4] << 6)).astype(np.uint8)
        bed.write(packed.tobytes(order="C"))

with Path(str(prefix) + ".long.pheno").open("w", encoding="utf-8") as f:
    f.write("FID IID TIME Y\n")
    for (fid, iid), row in zip(samples, Y):
        for t, y in zip(times, row):
            f.write(f"{fid}\t{iid}\t{t:.8f}\t{y:.8f}\n")

with Path(str(prefix) + ".causal.tsv").open("w", encoding="utf-8") as f:
    f.write("SNP\tPattern\n")
    f.write("rrmb_dyn_causal\tsign_switch\n")
    f.write("rrmb_stable_causal\tstable\n")
    f.write("rrmb_early_causal\tearly_specific\n")
    f.write("rrmb_late_causal\tlate_specific\n")
    f.write("rrmb_transient_causal\ttransient\n")

print(f"Generated RRM dynamic benchmark data: n={n}, m={m}, times={time_count}, observations={n * time_count}")
PY

printf 'ElapsedSeconds\tMaxRSSKB\n' > "${TIME_LOG}"
/usr/bin/time -f '%e\t%M' -o "${TIME_LOG}.tmp" \
  "${BIN}" \
    --bfile "${DATA_PREFIX}" \
    --long-pheno "${LONG_PHENO}" \
    --tpheno 1 \
    --mpheno 2 \
    --rrm-order "${RRM_ORDER}" \
    --rrm-basis "${RRM_BASIS}" \
    --rrm-curve-grid "${RRM_GRID}" \
    --rrm-top "${RRM_TOP}" \
    --rrm-sig-threshold "${RRM_SIG_THRESHOLD}" \
    --maf 0 \
    --gwas \
    --no-loco \
    --out "${OUT_PREFIX}" \
    > "${RESULT_DIR}/logs/cosmicscanner.log" 2>&1
cat "${TIME_LOG}.tmp" >> "${TIME_LOG}"
rm -f "${TIME_LOG}.tmp"

"${PYTHON_BIN}" - "${RESULT_DIR}" "${OUT_PREFIX}" "${N}" "${M}" "${TIMES}" "${RRM_ORDER}" "${RRM_GRID}" "${RRM_TOP}" "${RRM_BASIS}" <<'PY'
import csv
import sys
from pathlib import Path

result_dir = Path(sys.argv[1])
out_prefix = Path(sys.argv[2])
n, m, time_count = map(int, sys.argv[3:6])
rrm_order, rrm_grid, rrm_top = map(int, sys.argv[6:9])
rrm_basis = sys.argv[9]

summary_path = Path(str(out_prefix) + ".rrm.dynamic_summary.txt")
curve_path = Path(str(out_prefix) + ".rrm.dynamic_curve.txt")
joint_path = Path(str(out_prefix) + ".rrm_joint.assoc.txt")
null_path = Path(str(out_prefix) + ".rrm_joint.null.txt")
runtime_path = result_dir / "runtime.tsv"
summary_tsv = result_dir / "summary.tsv"
summary_md = result_dir / "summary.md"

summary = list(csv.DictReader(summary_path.open("r", encoding="utf-8"), delimiter="\t"))
curve = list(csv.DictReader(curve_path.open("r", encoding="utf-8"), delimiter="\t"))
joint = list(csv.DictReader(joint_path.open("r", encoding="utf-8"), delimiter="\t"))
null_rows = list(csv.DictReader(null_path.open("r", encoding="utf-8"), delimiter="\t"))
runtime = list(csv.DictReader(runtime_path.open("r", encoding="utf-8"), delimiter="\t"))[0]
null_mode = null_rows[0].get("NullMode", "NA") if null_rows else "NA"
ai_iter = null_rows[0].get("AIIter", "NA") if null_rows else "NA"

if len(summary) != rrm_top:
    raise SystemExit(f"Expected {rrm_top} summary rows, got {len(summary)}")
actual_grid = int(summary[0].get("GridPoints", rrm_grid)) if summary else rrm_grid
if len(curve) != rrm_top * actual_grid:
    raise SystemExit(f"Expected {rrm_top * actual_grid} curve rows, got {len(curve)}")
if len(joint) != m:
    raise SystemExit(f"Expected {m} joint rows, got {len(joint)}")

causal_snps = ["rrmb_dyn_causal", "rrmb_stable_causal", "rrmb_early_causal", "rrmb_late_causal", "rrmb_transient_causal"]
rows_by_snp = {row["SNP"]: row for row in summary}
dyn = rows_by_snp.get("rrmb_dyn_causal")
if dyn is None:
    raise SystemExit("Dynamic causal SNP missing from benchmark top summary")

with summary_tsv.open("w", encoding="utf-8") as out:
    out.write("Metric\tValue\n")
    out.write(f"N\t{n}\n")
    out.write(f"M\t{m}\n")
    out.write(f"TimePoints\t{time_count}\n")
    out.write(f"Observations\t{n * time_count}\n")
    out.write(f"RrmOrder\t{rrm_order}\n")
    out.write(f"RrmBasis\t{rrm_basis}\n")
    out.write(f"RequestedRrmGrid\t{rrm_grid}\n")
    out.write(f"ActualCurveGrid\t{actual_grid}\n")
    out.write(f"RrmTop\t{rrm_top}\n")
    out.write(f"NullMode\t{null_mode}\n")
    out.write(f"AIIter\t{ai_iter}\n")
    out.write(f"ElapsedSeconds\t{runtime['ElapsedSeconds']}\n")
    out.write(f"MaxRSSKB\t{runtime['MaxRSSKB']}\n")
    out.write(f"JointRows\t{len(joint)}\n")
    out.write(f"SummaryRows\t{len(summary)}\n")
    out.write(f"CurveRows\t{len(curve)}\n")
    for snp in causal_snps:
        row = rows_by_snp.get(snp)
        out.write(f"{snp}.Rank\t{row['Rank'] if row else 'NA'}\n")
        out.write(f"{snp}.JointP\t{row['JointP'] if row else 'NA'}\n")
        out.write(f"{snp}.Pattern\t{row['Pattern'] if row else 'NA'}\n")

with summary_md.open("w", encoding="utf-8") as out:
    out.write("# RRM Dynamic Benchmark Summary\n\n")
    out.write(f"- N: {n}\n")
    out.write(f"- M: {m}\n")
    out.write(f"- Time points: {time_count}\n")
    out.write(f"- Observations: {n * time_count}\n")
    out.write(f"- RRM order: {rrm_order}\n")
    out.write(f"- RRM basis: {rrm_basis}\n")
    out.write(f"- Null mode: {null_mode}\n")
    out.write(f"- AI iterations: {ai_iter}\n")
    out.write(f"- Requested curve grid: {rrm_grid}\n")
    out.write(f"- Actual curve grid: {actual_grid}\n")
    out.write(f"- Runtime: {runtime['ElapsedSeconds']} seconds\n")
    out.write(f"- Max RSS: {runtime['MaxRSSKB']} KB\n")
    out.write(f"- Dynamic causal rank: {dyn['Rank']}\n")
    out.write(f"- Dynamic causal pattern: {dyn['Pattern']}\n")

print(f"RRM dynamic benchmark completed: {summary_tsv}")
PY

echo "Done. RRM dynamic benchmark results: ${RESULT_DIR}"
