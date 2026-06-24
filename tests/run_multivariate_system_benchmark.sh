#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
N="${N:-240}"
M="${M:-1200}"
TRAITS="${TRAITS:-2}"
CHROMS="${CHROMS:-8}"
SEED="${SEED:-20260608}"
REFINE_TOP="${REFINE_TOP:-50}"
REFINE_P="${REFINE_P:-1}"
REFINE_DISCREPANCY_TOP="${REFINE_DISCREPANCY_TOP:-25}"
REFINE_DISCREPANCY="${REFINE_DISCREPANCY:-1}"
REFINE_SCORE_PREPASS="${REFINE_SCORE_PREPASS:-200}"
RUN_EMMVGWAS="${RUN_EMMVGWAS:-1}"
RUN_GEMMA="${RUN_GEMMA:-1}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/tests/results/multivariate_system_benchmark_n${N}_m${M}_t${TRAITS}}"
DATA_PREFIX="${RESULT_DIR}/data/sim_mv_bench"
GEMMA_WORK="${RESULT_DIR}/gemma_work"
SUMMARY_TSV="${RESULT_DIR}/summary.tsv"
METHOD_TSV="${RESULT_DIR}/method_metrics.tsv"
SUMMARY_MD="${RESULT_DIR}/summary.md"

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

mkdir -p "${RESULT_DIR}/data" "${GEMMA_WORK}/output" "${RESULT_DIR}/logs" "${RESULT_DIR}/emmvgwas"
rm -rf "${GEMMA_WORK}/output"
mkdir -p "${GEMMA_WORK}/output"
rm -f "${SUMMARY_TSV}" "${METHOD_TSV}" "${SUMMARY_MD}"

"${PYTHON_BIN}" - "${DATA_PREFIX}" "${N}" "${M}" "${TRAITS}" "${CHROMS}" "${SEED}" <<'PY'
import sys
from pathlib import Path
import numpy as np

prefix = Path(sys.argv[1])
n = int(sys.argv[2])
m = int(sys.argv[3])
t = int(sys.argv[4])
chroms = max(1, int(sys.argv[5]))
seed = int(sys.argv[6])
rng = np.random.default_rng(seed)
prefix.parent.mkdir(parents=True, exist_ok=True)

pop = np.arange(n, dtype=np.int16) % 4
samples = [(f"F{i+1:07d}", f"I{i+1:07d}") for i in range(n)]
bytes_per_snp = (n + 3) // 4
pad_n = bytes_per_snp * 4

causal_idx = sorted(set([max(0, min(m - 1, m // 17)), max(0, min(m - 1, m // 5)),
                         max(0, min(m - 1, m // 2)), max(0, min(m - 1, (4 * m) // 5))]))
effects = {}
for k, idx in enumerate(causal_idx):
    signs = np.array([1 if ((k + j) % 2 == 0) else -1 for j in range(t)], dtype=float)
    effects[idx] = signs * (0.18 + 0.08 * ((k + np.arange(t)) % 3))

shared = rng.normal(0.0, 0.45, size=(n, 1))
Y = 0.18 * pop[:, None] * (np.arange(t)[None, :] + 1) / t
Y = Y + shared * (0.8 - 0.1 * np.arange(t)[None, :]) + rng.normal(0.0, 0.65, size=(n, t))

with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for j in range(m):
        chrom = (j % chroms) + 1
        f.write(f"{chrom}\tmvb{j+1:09d}\t0\t{1000 + j * 13}\tA\tC\n")

with prefix.with_suffix(".bed").open("wb") as bed, prefix.with_suffix(".emmv.Z.csv").open("w", encoding="utf-8") as zf:
    bed.write(b"\x6c\x1b\x01")
    zf.write("SNP," + ",".join(iid for _, iid in samples) + "\n")
    chunk = 1024
    pop_shift = (pop.astype(np.float32) - 1.5)[None, :]
    for start in range(0, m, chunk):
        end = min(m, start + chunk)
        idx = np.arange(start, end, dtype=np.float64)
        base_p = 0.06 + 0.44 * ((idx.astype(np.int64) * 6151) % 1000) / 999.0
        shift = 0.035 * np.sin(idx * 0.017)
        p = np.clip(base_p[:, None].astype(np.float32) + shift[:, None].astype(np.float32) * pop_shift, 0.01, 0.95)
        g = (rng.random((end - start, n), dtype=np.float32) < p).astype(np.uint8)
        g += (rng.random((end - start, n), dtype=np.float32) < p).astype(np.uint8)

        for snp_idx, beta in effects.items():
            if start <= snp_idx < end:
                local = snp_idx - start
                centered = g[local].astype(np.float64) - 2.0 * p[local].astype(np.float64)
                Y[:] += centered[:, None] * beta[None, :]

        codes = np.empty((end - start, pad_n), dtype=np.uint8)
        codes[:, :n] = np.where(g == 2, 0, np.where(g == 1, 2, 3)).astype(np.uint8)
        if pad_n > n:
            codes[:, n:] = 0
        packed = (codes[:, 0::4] | (codes[:, 1::4] << 2) | (codes[:, 2::4] << 4) | (codes[:, 3::4] << 6)).astype(np.uint8)
        bed.write(packed.tobytes(order="C"))
        for local, row in enumerate(g):
            zf.write(f"mvb{start + local + 1:09d}," + ",".join(str(int(x) - 1) for x in row) + "\n")

with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for (fid, iid), y in zip(samples, Y[:, 0]):
        f.write(f"{fid}\t{iid}\t0\t0\t1\t{float(y):.8f}\n")

with prefix.with_suffix(".pheno").open("w", encoding="utf-8") as f:
    f.write("FID\tIID\t" + "\t".join(f"Trait{i+1}" for i in range(t)) + "\n")
    for (fid, iid), vals in zip(samples, Y):
        f.write(f"{fid}\t{iid}\t" + "\t".join(f"{float(v):.8f}" for v in vals) + "\n")

with prefix.with_suffix(".pheno.no_header").open("w", encoding="utf-8") as f:
    for vals in Y:
        f.write("\t".join(f"{float(v):.8f}" for v in vals) + "\n")

with prefix.with_suffix(".emmv.P.csv").open("w", encoding="utf-8") as f:
    f.write(",".join(f"Trait{i+1}" for i in range(t)) + "\n")
    for vals in Y:
        f.write(",".join(f"{float(v):.8f}" for v in vals) + "\n")

with prefix.with_suffix(".truth.tsv").open("w", encoding="utf-8") as f:
    f.write("SNP\t" + "\t".join(f"BETA_{i+1}" for i in range(t)) + "\n")
    for idx, beta in effects.items():
        f.write(f"mvb{idx+1:09d}\t" + "\t".join(f"{float(v):.8f}" for v in beta) + "\n")
PY

trait_list="$("${PYTHON_BIN}" - "${TRAITS}" <<'PY'
import sys
t = int(sys.argv[1])
print(','.join(str(i) for i in range(1, t + 1)))
PY
)"
IFS=',' read -r -a trait_args <<< "${trait_list}"

run_timed() {
  local method="$1"
  local prefix="$2"
  shift 2
  local log="${RESULT_DIR}/logs/${method}.log"
  local metric="${RESULT_DIR}/logs/${method}.time"
  /usr/bin/time -f 'seconds\t%e\nmax_rss_kb\t%M' -o "${metric}" "$@" >"${log}" 2>&1
  "${PYTHON_BIN}" - "${SUMMARY_TSV}" "${method}" "${prefix}" "${log}" "${metric}" <<'PY'
import sys
from pathlib import Path
summary, method, prefix, log, metric = map(Path, sys.argv[1:])
vals = {"seconds": "nan", "max_rss_kb": "nan"}
for line in metric.read_text(encoding="utf-8").splitlines():
    parts = line.split()
    if len(parts) == 2:
        vals[parts[0]] = parts[1]
if not summary.exists():
    summary.write_text("method\tseconds\tmax_rss_kb\tprefix\tlog\n", encoding="utf-8")
with summary.open("a", encoding="utf-8") as f:
    f.write(f"{method}\t{vals['seconds']}\t{vals['max_rss_kb']}\t{prefix}\t{log}\n")
PY
}

echo "[mv-system-benchmark] Running CS-fast..."
run_timed "CS-fast" "${RESULT_DIR}/cs_fast" "${BIN}" \
  --bfile "${DATA_PREFIX}" --pheno "${DATA_PREFIX}.pheno" --mpheno-list "${trait_list}" \
  --multivariate --mv-method emmv --mv-compatible emmvgwas --gwas --no-loco --maf 0.0 \
  --out "${RESULT_DIR}/cs_fast"

echo "[mv-system-benchmark] Running CS-hybrid..."
run_timed "CS-hybrid" "${RESULT_DIR}/cs_hybrid" "${BIN}" \
  --bfile "${DATA_PREFIX}" --pheno "${DATA_PREFIX}.pheno" --mpheno-list "${trait_list}" \
  --multivariate --mv-method hybrid --mv-test wald,lrt,score --mv-refine-p "${REFINE_P}" --mv-refine-top "${REFINE_TOP}" \
  --mv-refine-discrepancy-top "${REFINE_DISCREPANCY_TOP}" --mv-refine-discrepancy "${REFINE_DISCREPANCY}" \
  --mv-refine-score-prepass "${REFINE_SCORE_PREPASS}" \
  --gwas --no-loco --maf 0.0 --out "${RESULT_DIR}/cs_hybrid"

echo "[mv-system-benchmark] Running CS-LOCO-hybrid..."
run_timed "CS-LOCO-hybrid" "${RESULT_DIR}/cs_loco_hybrid" "${BIN}" \
  --bfile "${DATA_PREFIX}" --pheno "${DATA_PREFIX}.pheno" --mpheno-list "${trait_list}" \
  --multivariate --mv-method hybrid --mv-test wald,lrt,score --mv-refine-p "${REFINE_P}" --mv-refine-top "${REFINE_TOP}" \
  --mv-refine-discrepancy-top "${REFINE_DISCREPANCY_TOP}" --mv-refine-discrepancy "${REFINE_DISCREPANCY}" \
  --mv-refine-score-prepass "${REFINE_SCORE_PREPASS}" \
  --gwas --maf 0.0 --out "${RESULT_DIR}/cs_loco_hybrid"

if [[ "${RUN_EMMVGWAS}" == "1" ]]; then
  if RSCRIPT_EXE="$(resolve_rscript)"; then
    echo "[mv-system-benchmark] Running EMmvGWAS R..."
    run_timed "EMmvGWAS-R" "${RESULT_DIR}/emmvgwas/emmv" "${RSCRIPT_EXE}" \
      "${ROOT}/tests/run_emmvgwas_validation.R" "${DATA_PREFIX}.emmv.P.csv" "${DATA_PREFIX}.emmv.Z.csv" "${RESULT_DIR}/emmvgwas/emmv"
  else
    printf 'EMmvGWAS-R\tnan\tnan\t%s\tblocked_no_Rscript\n' "${RESULT_DIR}/emmvgwas/emmv" >> "${SUMMARY_TSV}"
  fi
fi

if [[ "${RUN_GEMMA}" == "1" ]]; then
  if GEMMA_EXE="$(resolve_gemma)"; then
    echo "[mv-system-benchmark] Running GEMMA kinship..."
    run_timed "GEMMA-kinship" "${GEMMA_WORK}/output/gemma_mv_kinship" bash -c \
      'cd "$1" && "$2" -bfile "$3" -gk 1 -o gemma_mv_kinship' _ "${GEMMA_WORK}" "${GEMMA_EXE}" "${DATA_PREFIX}"
    echo "[mv-system-benchmark] Running GEMMA multivariate LMM..."
    KINSHIP="${GEMMA_WORK}/output/gemma_mv_kinship.cXX.txt"
    run_timed "GEMMA-LMM" "${GEMMA_WORK}/output/gemma_mv_lmm" bash -c \
      'work="$1"; exe="$2"; bfile="$3"; pheno="$4"; kinship="$5"; shift 5; cd "${work}" && "${exe}" -bfile "${bfile}" -p "${pheno}" -n "$@" -k "${kinship}" -lmm -o gemma_mv_lmm' _ \
      "${GEMMA_WORK}" "${GEMMA_EXE}" "${DATA_PREFIX}" "${DATA_PREFIX}.pheno.no_header" "${KINSHIP}" "${trait_args[@]}"
  else
    printf 'GEMMA\tnan\tnan\t%s\tblocked_no_GEMMA\n' "${GEMMA_WORK}/output/gemma_mv_lmm" >> "${SUMMARY_TSV}"
  fi
fi

"${PYTHON_BIN}" - "${RESULT_DIR}" "${DATA_PREFIX}" "${N}" "${M}" "${TRAITS}" "${REFINE_TOP}" <<'PY'
import csv
import math
import sys
from pathlib import Path

result_dir = Path(sys.argv[1])
data_prefix = Path(sys.argv[2])
n, m, traits, refine_top = map(int, sys.argv[3:7])
summary_tsv = result_dir / "summary.tsv"
method_tsv = result_dir / "method_metrics.tsv"
summary_md = result_dir / "summary.md"

def read_tsv(path):
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def read_csv(path):
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def read_kv(path):
    out = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t", 1)
        if len(parts) == 2:
            out[parts[0]] = parts[1]
    return out

def fval(row, *names):
    lower = {k.lower(): k for k in row}
    for name in names:
        key = lower.get(name.lower())
        if key and row.get(key, "") not in ("", "NA", "nan", "-nan"):
            try:
                return float(row[key])
            except ValueError:
                pass
    return math.nan

def corr(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 3:
        return math.nan
    mx = sum(x for x, _ in pairs) / len(pairs)
    my = sum(y for _, y in pairs) / len(pairs)
    vx = sum((x - mx) ** 2 for x, _ in pairs)
    vy = sum((y - my) ** 2 for _, y in pairs)
    if vx <= 0 or vy <= 0:
        return math.nan
    return sum((x - mx) * (y - my) for x, y in pairs) / math.sqrt(vx * vy)

def logp(p):
    return -math.log10(max(p, 1e-300)) if math.isfinite(p) and p > 0 else math.nan

runtime_rows = read_tsv(summary_tsv)
runtime = {r["method"]: r for r in runtime_rows}

gemma_rows = read_tsv(result_dir / "gemma_work" / "output" / "gemma_mv_lmm.assoc.txt")
gemma = {r.get("rs", r.get("SNP")): r for r in gemma_rows}
emmv_rows = read_csv(result_dir / "emmvgwas" / "emmv_association_with_ids.csv")
emmv = {r.get("SNP_ID"): r for r in emmv_rows}

methods = [
    ("CS-fast", result_dir / "cs_fast.assoc.mv.txt", "P"),
    ("CS-hybrid", result_dir / "cs_hybrid.assoc.mv.txt", "P"),
    ("CS-LOCO-hybrid", result_dir / "cs_loco_hybrid.assoc.mv.txt", "P"),
    ("GEMMA", result_dir / "gemma_work" / "output" / "gemma_mv_lmm.assoc.txt", "p_wald"),
    ("EMmvGWAS-R", result_dir / "emmvgwas" / "emmv_association_with_ids.csv", "p"),
]

truth = {r["SNP"] for r in read_tsv(data_prefix.with_suffix(".truth.tsv"))}

metric_rows = []
for name, path, pcol in methods:
    rows = read_csv(path) if path.suffix == ".csv" else read_tsv(path)
    if name == "GEMMA":
        snp_get = lambda r: r.get("rs", r.get("SNP"))
    elif name == "EMmvGWAS-R":
        snp_get = lambda r: r.get("SNP_ID")
    else:
        snp_get = lambda r: r.get("SNP")
    p_by_snp = {snp_get(r): fval(r, pcol) for r in rows if snp_get(r)}
    ordered = sorted(((logp(p), snp) for snp, p in p_by_snp.items() if math.isfinite(logp(p))), reverse=True)
    top10 = {snp for _, snp in ordered[:10]}
    top50 = {snp for _, snp in ordered[:50]}
    gemma_corr = math.nan
    gemma_top10 = math.nan
    emmv_corr = math.nan
    if gemma and name != "GEMMA":
        xs, ys = [], []
        for snp, p in p_by_snp.items():
            if snp in gemma:
                xs.append(logp(p))
                ys.append(logp(fval(gemma[snp], "p_wald")))
        gemma_corr = corr(xs, ys)
        gemma_top = sorted(((logp(fval(r, "p_wald")), snp) for snp, r in gemma.items()), reverse=True)[:10]
        gemma_top10 = len(top10 & {snp for _, snp in gemma_top})
    if emmv and name != "EMmvGWAS-R":
        xs, ys = [], []
        for snp, p in p_by_snp.items():
            if snp in emmv:
                xs.append(logp(p))
                ys.append(logp(fval(emmv[snp], "p")))
        emmv_corr = corr(xs, ys)
    rt = runtime.get(name, {})
    seconds = fval(rt, "seconds")
    rss = fval(rt, "max_rss_kb")
    if name == "GEMMA":
        kin = runtime.get("GEMMA-kinship", {})
        lmm = runtime.get("GEMMA-LMM", {})
        seconds = fval(kin, "seconds") + fval(lmm, "seconds")
        rss = max(fval(kin, "max_rss_kb"), fval(lmm, "max_rss_kb"))
    refined_rows = 0
    stage1_mode = "NA"
    candidates_by_p = 0
    candidates_by_top = 0
    candidates_by_discrepancy = 0
    discrepancy_available = "NA"
    score_prepass_scores = 0
    if name.startswith("CS-hybrid") or name.startswith("CS-LOCO"):
        refined_rows = max(0, len(read_tsv(path.with_name(path.name.replace(".assoc.mv.txt", ".refined.assoc.mv.txt")))))
        hybrid_summary = read_kv(path.with_name(path.name.replace(".assoc.mv.txt", ".hybrid.summary.tsv")))
        stage1_mode = hybrid_summary.get("stage1_mode", "NA")
        candidates_by_p = int(float(hybrid_summary.get("candidates_by_p", "0")))
        candidates_by_top = int(float(hybrid_summary.get("candidates_by_top", "0")))
        candidates_by_discrepancy = int(float(hybrid_summary.get("candidates_by_discrepancy", "0")))
        discrepancy_available = hybrid_summary.get("discrepancy_available", "NA")
        score_prepass_scores = int(float(hybrid_summary.get("score_prepass_scores", "0")))
    metric_rows.append({
        "method": name, "n": n, "m": m, "traits": traits,
        "seconds": seconds, "max_rss_mb": rss / 1024 if math.isfinite(rss) else math.nan,
        "assoc_rows": len(rows), "refined_rows": refined_rows, "stage1_mode": stage1_mode,
        "candidates_by_p": candidates_by_p, "candidates_by_top": candidates_by_top,
        "candidates_by_discrepancy": candidates_by_discrepancy, "discrepancy_available": discrepancy_available,
        "score_prepass_scores": score_prepass_scores,
        "log10p_corr_vs_gemma": gemma_corr, "top10_overlap_vs_gemma": gemma_top10,
        "log10p_corr_vs_emmvgwas": emmv_corr,
        "causal_in_top10": len(top10 & truth), "causal_in_top50": len(top50 & truth),
    })

fields = list(metric_rows[0].keys())
with method_tsv.open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(metric_rows)

def fmt(x):
    if isinstance(x, float):
        if math.isnan(x): return "NA"
        if abs(x) >= 100: return f"{x:.1f}"
        return f"{x:.4g}"
    return str(x)

with summary_md.open("w", encoding="utf-8") as f:
    f.write(f"# Multivariate System Benchmark\n\nN={n}, M={m}, traits={traits}, refine_top={refine_top}\n\n")
    f.write("| Method | Seconds | Max RSS MB | Rows | Refined | Score audit | Disc cand | P corr vs GEMMA | Top10 vs GEMMA | P corr vs EMmvGWAS | Causal Top10 | Causal Top50 |\n")
    f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    for r in metric_rows:
        f.write("| {method} | {seconds} | {max_rss_mb} | {assoc_rows} | {refined_rows} | {score_prepass_scores} | {candidates_by_discrepancy} | {log10p_corr_vs_gemma} | {top10_overlap_vs_gemma} | {log10p_corr_vs_emmvgwas} | {causal_in_top10} | {causal_in_top50} |\n".format(**{k: fmt(v) for k, v in r.items()}))

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plot_rows = [r for r in metric_rows if math.isfinite(r["seconds"])]
    plt.figure(figsize=(8, 4.5))
    plt.bar([r["method"] for r in plot_rows], [r["seconds"] for r in plot_rows], color=["#3b82f6", "#14b8a6", "#f97316", "#6366f1", "#64748b"][:len(plot_rows)])
    plt.yscale("log")
    plt.ylabel("Seconds (log scale)")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(result_dir / "runtime.png", dpi=180)
    plt.savefig(result_dir / "runtime.pdf")
    plt.close()

    plt.figure(figsize=(8, 4.5))
    xs = [r["method"] for r in metric_rows if math.isfinite(r["log10p_corr_vs_gemma"])]
    ys = [r["log10p_corr_vs_gemma"] for r in metric_rows if math.isfinite(r["log10p_corr_vs_gemma"])]
    plt.bar(xs, ys, color="#0f766e")
    plt.ylim(max(0.0, min(ys) - 0.05) if ys else 0, 1.01)
    plt.ylabel("Genome-wide -log10(P) corr vs GEMMA")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(result_dir / "accuracy.png", dpi=180)
    plt.savefig(result_dir / "accuracy.pdf")
    plt.close()
except Exception as exc:
    with (result_dir / "plot_error.txt").open("w", encoding="utf-8") as f:
        f.write(str(exc) + "\n")

print(summary_md.read_text(encoding="utf-8"))
PY

echo "[mv-system-benchmark] Outputs:"
echo "  ${METHOD_TSV}"
echo "  ${SUMMARY_MD}"
echo "  ${RESULT_DIR}/runtime.png"
echo "  ${RESULT_DIR}/accuracy.png"
