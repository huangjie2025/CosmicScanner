#!/bin/bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

ROOT="/disk228/huangjie/PanOmics/CosmicScanner"
BENCH_SCRIPT="${ROOT}/tests/run_repeatability_performance_benchmark.sh"
DATA_ROOT=${ACC_DATA_ROOT:-/tmp/opencode/repeatability_repeatabel_accuracy_data}
RESULT_ROOT=${ACC_RESULT_ROOT:-/tmp/opencode/repeatability_repeatabel_accuracy_results}
R_LIB=${ACC_REPEATABL_R_LIB:-/tmp/opencode/Rlib_repeatabel}
N=${ACC_N:-80}
M=${ACC_M:-300}
R=${ACC_R:-3}
TIMEOUT=${ACC_TIMEOUT_SEC:-180}

mkdir -p "${DATA_ROOT}" "${RESULT_ROOT}"

BENCH_DATA_ROOT="${DATA_ROOT}" \
BENCH_RESULT_ROOT="${RESULT_ROOT}" \
BENCH_N_LIST="${N}" \
BENCH_M_LIST="${M}" \
BENCH_R_LIST="${R}" \
BENCH_METHODS=exact,repeatabel \
BENCH_TIMEOUT_SEC="${TIMEOUT}" \
BENCH_REPEATABL_R_LIB="${R_LIB}" \
bash "${BENCH_SCRIPT}"

BASE="${RESULT_ROOT}/n${N}_m${M}_r${R}"
RB_DIR="${BASE}.repeatabel.repeatabel_input"
if [[ ! -d "${RB_DIR}" ]]; then
  matches=("${BASE}".repeatabel*.repeatabel_input)
  if [[ -d "${matches[0]}" ]]; then
    RB_DIR="${matches[0]}"
  fi
fi
RB_OUT="${BASE}.repeatabel.full.tsv"
RB_R_SCRIPT="${BASE}.repeatabel.extract_full.R"

cat > "${RB_R_SCRIPT}" <<'RS'
.libPaths(c(Sys.getenv("REPEATABL_R_LIB"), .libPaths()))
suppressPackageStartupMessages(library(RepeatABEL))
suppressPackageStartupMessages(library(GenABEL))
rb_dir <- Sys.getenv("REPEATABL_INPUT_DIR")
out_file <- Sys.getenv("REPEATABL_OUT_FILE")
geno_raw <- file.path(rb_dir, "genos.raw")
if (!file.exists(geno_raw)) convert.snp.text(file.path(rb_dir, "genos.txt"), geno_raw)
gd <- load.gwaa.data(
  phenofile=file.path(rb_dir, "pheno_ind.dat"),
  genofile=geno_raw,
  force=TRUE,
  sort=FALSE,
  id="id"
)
phen <- read.table(file.path(rb_dir, "phen_repeat.dat"), header=TRUE)
fit <- rGLS(y ~ batch, genabel.data=gd, phenotype.data=phen, id.name="id")
res <- as.data.frame(fit@results)
res <- cbind(SNP=rownames(res), res)
write.table(res, out_file, sep="\t", quote=FALSE, row.names=FALSE)
RS

REPEATABL_R_LIB="${R_LIB}" \
REPEATABL_INPUT_DIR="${RB_DIR}" \
REPEATABL_OUT_FILE="${RB_OUT}" \
conda run -n R_GS Rscript "${RB_R_SCRIPT}"

python3 - <<PY
from pathlib import Path
import csv, math, statistics

base = Path("${BASE}")
cs_file = Path(str(base) + ".exact.assoc.txt")
rb_file = Path("${RB_OUT}")
summary_file = Path(str(base) + ".accuracy.tsv")

def read_cs(path):
    out = {}
    with path.open("r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                p = float(r["P"])
                beta = float(r["BETA"])
                se = float(r["SE"])
            except Exception:
                continue
            if math.isfinite(p) and p > 0 and math.isfinite(beta):
                out[r["SNP"]] = dict(beta=beta, se=se, p=p, logp=-math.log10(p))
    return out

def read_rb(path):
    out = {}
    with path.open("r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                p = float(r["P1df"])
                beta = float(r["effB"])
                se = float(r["se_effB"])
            except Exception:
                continue
            if math.isfinite(p) and p > 0 and math.isfinite(beta):
                out[r["SNP"]] = dict(beta=beta, se=se, p=p, logp=-math.log10(p))
    return out

def corr(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    vx = sum((v - mx) ** 2 for v in x)
    vy = sum((v - my) ** 2 for v in y)
    if vx <= 0 or vy <= 0:
        return float("nan")
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(vx * vy)

def top_overlap(a, b, k):
    ka = [s for s, _ in sorted(a.items(), key=lambda kv: kv[1]["p"])[:k]]
    kb = [s for s, _ in sorted(b.items(), key=lambda kv: kv[1]["p"])[:k]]
    return len(set(ka).intersection(kb))

cs = read_cs(cs_file)
rb = read_rb(rb_file)
common = sorted(set(cs).intersection(rb))
if len(common) < 10:
    raise SystemExit(f"Too few overlapping SNPs: {len(common)}")

cs_beta = [cs[s]["beta"] for s in common]
rb_beta = [rb[s]["beta"] for s in common]
cs_abs_beta = [abs(v) for v in cs_beta]
rb_abs_beta = [abs(v) for v in rb_beta]
cs_se = [cs[s]["se"] for s in common]
rb_se = [rb[s]["se"] for s in common]
cs_logp = [cs[s]["logp"] for s in common]
rb_logp = [rb[s]["logp"] for s in common]

beta_corr = corr(cs_beta, rb_beta)
beta_corr_flip = corr([-v for v in cs_beta], rb_beta)
rows = {
    "n_subjects": "${N}",
    "records_per_subject": "${R}",
    "n_snps_requested": "${M}",
    "cs_snps": len(cs),
    "repeatabel_snps": len(rb),
    "overlap": len(common),
    "beta_corr_signed": beta_corr,
    "beta_corr_best_orientation": max(beta_corr, beta_corr_flip),
    "abs_beta_corr": corr(cs_abs_beta, rb_abs_beta),
    "se_corr": corr(cs_se, rb_se),
    "log10p_corr": corr(cs_logp, rb_logp),
    "top10_overlap": top_overlap(cs, rb, min(10, len(common))),
    "top20_overlap": top_overlap(cs, rb, min(20, len(common))),
    "max_abs_log10p_diff": max(abs(cs[s]["logp"] - rb[s]["logp"]) for s in common),
    "median_abs_log10p_diff": statistics.median(abs(cs[s]["logp"] - rb[s]["logp"]) for s in common),
}
with summary_file.open("w", encoding="utf-8") as f:
    f.write("metric\tvalue\n")
    for k, v in rows.items():
        f.write(f"{k}\t{v}\n")
print(summary_file.read_text(encoding="utf-8"))
PY
