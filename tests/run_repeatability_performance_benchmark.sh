#!/bin/bash
set -euo pipefail

export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

ROOT="/disk228/huangjie/PanOmics/CosmicScanner"
BIN="${ROOT}/build/cosmicscanner"
DATA_ROOT=${BENCH_DATA_ROOT:-${ROOT}/tests/data/repeatability_perf_bench}
RESULT_ROOT=${BENCH_RESULT_ROOT:-${ROOT}/tests/results/repeatability_perf_bench}
SUMMARY="${RESULT_ROOT}/benchmark.summary.tsv"

BENCH_N_LIST=${BENCH_N_LIST:-200,500}
BENCH_M_LIST=${BENCH_M_LIST:-1000,5000}
BENCH_R_LIST=${BENCH_R_LIST:-2,5}
BENCH_METHODS=${BENCH_METHODS:-exact,fast,mean,repeatabel}
BENCH_REPEATABL_RSCRIPT=${BENCH_REPEATABL_RSCRIPT:-}
BENCH_REPEATABL_R_LIB=${BENCH_REPEATABL_R_LIB:-/tmp/opencode/Rlib_repeatabel}
BENCH_TIMEOUT_SEC=${BENCH_TIMEOUT_SEC:-0}
BENCH_REML_MAXIT=${BENCH_REML_MAXIT:-0}

mkdir -p "${DATA_ROOT}" "${RESULT_ROOT}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing cosmicscanner binary: ${BIN}" >&2
  exit 1
fi

detect_repeatabel() {
  if [[ -n "${BENCH_REPEATABL_RSCRIPT}" ]]; then
    "${BENCH_REPEATABL_RSCRIPT}" -e ".libPaths(c('${BENCH_REPEATABL_R_LIB}', .libPaths())); quit(status=as.integer(!requireNamespace('RepeatABEL', quietly=TRUE)))" >/dev/null 2>&1 && return 0
    return 1
  fi
  if command -v Rscript >/dev/null 2>&1; then
    Rscript -e ".libPaths(c('${BENCH_REPEATABL_R_LIB}', .libPaths())); quit(status=as.integer(!requireNamespace('RepeatABEL', quietly=TRUE)))" >/dev/null 2>&1 && return 0
  fi
  if command -v conda >/dev/null 2>&1; then
    conda run -n R_GS Rscript -e ".libPaths(c('${BENCH_REPEATABL_R_LIB}', .libPaths())); quit(status=as.integer(!requireNamespace('RepeatABEL', quietly=TRUE)))" >/dev/null 2>&1 && return 0
  fi
  return 1
}

REPEATABL_AVAILABLE=0
if detect_repeatabel; then
  REPEATABL_AVAILABLE=1
fi

if [[ ! -s "${SUMMARY}" ]]; then
  printf "n_subjects\tn_records_per_subject\tn_records\tn_snps\tmethod\tstatus\truntime_sec\tmax_rss_mb\tscan_mode\tn_output_snps\tnote\n" > "${SUMMARY}"
fi

contains_method() {
  local needle=$1
  IFS=',' read -ra methods <<< "${BENCH_METHODS}"
  for method in "${methods[@]}"; do
    [[ "${method}" == "${needle}" ]] && return 0
  done
  return 1
}

generate_dataset() {
  local n=$1 m=$2 r=$3 prefix=$4
  mkdir -p "$(dirname "${prefix}")"
  python3 - <<PY
from pathlib import Path
import numpy as np

n = int("${n}")
m = int("${m}")
r = int("${r}")
prefix = Path("${prefix}")
rng = np.random.default_rng(1000003 + n * 17 + m * 3 + r)

maf = rng.uniform(0.05, 0.45, size=m)
G = rng.binomial(2, maf, size=(n, m)).astype(float)
missing = rng.random((n, m)) < 0.005
G[missing] = np.nan
G_imp = G.copy()
for j in range(m):
    col = G_imp[:, j]
    mean = np.nanmean(col)
    if not np.isfinite(mean):
        mean = 2 * maf[j]
    col[np.isnan(col)] = mean
    G_imp[:, j] = col

with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
    for i in range(n):
        fid = iid = f"ID{i+1:06d}"
        f.write(f"{fid}\t{iid}\t0\t0\t0\t-9\n")

with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
    for j in range(m):
        chrom = 1 + (j % 5)
        f.write(f"{chrom}\trs{j+1:08d}\t0\t{1000 + 10*j}\tA\tG\n")

bed = bytearray(b"\x6c\x1b\x01")
bytes_per_snp = (n + 3) // 4
for j in range(m):
    block = bytearray(bytes_per_snp)
    for i in range(n):
        val = G[i, j]
        if np.isnan(val):
            code = 0b01
        elif val <= 0.5:
            code = 0b00
        elif val <= 1.5:
            code = 0b10
        else:
            code = 0b11
        block[i // 4] |= code << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix(".bed").write_bytes(bed)

causal = [min(7, m - 1), min(113, m - 1), min(503, m - 1)]
causal = sorted(set(causal))
effects = [0.45, -0.35, 0.30]
additive = np.zeros(n)
for eff, j in zip(effects, causal):
    additive += eff * (G_imp[:, j] - np.mean(G_imp[:, j]))
polygenic = rng.normal(0, 0.45, size=n)
pe = rng.normal(0, 0.50, size=n)

with (prefix.parent / "long.pheno").open("w", encoding="utf-8") as pf, \
     (prefix.parent / "record.qcovar").open("w", encoding="utf-8") as cf, \
     (prefix.parent / "mean.pheno").open("w", encoding="utf-8") as mf:
    for i in range(n):
        fid = iid = f"ID{i+1:06d}"
        ys = []
        for t in range(1, r + 1):
            batch = 1.0 if t % 2 == 0 else 0.0
            y = additive[i] + polygenic[i] + pe[i] + 0.25 * batch + rng.normal(0, 0.75)
            pf.write(f"{fid}\t{iid}\t{t}\t{y:.8f}\n")
            cf.write(f"{fid}\t{iid}\t{t}\t{batch:.1f}\n")
            ys.append(y)
        mf.write(f"{fid}\t{iid}\t{float(np.mean(ys)):.8f}\n")
PY
}

parse_seconds() {
  local time_file=$1
  python3 - <<PY
from pathlib import Path
import re
text = Path("${time_file}").read_text(errors="replace")
m = re.search(r"Elapsed \(wall clock\) time.*?: (?:(\d+):)?(\d+):(\d+(?:\.\d+)?)", text)
if not m:
    print("NA")
else:
    h = int(m.group(1) or 0)
    mi = int(m.group(2))
    sec = float(m.group(3))
    print(f"{h*3600 + mi*60 + sec:.3f}")
PY
}

parse_rss_mb() {
  local time_file=$1
  python3 - <<PY
from pathlib import Path
import re
text = Path("${time_file}").read_text(errors="replace")
m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
print(f"{int(m.group(1))/1024:.2f}" if m else "NA")
PY
}

scan_mode() {
  local null_file=$1
  python3 - <<PY
from pathlib import Path
import csv
p = Path("${null_file}")
if not p.exists():
    print("NA")
else:
    out = "NA"
    with p.open("r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("Item") == "scan_mode":
                out = row.get("Value", "NA")
                break
    print(out)
PY
}

count_assoc_rows() {
  local assoc_file=$1
  python3 - <<PY
from pathlib import Path
p = Path("${assoc_file}")
if not p.exists():
    print(0)
else:
    with p.open("r", encoding="utf-8") as f:
        print(max(0, sum(1 for line in f if line.strip()) - 1))
PY
}

append_summary() {
  local n=$1 r=$2 m=$3 method=$4 status=$5 runtime=$6 rss=$7 mode=$8 rows=$9 note=${10}
  local nrec=$((n * r))
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${n}" "${r}" "${nrec}" "${m}" "${method}" "${status}" "${runtime}" "${rss}" "${mode}" "${rows}" "${note}" >> "${SUMMARY}"
}

run_cs() {
  local n=$1 m=$2 r=$3 method=$4 data_dir=$5 out=$6
  local time_file="${out}.time.txt"
  local status="ok"
  local rc=0
  local reml_args=()
  local note="NA"
  if [[ "${BENCH_REML_MAXIT}" -gt 0 ]]; then
    reml_args=(--reml-maxit "${BENCH_REML_MAXIT}")
    note="reml_maxit=${BENCH_REML_MAXIT}"
  fi
  if [[ "${method}" == "exact" ]]; then
    if [[ "${BENCH_TIMEOUT_SEC}" -gt 0 ]]; then
      /usr/bin/time -v -o "${time_file}" timeout "${BENCH_TIMEOUT_SEC}" "${BIN}" --bfile "${data_dir}/sim" \
        --pheno "${data_dir}/long.pheno" --tpheno 1 --mpheno 2 --qcovar "${data_dir}/record.qcovar" \
        --repeatability-model --repeatability-mode exact --gwas --no-loco --maf 0.01 "${reml_args[@]}" \
        --out "${out}" >"${out}.log" 2>&1 || rc=$?
    else
      /usr/bin/time -v -o "${time_file}" "${BIN}" --bfile "${data_dir}/sim" \
        --pheno "${data_dir}/long.pheno" --tpheno 1 --mpheno 2 --qcovar "${data_dir}/record.qcovar" \
        --repeatability-model --repeatability-mode exact --gwas --no-loco --maf 0.01 "${reml_args[@]}" \
        --out "${out}" >"${out}.log" 2>&1 || rc=$?
    fi
  elif [[ "${method}" == "fast" ]]; then
    if [[ "${BENCH_TIMEOUT_SEC}" -gt 0 ]]; then
      /usr/bin/time -v -o "${time_file}" timeout "${BENCH_TIMEOUT_SEC}" "${BIN}" --bfile "${data_dir}/sim" \
        --pheno "${data_dir}/long.pheno" --tpheno 1 --mpheno 2 --qcovar "${data_dir}/record.qcovar" \
        --repeatability-model --repeatability-mode fast --gwas --no-loco --maf 0.01 "${reml_args[@]}" \
        --out "${out}" >"${out}.log" 2>&1 || rc=$?
    else
      /usr/bin/time -v -o "${time_file}" "${BIN}" --bfile "${data_dir}/sim" \
        --pheno "${data_dir}/long.pheno" --tpheno 1 --mpheno 2 --qcovar "${data_dir}/record.qcovar" \
        --repeatability-model --repeatability-mode fast --gwas --no-loco --maf 0.01 "${reml_args[@]}" \
        --out "${out}" >"${out}.log" 2>&1 || rc=$?
    fi
  else
    if [[ "${BENCH_TIMEOUT_SEC}" -gt 0 ]]; then
      /usr/bin/time -v -o "${time_file}" timeout "${BENCH_TIMEOUT_SEC}" "${BIN}" --bfile "${data_dir}/sim" \
        --pheno "${data_dir}/mean.pheno" --mpheno 1 --gwas --no-loco --maf 0.01 \
        --out "${out}" >"${out}.log" 2>&1 || rc=$?
    else
      /usr/bin/time -v -o "${time_file}" "${BIN}" --bfile "${data_dir}/sim" \
        --pheno "${data_dir}/mean.pheno" --mpheno 1 --gwas --no-loco --maf 0.01 \
        --out "${out}" >"${out}.log" 2>&1 || rc=$?
    fi
  fi
  if [[ "${rc}" -eq 124 ]]; then
    status="timeout"
  elif [[ "${rc}" -ne 0 ]]; then
    status="failed"
  fi
  local runtime rss mode rows
  runtime=$(parse_seconds "${time_file}")
  rss=$(parse_rss_mb "${time_file}")
  if [[ "${method}" == "mean" ]]; then
    mode="single_trait_lmm_mean"
  else
    mode=$(scan_mode "${out}.repeatability.null.txt")
  fi
  rows=$(count_assoc_rows "${out}.assoc.txt")
  append_summary "${n}" "${r}" "${m}" "${method}" "${status}" "${runtime}" "${rss}" "${mode}" "${rows}" "${note}"
}

run_repeatabel() {
  local n=$1 m=$2 r=$3 data_dir=$4 out=$5
  local rb_dir="${out}.repeatabel_input"
  local r_script="${out}.repeatabel.R"
  local time_file="${out}.time.txt"
  mkdir -p "${rb_dir}"
  python3 - <<PY
from pathlib import Path

data = Path("${data_dir}")
rb = Path("${rb_dir}")
prefix = data / "sim"
fam = [line.split() for line in prefix.with_suffix(".fam").read_text().splitlines() if line.strip()]
bim = [line.split() for line in prefix.with_suffix(".bim").read_text().splitlines() if line.strip()]
n = len(fam)
m = len(bim)
bed = prefix.with_suffix(".bed").read_bytes()
if bed[:3] != b"\x6c\x1b\x01":
    raise SystemExit("Unexpected PLINK bed header")
bytes_per_snp = (n + 3) // 4
with (rb / "genos.txt").open("w", encoding="utf-8") as f:
    f.write(" ".join(row[1] for row in fam) + "\n")
    f.write(" ".join(row[1] for row in bim) + "\n")
    f.write(" ".join(row[0] for row in bim) + "\n")
    f.write(" ".join(row[3] for row in bim) + "\n")
    for j in range(m):
        block = bed[3 + j * bytes_per_snp: 3 + (j + 1) * bytes_per_snp]
        vals = []
        for i in range(n):
            code = (block[i // 4] >> ((i % 4) * 2)) & 0b11
            vals.append({0b00: "1", 0b10: "2", 0b11: "3", 0b01: "0"}[code])
        f.write(" ".join(vals) + "\n")
with (rb / "pheno_ind.dat").open("w", encoding="utf-8") as f:
    f.write("id sex\n")
    for row in fam:
        f.write(f'"{row[1]}" 0\n')
qcov = {}
with (data / "record.qcovar").open("r", encoding="utf-8") as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 4:
            qcov[(parts[1], parts[2])] = parts[3]
with (data / "long.pheno").open("r", encoding="utf-8") as inp, (rb / "phen_repeat.dat").open("w", encoding="utf-8") as outp:
    outp.write("id time y batch\n")
    for line in inp:
        parts = line.split()
        if len(parts) >= 4:
            batch = qcov.get((parts[1], parts[2]), "0")
            outp.write(f'"{parts[1]}" {parts[2]} {parts[3]} {batch}\n')
PY
  python3 - <<PY
from pathlib import Path
p = Path("${r_script}")
p.write_text(r'''
.libPaths(c(Sys.getenv("REPEATABL_R_LIB"), .libPaths()))
suppressPackageStartupMessages(library(RepeatABEL))
suppressPackageStartupMessages(library(GenABEL))
rb_dir <- Sys.getenv("REPEATABL_INPUT_DIR")
out_prefix <- Sys.getenv("REPEATABL_OUT_PREFIX")
geno_txt <- file.path(rb_dir, "genos.txt")
geno_raw <- file.path(rb_dir, "genos.raw")
if (!file.exists(geno_raw)) convert.snp.text(geno_txt, geno_raw)
gd <- load.gwaa.data(
  phenofile=file.path(rb_dir, "pheno_ind.dat"),
  genofile=geno_raw,
  force=TRUE,
  sort=FALSE,
  id="id"
)
phen <- read.table(file.path(rb_dir, "phen_repeat.dat"), header=TRUE)
fit <- rGLS(y ~ batch, genabel.data=gd, phenotype.data=phen, id.name="id")
sink(paste0(out_prefix, ".repeatabel.summary.txt"))
print(summary(fit))
if (!is.null(fit@call\$hglm)) print(summary(fit@call\$hglm))
sink()
write.table(data.frame(status="ok", n_snps=ncol(gtdata(gd))), paste0(out_prefix, ".repeatabel.status.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
''', encoding="utf-8")
PY
  local rc=0
  if [[ -n "${BENCH_REPEATABL_RSCRIPT}" ]]; then
    REPEATABL_R_LIB="${BENCH_REPEATABL_R_LIB}" REPEATABL_INPUT_DIR="${rb_dir}" REPEATABL_OUT_PREFIX="${out}" \
      /usr/bin/time -v -o "${time_file}" timeout "${BENCH_TIMEOUT_SEC:-0}" "${BENCH_REPEATABL_RSCRIPT}" "${r_script}" >"${out}.log" 2>&1 || rc=$?
  elif command -v Rscript >/dev/null 2>&1; then
    REPEATABL_R_LIB="${BENCH_REPEATABL_R_LIB}" REPEATABL_INPUT_DIR="${rb_dir}" REPEATABL_OUT_PREFIX="${out}" \
      /usr/bin/time -v -o "${time_file}" timeout "${BENCH_TIMEOUT_SEC:-0}" Rscript "${r_script}" >"${out}.log" 2>&1 || rc=$?
  else
    REPEATABL_R_LIB="${BENCH_REPEATABL_R_LIB}" REPEATABL_INPUT_DIR="${rb_dir}" REPEATABL_OUT_PREFIX="${out}" \
      /usr/bin/time -v -o "${time_file}" timeout "${BENCH_TIMEOUT_SEC:-0}" conda run -n R_GS Rscript "${r_script}" >"${out}.log" 2>&1 || rc=$?
  fi
  local status="ok"
  [[ "${rc}" -eq 124 ]] && status="timeout"
  [[ "${rc}" -ne 0 && "${rc}" -ne 124 ]] && status="failed"
  local runtime rss rows
  runtime=$(parse_seconds "${time_file}")
  rss=$(parse_rss_mb "${time_file}")
  rows=0
  [[ -s "${out}.repeatabel.status.tsv" ]] && rows="${m}"
  append_summary "${n}" "${r}" "${m}" "repeatabel" "${status}" "${runtime}" "${rss}" "RepeatABEL_rGLS" "${rows}" "Rlib=${BENCH_REPEATABL_R_LIB}"
}

IFS=',' read -ra n_values <<< "${BENCH_N_LIST}"
IFS=',' read -ra m_values <<< "${BENCH_M_LIST}"
IFS=',' read -ra r_values <<< "${BENCH_R_LIST}"

for n in "${n_values[@]}"; do
  for m in "${m_values[@]}"; do
    for r in "${r_values[@]}"; do
      case_dir="${DATA_ROOT}/n${n}_m${m}_r${r}"
      prefix="${case_dir}/sim"
      if [[ ! -s "${prefix}.bed" || ! -s "${prefix}.bim" || ! -s "${prefix}.fam" ]]; then
        generate_dataset "${n}" "${m}" "${r}" "${prefix}"
      fi
      if contains_method exact; then
        run_cs "${n}" "${m}" "${r}" exact "${case_dir}" "${RESULT_ROOT}/n${n}_m${m}_r${r}.exact"
      fi
      if contains_method fast; then
        run_cs "${n}" "${m}" "${r}" fast "${case_dir}" "${RESULT_ROOT}/n${n}_m${m}_r${r}.fast"
      fi
      if contains_method mean; then
        run_cs "${n}" "${m}" "${r}" mean "${case_dir}" "${RESULT_ROOT}/n${n}_m${m}_r${r}.mean"
      fi
      if contains_method repeatabel; then
        if [[ "${REPEATABL_AVAILABLE}" -eq 1 ]]; then
          run_repeatabel "${n}" "${m}" "${r}" "${case_dir}" "${RESULT_ROOT}/n${n}_m${m}_r${r}.repeatabel"
        else
          append_summary "${n}" "${r}" "${m}" "repeatabel" "skipped" "NA" "NA" "NA" "0" "RepeatABEL package unavailable"
        fi
      fi
    done
  done
done

cat "${SUMMARY}"
