#!/usr/bin/env bash
# Large-scale biobank-mode stress/accuracy test for CosmicScanner.
#
# Purpose:
#   1. Verify that large-sample optimizations still produce consistent results.
#   2. Measure runtime and peak memory for PLINK/PGEN/BGEN biobank-mode paths.
#   3. Optionally compare against REGENIE on tractable settings.
#
# Recommended tiers:
#   Quick correctness:
#     N=10000 M=50000 MODEL_M=5000 RUN_PGEN=1 RUN_BGEN=1 RUN_REGENIE=1 bash tests/run_biobank_stress.sh
#   Medium stress:
#     N=50000 M=200000 MODEL_M=10000 RUN_PGEN=1 RUN_BGEN=0 RUN_REGENIE=0 bash tests/run_biobank_stress.sh
#   Large stress:
#     N=100000 M=1000000 MODEL_M=20000 RUN_PGEN=1 RUN_BGEN=0 RUN_REGENIE=0 bash tests/run_biobank_stress.sh
#
# Notes:
#   - BGEN conversion can be slow and space-heavy at large M; keep RUN_BGEN=0 for first large tests.
#   - REGENIE comparison is intended for small/medium tiers only.
#   - This script writes under tests/data/biobank_stress_* and tests/results/biobank_stress_*.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
PLINK2="${PLINK2_BIN:-plink2}"
REGENIE="${REGENIE_BIN:-regenie}"
PLINK2="$(command -v "${PLINK2}" || printf '%s' "${PLINK2}")"
REGENIE="$(command -v "${REGENIE}" || printf '%s' "${REGENIE}")"

N="${N:-10000}"
M="${M:-50000}"
MODEL_M="${MODEL_M:-5000}"
SEED="${SEED:-20260601}"
CHROMS="${CHROMS:-1}"
BLOCK_SIZE="${BLOCK_SIZE:-1000}"
RUN_PLINK="${RUN_PLINK:-1}"
RUN_PGEN="${RUN_PGEN:-1}"
RUN_BGEN="${RUN_BGEN:-0}"
RUN_REGENIE="${RUN_REGENIE:-0}"
GENERATE_DATA="${GENERATE_DATA:-1}"
CONVERT_FORMATS="${CONVERT_FORMATS:-1}"
KEEP_DATA="${KEEP_DATA:-1}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

TAG="n${N}_m${M}_model${MODEL_M}_b${BLOCK_SIZE}"
DATA_DIR="${ROOT}/tests/data/biobank_stress_${TAG}"
RESULT_DIR="${ROOT}/tests/results/biobank_stress_${TAG}"
DATA_PREFIX="${DATA_DIR}/stress"
MODEL_SNPS="${DATA_DIR}/model_snps.txt"
COVAR_FILE="${DATA_DIR}/covar.txt"
SUMMARY="${RESULT_DIR}/stress_summary.tsv"
COMPARE="${RESULT_DIR}/cross_format_compare.tsv"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

for cmd_var in BIN PLINK2; do
  cmd_path="${!cmd_var}"
  if [[ ! -x "${cmd_path}" ]]; then
    echo "Missing ${cmd_var}: ${cmd_path}" >&2
    exit 2
  fi
done
if [[ "${RUN_REGENIE}" == "1" && ! -x "${REGENIE}" ]]; then
  echo "Missing REGENIE: ${REGENIE}" >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}" "${RESULT_DIR}/cosmic_plink" \
         "${RESULT_DIR}/cosmic_pgen" "${RESULT_DIR}/cosmic_bgen" "${RESULT_DIR}/regenie"

write_summary_header() {
  printf 'path\tn\tm\tmodel_m\tblock_size\tseconds\tmax_rss_kb\tassoc_rows\tlog\n' > "${SUMMARY}"
}

count_assoc_rows() {
  local file="$1"
  if [[ ! -f "${file}" ]]; then
    printf '0'
    return
  fi
  "${PYTHON_BIN}" - "${file}" <<'PY'
import sys
from pathlib import Path
p = Path(sys.argv[1])
print(max(0, sum(1 for _ in p.open()) - 1))
PY
}

append_summary_row() {
  local label="$1"
  local seconds="$2"
  local rss="$3"
  local rows="$4"
  local log="$5"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${label}" "${N}" "${M}" "${MODEL_M}" "${BLOCK_SIZE}" \
    "${seconds}" "${rss}" "${rows}" "${log}" >> "${SUMMARY}"
}

run_timed() {
  local label="$1"
  local out_prefix="$2"
  local log="$3"
  shift 3
  local metric="${log}.time"
  echo "[biobank-stress] Running ${label}..."
  /usr/bin/time -f 'seconds\t%e\nmax_rss_kb\t%M' -o "${metric}" "$@" >"${log}" 2>&1
  local seconds rss rows
  seconds="$(awk '$1=="seconds"{print $2}' "${metric}")"
  rss="$(awk '$1=="max_rss_kb"{print $2}' "${metric}")"
  rows="$(count_assoc_rows "${out_prefix}.biobank.assoc.txt")"
  append_summary_row "${label}" "${seconds}" "${rss}" "${rows}" "${log}"
}

if [[ "${GENERATE_DATA}" == "1" ]]; then
  echo "[biobank-stress] Generating PLINK test data: N=${N}, M=${M}, MODEL_M=${MODEL_M}, seed=${SEED}"
  "${PYTHON_BIN}" - "${DATA_PREFIX}" "${N}" "${M}" "${MODEL_M}" "${SEED}" "${CHROMS}" <<'PY'
import sys
from pathlib import Path

try:
    import numpy as np
except Exception as exc:
    raise SystemExit(f"Need numpy for stress data generation. Set PYTHON_BIN to a Python with numpy. Error: {exc}")

prefix = Path(sys.argv[1])
n = int(sys.argv[2])
m = int(sys.argv[3])
model_m = int(sys.argv[4])
seed = int(sys.argv[5])
chroms = max(1, int(sys.argv[6]))
rng = np.random.default_rng(seed)

prefix.parent.mkdir(parents=True, exist_ok=True)
bytes_per_snp = (n + 3) // 4
pad_n = bytes_per_snp * 4
chunk = 512

pop = np.arange(n, dtype=np.int32) % 6
pop_center = (pop.astype(np.float32) - pop.mean()).reshape(1, -1)
age = rng.normal(0, 1, size=n).astype(np.float64)
sex = (np.arange(n) % 2).astype(np.float64)
phen = 0.08 * age + 0.12 * sex + 0.05 * pop.astype(np.float64) + rng.normal(0, 1.0, size=n)

causal = {}
for frac, eff in [(0.03, 0.18), (0.11, -0.16), (0.29, 0.14), (0.51, -0.12), (0.77, 0.10)]:
    idx = min(m - 1, max(0, int(frac * m)))
    causal[idx] = eff

with prefix.with_suffix('.bim').open('w', encoding='utf-8') as f:
    for j in range(m):
        chrom = (j % chroms) + 1
        f.write(f"{chrom}\trs{j+1:09d}\t0\t{1000 + j * 7}\tA\tC\n")

with prefix.with_suffix('.bed').open('wb') as bed:
    bed.write(b'\x6c\x1b\x01')
    for start in range(0, m, chunk):
        end = min(m, start + chunk)
        idx = np.arange(start, end, dtype=np.float64)
        base_p = 0.03 + 0.47 * ((idx.astype(np.int64) * 7919) % 1000) / 999.0
        # Mild population structure + sinusoidal genome-wide drift.
        drift = (0.015 * np.sin(idx * 0.017)).astype(np.float32).reshape(-1, 1)
        p = np.clip(base_p.astype(np.float32).reshape(-1, 1) + drift * pop_center, 0.005, 0.995)
        g = (rng.random((end - start, n), dtype=np.float32) < p).astype(np.uint8)
        g += (rng.random((end - start, n), dtype=np.float32) < p).astype(np.uint8)

        for snp_idx, eff in causal.items():
            if start <= snp_idx < end:
                local = snp_idx - start
                phen += eff * (g[local].astype(np.float64) - 2.0 * p[local].astype(np.float64))

        codes = np.empty((end - start, pad_n), dtype=np.uint8)
        codes[:, :n] = np.where(g == 2, 0, np.where(g == 1, 2, 3)).astype(np.uint8)
        if pad_n > n:
            codes[:, n:] = 0
        packed = (codes[:, 0::4]
                  | (codes[:, 1::4] << 2)
                  | (codes[:, 2::4] << 4)
                  | (codes[:, 3::4] << 6)).astype(np.uint8)
        bed.write(packed.tobytes(order='C'))
        if (start // chunk + 1) % 50 == 0:
            print(f"  generated {end}/{m} SNPs", flush=True)

with prefix.with_suffix('.fam').open('w', encoding='utf-8') as f:
    for i, y in enumerate(phen):
        f.write(f"F{i+1:08d}\tI{i+1:08d}\t0\t0\t{int(sex[i]) + 1}\t{float(y):.8f}\n")

with prefix.with_suffix('.pheno').open('w', encoding='utf-8') as f:
    f.write('FID\tIID\tPHENO\n')
    for i, y in enumerate(phen):
        f.write(f"F{i+1:08d}\tI{i+1:08d}\t{float(y):.8f}\n")

with (prefix.parent / 'covar.txt').open('w', encoding='utf-8') as f:
    f.write('FID\tIID\tAGE\tSEX\tPOP1\tPOP2\tPOP3\tPOP4\tPOP5\n')
    for i in range(n):
        pcs = [1 if pop[i] == k else 0 for k in range(1, 6)]
        f.write(f"F{i+1:08d}\tI{i+1:08d}\t{age[i]:.6f}\t{int(sex[i])}\t" + '\t'.join(map(str, pcs)) + '\n')

with (prefix.parent / 'model_snps.txt').open('w', encoding='utf-8') as f:
    f.write('SNP\n')
    for j in range(min(model_m, m)):
        f.write(f"rs{j+1:09d}\n")

with prefix.with_suffix('.truth.txt').open('w', encoding='utf-8') as f:
    f.write('SNP\tTRUE_BETA\n')
    for idx, eff in sorted(causal.items()):
        f.write(f"rs{idx+1:09d}\t{eff}\n")
PY
fi

if [[ "${CONVERT_FORMATS}" == "1" ]]; then
  if [[ "${RUN_PGEN}" == "1" ]]; then
    echo "[biobank-stress] Converting PLINK -> PGEN..."
    "${PLINK2}" --bfile "${DATA_PREFIX}" --make-pgen --out "${RESULT_DIR}/cosmic_pgen/stress" \
      >"${RESULT_DIR}/plink2_pgen.log" 2>&1
  fi
  if [[ "${RUN_BGEN}" == "1" ]]; then
    echo "[biobank-stress] Converting PLINK -> BGEN..."
    "${PLINK2}" --bfile "${DATA_PREFIX}" --export bgen-1.2 --out "${RESULT_DIR}/cosmic_bgen/stress" \
      >"${RESULT_DIR}/plink2_bgen.log" 2>&1
  fi
fi

write_summary_header

if [[ "${RUN_REGENIE}" == "1" ]]; then
  echo "[biobank-stress] Running REGENIE baseline..."
  /usr/bin/time -f 'seconds\t%e\nmax_rss_kb\t%M' -o "${RESULT_DIR}/regenie/step1.log.time" \
    "${REGENIE}" --step 1 --bed "${DATA_PREFIX}" --phenoFile "${DATA_PREFIX}.pheno" \
    --covarFile "${COVAR_FILE}" --bsize "${BLOCK_SIZE}" --extract "${MODEL_SNPS}" \
    --pred "${RESULT_DIR}/regenie/step1_pred" --out "${RESULT_DIR}/regenie/step1" \
    >"${RESULT_DIR}/regenie/step1.log" 2>&1
  /usr/bin/time -f 'seconds\t%e\nmax_rss_kb\t%M' -o "${RESULT_DIR}/regenie/step2.log.time" \
    "${REGENIE}" --step 2 --bed "${DATA_PREFIX}" --phenoFile "${DATA_PREFIX}.pheno" \
    --covarFile "${COVAR_FILE}" --bsize "${BLOCK_SIZE}" \
    --pred "${RESULT_DIR}/regenie/step1_pred.list" --out "${RESULT_DIR}/regenie/step2" \
    >"${RESULT_DIR}/regenie/step2.log" 2>&1
  regenie_step1_s="$(awk '$1=="seconds"{print $2}' "${RESULT_DIR}/regenie/step1.log.time")"
  regenie_step1_rss="$(awk '$1=="max_rss_kb"{print $2}' "${RESULT_DIR}/regenie/step1.log.time")"
  regenie_step2_s="$(awk '$1=="seconds"{print $2}' "${RESULT_DIR}/regenie/step2.log.time")"
  regenie_step2_rss="$(awk '$1=="max_rss_kb"{print $2}' "${RESULT_DIR}/regenie/step2.log.time")"
  regenie_files=("${RESULT_DIR}"/regenie/step2_*.regenie)
  regenie_rows=0
  if [[ -f "${regenie_files[0]}" ]]; then
    regenie_rows="$(count_assoc_rows "${regenie_files[0]}")"
  fi
  regenie_total_s="$(awk -v a="${regenie_step1_s}" -v b="${regenie_step2_s}" 'BEGIN{printf "%.2f", a+b}')"
  regenie_total_rss="$(awk -v a="${regenie_step1_rss}" -v b="${regenie_step2_rss}" 'BEGIN{print (a>b)?a:b}')"
  append_summary_row "REGENIE_STEP1" "${regenie_step1_s}" "${regenie_step1_rss}" 0 "${RESULT_DIR}/regenie/step1.log"
  append_summary_row "REGENIE_STEP2" "${regenie_step2_s}" "${regenie_step2_rss}" "${regenie_rows}" "${RESULT_DIR}/regenie/step2.log"
  append_summary_row "REGENIE_TOTAL" "${regenie_total_s}" "${regenie_total_rss}" "${regenie_rows}" "${RESULT_DIR}/regenie"
fi

COVAR_ARGS="--qcovar ${COVAR_FILE}"

if [[ "${RUN_PLINK}" == "1" ]]; then
  run_timed "Cosmic_PLINK" "${RESULT_DIR}/cosmic_plink/stress" "${RESULT_DIR}/cosmic_plink/stress.log" \
    "${BIN}" --bfile "${DATA_PREFIX}" --pheno "${DATA_PREFIX}.pheno" --mpheno 1 --gwas \
    --biobank-mode --block-size "${BLOCK_SIZE}" --regenie-model-snps "${MODEL_SNPS}" \
    ${COVAR_ARGS} --maf 0.0 --out "${RESULT_DIR}/cosmic_plink/stress"
fi

if [[ "${RUN_PGEN}" == "1" ]]; then
  run_timed "Cosmic_PGEN" "${RESULT_DIR}/cosmic_pgen/stress" "${RESULT_DIR}/cosmic_pgen/stress.log" \
    "${BIN}" --pfile "${RESULT_DIR}/cosmic_pgen/stress" --pheno "${DATA_PREFIX}.pheno" --mpheno 1 --gwas \
    --biobank-mode --block-size "${BLOCK_SIZE}" --regenie-model-snps "${MODEL_SNPS}" \
    ${COVAR_ARGS} --maf 0.0 --out "${RESULT_DIR}/cosmic_pgen/stress"
fi

if [[ "${RUN_BGEN}" == "1" ]]; then
  run_timed "Cosmic_BGEN" "${RESULT_DIR}/cosmic_bgen/stress" "${RESULT_DIR}/cosmic_bgen/stress.log" \
    "${BIN}" --bgen "${RESULT_DIR}/cosmic_bgen/stress.bgen" --sample "${RESULT_DIR}/cosmic_bgen/stress.sample" \
    --pheno "${DATA_PREFIX}.pheno" --mpheno 1 --gwas --biobank-mode --block-size "${BLOCK_SIZE}" \
    --regenie-model-snps "${MODEL_SNPS}" ${COVAR_ARGS} \
    --maf 0.0 --out "${RESULT_DIR}/cosmic_bgen/stress"
fi

echo "[biobank-stress] Comparing available outputs..."
"${PYTHON_BIN}" - "${RESULT_DIR}" "${COMPARE}" <<'PY'
import csv, glob, math, sys
from pathlib import Path

res_dir = Path(sys.argv[1])
out_path = Path(sys.argv[2])

paths = {
    'Cosmic_PLINK': res_dir / 'cosmic_plink' / 'stress.biobank.assoc.txt',
    'Cosmic_PGEN': res_dir / 'cosmic_pgen' / 'stress.biobank.assoc.txt',
    'Cosmic_BGEN': res_dir / 'cosmic_bgen' / 'stress.biobank.assoc.txt',
}
regenie = sorted((res_dir / 'regenie').glob('step2_*.regenie'))
if regenie:
    paths['REGENIE'] = regenie[0]

def read_assoc(path):
    if not path.exists():
        return {}
    first = path.open().readline().rstrip('\n')
    rows = {}
    if '\t' in first:
        with path.open() as f:
            reader = csv.DictReader(f, delimiter='\t')
            for r in reader:
                snp = r.get('SNP') or r.get('ID') or ''
                if not snp:
                    continue
                rows[snp] = {
                    'BETA': float(r.get('BETA', 'nan') or 'nan'),
                    'SE': float(r.get('SE', 'nan') or 'nan'),
                    'P': float(r.get('P', 'nan') or 'nan'),
                }
    else:
        with path.open() as f:
            header = f.readline().split()
            for line in f:
                parts = line.split()
                if len(parts) != len(header):
                    continue
                r = dict(zip(header, parts))
                snp = r.get('ID') or r.get('SNP') or ''
                if not snp:
                    continue
                log10p = float(r.get('LOG10P', 'nan') or 'nan')
                p = 10.0 ** (-log10p) if math.isfinite(log10p) else float('nan')
                rows[snp] = {
                    'BETA': float(r.get('BETA', 'nan') or 'nan'),
                    'SE': float(r.get('SE', 'nan') or 'nan'),
                    'P': p,
                }
    return rows

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

tables = {k: read_assoc(v) for k, v in paths.items()}
labels = [k for k, v in tables.items() if v]

with out_path.open('w', encoding='utf-8') as out:
    out.write('a\tb\toverlap\tbeta_corr\tse_corr\tlog10p_corr\tmean_abs_beta\tmax_abs_beta\tmean_abs_se\tmax_abs_se\tmean_abs_p\tmax_abs_p\n')
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            common = sorted(set(tables[a]) & set(tables[b]))
            bet_a = [tables[a][s]['BETA'] for s in common]
            bet_b = [tables[b][s]['BETA'] for s in common]
            se_a = [tables[a][s]['SE'] for s in common]
            se_b = [tables[b][s]['SE'] for s in common]
            p_a = [tables[a][s]['P'] for s in common]
            p_b = [tables[b][s]['P'] for s in common]
            lp_a = [-math.log10(max(p, 1e-300)) if math.isfinite(p) else math.nan for p in p_a]
            lp_b = [-math.log10(max(p, 1e-300)) if math.isfinite(p) else math.nan for p in p_b]
            def mad(xs, ys):
                vals = [abs(x - y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
                return (sum(vals) / len(vals), max(vals)) if vals else (math.nan, math.nan)
            mab, xab = mad(bet_a, bet_b)
            mas, xas = mad(se_a, se_b)
            mapv, xapv = mad(p_a, p_b)
            out.write(f"{a}\t{b}\t{len(common)}\t{corr(bet_a, bet_b):.12g}\t{corr(se_a, se_b):.12g}\t{corr(lp_a, lp_b):.12g}\t{mab:.12g}\t{xab:.12g}\t{mas:.12g}\t{xas:.12g}\t{mapv:.12g}\t{xapv:.12g}\n")

print(f"Wrote comparison: {out_path}")
PY

echo "[biobank-stress] Summary: ${SUMMARY}"
cat "${SUMMARY}"
echo "[biobank-stress] Comparison: ${COMPARE}"
cat "${COMPARE}"

if [[ "${KEEP_DATA}" != "1" ]]; then
  echo "[biobank-stress] Removing generated data directory: ${DATA_DIR}"
  rm -rf "${DATA_DIR}"
fi
