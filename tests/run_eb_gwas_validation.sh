#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/eb_gwas_validation"
RESULT_DIR="${ROOT}/tests/results/eb_gwas_validation"
PREFIX="${DATA_DIR}/ebval"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}"

python3 - "${PREFIX}" <<'PY'
import math, random, sys
from pathlib import Path

prefix = Path(sys.argv[1])
rng = random.Random(20260624)
n, m = 84, 260
causal = 57
samples = [(f"F{i+1:04d}", f"I{i+1:04d}", 0 if i < n//2 else 1) for i in range(n)]
genos, variants = [], []
for j in range(m):
    chrom = str((j % 4) + 1)
    base = 0.08 + 0.42 * ((j * 37) % 101) / 100.0
    shift = 0.06 * math.sin(j * 0.19)
    row = []
    for _, _, pop in samples:
        p = min(0.95, max(0.02, base + (shift if pop else -shift)))
        row.append((1 if rng.random() < p else 0) + (1 if rng.random() < p else 0))
    genos.append(row)
    variants.append((chrom, f"rs{j+1:06d}", 1000 + 10 * j, "A", "C"))

prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix('.fam').open('w') as f:
    for fid, iid, _ in samples:
        f.write(f"{fid}\t{iid}\t0\t0\t1\t-9\n")
with prefix.with_suffix('.bim').open('w') as f:
    for v in variants:
        f.write('\t'.join(map(str, v)) + '\n')

def code(g):
    return 0b00 if g == 2 else (0b10 if g == 1 else 0b11)

bps = (n + 3) // 4
bed = bytearray(b"\x6c\x1b\x01")
for row in genos:
    block = bytearray(bps)
    for i, g in enumerate(row):
        block[i // 4] |= code(g) << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix('.bed').write_bytes(bed)

u = [rng.gauss(0, 0.20) for _ in range(n)]
with prefix.with_suffix('.single_null.pheno').open('w') as fnull, prefix.with_suffix('.single_power.pheno').open('w') as fpow:
    for i, (fid, iid, pop) in enumerate(samples):
        g = genos[causal][i] - 1.0
        y_null = 0.15 * pop + rng.gauss(0, 0.80)
        y_power = 0.15 * pop + 0.95 * g + rng.gauss(0, 0.60)
        fnull.write(f"{fid}\t{iid}\t{y_null:.8f}\n")
        fpow.write(f"{fid}\t{iid}\t{y_power:.8f}\n")

with prefix.with_suffix('.rep_null.pheno').open('w') as fnull, prefix.with_suffix('.rep_power.pheno').open('w') as fpow:
    for i, (fid, iid, pop) in enumerate(samples):
        g = genos[causal][i] - 1.0
        for r in range(3):
            y_null = 0.12 * pop + u[i] + 0.03 * r + rng.gauss(0, 0.45)
            y_power = 0.12 * pop + u[i] + 0.03 * r + 0.80 * g + rng.gauss(0, 0.35)
            fnull.write(f"{fid}\t{iid}\t{y_null:.8f}\n")
            fpow.write(f"{fid}\t{iid}\t{y_power:.8f}\n")

with prefix.with_suffix('.mv_null.pheno').open('w') as fnull, prefix.with_suffix('.mv_power.pheno').open('w') as fpow:
    for i, (fid, iid, pop) in enumerate(samples):
        g = genos[causal][i] - 1.0
        e1 = rng.gauss(0, 0.65); e2 = 0.35 * e1 + rng.gauss(0, 0.55)
        fnull.write(f"{fid}\t{iid}\t{0.10*pop + e1:.8f}\t{0.05*pop + e2:.8f}\n")
        e1 = rng.gauss(0, 0.55); e2 = 0.35 * e1 + rng.gauss(0, 0.50)
        fpow.write(f"{fid}\t{iid}\t{0.10*pop + 0.80*g + e1:.8f}\t{0.05*pop - 0.60*g + e2:.8f}\n")

with prefix.with_suffix('.truth.txt').open('w') as f:
    f.write(f"causal_snp\trs{causal+1:06d}\n")
PY

run_single() {
  local pheno="$1" out="$2"
  "${BIN}" --bfile "${PREFIX}" --pheno "${pheno}" --mpheno 1 --gwas --model lmm-eb --adjust eb-mbonf --no-loco --maf 0.0 --out "${out}"
}

run_repeatability() {
  local pheno="$1" out="$2"
  "${BIN}" --bfile "${PREFIX}" --pheno "${pheno}" --mpheno 1 --repeatability-model --gwas --model lmm-eb --adjust eb-mbonf --no-loco --maf 0.0 --out "${out}"
}

run_mv() {
  local pheno="$1" out="$2"
  "${BIN}" --bfile "${PREFIX}" --pheno "${pheno}" --mpheno-list 1,2 --multivariate --mv-method eb --gwas --no-loco --maf 0.0 --out "${out}"
}

run_single "${PREFIX}.single_null.pheno" "${RESULT_DIR}/single_null"
run_single "${PREFIX}.single_power.pheno" "${RESULT_DIR}/single_power"
run_repeatability "${PREFIX}.rep_null.pheno" "${RESULT_DIR}/rep_null"
run_repeatability "${PREFIX}.rep_power.pheno" "${RESULT_DIR}/rep_power"
run_mv "${PREFIX}.mv_null.pheno" "${RESULT_DIR}/mv_null"
run_mv "${PREFIX}.mv_power.pheno" "${RESULT_DIR}/mv_power"

python3 - "${RESULT_DIR}" <<'PY'
import csv, math, sys
from pathlib import Path

result_dir = Path(sys.argv[1])
causal = 'rs000058'

cases = [
    ('single_null', 'single', False, result_dir/'single_null.assoc.txt', result_dir/'single_null.eb.summary.txt'),
    ('single_power', 'single', True, result_dir/'single_power.assoc.txt', result_dir/'single_power.eb.summary.txt'),
    ('rep_null', 'repeatability', False, result_dir/'rep_null.assoc.txt', result_dir/'rep_null.repeatability.eb.summary.txt'),
    ('rep_power', 'repeatability', True, result_dir/'rep_power.assoc.txt', result_dir/'rep_power.repeatability.eb.summary.txt'),
    ('mv_null', 'multivariate', False, result_dir/'mv_null.assoc.mv.txt', result_dir/'mv_null.eb.mv.summary.txt'),
    ('mv_power', 'multivariate', True, result_dir/'mv_power.assoc.mv.txt', result_dir/'mv_power.eb.mv.summary.txt'),
]

def read_summary(path):
    vals = {}
    for line in open(path):
        k, _, v = line.rstrip('\n').partition('\t')
        vals[k] = v
    return vals

def chisq_from_p1(p):
    # monotonic proxy sufficient for lambda median: use -log10(P) scale when scipy is unavailable
    return -math.log10(max(p, 1e-300))

summary_rows = []
for name, model, is_power, assoc_path, summary_path in cases:
    rows = list(csv.DictReader(open(assoc_path), delimiter='\t'))
    if len(rows) != 260:
        raise SystemExit(f"{name}: expected 260 rows, got {len(rows)}")
    pvals = []
    bad = 0
    conf_sum = 0.0
    conf_upper = 2.0 if model == 'multivariate' else 1.0
    for r in rows:
        for col in ('P', 'LOG10P') if model != 'multivariate' else ('P',):
            if col in r and not math.isfinite(float(r[col])): bad += 1
        for col in ('MARKER_LAMBDA', 'EB_CONFIDENCE'):
            if col not in r:
                raise SystemExit(f"{name}: missing EB column {col}")
            val = float(r[col])
            if not math.isfinite(val):
                raise SystemExit(f"{name}: non-finite {col}")
            if col == 'MARKER_LAMBDA' and val < -1e-12:
                raise SystemExit(f"{name}: negative MARKER_LAMBDA={val}")
            if col == 'EB_CONFIDENCE':
                if val < -1e-12 or val > conf_upper + 1e-8:
                    raise SystemExit(f"{name}: EB_CONFIDENCE out of range ({val}, upper={conf_upper})")
                conf_sum += max(0.0, min(conf_upper, val))
        pvals.append(float(r['P']))
    if bad:
        raise SystemExit(f"{name}: non-finite association values count={bad}")
    ordered = sorted(rows, key=lambda r: float(r['P']))
    top = ordered[0]['SNP']
    rank = next((i+1 for i, r in enumerate(ordered) if r['SNP'] == causal), None)
    vals = read_summary(summary_path)
    meff = float(vals['effective_tests'])
    tested = int(float(vals['tested_markers']))
    eb_thr = float(vals['eb_mbonf_0.05'])
    bonf_thr = float(vals['bonferroni_0.05'])
    if tested != 260:
        raise SystemExit(f"{name}: expected tested_markers=260, got {tested}")
    if abs(conf_sum - meff) > max(1e-6, 1e-6 * max(1.0, meff)):
        raise SystemExit(f"{name}: effective_tests mismatch summary={meff} assoc_sum={conf_sum}")
    upper = 520 if model == 'multivariate' else 260
    if not (0.0 <= meff <= upper):
        raise SystemExit(f"{name}: invalid effective_tests={meff}")
    n_bonf = sum(p < bonf_thr for p in pvals)
    n_eb = sum(p < eb_thr for p in pvals)
    if not is_power and n_eb > 5:
        raise SystemExit(f"{name}: null EB-mBonf hits too high ({n_eb})")
    if is_power and (rank is None or rank > 10):
        raise SystemExit(f"{name}: causal {causal} rank too low ({rank}) top={top}")
    median_logp = sorted(chisq_from_p1(p) for p in pvals)[len(pvals)//2]
    summary_rows.append({
        'case': name, 'model': model, 'scenario': 'power' if is_power else 'null',
        'rows': len(rows), 'effective_tests': meff, 'bonf_hits': n_bonf,
        'eb_mbonf_hits': n_eb, 'causal_rank': rank or '', 'top_snp': top,
        'median_log10p': median_logp,
    })

out_tsv = result_dir/'eb_gwas_validation.summary.tsv'
with out_tsv.open('w', newline='') as f:
    fieldnames = ['case','model','scenario','rows','effective_tests','bonf_hits','eb_mbonf_hits','causal_rank','top_snp','median_log10p']
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    w.writeheader(); w.writerows(summary_rows)

out_md = result_dir/'eb_gwas_validation.summary.md'
with out_md.open('w') as f:
    f.write('# EB GWAS Validation Summary\n\n')
    f.write('| Case | Model | Scenario | Rows | Effective tests | Bonf hits | EB-mBonf hits | Causal rank | Top SNP |\n')
    f.write('|---|---|---:|---:|---:|---:|---:|---:|---|\n')
    for r in summary_rows:
        f.write(f"| {r['case']} | {r['model']} | {r['scenario']} | {r['rows']} | {float(r['effective_tests']):.4g} | {r['bonf_hits']} | {r['eb_mbonf_hits']} | {r['causal_rank']} | {r['top_snp']} |\n")
print(f"EB GWAS validation passed: {out_tsv}")
PY
