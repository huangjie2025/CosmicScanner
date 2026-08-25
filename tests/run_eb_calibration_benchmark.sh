#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
DATA_DIR="${ROOT}/tests/data/eb_calibration_benchmark"
CAL_LABEL="${CAL_LABEL:-}"
if [[ -n "${CAL_LABEL}" ]]; then
  RESULT_DIR="${ROOT}/tests/results/eb_calibration_benchmark_${CAL_LABEL}"
else
  RESULT_DIR="${ROOT}/tests/results/eb_calibration_benchmark"
fi
PREFIX="${DATA_DIR}/ebcal"
REPS="${REPS:-4}"
N="${N:-96}"
M="${M:-320}"
EFFECT="${EFFECT:-0.80}"
NOISE="${NOISE:-0.70}"
POP_SHIFT="${POP_SHIFT:-0.08}"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}"

python3 - "${PREFIX}" "${REPS}" "${N}" "${M}" "${EFFECT}" "${NOISE}" "${POP_SHIFT}" <<'PY'
import math, random, sys
from pathlib import Path

prefix = Path(sys.argv[1])
reps = int(sys.argv[2])
n = int(sys.argv[3])
m = int(sys.argv[4])
effect_power = float(sys.argv[5])
noise = float(sys.argv[6])
pop_shift = float(sys.argv[7])
rng = random.Random(20260626)
causal = min(m - 1, max(0, int(m * 0.285)))
samples = [(f'F{i+1:04d}', f'I{i+1:04d}', 0 if i < n//2 else 1) for i in range(n)]
genos, variants = [], []
for j in range(m):
    chrom = str((j % 4) + 1)
    base = 0.08 + 0.42 * ((j * 43) % 101) / 100.0
    shift = pop_shift * math.sin(j * 0.13)
    row = []
    for _, _, pop in samples:
        p = min(0.95, max(0.02, base + (shift if pop else -shift)))
        row.append((1 if rng.random() < p else 0) + (1 if rng.random() < p else 0))
    genos.append(row)
    variants.append((chrom, f'rs{j+1:06d}', 1000 + 10*j, 'A', 'C'))

prefix.parent.mkdir(parents=True, exist_ok=True)
with prefix.with_suffix('.fam').open('w') as f:
    for fid, iid, _ in samples:
        f.write(f'{fid}\t{iid}\t0\t0\t1\t-9\n')
with prefix.with_suffix('.bim').open('w') as f:
    for chrom, snp, bp, a1, a2 in variants:
        f.write(f'{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n')

def code(g):
    return 0b00 if g == 2 else (0b10 if g == 1 else 0b11)
bps = (n + 3) // 4
bed = bytearray(b'\x6c\x1b\x01')
for row in genos:
    block = bytearray(bps)
    for i, g in enumerate(row):
        block[i // 4] |= code(g) << ((i % 4) * 2)
    bed.extend(block)
prefix.with_suffix('.bed').write_bytes(bed)

for rep in range(1, reps + 1):
    local = random.Random(20260626 + rep)
    for scenario, effect in [('null', 0.0), ('power', effect_power)]:
        with prefix.with_suffix(f'.{scenario}.{rep}.pheno').open('w') as f:
            for i, (fid, iid, pop) in enumerate(samples):
                g = genos[causal][i] - 1.0
                y = 0.18 * pop + effect * g + local.gauss(0, noise)
                f.write(f'{fid}\t{iid}\t{y:.8f}\n')
with prefix.with_suffix('.truth.txt').open('w') as f:
    f.write(f'causal_snp\trs{causal+1:06d}\n')
PY

for rep in $(seq 1 "${REPS}"); do
  for scenario in null power; do
    "${BIN}" --bfile "${PREFIX}" --pheno "${PREFIX}.${scenario}.${rep}.pheno" --mpheno 1 --gwas --model lmm --no-loco --maf 0.0 --out "${RESULT_DIR}/${scenario}_${rep}.standard" >/dev/null
    "${BIN}" --bfile "${PREFIX}" --pheno "${PREFIX}.${scenario}.${rep}.pheno" --mpheno 1 --gwas --model lmm-eb --adjust eb-mbonf --no-loco --maf 0.0 --out "${RESULT_DIR}/${scenario}_${rep}.eb" >/dev/null
  done
done

python3 - "${RESULT_DIR}" "${REPS}" "${M}" <<'PY'
import csv, math, sys
from pathlib import Path

result_dir = Path(sys.argv[1])
reps = int(sys.argv[2])
m = int(sys.argv[3])
causal = f'rs{min(m - 1, max(0, int(m * 0.285))) + 1:06d}'
rows = []
for scenario in ['null', 'power']:
    for rep in range(1, reps + 1):
        std = list(csv.DictReader(open(result_dir/f'{scenario}_{rep}.standard.assoc.txt'), delimiter='\t'))
        eb = list(csv.DictReader(open(result_dir/f'{scenario}_{rep}.eb.assoc.txt'), delimiter='\t'))
        vals = {}
        for line in open(result_dir/f'{scenario}_{rep}.eb.eb.summary.txt'):
            k, _, v = line.rstrip('\n').partition('\t')
            vals[k] = v
        std_p = [float(r['P']) for r in std]
        eb_p = [float(r['P']) for r in eb]
        eb_thr = float(vals['eb_mbonf_0.05'])
        bonf_thr = float(vals['bonferroni_0.05'])
        std_top = sorted(std, key=lambda r: float(r['P']))
        eb_top = sorted(eb, key=lambda r: float(r['P']))
        causal_rank = next((i + 1 for i, r in enumerate(eb_top) if r['SNP'] == causal), '')
        rows.append({
            'scenario': scenario,
            'rep': rep,
            'markers': len(eb),
            'standard_bonf_hits': sum(p < bonf_thr for p in std_p),
            'eb_bonf_hits': sum(p < bonf_thr for p in eb_p),
            'eb_mbonf_hits': sum(p < eb_thr for p in eb_p),
            'effective_tests': float(vals['effective_tests']),
            'effective_tests_per_marker': float(vals['effective_tests_per_marker']),
            'standard_top1': std_top[0]['SNP'],
            'eb_top1': eb_top[0]['SNP'],
            'causal_rank': causal_rank,
        })

with (result_dir/'eb_calibration.summary.tsv').open('w', newline='') as f:
    fields = ['scenario','rep','markers','standard_bonf_hits','eb_bonf_hits','eb_mbonf_hits','effective_tests','effective_tests_per_marker','standard_top1','eb_top1','causal_rank']
    w = csv.DictWriter(f, fields, delimiter='\t')
    w.writeheader(); w.writerows(rows)

agg = []
for scenario in ['null', 'power']:
    subset = [r for r in rows if r['scenario'] == scenario]
    agg.append({
        'scenario': scenario,
        'reps': len(subset),
        'standard_bonf_reps_with_hits': sum(r['standard_bonf_hits'] > 0 for r in subset),
        'eb_bonf_reps_with_hits': sum(r['eb_bonf_hits'] > 0 for r in subset),
        'eb_mbonf_reps_with_hits': sum(r['eb_mbonf_hits'] > 0 for r in subset),
        'mean_effective_tests': sum(r['effective_tests'] for r in subset) / len(subset),
        'power_causal_top10': sum((r['causal_rank'] != '' and int(r['causal_rank']) <= 10) for r in subset),
    })
with (result_dir/'eb_calibration.aggregate.tsv').open('w', newline='') as f:
    fields = ['scenario','reps','standard_bonf_reps_with_hits','eb_bonf_reps_with_hits','eb_mbonf_reps_with_hits','mean_effective_tests','power_causal_top10']
    w = csv.DictWriter(f, fields, delimiter='\t')
    w.writeheader(); w.writerows(agg)
with (result_dir/'eb_calibration.summary.md').open('w') as f:
    f.write('# EB Calibration Benchmark\n\n')
    f.write('| Scenario | Reps | Standard Bonf hit reps | EB Bonf hit reps | EB-mBonf hit reps | Mean effective tests | Causal top10 reps |\n')
    f.write('|---|---:|---:|---:|---:|---:|---:|\n')
    for r in agg:
        f.write(f"| {r['scenario']} | {r['reps']} | {r['standard_bonf_reps_with_hits']} | {r['eb_bonf_reps_with_hits']} | {r['eb_mbonf_reps_with_hits']} | {r['mean_effective_tests']:.4g} | {r['power_causal_top10']} |\n")

null_eb_hits = next(r for r in agg if r['scenario'] == 'null')['eb_mbonf_reps_with_hits']
power_top10 = next(r for r in agg if r['scenario'] == 'power')['power_causal_top10']
if null_eb_hits > reps:
    raise SystemExit('invalid null hit count')
if power_top10 < max(1, reps // 2):
    raise SystemExit(f'weak power benchmark: causal top10 reps={power_top10}/{reps}')
print(f'EB calibration benchmark passed: {result_dir / "eb_calibration.aggregate.tsv"}')
PY
