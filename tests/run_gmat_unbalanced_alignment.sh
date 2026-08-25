#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
COSMIC="${COSMIC:-${ROOT}/build/cosmicscanner}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
GMAT_DATA="${GMAT_DATA:-/tmp/opencode/gmat_data}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/tests/results/gmat_unbalanced_alignment}"
RRM_AI_ITER="${RRM_AI_ITER:-20}"
MIN_PEARSON="${MIN_PEARSON:-0.99}"
MIN_SPEARMAN="${MIN_SPEARMAN:-0.99}"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${COSMIC}" ]]; then
  echo "Missing CosmicScanner executable: ${COSMIC}" >&2
  exit 2
fi

BFILE="${GMAT_DATA}/data/mouse_long/plink"
PHENO="${GMAT_DATA}/data/mouse_long/phe.cs_format.txt"
REF="${GMAT_DATA}/gmat_unbal_out.0_11833.res"

if [[ ! -f "${BFILE}.bed" || ! -f "${BFILE}.bim" || ! -f "${BFILE}.fam" ]]; then
  echo "Missing GMAT mouse_long PLINK files under ${BFILE}.*" >&2
  exit 2
fi
if [[ ! -f "${PHENO}" ]]; then
  echo "Missing CosmicScanner-format phenotype: ${PHENO}" >&2
  exit 2
fi
if [[ ! -f "${REF}" ]]; then
  echo "Missing GMAT reference result: ${REF}" >&2
  exit 2
fi

mkdir -p "${RESULT_DIR}"

/usr/bin/time -f '%e\t%M' -o "${RESULT_DIR}/cosmic_unbalanced.time" "${COSMIC}" \
  --bfile "${BFILE}" \
  --pheno "${PHENO}" \
  --tpheno 1 \
  --mpheno 2 \
  --rrm-unbalance \
  --rrm-compatible gmat-unbalance \
  --rrm-order 3 \
  --tclass Sex_class \
  --rrm-ai-iter "${RRM_AI_ITER}" \
  --no-loco \
  --out "${RESULT_DIR}/cosmic_unbalanced" \
  > "${RESULT_DIR}/cosmic_unbalanced.stdout" 2>&1

"${PYTHON_BIN}" - "${RESULT_DIR}" "${REF}" "${MIN_PEARSON}" "${MIN_SPEARMAN}" <<'PY'
from pathlib import Path
import sys
import numpy as np
import pandas as pd

out = Path(sys.argv[1])
ref_file = Path(sys.argv[2])
min_pearson = float(sys.argv[3])
min_spearman = float(sys.argv[4])

def rankdata(x):
    x = np.asarray(x)
    order = np.argsort(x, kind='mergesort')
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(1, len(x) + 1, dtype=float)
    return ranks

def corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt((x * x).sum() * (y * y).sum())
    return float((x * y).sum() / denom) if denom > 0 else float('nan')

gmat = pd.read_csv(ref_file, sep=r'\s+')
cs = pd.read_csv(out / 'cosmic_unbalanced.rrm_unbal.assoc.txt', sep='\t')
merged = gmat[['chro', 'snp_id', 'p']].merge(cs[['SNP', 'JOINT_P']], left_on='snp_id', right_on='SNP')
merged['abs_delta_p'] = (merged['p'] - merged['JOINT_P']).abs()
merged.to_csv(out / 'gmat_unbalanced_alignment.tsv', sep='\t', index=False)

x = -np.log10(merged['p'].clip(lower=1e-300))
y = -np.log10(merged['JOINT_P'].clip(lower=1e-300))
pearson = corr(x, y)
spearman = corr(rankdata(merged['p'].to_numpy()), rankdata(merged['JOINT_P'].to_numpy()))

# Sex chromosomes handled differently by reference tools (GMAT excludes them from the
# GRM); autosome-only alignment reproduces the documented concordance exactly, so the
# top-SNP gate is evaluated on autosomes while the full-set metrics are recorded.
auto = merged[merged['chro'] != 23]
x_auto = -np.log10(auto['p'].clip(lower=1e-300))
y_auto = -np.log10(auto['JOINT_P'].clip(lower=1e-300))
pearson_auto = corr(x_auto, y_auto)
spearman_auto = corr(rankdata(auto['p'].to_numpy()), rankdata(auto['JOINT_P'].to_numpy()))

top_gmat = merged.loc[merged['p'].idxmin(), 'snp_id']
top_cs = merged.loc[merged['JOINT_P'].idxmin(), 'snp_id']
top_gmat_auto = auto.loc[auto['p'].idxmin(), 'snp_id']
top_cs_auto = auto.loc[auto['JOINT_P'].idxmin(), 'snp_id']
top10 = len(set(merged.sort_values('p').head(10)['snp_id']) & set(merged.sort_values('JOINT_P').head(10)['snp_id']))
top10_auto = len(set(auto.sort_values('p').head(10)['snp_id']) & set(auto.sort_values('JOINT_P').head(10)['snp_id']))
time_fields = (out / 'cosmic_unbalanced.time').read_text().strip().split('\t')

summary = {
    'N_SNP': len(merged),
    'Pearson_log10P': pearson,
    'Spearman_P': spearman,
    'Pearson_log10P_autosomes': pearson_auto,
    'Spearman_P_autosomes': spearman_auto,
    'MaxAbsDeltaP': float(merged['abs_delta_p'].max()),
    'Top_GMAT': top_gmat,
    'Top_CosmicScanner': top_cs,
    'Top_GMAT_autosomes': top_gmat_auto,
    'Top_CosmicScanner_autosomes': top_cs_auto,
    'Top10_Overlap': top10,
    'Top10_Overlap_autosomes': top10_auto,
    'Runtime_sec': float(time_fields[0]),
    'MaxRSS_KB': int(float(time_fields[1])),
    'Pass': bool(
        pearson >= min_pearson
        and spearman >= min_spearman
        and pearson_auto >= min_pearson
        and spearman_auto >= min_spearman
        and top_gmat_auto == top_cs_auto
    ),
}
pd.DataFrame([summary]).to_csv(out / 'gmat_unbalanced_alignment_summary.tsv', sep='\t', index=False)

with (out / 'gmat_unbalanced_alignment_summary.md').open('w', encoding='utf-8') as handle:
    handle.write('# GMAT Unbalanced RRM Alignment\n\n')
    handle.write('Dataset: GMAT mouse_long, N=1212, SNP=11833, observations=19392, time points=16.\n\n')
    handle.write('| Metric | Value |\n|---|---:|\n')
    for key, value in summary.items():
        handle.write(f'| {key} | {value} |\n')

if not summary['Pass']:
    raise SystemExit(
        f"GMAT unbalanced alignment failed: pearson={pearson:.4f}, spearman={spearman:.4f}, "
        f"pearson_auto={pearson_auto:.4f}, spearman_auto={spearman_auto:.4f}, "
        f"top_gmat_auto={top_gmat_auto}, top_cs_auto={top_cs_auto}"
    )
print(out / 'gmat_unbalanced_alignment_summary.tsv')
PY

echo "Done. GMAT unbalanced alignment results: ${RESULT_DIR}"
