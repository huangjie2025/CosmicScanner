#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GMA_SRC="${GMA_SRC:-/tmp/opencode/gma_py3/tests}"
GMA_PYTHONPATH="${GMA_PYTHONPATH:-/tmp/opencode/gma_py3_deps:/tmp/opencode/gma_py3_install}"
PLINK2="${PLINK2:-${PLINK2_BIN:-plink2}}"
PLINK2="$(command -v "${PLINK2}" || printf '%s' "${PLINK2}")"
COSMIC="${COSMIC:-${ROOT}/build/cosmicscanner}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
N="${N:-120}"
M="${M:-20}"
RESULT_DIR="${RESULT_DIR:-${ROOT}/tests/results/gma_fixreg_alignment_n${N}_m${M}}"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${PLINK2}" ]]; then
  echo "Missing plink2 executable: ${PLINK2}" >&2
  exit 2
fi
if [[ ! -d "${GMA_SRC}" ]]; then
  echo "Missing converted GMA test directory: ${GMA_SRC}" >&2
  exit 2
fi
if [[ ! -x "${COSMIC}" ]]; then
  echo "Missing CosmicScanner executable: ${COSMIC}" >&2
  exit 2
fi

mkdir -p "${RESULT_DIR}"

PYTHONPATH="${GMA_PYTHONPATH}" "${PYTHON_BIN}" - "${GMA_SRC}" "${RESULT_DIR}" "${N}" "${M}" <<'PY'
from pathlib import Path
import sys
import pandas as pd

gma_src = Path(sys.argv[1])
out = Path(sys.argv[2])
n = int(sys.argv[3])
m = int(sys.argv[4])
fam = pd.read_csv(gma_src / 'plink.fam', sep=r'\s+', header=None,
                  names=['FID', 'IID', 'PID', 'MID', 'SEX0', 'PHENO'], dtype={'FID': str, 'IID': str})
fam[['FID', 'IID']].iloc[:n].to_csv(out / 'keep.txt', sep='\t', index=False, header=False)
bim = pd.read_csv(gma_src / 'plink.bim', sep=r'\s+', header=None)
bim.iloc[:m, 1].to_csv(out / 'extract.txt', index=False, header=False)
PY

"${PLINK2}" \
  --bfile "${GMA_SRC}/plink" \
  --keep "${RESULT_DIR}/keep.txt" \
  --extract "${RESULT_DIR}/extract.txt" \
  --make-bed \
  --out "${RESULT_DIR}/plink_subset" \
  > "${RESULT_DIR}/plink_subset.stdout" 2>&1

PYTHONPATH="${GMA_PYTHONPATH}" "${PYTHON_BIN}" - "${GMA_SRC}" "${RESULT_DIR}" <<'PY'
from pathlib import Path
import sys
import pandas as pd

gma_src = Path(sys.argv[1])
out = Path(sys.argv[2])
fam = pd.read_csv(out / 'plink_subset.fam', sep=r'\s+', header=None,
                  names=['FID', 'IID', 'PID', 'MID', 'SEX0', 'PHENO'], dtype={'FID': str, 'IID': str})
ph = pd.read_csv(gma_src / 'phe.balance.txt', sep=r'\s+', dtype={'ID': str})
merged = fam[['FID', 'IID']].merge(ph, left_on='IID', right_on='ID', how='inner')
balance = merged[['ID', 'Sex'] + [f'trait{i}' for i in range(1, 17)]]
balance.to_csv(out / 'phe.balance.subset.txt', sep='\t', index=False)
cs_pheno = merged[['FID', 'IID'] + [f'trait{i}' for i in range(1, 17)]]
cs_pheno.to_csv(out / 'phe.balance.cosmic.txt', sep='\t', index=False, header=False)
PY

PYTHONPATH="${GMA_PYTHONPATH}" /usr/bin/time -f '%e\t%M' -o "${RESULT_DIR}/gma_reference.time" "${PYTHON_BIN}" - "${RESULT_DIR}" <<'PY'
from pathlib import Path
import sys
import numpy as np
from gma import cal_kin, balance_varcom, balance_longwas_fixreg

out = Path(sys.argv[1])
bed = str(out / 'plink_subset')
cal_kin(bed)
tpoint = np.array(list(range(16))) + 1.0
trait = list(range(2, 18))
res_var = balance_varcom(str(out / 'phe.balance.subset.txt'), 'ID', tpoint, trait, bed + '.grm',
                         prefix_outfile=str(out / 'gma_reference_varcom'))
balance_longwas_fixreg(str(out / 'phe.balance.subset.txt'), 'ID', tpoint, trait, bed + '.grm', bed,
                       res_var['variances'], snp_lst=None, prefix_outfile=str(out / 'gma_reference_fixreg'))
PY

/usr/bin/time -f '%e\t%M' -o "${RESULT_DIR}/cs_compatible.time" "${COSMIC}" \
  --bfile "${RESULT_DIR}/plink_subset" \
  --gwas \
  --rrm-compatible gma-fixreg \
  --rrm-compatible-varcom "${RESULT_DIR}/gma_reference_varcom.var" \
  --pheno "${RESULT_DIR}/phe.balance.cosmic.txt" \
  --mpheno-list 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 \
  --rrm-order 3 \
  --maf 0 \
  --no-loco \
  --out "${RESULT_DIR}/cs_compatible" \
  > "${RESULT_DIR}/cs_compatible.stdout" 2>&1

PYTHONPATH="${GMA_PYTHONPATH}" "${PYTHON_BIN}" - "${RESULT_DIR}" <<'PY'
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

out = Path(sys.argv[1])
gma = pd.read_csv(out / 'gma_reference_fixreg.res', sep=r'\s+')
cs = pd.read_csv(out / 'cs_compatible.rrm.gma_fixreg.assoc.txt', sep=r'\s+')
cs = cs.rename(columns={'SNP': 'snp_ID'})
merged = gma[['snp_ID', 'chi_val', 'p_val']].merge(cs[['snp_ID', 'chi_val', 'p_val']], on='snp_ID', suffixes=('_gma', '_cs'))
merged['abs_delta_p'] = (merged['p_val_gma'] - merged['p_val_cs']).abs()
merged['abs_delta_chi'] = (merged['chi_val_gma'] - merged['chi_val_cs']).abs()
merged.to_csv(out / 'gma_fixreg_alignment.tsv', sep='\t', index=False)
x = -np.log10(merged['p_val_gma'].astype(float).clip(lower=1e-300))
y = -np.log10(merged['p_val_cs'].astype(float).clip(lower=1e-300))
summary = {
    'N_SNP': len(merged),
    'Pearson_log10P': pearsonr(x, y).statistic,
    'Spearman_P': spearmanr(merged['p_val_gma'], merged['p_val_cs']).statistic,
    'MaxAbsDeltaP': merged['abs_delta_p'].max(),
    'MaxAbsDeltaChi': merged['abs_delta_chi'].max(),
    'Top_GMA': merged.loc[merged['p_val_gma'].idxmin(), 'snp_ID'],
    'Top_CS_Compatible': merged.loc[merged['p_val_cs'].idxmin(), 'snp_ID'],
}
pd.DataFrame([summary]).to_csv(out / 'gma_fixreg_alignment_summary.tsv', sep='\t', index=False)
gma_time = (out / 'gma_reference.time').read_text().strip().split('\t')
cs_time = (out / 'cs_compatible.time').read_text().strip().split('\t')
sample_n = pd.read_csv(out / 'phe.balance.subset.txt', sep=r'\s+').shape[0]
with (out / 'gma_fixreg_alignment_summary.md').open('w', encoding='utf-8') as handle:
    handle.write('# GMA Fixreg Alignment\n\n')
    handle.write(f"Dataset: official GMA subset, N={sample_n}, SNP={len(merged)}, TIME=16.\n\n")
    handle.write('| Metric | Value |\n|---|---:|\n')
    for key, value in summary.items():
        handle.write(f'| {key} | {value} |\n')
    handle.write('\n| Method | Runtime(s) | MaxRSS(KB) |\n|---|---:|---:|\n')
    handle.write(f'| GMA reference wrapper | {gma_time[0]} | {gma_time[1]} |\n')
    handle.write(f'| CosmicScanner gma-fixreg compatible | {cs_time[0]} | {cs_time[1]} |\n')
print(out / 'gma_fixreg_alignment_summary.tsv')
PY

echo "Done. Alignment results: ${RESULT_DIR}"
