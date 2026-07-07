#!/usr/bin/env bash
set -euo pipefail

COSMIC="/disk228/huangjie/PanOmics/CosmicScanner/build/cosmicscanner"
BASE="/disk228/huangjie/PanOmics/CosmicScanner/results/human_ukb_benchmark"
PFILE="${BASE}/prepared/human_chr1_22_merged"
PHENO="${BASE}/prepared/human_ukb_10traits.standardized.pheno"
OUTDIR="${BASE}/cosmicscanner_full"
THREADS="${THREADS:-10}"
EXPECTED_LINES=2026363

mkdir -p "${OUTDIR}"

for traits in 2 4 6 8 10; do
  prefix="${OUTDIR}/human_t${traits}_full"
  assoc="${prefix}.assoc.mv.txt"
  log="${prefix}.run.log"
  cols="$(seq -s, 1 "${traits}")"

  if [[ -s "${assoc}" ]] && [[ "$(wc -l < "${assoc}")" -eq "${EXPECTED_LINES}" ]]; then
    echo "SKIP t=${traits}: existing complete output ${assoc}"
    continue
  fi

  echo "RUN t=${traits}: $(date)"
  /usr/bin/time -v "${COSMIC}" \
    --pfile "${PFILE}" \
    --pheno "${PHENO}" \
    --mpheno-list "${cols}" \
    --multivariate --mv-method emmv --mv-compatible emmvgwas \
    --gwas --no-loco --threads "${THREADS}" \
    --out "${prefix}" 2>&1 | tee "${log}"
  echo "DONE t=${traits}: $(date)"

  if [[ ! -s "${assoc}" ]]; then
    echo "ERROR t=${traits}: missing output ${assoc}" >&2
    exit 1
  fi
  lines="$(wc -l < "${assoc}")"
  if [[ "${lines}" -ne "${EXPECTED_LINES}" ]]; then
    echo "ERROR t=${traits}: expected ${EXPECTED_LINES} lines, observed ${lines}" >&2
    exit 1
  fi
done
