#!/usr/bin/env bash
set -euo pipefail

ROOT="${COSMICSCANNER_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
BFILE="${ROOT}/tests/data/gemma_accuracy/gemma_accuracy"
WORK="${INPUT_ADAPT_WORK:-/tmp/opencode/cosmicscanner_input_adaptation}"

mkdir -p "${WORK}"
rm -rf "${WORK}/long" "${WORK}/wide" "${WORK}/preview" "${WORK}/strict" "${WORK}/config" "${WORK}/infer_only" "${WORK}/rrm" "${WORK}/rrm_smoke"* "${WORK}/lmm" "${WORK}/lmm_smoke"* "${WORK}/mv" "${WORK}/mv_smoke"* "${WORK}/realdata_cpp" "${WORK}/realdata_cpp_smoke" "${WORK}/realdata_vol_smoke"* "${WORK}/extract_100_snps.txt" "${WORK}/extract_50_snps.txt"

if [[ ! -x "${BIN}" ]]; then
  echo "ERROR: cosmicscanner binary not found: ${BIN}" >&2
  exit 1
fi

cat > "${WORK}/raw_long.csv" <<'CSV'
pig_id,collect_date,volume,concentration,age,temp,farm,line
I0001,2024-01-01,10.5,210,300,37.1,A,L1
I0001,2024-01-01,11.0,215,300,37.2,A,L1
I0002,2024-01-03,8.5,180,290,36.9,B,L2
I0003,2024-01-04,9.1,190,310,37.0,A,L1
NO_MATCH,2024-01-04,9.9,199,310,37.0,A,L1
CSV

"${BIN}" \
  --bfile "${BFILE}" \
  --input-file "${WORK}/raw_long.csv" \
  --input-model repeatability \
  --id-col pig_id \
  --time-col collect_date \
  --trait-cols volume,concentration \
  --record-qcovars age,temp \
  --subject-ccovars farm,line \
  --time-parser date \
  --prepare-input-only \
  --prepared-input-dir "${WORK}/long"

python3 - <<PY
from pathlib import Path
root = Path('${WORK}/long')
summary = dict(line.rstrip('\n').split('\t') for line in (root/'summary.tsv').read_text().splitlines()[1:])
assert summary['matched_rows'] == '4', summary
assert summary['matched_ids'] == '3', summary
assert summary['duplicate_time_keys'] == '0', summary
assert summary['time_multiplier'] == '10', summary
assert summary['nonmissing_volume'] == '4', summary
assert summary['nonmissing_concentration'] == '4', summary
audit = (root/'input_audit.json').read_text()
assert '"time_index_multiplier": 10' in audit
record = (root/'record_manifest.tsv').read_text().splitlines()
assert any('\t1\tR2\t2024-01-01\t0\t1\t3' in line for line in record), record
PY

set +e
"${BIN}" \
  --bfile "${BFILE}" \
  --input-file "${WORK}/raw_long.csv" \
  --input-model repeatability \
  --id-col pig_id \
  --time-col collect_date \
  --trait-cols volume \
  --time-parser date \
  --time-unique-policy strict \
  --prepare-input-only \
  --prepared-input-dir "${WORK}/strict" > "${WORK}/strict.log" 2>&1
STRICT_RC=$?
set -e
if [[ "${STRICT_RC}" -eq 0 ]]; then
  echo "ERROR: strict duplicate-time policy should have failed" >&2
  exit 1
fi
if ! grep -q "DUPLICATE_RECORD_KEY" "${WORK}/strict.log"; then
  echo "ERROR: strict duplicate-time failure did not include DUPLICATE_RECORD_KEY" >&2
  exit 1
fi

"${BIN}" \
  --bfile "${BFILE}" \
  --input-file "${WORK}/raw_long.csv" \
  --input-model repeatability \
  --id-col pig_id \
  --time-col collect_date \
  --trait-cols volume,concentration \
  --record-qcovars age,temp \
  --subject-ccovars farm,line \
  --time-parser date \
  --preview-schema \
  --prepared-input-dir "${WORK}/preview"

test -s "${WORK}/preview/schema_preview.tsv"
test -s "${WORK}/long/covariate_audit.tsv"
test -s "${WORK}/long/subject_covar_check.tsv"

"${BIN}" --input-file "${WORK}/raw_long.csv" --write-template-config "${WORK}/template_input_prepare.json"
test -s "${WORK}/template_input_prepare.json"

cat > "${WORK}/input_prepare.json" <<JSON
{
  "input_file": "${WORK}/raw_long.csv",
  "input_format": "csv",
  "input_model": "repeatability",
  "prepared_input_dir": "${WORK}/config",
  "id_col": "pig_id",
  "time_col": "collect_date",
  "trait_cols": ["volume", "concentration"],
  "record_qcovars": ["age", "temp"],
  "subject_ccovars": ["farm", "line"],
  "time_parser": "date",
  "time_unique_policy": "within-time-index",
  "time_index_multiplier": "auto"
}
JSON

"${BIN}" \
  --bfile "${BFILE}" \
  --input-config "${WORK}/input_prepare.json" \
  --prepare-input-only

python3 - <<PY
from pathlib import Path
root = Path('${WORK}/config')
summary = dict(line.rstrip('\n').split('\t') for line in (root/'summary.tsv').read_text().splitlines()[1:])
assert summary['matched_rows'] == '4', summary
assert (root/'covariate_audit.tsv').exists()
assert (root/'unresolved_columns.tsv').exists()
audit = (root/'input_audit.json').read_text()
assert '"mixed_subject_covariate_warnings"' in audit
PY

"${BIN}" \
  --bfile "${BFILE}" \
  --input-file "${WORK}/raw_long.csv" \
  --input-model repeatability \
  --infer-columns \
  --preview-schema \
  --prepared-input-dir "${WORK}/infer_only"
test -s "${WORK}/infer_only/column_inference.tsv"
test -s "${WORK}/infer_only/inferred_input_config.json"
python3 - <<PY
from pathlib import Path
text = Path('${WORK}/infer_only/column_inference.tsv').read_text()
assert 'pig_id\tid_col' in text, text
assert 'collect_date\ttime_col' in text, text
PY

python3 - <<PY
from pathlib import Path
work = Path('${WORK}')
with (work/'raw_rrm.csv').open('w') as out:
    out.write('pig_id,day,volume,temp,season\n')
    for i in range(1, 31):
        for day in range(1, 5):
            season = 'S1' if day <= 2 else 'S2'
            out.write(f'I{i:04d},{day},{10 + i * 0.1 + day * 0.2:.4f},{36 + day * 0.1:.2f},{season}\n')
with Path('${BFILE}.bim').open() as inp, (work/'extract_50_snps.txt').open('w') as out:
    for i, line in enumerate(inp):
        if i >= 50:
            break
        out.write(line.split()[1] + '\n')
PY

"${BIN}" \
  --bfile "${BFILE}" \
  --input-file "${WORK}/raw_rrm.csv" \
  --input-model rrm \
  --id-col pig_id \
  --time-col day \
  --trait-cols volume \
  --record-qcovars temp \
  --record-ccovars season \
  --time-parser integer \
  --prepare-input-only \
  --prepared-input-dir "${WORK}/rrm"

python3 - <<PY
from pathlib import Path
rrm = Path('${WORK}/rrm/rrm/volume.pheno.tsv').read_text().splitlines()
rep = Path('${WORK}/rrm/traits/volume.pheno.tsv').read_text().splitlines()
assert rrm[0] == 'FID\tIID\tTIME\tVALUE\ttemp\tseason'
assert len(rrm) == 121, len(rrm)
assert rrm[1].split('\t')[2] == '1', rrm[1]
assert rep[1].split('\t')[2] == '10', rep[1]
PY

"${BIN}" \
  --bfile "${BFILE}" \
  --extract "${WORK}/extract_50_snps.txt" \
  --prepared-input-dir "${WORK}/rrm" \
  --input-model rrm \
  --input-trait volume \
  --rrm-order 1 \
  --gwas \
  --no-loco \
  --out "${WORK}/rrm_smoke"

test -s "${WORK}/rrm_smoke.rrm_unbal.assoc.txt"
test -s "${WORK}/rrm_smoke.rrm_unbal.null.txt"

python3 - <<PY
from pathlib import Path
work = Path('${WORK}')
with (work/'raw_lmm.csv').open('w') as out:
    out.write('pig_id,day,volume,concentration,age,farm\n')
    for i in range(1, 31):
        out.write(f'I{i:04d},1,{10 + i * 0.1:.4f},{100 + i * 0.2:.4f},{200+i},F{i%3}\n')
PY

"${BIN}" \
  --bfile "${BFILE}" \
  --extract "${WORK}/extract_50_snps.txt" \
  --input-file "${WORK}/raw_lmm.csv" \
  --input-model lmm \
  --id-col pig_id \
  --time-col day \
  --trait-cols volume,concentration \
  --input-trait volume \
  --subject-qcovars age \
  --subject-ccovars farm \
  --time-parser integer \
  --gwas \
  --no-loco \
  --prepared-input-dir "${WORK}/lmm" \
  --out "${WORK}/lmm_smoke"
test -s "${WORK}/lmm_smoke.assoc.txt"

"${BIN}" \
  --bfile "${BFILE}" \
  --extract "${WORK}/extract_50_snps.txt" \
  --input-file "${WORK}/raw_lmm.csv" \
  --input-model multivariate \
  --id-col pig_id \
  --time-col day \
  --trait-cols volume,concentration \
  --subject-qcovars age \
  --subject-ccovars farm \
  --time-parser integer \
  --gwas \
  --no-loco \
  --prepared-input-dir "${WORK}/mv" \
  --out "${WORK}/mv_smoke"
test -s "${WORK}/mv_smoke.assoc.mv.txt"

cat > "${WORK}/raw_wide.csv" <<'CSV'
pig_id,age,farm,volume_d1,day2_volume,concentration_1,concentration_T2
I0001,300,A,10,11,210,215
I0002,290,B,8,9,180,185
CSV

"${BIN}" \
  --bfile "${BFILE}" \
  --input-file "${WORK}/raw_wide.csv" \
  --input-model repeatability \
  --id-col pig_id \
  --trait-cols volume_d1,day2_volume,concentration_1,concentration_T2 \
  --subject-qcovars age \
  --subject-ccovars farm \
  --time-parser integer \
  --wide-to-long \
  --prepare-input-only \
  --prepared-input-dir "${WORK}/wide"

python3 - <<PY
from pathlib import Path
root = Path('${WORK}/wide')
summary = dict(line.rstrip('\n').split('\t') for line in (root/'summary.tsv').read_text().splitlines()[1:])
assert summary['nonmissing_volume'] == '4', summary
assert summary['nonmissing_concentration'] == '4', summary
assert (root/'traits/volume.pheno.tsv').exists()
assert (root/'traits/concentration.pheno.tsv').exists()
assert sum(1 for _ in (root/'traits/volume.pheno.tsv').open()) == 5
PY

REAL_CSV="${COSMIC_REAL_CSV:-}"
REAL_BFILE="${COSMIC_REAL_BFILE:-${ROOT}/tests/repeatability_pig_semen_realdata_20260612/geno_subsets/4683_5k}"
REAL_BASE="${COSMIC_REAL_BASE:-${ROOT}/tests/repeatability_pig_semen_realdata_20260612/summary.tsv}"
if [[ -n "${COSMIC_REAL_CSV:-}${COSMIC_REAL_BFILE:-}${COSMIC_REAL_BASE:-}" || ( -f "${REAL_BFILE}.fam" && -f "${REAL_BASE}" ) ]]; then
  if [[ ! -f "${REAL_CSV}" || ! -f "${REAL_BFILE}.fam" || ! -f "${REAL_BASE}" ]]; then
    echo "Real-data validation requires readable COSMIC_REAL_CSV, COSMIC_REAL_BFILE (.fam), and COSMIC_REAL_BASE; refusing to silently skip a configured case." >&2
    exit 2
  fi
fi
if [[ -f "${REAL_CSV}" && -f "${REAL_BFILE}.fam" && -f "${REAL_BASE}" ]]; then
  "${BIN}" \
    --bfile "${REAL_BFILE}" \
    --input-file "${REAL_CSV}" \
    --input-model repeatability \
    --id-col id \
    --time-col date \
    --trait-cols vol,con,mot,Progressive,head_ABR,tail_ABR,tail,ABR,N_total,SUM \
    --record-qcovars age,age_month,semen_interval \
    --record-ccovars farm,breed,Season,Year,Month \
    --time-parser date \
    --prepare-input-only \
    --prepared-input-dir "${WORK}/realdata_cpp"

  python3 - <<PY
from pathlib import Path
base = dict(line.rstrip('\n').split('\t') for line in Path('${REAL_BASE}').read_text().splitlines()[1:])
new = dict(line.rstrip('\n').split('\t') for line in Path('${WORK}/realdata_cpp/summary.tsv').read_text().splitlines()[1:])
keys = ['matched_ids','matched_rows','duplicate_time_keys','time_multiplier','nonmissing_vol','nonmissing_con','nonmissing_mot','nonmissing_Progressive','nonmissing_head_ABR','nonmissing_tail_ABR','nonmissing_tail','nonmissing_ABR','nonmissing_N_total','nonmissing_SUM']
diff = {k:(base.get(k), new.get(k)) for k in keys if base.get(k) != new.get(k)}
assert not diff, diff
PY

  python3 - <<PY
from pathlib import Path
src = Path('${WORK}/realdata_cpp')
dst = Path('${WORK}/realdata_cpp_smoke')
(dst/'traits').mkdir(parents=True, exist_ok=True)
ids = []
with (src/'traits/vol.pheno.tsv').open() as handle:
    next(handle)
    for line in handle:
        iid = line.split('\t')[1]
        if iid not in ids:
            ids.append(iid)
        if len(ids) >= 50:
            break
keep = set(ids)
def subset_file(infile, outfile):
    with infile.open() as inp, outfile.open('w') as out:
        header = next(inp)
        out.write(header)
        for line in inp:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 2 and fields[1] in keep:
                out.write(line)
subset_file(src/'traits/vol.pheno.tsv', dst/'traits/vol.pheno.tsv')
subset_file(src/'model_qcovar.tsv', dst/'model_qcovar.tsv')
subset_file(src/'model_covar.tsv', dst/'model_covar.tsv')
subset_file(src/'record_manifest.tsv', dst/'record_manifest.tsv')
with Path('${REAL_BFILE}.bim').open() as inp, Path('${WORK}/extract_100_snps.txt').open('w') as out:
    for i, line in enumerate(inp):
        if i >= 100:
            break
        out.write(line.split()[1] + '\n')
PY

  "${BIN}" \
    --bfile "${REAL_BFILE}" \
    --extract "${WORK}/extract_100_snps.txt" \
    --prepared-input-dir "${WORK}/realdata_cpp_smoke" \
    --input-trait vol \
    --repeatability-model \
    --repeatability-mode fast \
    --gwas \
    --no-loco \
    --out "${WORK}/realdata_vol_smoke"

  python3 - <<PY
from pathlib import Path
import math
null = Path('${WORK}/realdata_vol_smoke.repeatability.null.txt')
assoc = Path('${WORK}/realdata_vol_smoke.assoc.txt')
assert null.exists() and null.stat().st_size > 0, null
assert assoc.exists() and assoc.stat().st_size > 0, assoc
lines = assoc.read_text().splitlines()
assert len(lines) > 1, 'assoc has no result rows'
header = lines[0].split('\t')
pidx = header.index('P') if 'P' in header else -1
assert pidx >= 0, header
finite = 0
for line in lines[1:200]:
    fields = line.split('\t')
    try:
        p = float(fields[pidx])
        finite += math.isfinite(p) and 0.0 <= p <= 1.0
    except Exception:
        pass
assert finite > 0, 'no finite p-values in first assoc rows'
PY
else
  echo "Skipping real pig semen regression: source CSV or genotype subset is unavailable."
fi

echo "Input adaptation tests passed. Results: ${WORK}"
