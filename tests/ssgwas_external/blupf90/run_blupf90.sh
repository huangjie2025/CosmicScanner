#!/usr/bin/env bash
# Phase C external comparison runner (BLUPF90 + postGSf90).
# SKIP-safe: exits 0 with a notice when renf90/postGSf90 are absent.
#   bash run_blupf90.sh --bed PREFIX --ped PED [--pheno PHE] --out DIR
set -euo pipefail
cd "$(dirname "$0")"

VA=0.616532
VE=1.85018

if ! command -v renf90 >/dev/null 2>&1 || ! command -v postGSf90 >/dev/null 2>&1; then
  echo "SKIP: renf90/postGSf90 not found; BLUPF90 comparison not run (see README.md)"
  exit 0
fi

BED=${BED:-}; PED=${PED:-}; PHENO=${PHENO:-}; OUT=${OUT:-}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bed) BED=$2; shift 2;;
    --ped) PED=$2; shift 2;;
    --pheno) PHENO=$2; shift 2;;
    --out) OUT=$2; shift 2;;
    *) shift;;
  esac
done
mkdir -p "$OUT"
python3 renum_inputs.py --bed "$BED" --ped "$PED" --pheno "$PHENO" --out "$OUT"
N_SNP=$(wc -l < "$BED.bim")

cat > "$OUT/renf90.par" <<EOF
DATAFILE: input.dat
NUMBER_OF_TRAITS: 1
PEDIGREE: input.ped
PED_DEPTH: 3
MODEL: y = mu + animal
RANDOM: animal
OPTION solution animal
OPTION SNP_FILE 1 $N_SNP
OPTION fix_eff_alpha $VA
OPTION residual $VE
OPTION window 10
EOF

(cd "$OUT" && renf90 > renf90.log && postGSf90 > postGSf90.log)
echo "run ok: $OUT/solutions, $OUT/postGSf90.blp (see README.md for comparison protocol)"
