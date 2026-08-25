#!/usr/bin/env bash
# CosmicScanner --maternal-snp-test (direct|maternal|joint) simulation validation.
# Runs the full scenario suite in tests/run_maternal_snp_effect_simulation.py:
#   A. Null          - no system inflation for maternal/joint (with maternal random effects)
#   B. Direct-only   - direct finds causal, maternal weaker, joint finds it
#   C. Maternal-only - maternal finds causal, direct weaker, joint finds it
#   D. Same-direction direct+maternal
#   E. Opposite-sign effects (joint recovers association hidden in marginals)
#   F. Direct-maternal covariance in data produces no spurious SNP signal
#   G. Missing dam genotype (10/30/50%) safely handled with usable-N reporting
#   H. Near-identical offspring/dam vectors -> joint COLLINEAR/rank diagnostics
# Usage:
#   bash tests/run_test_maternal_snp_effects.sh [--n N] [--m M] [--nf NF] [--seed S] [--out DIR]
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${ROOT}/build/cosmicscanner"
[ -x "${BIN}" ] || { echo "FAIL: ${BIN} not found; build first"; exit 2; }

N=250; M=1200; NF=120; SEED=7
OUT="${ROOT}/tests/results/maternal_snp_effects_sim"
EXTRA=()
while [ $# -gt 0 ]; do
  case "$1" in
    --n) N="$2"; shift 2;;
    --m) M="$2"; shift 2;;
    --nf) NF="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    *) EXTRA+=("$1"); shift;;
  esac
done

echo "=== Maternal SNP-effect simulation suite (n=$N m=$M nf=$NF seed=$SEED) ==="
python3 "${ROOT}/tests/run_maternal_snp_effect_simulation.py" \
  --bin "${BIN}" --out "${OUT}" --n "${N}" --m "${M}" --nf "${NF}" --seed "${SEED}" "${EXTRA[@]}"
status=$?
exit $status
