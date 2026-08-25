#!/usr/bin/env bash
# ssGBLUP back-solving (Phase 6) validation orchestrator.
#   bash tests/run_ssgwas_backsolve.sh [--build]
set -euo pipefail
cd "$(dirname "$0")/.."

BIN=${BIN:-build/cosmicscanner}
OUT=${OUT:-${COSMIC_BACKSOLVE_OUT:-tests/results/ssgwas_backsolve}}
BUILD=0
if [[ "${1:-}" == "--build" ]]; then
  BUILD=1
fi

if [[ $BUILD -eq 1 ]]; then
  make -C build cosmicscanner -j4
fi

mkdir -p "$OUT"
run_py() {
  local name=$1; shift
  printf '\n=== %s ===\n' "$name"
  python3 "$@"
}

run_py "analytic (exact linear-algebra gate vs brute force)" \
  tests/run_ssgwas_backsolve_analytic.py --bin "$BIN" --outdir "$OUT/analytic"
run_py "simulation (causal recovery / typed+untyped / shuffle invariance)" \
  tests/run_ssgwas_backsolve_simulation.py --bin "$BIN" --outdir "$OUT/simulation"

printf '\nPASS: ssGBLUP back-solving full validation (analytic + simulation)\n'