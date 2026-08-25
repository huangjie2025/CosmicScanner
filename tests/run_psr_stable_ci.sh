#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "${1:-}" == "--build" ]]; then
  make -C "${ROOT}/build" cosmicscanner -j"${JOBS:-4}"
fi

echo "[PSR-GWAS] Running stable CI regression via core GWAS suite..."
bash "${ROOT}/tests/run_gwas_core_regression.sh" --fast --with-psr

echo "[PSR-GWAS] Stable CI regression passed."
