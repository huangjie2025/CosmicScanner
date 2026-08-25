#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${ROOT_DIR}/bin"

echo "CosmicGWAS installation check"
echo "Checking executables..."

if [[ ! -d "${BIN_DIR}" ]]; then
  echo "ERROR: bin/ directory not found under ${ROOT_DIR}" >&2
  exit 1
fi

if [[ ! -f "${BIN_DIR}/cosmicscanner" ]]; then
  echo "ERROR: bin/cosmicscanner not found under ${ROOT_DIR}" >&2
  exit 1
fi

chmod +x "${BIN_DIR}/cosmicscanner"

echo "Checking shared libraries..."
if [[ -x "${ROOT_DIR}/scripts/check_dependencies.sh" ]]; then
  "${ROOT_DIR}/scripts/check_dependencies.sh" "${ROOT_DIR}"
else
  echo "WARNING: scripts/check_dependencies.sh not found; skipping ldd checks." >&2
fi

echo "All required binaries are executable."
echo "To use CosmicGWAS, add the following line to your shell profile:"
echo "export PATH=${BIN_DIR}:\$PATH"
