#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BIN_DIR="${ROOT_DIR}/bin"
LIB_DIR="${ROOT_DIR}/lib"

if [[ ! -d "${BIN_DIR}" ]]; then
  echo "ERROR: bin/ directory not found under ${ROOT_DIR}" >&2
  exit 1
fi

if ! command -v ldd >/dev/null 2>&1; then
  echo "ERROR: ldd not found; cannot check runtime dependencies." >&2
  exit 1
fi

targets=()
if [[ -x "${BIN_DIR}/cosmicscanner" ]]; then
  targets+=("${BIN_DIR}/cosmicscanner")
fi

if [[ "${#targets[@]}" -eq 0 ]]; then
  echo "ERROR: no executable targets found under ${BIN_DIR}" >&2
  exit 1
fi

status=0
for exe in "${targets[@]}"; do
  rel="${exe#${ROOT_DIR}/}"
  echo "===== ${rel} ====="
  output="$(LD_LIBRARY_PATH="${LIB_DIR}:${LD_LIBRARY_PATH:-}" ldd "${exe}" || true)"
  echo "${output}"
  if grep -q "not found" <<<"${output}"; then
    status=1
  fi
  echo
done

if [[ "${status}" -ne 0 ]]; then
  echo "ERROR: missing runtime dependencies detected." >&2
  exit 1
fi

echo "Dependency check passed."
