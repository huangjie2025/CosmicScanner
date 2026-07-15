#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${COSMICSCANNER_BUILD_DIR:-${ROOT_DIR}/build}"
BUILD_TYPE="${COSMICSCANNER_BUILD_TYPE:-Release}"
JOBS="${COSMICSCANNER_BUILD_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
ENABLE_BLAS="${COSMICSCANNER_ENABLE_BLAS:-OFF}"

CMAKE_ARGS=(
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
  -DCOSMICSCANNER_ENABLE_BLAS="${ENABLE_BLAS}"
)
if [[ -n "${CMAKE_PREFIX_PATH:-}" ]]; then
  CMAKE_ARGS+=( -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" )
elif [[ -n "${CONDA_PREFIX:-}" ]]; then
  CMAKE_ARGS+=( -DCMAKE_PREFIX_PATH="${CONDA_PREFIX}" )
fi

cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" "${CMAKE_ARGS[@]}"
cmake --build "${BUILD_DIR}" --target cosmicscanner --parallel "${JOBS}"
EXE_SUFFIX=""
case "${OSTYPE:-}" in
  msys*|cygwin*) EXE_SUFFIX=".exe" ;;
esac
"${BUILD_DIR}/cosmicscanner${EXE_SUFFIX}" --help >/dev/null

echo "Build passed: ${BUILD_DIR}/cosmicscanner${EXE_SUFFIX}"
