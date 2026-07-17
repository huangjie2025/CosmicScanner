#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION="$(tr -d '[:space:]' < "${ROOT_DIR}/VERSION")"
COMMIT="$(git -C "${ROOT_DIR}" rev-parse --short=12 HEAD)"
OUT_DIR="${COSMICSCANNER_DIST_DIR:-${ROOT_DIR}/dist}"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/cosmicscanner_source_pkg.XXXXXX")"
PKG_NAME="CosmicScanner-${VERSION}-source"

mkdir -p "${OUT_DIR}" "${TMP_DIR}/${PKG_NAME}"
git -C "${ROOT_DIR}" archive --format=tar HEAD | tar -x -C "${TMP_DIR}/${PKG_NAME}"
printf '%s\n' "${COMMIT}" >"${TMP_DIR}/${PKG_NAME}/COSMICSCANNER_GIT_COMMIT"
tar -czf "${OUT_DIR}/${PKG_NAME}.tar.gz" -C "${TMP_DIR}" "${PKG_NAME}"
echo "${OUT_DIR}/${PKG_NAME}.tar.gz"
