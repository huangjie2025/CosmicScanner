#!/usr/bin/env bash
# Build a portable Linux release tarball for CosmicGWAS/CosmicScanner.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -f "${ROOT_DIR}/VERSION" ]]; then
  DEFAULT_VERSION="$(tr -d '[:space:]' < "${ROOT_DIR}/VERSION")"
else
  DEFAULT_VERSION="v1.1.0"
fi

PACKAGE_VERSION="${COSMICGWAS_PACKAGE_VERSION:-${DEFAULT_VERSION}}"
PACKAGE_NAME="CosmicGWAS-${PACKAGE_VERSION}-linux-x86_64"
BUILD_DIR="${COSMICGWAS_BUILD_DIR:-${ROOT_DIR}/build-release}"
PACKAGE_ROOT="${COSMICGWAS_PACKAGE_ROOT:-${ROOT_DIR}/dist}"
PREFIX="${PACKAGE_ROOT}/${PACKAGE_NAME}"
TARBALL="${PACKAGE_ROOT}/${PACKAGE_NAME}.tar.gz"
CHECKSUM="${TARBALL}.sha256"
JOBS="${COSMICGWAS_BUILD_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

if [[ -z "${CMAKE_BIN:-}" ]]; then
  if command -v cmake >/dev/null 2>&1; then
    CMAKE_BIN="$(command -v cmake)"
  else
    echo "ERROR: no CMake found. Set CMAKE_BIN=/path/to/cmake and rerun." >&2
    exit 1
  fi
fi

CXX_BIN="${CXX:-}"
if [[ -z "${CXX_BIN}" ]]; then
  if command -v c++ >/dev/null 2>&1; then
    CXX_BIN="$(command -v c++)"
  else
    echo "ERROR: no C++ compiler found. Set CXX=/path/to/g++ and rerun." >&2
    exit 1
  fi
fi

rm -rf "${BUILD_DIR}" "${PREFIX}" "${TARBALL}" "${CHECKSUM}"
mkdir -p "${PACKAGE_ROOT}"

# Portable release build: disable BLAS/LAPACK by default so the binary does not
# accidentally depend on site-specific MKL/OpenBLAS libraries. Set
# COSMICGWAS_ENABLE_BLAS=ON to opt into accelerated BLAS/LAPACK and bundle the
# resolved runtime libraries when ldd can locate them.
CMAKE_BLAS_ARGS=()
if [[ "${COSMICGWAS_ENABLE_BLAS:-OFF}" != "ON" ]]; then
  CMAKE_BLAS_ARGS+=( -DCMAKE_DISABLE_FIND_PACKAGE_BLAS=ON
                     -DCMAKE_DISABLE_FIND_PACKAGE_LAPACK=ON )
  echo "Portable build: BLAS/LAPACK disabled (set COSMICGWAS_ENABLE_BLAS=ON to enable)."
else
  echo "BLAS/LAPACK enabled: MKL/OpenBLAS runtime libs will be bundled when present."
fi

CMAKE_PREFIX_ARGS=()
PREFIX_CANDIDATES=()
[[ -n "${CMAKE_PREFIX_PATH:-}" ]] && PREFIX_CANDIDATES+=("${CMAKE_PREFIX_PATH}")
[[ -n "${CONDA_PREFIX:-}" ]] && PREFIX_CANDIDATES+=("${CONDA_PREFIX}")
[[ -f "/disk228/huangjie/software/miniconda3/share/eigen3/cmake/Eigen3Config.cmake" ]] && PREFIX_CANDIDATES+=("/disk228/huangjie/software/miniconda3")
if [[ "${#PREFIX_CANDIDATES[@]}" -gt 0 ]]; then
  IFS=';'; CMAKE_PREFIX_JOINED="${PREFIX_CANDIDATES[*]}"; unset IFS
  CMAKE_PREFIX_ARGS+=( -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_JOINED}" )
fi

CXX="${CXX_BIN}" "${CMAKE_BIN}" -S "${ROOT_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release "${CMAKE_BLAS_ARGS[@]}" "${CMAKE_PREFIX_ARGS[@]}"

"${CMAKE_BIN}" --build "${BUILD_DIR}" --target cosmicscanner --parallel "${JOBS}"
"${CMAKE_BIN}" --install "${BUILD_DIR}" --prefix "${PREFIX}"

mkdir -p "${PREFIX}/lib" "${PREFIX}/docs" "${PREFIX}/examples" "${PREFIX}/scripts" "${PREFIX}/test"

for file in README.md CHANGELOG.md LICENSE VERSION install.sh; do
  if [[ -f "${ROOT_DIR}/${file}" ]]; then
    cp "${ROOT_DIR}/${file}" "${PREFIX}/${file}"
  fi
done
chmod +x "${PREFIX}/install.sh"

for script in check_dependencies.sh package_linux.sh; do
  if [[ -f "${ROOT_DIR}/scripts/${script}" ]]; then
    cp "${ROOT_DIR}/scripts/${script}" "${PREFIX}/scripts/${script}"
    chmod +x "${PREFIX}/scripts/${script}"
  fi
done

[[ -d "${ROOT_DIR}/docs" ]] && cp -a "${ROOT_DIR}/docs/." "${PREFIX}/docs/"

if [[ -d "${ROOT_DIR}/examples" ]]; then
  cp -a "${ROOT_DIR}/examples/." "${PREFIX}/examples/"
elif [[ -d "${ROOT_DIR}/../CosmicScanner_Release/example" ]]; then
  cp -a "${ROOT_DIR}/../CosmicScanner_Release/example/." "${PREFIX}/examples/"
fi

cat > "${PREFIX}/test/run_demo.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${ROOT_DIR}/bin/cosmicscanner"
OUT_DIR="${ROOT_DIR}/test/results"
mkdir -p "${OUT_DIR}"

"${BIN}" --help >/dev/null

if [[ -f "${ROOT_DIR}/examples/test_data.bed" ]]; then
  "${BIN}" \
    --bfile "${ROOT_DIR}/examples/test_data" \
    --pheno "${ROOT_DIR}/examples/test_data.pheno" \
    --gwas --no-loco \
    --out "${OUT_DIR}/demo"
  test -s "${OUT_DIR}/demo.assoc.txt"
fi

echo "CosmicGWAS demo passed."
EOF
chmod +x "${PREFIX}/test/run_demo.sh"

cat > "${PREFIX}/test/README.md" <<'EOF'
# CosmicGWAS Release Test

Run the packaged smoke test from the unpacked release directory:

```bash
bash test/run_demo.sh
```

The script checks `bin/cosmicscanner --help` and, when the bundled PLINK demo
data are present, runs a minimal no-LOCO GWAS into `test/results/`.
EOF

LDD_REPORT="${PREFIX}/docs/ldd_report.txt"

run_ldd_report() {
  local report_path="$1"
  : > "${report_path}"
  while IFS= read -r -d '' exe; do
    {
      echo "===== ${exe#${PREFIX}/} ====="
      LD_LIBRARY_PATH="${PREFIX}/lib:${LD_LIBRARY_PATH:-}" ldd "${exe}"
      echo
    } >> "${report_path}"
  done < <(find "${PREFIX}/bin" -maxdepth 1 -type f -perm -111 -print0 | sort -z)
}

run_ldd_report "${LDD_REPORT}"
if grep -q "not found" "${LDD_REPORT}"; then
  echo "ERROR: ldd found missing runtime libraries before bundling" >&2
  cat "${LDD_REPORT}" >&2
  exit 1
fi

copy_runtime_from_ldd() {
  local soname="$1"
  local lib_path
  lib_path="$(awk -v soname="${soname}" '$1 == soname { for (i=1; i<=NF; ++i) if ($i ~ /^\//) { print $i; exit } }' "${LDD_REPORT}")"
  if [[ -n "${lib_path}" && -f "${lib_path}" ]]; then
    cp -L "${lib_path}" "${PREFIX}/lib/${soname}"
  elif grep -q "${soname}" "${LDD_REPORT}"; then
    echo "WARNING: ${soname} is required but could not be copied; relying on target system runtime." >&2
  fi
}

copy_runtime_from_ldd "libgomp.so.1"
copy_runtime_from_ldd "libstdc++.so.6"
copy_runtime_from_ldd "libgcc_s.so.1"
copy_runtime_from_ldd "libzstd.so.1"
copy_runtime_from_ldd "libz.so.1"

if [[ "${COSMICGWAS_ENABLE_BLAS:-OFF}" == "ON" ]]; then
  copy_runtime_from_ldd "libmkl_intel_lp64.so.2"
  copy_runtime_from_ldd "libmkl_intel_thread.so.2"
  copy_runtime_from_ldd "libmkl_core.so.2"
  copy_runtime_from_ldd "libiomp5.so"
  copy_runtime_from_ldd "libopenblas.so.0"
  copy_runtime_from_ldd "libblas.so.3"
  copy_runtime_from_ldd "liblapack.so.3"
fi

run_ldd_report "${LDD_REPORT}"
if grep -q "not found" "${LDD_REPORT}"; then
  echo "ERROR: ldd found missing runtime libraries after bundling" >&2
  cat "${LDD_REPORT}" >&2
  exit 1
fi

tar -C "${PACKAGE_ROOT}" -czf "${TARBALL}" "${PACKAGE_NAME}"
(
  cd "${PACKAGE_ROOT}"
  sha256sum "${PACKAGE_NAME}.tar.gz" > "${PACKAGE_NAME}.tar.gz.sha256"
)

echo "Package directory: ${PREFIX}"
echo "Tarball: ${TARBALL}"
echo "SHA256: ${CHECKSUM}"
