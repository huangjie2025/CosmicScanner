#!/usr/bin/env bash
# Build a portable Linux release tarball for CosmicScanner/CosmicScanner.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -f "${ROOT_DIR}/VERSION" ]]; then
  DEFAULT_VERSION="$(tr -d '[:space:]' < "${ROOT_DIR}/VERSION")"
else
  DEFAULT_VERSION="v1.1.0"
fi

PACKAGE_VERSION="${COSMICSCANNER_PACKAGE_VERSION:-${DEFAULT_VERSION}}"
PACKAGE_NAME="CosmicScanner-${PACKAGE_VERSION}-linux-x86_64"
BUILD_DIR="${COSMICSCANNER_BUILD_DIR:-${ROOT_DIR}/build-release}"
PACKAGE_ROOT="${COSMICSCANNER_PACKAGE_ROOT:-${ROOT_DIR}/dist}"
SOURCE_PACKAGE_DIR="${COSMICSCANNER_SOURCE_DIR:-${ROOT_DIR}/dist/CosmicScanner-${PACKAGE_VERSION}-source}"
PREFIX="${PACKAGE_ROOT}/${PACKAGE_NAME}"
TARBALL="${PACKAGE_ROOT}/${PACKAGE_NAME}.tar.gz"
CHECKSUM="${TARBALL}.sha256"
JOBS="${COSMICSCANNER_BUILD_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

if [[ ! -f "${SOURCE_PACKAGE_DIR}/CMakeLists.txt" ]]; then
  echo "ERROR: focused source package not found at ${SOURCE_PACKAGE_DIR}" >&2
  echo "Run scripts/package_source.sh before scripts/package_linux.sh." >&2
  exit 1
fi

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
# COSMICSCANNER_ENABLE_BLAS=ON to opt into accelerated BLAS/LAPACK and bundle the
# resolved runtime libraries when ldd can locate them.
CMAKE_BLAS_ARGS=()
if [[ "${COSMICSCANNER_ENABLE_BLAS:-OFF}" != "ON" ]]; then
  CMAKE_BLAS_ARGS+=( -DCOSMICSCANNER_ENABLE_BLAS=OFF )
  echo "Portable build: BLAS/LAPACK disabled (set COSMICSCANNER_ENABLE_BLAS=ON to enable)."
else
  CMAKE_BLAS_ARGS+=( -DCOSMICSCANNER_ENABLE_BLAS=ON )
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

CXX="${CXX_BIN}" "${CMAKE_BIN}" -S "${SOURCE_PACKAGE_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release "${CMAKE_BLAS_ARGS[@]}" "${CMAKE_PREFIX_ARGS[@]}"

"${CMAKE_BIN}" --build "${BUILD_DIR}" --target cosmicscanner --parallel "${JOBS}"
"${CMAKE_BIN}" --install "${BUILD_DIR}" --prefix "${PREFIX}"

mkdir -p "${PREFIX}/lib" "${PREFIX}/docs" "${PREFIX}/examples" "${PREFIX}/scripts" "${PREFIX}/test"

for file in CHANGELOG.md LICENSE VERSION; do
  if [[ -f "${SOURCE_PACKAGE_DIR}/${file}" ]]; then
    cp "${SOURCE_PACKAGE_DIR}/${file}" "${PREFIX}/${file}"
  fi
done
cp "${ROOT_DIR}/install.sh" "${PREFIX}/install.sh"
chmod +x "${PREFIX}/install.sh"

for script in check_dependencies.sh; do
  if [[ -f "${ROOT_DIR}/scripts/${script}" ]]; then
    cp "${ROOT_DIR}/scripts/${script}" "${PREFIX}/scripts/${script}"
    chmod +x "${PREFIX}/scripts/${script}"
  fi
done

cat > "${PREFIX}/README.md" <<'EOF'
# CosmicScanner v1.1.0

CosmicScanner is the first packaged Linux release of the CosmicScanner GWAS module.
This release is intentionally scoped to three Gaussian LMM GWAS workflows and
one relationship-matrix utility:

- Single-trait GWAS: standard LM/LMM with PLINK/BGEN/PGEN input and optional LOCO.
- Multivariate GWAS: multi-trait LMM, including EMmvGWAS-compatible and GEMMA-style validation routes.
- Repeatability GWAS: repeated-record mixed model for subjects with multiple observations.
- A/G/H matrices: A-inverse, genomic G, and H-inverse construction.

Other experimental or broader CosmicScanner model entries may still exist inside
the executable, but they are not part of this binary package's supported release
scope. Do not treat RRM, PSR/ordinal/count, ssGWAS association scans, EB-GWAS, survival, maternal,
GxE, rare-variant, BSLMM/DAP, or summary-stat utilities as production-supported
features in this v1.1.0 package.

## Quick Start

```bash
bash install.sh
export PATH=$(pwd)/bin:$PATH
bash test/run_demo.sh
```

## Main Entry

```bash
cosmicscanner --bfile examples/test_data --pheno examples/test_data.pheno \
  --mpheno 1 --gwas --scan-tier semi --out results/single_trait
```

See `docs/USER_MANUAL_EN.md` (English) and `docs/USER_MANUAL_CN.md` (中文) for the complete user manual.
EOF

cp -a "${SOURCE_PACKAGE_DIR}/docs/." "${PREFIX}/docs/"
cp -a "${SOURCE_PACKAGE_DIR}/examples/." "${PREFIX}/examples/"

cat > "${PREFIX}/test/run_demo.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${ROOT_DIR}/bin/cosmicscanner"
OUT_DIR="${COSMICSCANNER_VALIDATION_OUT:-${ROOT_DIR}/test/results}"
rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

"${BIN}" --help >/dev/null

if [[ -f "${ROOT_DIR}/examples/test_data.bed" ]]; then
  "${BIN}" \
    --bfile "${ROOT_DIR}/examples/test_data" \
    --pheno "${ROOT_DIR}/examples/test_data.pheno" \
    --mpheno 1 --gwas --scan-tier semi --no-loco \
    --chr 1 --from-bp 10000 --to-bp 20000 \
    --out "${OUT_DIR}/single_trait"
  test -s "${OUT_DIR}/single_trait.assoc.txt"
  test "$(wc -l < "${OUT_DIR}/single_trait.assoc.txt")" -eq 12

  "${BIN}" \
    --bfile "${ROOT_DIR}/examples/test_data" \
    --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
    --mpheno-list 1,2 \
    --multivariate --gwas --scan-tier semi --no-loco \
    --out "${OUT_DIR}/multivariate"
  test -s "${OUT_DIR}/multivariate.assoc.mv.txt"

  "${BIN}" \
    --bfile "${ROOT_DIR}/examples/test_data" \
    --pheno "${ROOT_DIR}/examples/demo_repeatability.pheno" \
    --qcovar "${ROOT_DIR}/examples/demo_repeatability.qcovar" \
    --tpheno 1 --mpheno 2 \
    --repeatability-model --repeatability-mode exact --gwas --scan-tier semi --no-loco \
    --out "${OUT_DIR}/repeatability_semi"
  test -s "${OUT_DIR}/repeatability_semi.assoc.txt"
  test -s "${OUT_DIR}/repeatability_semi.repeatability.null.txt"
  test -s "${OUT_DIR}/repeatability_semi.scan.meta.txt"

  "${BIN}" \
    --bfile "${ROOT_DIR}/examples/test_data" \
    --pheno "${ROOT_DIR}/examples/demo_repeatability.pheno" \
    --qcovar "${ROOT_DIR}/examples/demo_repeatability.qcovar" \
    --tpheno 1 --mpheno 2 \
    --repeatability-model --repeatability-mode fast --no-loco --gwas \
    --out "${OUT_DIR}/repeatability_fast"
  test -s "${OUT_DIR}/repeatability_fast.assoc.txt"

  "${BIN}" \
    --bfile "${ROOT_DIR}/examples/test_data" \
    --make-agh \
    --ped "${ROOT_DIR}/examples/demo.ped" \
    --genotyped-list "${ROOT_DIR}/examples/demo.genotyped.txt" \
    --out "${OUT_DIR}/agh"
  test -s "${OUT_DIR}/agh.Ainv.bin"
  test -s "${OUT_DIR}/agh.Ainv.id"
  test -s "${OUT_DIR}/agh.G.bin"
  test -s "${OUT_DIR}/agh.G.id"
  test -s "${OUT_DIR}/agh.Hinv.bin"
  test -s "${OUT_DIR}/agh.Hinv.id"
  test -s "${OUT_DIR}/agh.agh.meta.txt"
fi

echo "CosmicScanner single-trait, multivariate, repeatability, and A/G/H demos passed."
EOF
chmod +x "${PREFIX}/test/run_demo.sh"

cat > "${PREFIX}/test/README.md" <<'EOF'
# CosmicScanner Release Test

Run the packaged smoke test from the unpacked release directory:

```bash
bash test/run_demo.sh
```

The script checks `bin/cosmicscanner --help` and runs four packaged smoke tests:
single-trait GWAS, multivariate GWAS, repeatability GWAS, and A/G/H construction.
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

if [[ "${COSMICSCANNER_ENABLE_BLAS:-OFF}" == "ON" ]]; then
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
