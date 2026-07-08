#!/usr/bin/env bash
# Build a curated source tarball for the focused CosmicGWAS release.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -f "${ROOT_DIR}/VERSION" ]]; then
  DEFAULT_VERSION="$(tr -d '[:space:]' < "${ROOT_DIR}/VERSION")"
else
  DEFAULT_VERSION="v1.1.0"
fi

PACKAGE_VERSION="${COSMICGWAS_PACKAGE_VERSION:-${DEFAULT_VERSION}}"
PACKAGE_NAME="CosmicGWAS-${PACKAGE_VERSION}-source"
PACKAGE_ROOT="${COSMICGWAS_SOURCE_PACKAGE_ROOT:-${ROOT_DIR}/dist}"
PREFIX="${PACKAGE_ROOT}/${PACKAGE_NAME}"
TARBALL="${PACKAGE_ROOT}/${PACKAGE_NAME}.tar.gz"
CHECKSUM="${TARBALL}.sha256"

rm -rf "${PREFIX}" "${TARBALL}" "${CHECKSUM}"
mkdir -p "${PREFIX}" "${PACKAGE_ROOT}"

copy_file() {
  local src="$1"
  local dst="$2"
  if [[ -f "${src}" ]]; then
    mkdir -p "$(dirname "${dst}")"
    cp "${src}" "${dst}"
  fi
}

copy_dir() {
  local src="$1"
  local dst="$2"
  if [[ -d "${src}" ]]; then
    mkdir -p "${dst}"
    cp -a "${src}/." "${dst}/"
  fi
}

copy_file "${ROOT_DIR}/CMakeLists.txt" "${PREFIX}/CMakeLists.txt"
copy_file "${ROOT_DIR}/VERSION" "${PREFIX}/VERSION"
copy_file "${ROOT_DIR}/LICENSE" "${PREFIX}/LICENSE"
copy_file "${ROOT_DIR}/CHANGELOG.md" "${PREFIX}/CHANGELOG.md"
copy_file "${ROOT_DIR}/install.sh" "${PREFIX}/install.sh"
chmod +x "${PREFIX}/install.sh"
copy_dir "${ROOT_DIR}/src" "${PREFIX}/src"

mkdir -p "${PREFIX}/libCosmicCore"
copy_file "${ROOT_DIR}/libCosmicCore/CMakeLists.txt" "${PREFIX}/libCosmicCore/CMakeLists.txt"
copy_dir "${ROOT_DIR}/libCosmicCore/src" "${PREFIX}/libCosmicCore/src"

mkdir -p "${PREFIX}/scripts" "${PREFIX}/docs" "${PREFIX}/examples" "${PREFIX}/test"
copy_file "${ROOT_DIR}/scripts/package_linux.sh" "${PREFIX}/scripts/package_linux.sh"
copy_file "${ROOT_DIR}/scripts/check_dependencies.sh" "${PREFIX}/scripts/check_dependencies.sh"
copy_file "${ROOT_DIR}/scripts/package_source.sh" "${PREFIX}/scripts/package_source.sh"
chmod +x "${PREFIX}/scripts/"*.sh

cat > "${PREFIX}/README.md" <<'EOF'
# CosmicGWAS Source Package

This source package is curated for the first CosmicGWAS Linux release. It keeps
only the source files and minimal release assets needed to build and evaluate
the three supported workflows:

- Single-trait GWAS
- Multivariate GWAS
- Repeatability GWAS

The `src/` tree is intentionally kept intact so the checked-in build remains
consistent with the development executable. Broader model entry points may be
present in source files, but this source package documents and validates only
the three supported workflows above.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCOSMICSCANNER_BUILD_TESTS=OFF
cmake --build build --target cosmicscanner -j
./build/cosmicscanner --help
```

If CMake cannot find Eigen3, set `CMAKE_PREFIX_PATH` to the prefix containing
`share/eigen3/cmake/Eigen3Config.cmake`.

## Smoke Test

```bash
bash test/run_source_smoke.sh
```

The smoke test builds `cosmicscanner` with tests disabled and runs the packaged
single-trait, multivariate, and repeatability demo analyses.
EOF

cat > "${PREFIX}/docs/USAGE_CN.md" <<'EOF'
# CosmicGWAS v1.1.0 源码包说明

本源码包只面向三个首发板块：单性状 GWAS、多性状 GWAS、重复力 GWAS。
源码中的其他模型入口不作为本源码包的正式交付范围。

## 构建

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCOSMICSCANNER_BUILD_TESTS=OFF
cmake --build build --target cosmicscanner -j
```

## 单性状 GWAS

```bash
./build/cosmicscanner --bfile examples/test_data --pheno examples/test_data.pheno --gwas --out results/single_trait
```

## 多性状 GWAS

```bash
./build/cosmicscanner \
  --bfile examples/test_data \
  --pheno examples/demo_multivariate.pheno \
  --mpheno-list 1,2 \
  --multivariate --mv-method emmv --no-loco --gwas \
  --out results/multivariate
```

## 重复力 GWAS

```bash
./build/cosmicscanner \
  --bfile examples/test_data \
  --pheno examples/demo_repeatability.pheno \
  --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode exact --gwas \
  --out results/repeatability
```

包内快速 smoke 使用 `--repeatability-mode fast`，正式主分析建议使用 `exact`。
EOF

if [[ -f "${ROOT_DIR}/docs/REPEATABILITY_GWAS_RELEASE_CN.md" ]]; then
  cp "${ROOT_DIR}/docs/REPEATABILITY_GWAS_RELEASE_CN.md" "${PREFIX}/docs/REPEATABILITY_GWAS_RELEASE_CN.md"
fi

if [[ -d "${ROOT_DIR}/examples" ]]; then
  cp -a "${ROOT_DIR}/examples/." "${PREFIX}/examples/"
elif [[ -d "${ROOT_DIR}/../CosmicScanner_Release/example" ]]; then
  cp -a "${ROOT_DIR}/../CosmicScanner_Release/example/." "${PREFIX}/examples/"
fi

python3 - <<PY
from pathlib import Path
root = Path("${PREFIX}")
src = root / "examples" / "test_data.pheno"
mv = root / "examples" / "demo_multivariate.pheno"
rep = root / "examples" / "demo_repeatability.pheno"
rows = []
if src.exists():
    for line in src.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[2] not in {"NA", "-9", "nan", "NaN"}:
            try:
                y = float(parts[2])
            except ValueError:
                continue
            rows.append((parts[0], parts[1], y))
    mv.write_text("".join(f"{fid}\t{iid}\t{y:.8g}\t{0.5 * y + 0.1:.8g}\n" for fid, iid, y in rows))
    rep.write_text("".join(
        f"{fid}\t{iid}\t1\t{y:.8g}\n{fid}\t{iid}\t2\t{y + 0.05:.8g}\n"
        for fid, iid, y in rows[:120]
    ))
PY

cat > "${PREFIX}/test/run_source_smoke.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${COSMICGWAS_SOURCE_BUILD_DIR:-${ROOT_DIR}/build}"
OUT_DIR="${ROOT_DIR}/test/results"
JOBS="${COSMICGWAS_BUILD_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

mkdir -p "${OUT_DIR}"

CMAKE_PREFIX_ARGS=()
PREFIX_CANDIDATES=()
[[ -n "${CMAKE_PREFIX_PATH:-}" ]] && PREFIX_CANDIDATES+=("${CMAKE_PREFIX_PATH}")
[[ -n "${CONDA_PREFIX:-}" ]] && PREFIX_CANDIDATES+=("${CONDA_PREFIX}")
[[ -f "/disk228/huangjie/software/miniconda3/share/eigen3/cmake/Eigen3Config.cmake" ]] && PREFIX_CANDIDATES+=("/disk228/huangjie/software/miniconda3")
if [[ "${#PREFIX_CANDIDATES[@]}" -gt 0 ]]; then
  IFS=';'; CMAKE_PREFIX_JOINED="${PREFIX_CANDIDATES[*]}"; unset IFS
  CMAKE_PREFIX_ARGS+=( -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_JOINED}" )
fi

cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCOSMICSCANNER_BUILD_TESTS=OFF \
  "${CMAKE_PREFIX_ARGS[@]}"
cmake --build "${BUILD_DIR}" --target cosmicscanner --parallel "${JOBS}"

BIN="${BUILD_DIR}/cosmicscanner"
"${BIN}" --help >/dev/null

"${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/test_data.pheno" \
  --gwas --no-loco \
  --out "${OUT_DIR}/single_trait"
test -s "${OUT_DIR}/single_trait.assoc.txt"

"${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 \
  --multivariate --mv-method emmv --no-loco --gwas \
  --out "${OUT_DIR}/multivariate"
test -s "${OUT_DIR}/multivariate.assoc.mv.txt"

"${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_repeatability.pheno" \
  --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode fast --gwas \
  --out "${OUT_DIR}/repeatability"
test -s "${OUT_DIR}/repeatability.assoc.txt"

echo "CosmicGWAS source smoke passed."
EOF
chmod +x "${PREFIX}/test/run_source_smoke.sh"

tar -C "${PACKAGE_ROOT}" -czf "${TARBALL}" "${PACKAGE_NAME}"
(
  cd "${PACKAGE_ROOT}"
  sha256sum "${PACKAGE_NAME}.tar.gz" > "${PACKAGE_NAME}.tar.gz.sha256"
)

echo "Source package directory: ${PREFIX}"
echo "Source tarball: ${TARBALL}"
echo "SHA256: ${CHECKSUM}"
