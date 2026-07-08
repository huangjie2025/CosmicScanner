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

for file in CHANGELOG.md LICENSE VERSION install.sh; do
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

cat > "${PREFIX}/README.md" <<'EOF'
# CosmicGWAS v1.1.0

CosmicGWAS is the first packaged Linux release of the CosmicScanner GWAS module.
This release is intentionally scoped to three production-ready GWAS workflows:

- Single-trait GWAS: standard LM/LMM with PLINK/BGEN/PGEN input and optional LOCO.
- Multivariate GWAS: multi-trait LMM, including EMmvGWAS-compatible and GEMMA-style validation routes.
- Repeatability GWAS: repeated-record mixed model for subjects with multiple observations.

Other experimental or broader CosmicScanner model entries may still exist inside
the executable, but they are not part of this binary package's supported release
scope. Do not treat RRM, PSR/ordinal/count, ssGWAS, EB-GWAS, survival, maternal,
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
cosmicscanner --bfile examples/test_data --pheno examples/test_data.pheno --gwas --out results/single_trait
```

See `docs/USAGE_CN.md` for the three supported workflows and release boundaries.
EOF

cat > "${PREFIX}/docs/USAGE_CN.md" <<'EOF'
# CosmicGWAS v1.1.0 三板块首发说明

本 Linux 二进制包只按三个板块交付和验收：单性状 GWAS、多性状 GWAS、重复力 GWAS。
源码中的其他模型入口不作为本首发包的正式支持范围。

## 1. 单性状 GWAS

适用场景：一个个体一条表型记录，连续性状，常规 LM/LMM GWAS。

```bash
bin/cosmicscanner \
  --bfile examples/test_data \
  --pheno examples/test_data.pheno \
  --gwas \
  --out results/single_trait
```

常用选项：

- `--no-loco`：关闭 LOCO，使用全局 GRM。
- `--lm 1` 或 `--lm-mode 4`：GEMMA-style 普通线性模型。
- `--lmm-gemma-mode 1|2|3|4`：GEMMA-style LMM Wald/LRT/Score/all。
- `--extract`、`--exclude`、`--chr`、`--from-bp`、`--to-bp`：限制 SNP 集合或基因组区间。

主要输出：`*.assoc.txt` 或 GEMMA-style 对应关联结果文件。

## 2. 多性状 GWAS

适用场景：同一批个体有多个连续表型，需要联合检验 SNP 对多个性状的效应。

表型文件格式：`FID IID TRAIT1 TRAIT2 ...`。

```bash
bin/cosmicscanner \
  --bfile examples/test_data \
  --pheno examples/demo_multivariate.pheno \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method emmv \
  --no-loco \
  --gwas \
  --out results/multivariate
```

推荐首发使用：

- `--mv-method emmv`：快速多性状主路线。
- `--emmvgwas`：EMmvGWAS-compatible 便捷入口，适合与 EMmvGWAS 路线对齐。
- `--mv-method gemma-full`：小数据或候选区间的 GEMMA-style 精确/参考路线。

主要输出：`*.assoc.mv.txt`。

## 3. 重复力 GWAS

适用场景：同一个体存在多条重复记录，例如多批次、多次采样、多次测定，但不强调动态曲线形状。

长格式表型文件：`FID IID TIME PHENO` 或 `FID IID PHENO`。

```bash
bin/cosmicscanner \
  --bfile examples/test_data \
  --pheno examples/demo_repeatability.pheno \
  --tpheno 1 \
  --mpheno 2 \
  --repeatability-model \
  --repeatability-mode exact \
  --gwas \
  --out results/repeatability
```

推荐首发使用：

- `--repeatability-mode exact`：正式主分析。
- `--repeatability-mode fast`：快速预筛，不能替代正式结果。
- `--repeatability-grm-cache <prefix>`：重复分析时缓存 subject-level GRM。

主要输出：

- `*.assoc.txt`：重复力 GWAS 关联结果。
- `*.repeatability.null.txt`：null model 方差组分、h2 和重复力摘要。

## 首发包不支持范围

以下模型/工具可能仍在可执行文件帮助信息中出现，但不作为本包正式交付内容：RRM/dynamic GWAS、PSR/ordinal/count、ssGWAS、EB-GWAS、survival、maternal、GxE、rare-variant、BSLMM/DAP、summary-stat 工具、biobank two-step 原型和输入适配平台。

这些功能留在源码开发线或后续专项发布包中，不建议在本 v1.1.0 包中作为生产结论使用。
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
fi

echo "CosmicGWAS single-trait, multivariate, and repeatability demos passed."
EOF
chmod +x "${PREFIX}/test/run_demo.sh"

cat > "${PREFIX}/test/README.md" <<'EOF'
# CosmicGWAS Release Test

Run the packaged smoke test from the unpacked release directory:

```bash
bash test/run_demo.sh
```

The script checks `bin/cosmicscanner --help` and runs three packaged smoke tests:
single-trait GWAS, multivariate GWAS, and repeatability GWAS.
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
