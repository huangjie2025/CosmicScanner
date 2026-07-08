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
cat > "${PREFIX}/CHANGELOG.md" <<'EOF'
# CosmicGWAS Source Package Changelog

## v1.1.0-source

- 首发源码包只保留单性状 GWAS、多性状 GWAS 和重复力 GWAS 三个板块。
- 源码目录统一为单一 `src/`，原核心库源码并入 `src/core/`，不再保留单独 `libCosmicCore/` 目录。
- 删除首发源码包无关模型和工具源码，包括 RRM、PSR/count/ordinal、ssGWAS、maternal、survival、GxE、rare-variant、biobank/REGENIE、BSLMM/DAP 和 summary-stat 工具。
- 新增 focused CLI 和 focused CMake，仅构建三板块所需入口、I/O、GRM、LMM/MvLMM、repeatability 和输出代码。
- 包内 smoke test 覆盖单性状、多性状和重复力 GWAS。
EOF
copy_file "${ROOT_DIR}/install.sh" "${PREFIX}/install.sh"
chmod +x "${PREFIX}/install.sh"
copy_dir "${ROOT_DIR}/src" "${PREFIX}/src"
copy_dir "${ROOT_DIR}/libCosmicCore/src" "${PREFIX}/src/core"

rm -rf "${PREFIX}/src/input"
rm -f \
  "${PREFIX}/src/gemma_calccor.cpp" "${PREFIX}/src/gemma_calccor.h" \
  "${PREFIX}/src/gemma_predictor.cpp" "${PREFIX}/src/gemma_predictor.h" \
  "${PREFIX}/src/gemma_summary_tools.cpp" "${PREFIX}/src/gemma_summary_tools.h" \
  "${PREFIX}/src/gwas_app.cpp" "${PREFIX}/src/gwas_app.h" \
  "${PREFIX}/src/gwas_biobank_kernels.cpp" "${PREFIX}/src/gwas_biobank_scan.cpp" \
  "${PREFIX}/src/gwas_rrm_scan.cpp" "${PREFIX}/src/gwas_scanner_rv.cpp" \
  "${PREFIX}/src/rare_variant_test.cpp" "${PREFIX}/src/rare_variant_test.h" \
  "${PREFIX}/src/ridge_regenie.cpp" "${PREFIX}/src/ridge_regenie.h" \
  "${PREFIX}/src/models/count_model.cpp" "${PREFIX}/src/models/count_model.h" \
  "${PREFIX}/src/models/gemma_bslmm_model.cpp" "${PREFIX}/src/models/gemma_bslmm_model.h" \
  "${PREFIX}/src/models/gemma_vc_model.cpp" "${PREFIX}/src/models/gemma_vc_model.h" \
  "${PREFIX}/src/models/gxe_model.cpp" "${PREFIX}/src/models/gxe_model.h" \
  "${PREFIX}/src/models/maternal_model.cpp" "${PREFIX}/src/models/maternal_model.h" \
  "${PREFIX}/src/models/ordinal_model.cpp" "${PREFIX}/src/models/ordinal_model.h" \
  "${PREFIX}/src/models/rrm_model.cpp" "${PREFIX}/src/models/rrm_model.h" \
  "${PREFIX}/src/models/single_step_model.cpp" "${PREFIX}/src/models/single_step_model.h" \
  "${PREFIX}/src/models/survival_model.cpp" "${PREFIX}/src/models/survival_model.h" \
  "${PREFIX}/src/models/threshold_model.cpp" "${PREFIX}/src/models/threshold_model.h" \
  "${PREFIX}/src/models/multivariate_gemma_gxe.cpp" \
  "${PREFIX}/src/models/multivariate_gemma_refit.cpp" \
  "${PREFIX}/src/models/multivariate_hybrid.cpp"

cat > "${PREFIX}/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.15)
project(CosmicGWAS VERSION 1.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(GNUInstallDirs)

find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)
find_package(zstd REQUIRED)
option(COSMICGWAS_ENABLE_BLAS "Enable BLAS/LAPACK acceleration" OFF)

file(GLOB CORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/core/*.cpp")
add_library(cosmic_core STATIC ${CORE_SOURCES})
target_include_directories(cosmic_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/core)
target_link_libraries(cosmic_core PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX ZLIB::ZLIB zstd::libzstd)

if(COSMICGWAS_ENABLE_BLAS)
  find_package(BLAS)
  if(BLAS_FOUND)
    target_link_libraries(cosmic_core PUBLIC ${BLAS_LIBRARIES})
    target_compile_definitions(cosmic_core PUBLIC EIGEN_USE_BLAS)
  endif()

  find_package(LAPACK)
  if(LAPACK_FOUND)
    target_link_libraries(cosmic_core PUBLIC ${LAPACK_LIBRARIES})
    target_compile_definitions(cosmic_core PUBLIC EIGEN_USE_LAPACKE)
  endif()
endif()

set(GWAS_SOURCES
  src/main.cpp
  src/gwas_output.cpp
  src/gwas_scanner.cpp
  src/lmm_engine.cpp
  src/mvlmm_engine.cpp
  src/spa_test.cpp
  src/models/univariate_model.cpp
  src/models/multivariate_model.cpp
  src/models/multivariate_output.cpp
  src/models/repeatability_cache.cpp
  src/models/repeatability_model.cpp
  src/models/repeatability_output.cpp
  src/models/repeatability_refit_validation.cpp
  src/models/repeatability_scan_exact.cpp
  src/models/repeatability_scan_fast.cpp
  src/models/repeatability_scan_stcg.cpp
  src/models/repeatability_utils.cpp
)

add_executable(cosmicscanner ${GWAS_SOURCES})
target_include_directories(cosmicscanner PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/models
  ${CMAKE_CURRENT_SOURCE_DIR}/src/core
)
target_link_libraries(cosmicscanner PRIVATE cosmic_core OpenMP::OpenMP_CXX)
target_compile_definitions(cosmicscanner PRIVATE COSMICSCANNER_VERSION="${PROJECT_VERSION}")

install(TARGETS cosmicscanner RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
EOF

cat > "${PREFIX}/src/config.h" <<'EOF'
#pragma once
#include <string>
#include <vector>

namespace cosmic {

struct GwasConfig {
    std::string bfile_prefix;
    std::string bgen_file;
    std::string sample_file;
    std::string pfile_prefix;
    std::string pgen_file;
    std::string pvar_file;
    std::string psam_file;
    std::string out_prefix = "cosmic_gwas";

    std::string pheno_file;
    std::string covar_file;
    std::string qcovar_file;
    int mpheno = 1;
    int tpheno = 2;
    bool tpheno_set = false;
    std::vector<int> mpheno_list;

    bool run_gwas = false;
    bool no_loco = false;
    bool fast_lmm = false;
    bool linear_model = false;
    int linear_model_mode = 4;
    int lmm_gemma_mode = 0;
    bool lmm_gemma_refit = false;
    bool repeatability_model = false;
    std::string repeatability_mode = "auto";
    std::string repeatability_grm_cache_prefix;
    std::string repeatability_weight_file;
    bool repeatability_mixed_precision = false;
    bool repeatability_refit_validation = false;
    int repeatability_refit_limit = 0;
    std::string vce_refine_mode = "auto";
    bool multivariate = false;
    std::string multivariate_method = "emmv";
    std::string multivariate_compatible;
    std::string multivariate_tests = "wald,lrt,score";
    int multivariate_gemma_mode = 4;
    double multivariate_pnr = 0.001;
    int multivariate_emi = 10000;
    int multivariate_nri = 100;
    double multivariate_emp = 0.0001;
    double multivariate_nrp = 0.0001;
    bool multivariate_crt = false;
    int multivariate_refine_top = 100;
    double multivariate_refine_p = 1e-3;
    int multivariate_refine_discrepancy_top = 0;
    double multivariate_refine_discrepancy = 1.0;
    int multivariate_refine_score_prepass = 0;
    double multivariate_refit_p = 1.0;
    bool multivariate_use_aireml = false;
    int multivariate_em_min_iter = 10;
    int multivariate_ai_max_iter = 50;
    double multivariate_ai_tol = 1e-8;
    bool multivariate_output_vbeta = false;
    bool multivariate_emmvgwas_output = false;
    bool multivariate_emmvgwas_vc_only = false;
    std::string multivariate_emmvgwas_G_file;
    std::string multivariate_emmvgwas_R_file;
    bool multivariate_eb = false;
    int env_covar_col = -1;

    int grm_alg = 1;
    double min_maf = 0.01;
    double min_hwe = 0.0;
    std::string lmm_null_solver = "brent";
    int reml_maxit = 50;
    double reml_tol = 1e-6;
    double reml_vg = -1.0;
    double reml_vpe = -1.0;
    double reml_ve = -1.0;
    bool use_covariate_projection = true;
    bool out_ssf = false;
    int num_threads = 0;
    std::string extract_file;
    std::string exclude_file;
    std::string scan_chr;
    long long scan_from_bp = -1;
    long long scan_to_bp = -1;

    // Unsupported fields kept only so retained implementation files compile.
    bool eb_lmm = false;
    bool use_spa = false;
    bool run_skat = false;
    bool run_burden = false;
    bool biobank_mode = false;
    int block_size = 1000;
    std::string set_file;
    std::string regenie_model_snps;
    std::string regenie_write_pred;
    std::string regenie_use_pred;
    int regenie_step = 0;
    bool make_grm_sparse = false;
    bool force_slow = false;
    double sparse_cutoff = 0.05;
    std::string snp_weight_file;
    std::string pop_class_file;
};

}  // namespace cosmic
EOF

cat > "${PREFIX}/src/main.cpp" <<'EOF'
#include <chrono>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/resource.h>
#include <omp.h>

#include "config.h"
#include "data_loader.h"
#include "models/multivariate_model.h"
#include "models/repeatability_model.h"
#include "models/univariate_model.h"

namespace {

std::vector<int> parse_int_list(const std::string& s) {
    std::vector<int> out;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (!token.empty()) out.push_back(std::stoi(token));
    }
    return out;
}

void print_usage(const char* prog) {
    std::cout << "CosmicGWAS focused source build\n"
              << "Usage: " << prog << " --bfile <prefix> --pheno <file> --gwas [options]\n\n"
              << "Input:\n"
              << "  --bfile <prefix>             PLINK BED/BIM/FAM prefix\n"
              << "  --pfile <prefix>             PLINK2 PGEN/PVAR/PSAM prefix\n"
              << "  --bgen <file> --sample <f>   BGEN input and sample file\n"
              << "  --pheno <file>               Phenotype file\n"
              << "  --qcovar <file>              Quantitative covariates\n"
              << "  --covar <file>               Categorical covariates\n"
              << "  --out <prefix>               Output prefix\n\n"
              << "Single-trait GWAS:\n"
              << "  --gwas [--no-loco] [--fast-lmm]\n"
              << "  --lm <1|2|3|4>               Linear model Wald/LRT/Score/all\n"
              << "  --lmm-gemma-mode <1|2|3|4>   LMM Wald/LRT/Score/all\n\n"
              << "Multivariate GWAS:\n"
              << "  --multivariate --mpheno-list 1,2[,3...] --mv-method emmv\n\n"
              << "Repeatability GWAS:\n"
              << "  --repeatability-model --tpheno <time_col> --mpheno <trait_col>\n"
              << "  --repeatability-mode exact|fast|auto\n\n"
              << "Common filters:\n"
              << "  --maf <x> --hwe <x> --extract <file> --exclude <file>\n"
              << "  --chr <chr> --from-bp <bp> --to-bp <bp> --threads <n>\n";
}

}  // namespace

int main(int argc, char** argv) {
    using namespace cosmic;
    GwasConfig cfg;

    try {
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            auto need = [&](const char* name) -> std::string {
                if (i + 1 >= argc) throw std::runtime_error(std::string("Missing value for ") + name);
                return argv[++i];
            };
            if (a == "--help" || a == "-h") { print_usage(argv[0]); return 0; }
            else if (a == "--bfile") cfg.bfile_prefix = need("--bfile");
            else if (a == "--pfile") cfg.pfile_prefix = need("--pfile");
            else if (a == "--pgen") cfg.pgen_file = need("--pgen");
            else if (a == "--pvar") cfg.pvar_file = need("--pvar");
            else if (a == "--psam") cfg.psam_file = need("--psam");
            else if (a == "--bgen") cfg.bgen_file = need("--bgen");
            else if (a == "--sample") cfg.sample_file = need("--sample");
            else if (a == "--pheno" || a == "-p") cfg.pheno_file = need("--pheno");
            else if (a == "--qcovar") cfg.qcovar_file = need("--qcovar");
            else if (a == "--covar") cfg.covar_file = need("--covar");
            else if (a == "--out" || a == "-o") cfg.out_prefix = need("--out");
            else if (a == "--mpheno" || a == "-n") cfg.mpheno = std::stoi(need("--mpheno"));
            else if (a == "--tpheno") { cfg.tpheno = std::stoi(need("--tpheno")); cfg.tpheno_set = true; }
            else if (a == "--mpheno-list") cfg.mpheno_list = parse_int_list(need("--mpheno-list"));
            else if (a == "--gwas") cfg.run_gwas = true;
            else if (a == "--no-loco") cfg.no_loco = true;
            else if (a == "--fast-lmm") cfg.fast_lmm = true;
            else if (a == "--lm" || a == "-lm" || a == "--lm-mode") { cfg.linear_model = true; cfg.linear_model_mode = std::stoi(need(a.c_str())); }
            else if (a == "--lmm-gemma-mode" || a == "-lmm") cfg.lmm_gemma_mode = std::stoi(need(a.c_str()));
            else if (a == "--lmm-gemma-refit") cfg.lmm_gemma_refit = true;
            else if (a == "--multivariate") cfg.multivariate = true;
            else if (a == "--mv-method") cfg.multivariate_method = need("--mv-method");
            else if (a == "--mv-output-vbeta") cfg.multivariate_output_vbeta = true;
            else if (a == "--emmvgwas") { cfg.multivariate = true; cfg.multivariate_method = "emmv"; cfg.multivariate_compatible = "emmvgwas"; cfg.no_loco = true; }
            else if (a == "--repeatability-model") cfg.repeatability_model = true;
            else if (a == "--repeatability-mode") cfg.repeatability_mode = need("--repeatability-mode");
            else if (a == "--repeatability-grm-cache") cfg.repeatability_grm_cache_prefix = need("--repeatability-grm-cache");
            else if (a == "--repeatability-weight") cfg.repeatability_weight_file = need("--repeatability-weight");
            else if (a == "--maf") cfg.min_maf = std::stod(need("--maf"));
            else if (a == "--hwe") cfg.min_hwe = std::stod(need("--hwe"));
            else if (a == "--extract") cfg.extract_file = need("--extract");
            else if (a == "--exclude") cfg.exclude_file = need("--exclude");
            else if (a == "--chr") cfg.scan_chr = need("--chr");
            else if (a == "--from-bp") cfg.scan_from_bp = std::stoll(need("--from-bp"));
            else if (a == "--to-bp") cfg.scan_to_bp = std::stoll(need("--to-bp"));
            else if (a == "--threads") cfg.num_threads = std::stoi(need("--threads"));
            else if (a == "--lmm-null-solver") cfg.lmm_null_solver = need("--lmm-null-solver");
            else if (a == "--reml-maxit") cfg.reml_maxit = std::stoi(need("--reml-maxit"));
            else if (a == "--reml-tol") cfg.reml_tol = std::stod(need("--reml-tol"));
            else throw std::runtime_error("Unsupported option in focused source build: " + a);
        }

        if (!cfg.run_gwas) throw std::runtime_error("This focused build requires --gwas.");
        if (cfg.num_threads > 0) omp_set_num_threads(cfg.num_threads);

        #ifndef COSMICSCANNER_VERSION
        #define COSMICSCANNER_VERSION "unknown"
        #endif
        std::cout << "CosmicGWAS focused source v" << COSMICSCANNER_VERSION << "\n";

        PlinkReader plink_reader;
        BgenReader bgen_reader;
        PgenReader pgen_reader;
        auto meta = DataLoader::loadGenotype(cfg.bfile_prefix, cfg.bgen_file, cfg.sample_file,
                                             cfg.pfile_prefix, cfg.pgen_file, cfg.pvar_file, cfg.psam_file,
                                             plink_reader, bgen_reader, pgen_reader);

        auto start = std::chrono::high_resolution_clock::now();
        if (cfg.repeatability_model) {
            auto opts = RepeatabilityModel::Options::fromConfig(cfg);
            RepeatabilityModel model(opts);
            model.prepare(meta.fids, meta.iids, meta.n_samples, plink_reader, bgen_reader, pgen_reader, meta.is_bgen, meta.is_pgen);
            model.run();
        } else if (cfg.multivariate) {
            if (cfg.multivariate_method != "emmv" && cfg.multivariate_method != "wald") {
                throw std::runtime_error("Focused source build supports --mv-method emmv or wald only.");
            }
            auto opts = MultivariateModel::Options::fromConfig(cfg);
            MultivariateModel model(opts);
            model.prepare(meta.fids, meta.iids, meta.n_samples, plink_reader, bgen_reader, pgen_reader, meta.is_bgen, meta.is_pgen);
            model.run();
        } else {
            auto opts = UnivariateModel::Options::fromConfig(cfg);
            UnivariateModel model(opts);
            model.prepare(meta.fids, meta.iids, meta.n_samples, plink_reader, bgen_reader, pgen_reader, meta.is_bgen, meta.is_pgen);
            model.run();
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        std::cout << "Total computation time: " << elapsed.count() << " seconds\n";
        std::cout << "Peak memory usage: " << (usage.ru_maxrss / 1024.0) << " MB\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
EOF

python3 - "${PREFIX}" <<'PY'
from pathlib import Path
root = Path(__import__('sys').argv[1])
uv = root / "src/models/univariate_model.cpp"
text = uv.read_text()
old = '''    if (options_.run_skat || options_.run_burden) {
        scanner.runRareVariantScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_, keep_);
    } else if (options_.biobank_mode) {
        scanner.runBiobankScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_, keep_);
    } else {
        scanner.runLocoScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_, keep_);
    }
'''
new = '''    if (options_.run_skat || options_.run_burden || options_.biobank_mode) {
        throw std::runtime_error("Focused source build supports standard single-trait LM/LMM only; rare-variant and biobank routes are not included.");
    }
    scanner.runLocoScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_, keep_);
'''
uv.write_text(text.replace(old, new))

mv = root / "src/models/multivariate_model.cpp"
text = mv.read_text()
old = '''    if (options_.use_gemma_gxe) {
        runGemmaGxe();
        return;
    }

    if (options_.use_hybrid) {
        runHybrid();
        return;
    }

    if (options_.use_gemma_refit) {
        runGemmaRefit();
        return;
    }
'''
new = '''    if (options_.use_gemma_gxe || options_.use_hybrid || options_.use_gemma_refit) {
        throw std::runtime_error("Focused source build supports multivariate emmv/wald routes only.");
    }
'''
mv.write_text(text.replace(old, new))

hdr = root / "src/gwas_scanner.h"
text = hdr.read_text()
text = text.replace('#include "rare_variant_test.h"\n', '')
text = text.replace('#include "ridge_regenie.h"\n', '')
for block in [
'''    // Run Biobank mode scan (Matrix-Free RHE + PCG)
    void runBiobankScan(
        const Eigen::VectorXd& y,
        const Eigen::MatrixXf& X,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        bool is_bgen,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );

''',
'''    // Run Rare Variant Association Test (SKAT / Burden)
    void runRareVariantScan(
        const Eigen::VectorXd& y,
        const Eigen::MatrixXf& X,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        bool is_bgen,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );

''',
'''    // Run Random Regression Model (RRM) scan for longitudinal data
    void runRrmScan(
        const std::vector<PhenoReader::PhenoRecord>& records,
        int rrm_order,
        const Eigen::MatrixXf& X_base,
        const std::vector<std::string>& fids,
        const std::vector<std::string>& iids,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        bool is_bgen,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );

''',
'''    // Run unbalanced longitudinal (RRM) scan with pre-built fixed effects per observation
    void runRrmUnbalancedScan(
        const std::vector<PhenoReader::PhenoRecord>& records,
        int rrm_order,
        const Eigen::MatrixXf& X_base,
        const Eigen::MatrixXf& X_tdep,         // time-dependent covariates (n_obs x n_tdep)
        const std::vector<std::string>& tdep_names,
        const std::vector<std::string>& fids,
        const std::vector<std::string>& iids,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        bool is_bgen,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );

''']:
    text = text.replace(block, '')
hdr.write_text(text)
PY

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
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
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
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
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
