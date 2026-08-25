#!/usr/bin/env bash
# Build a curated source tarball for the focused CosmicScanner release.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -f "${ROOT_DIR}/VERSION" ]]; then
  DEFAULT_VERSION="$(tr -d '[:space:]' < "${ROOT_DIR}/VERSION")"
else
  DEFAULT_VERSION="v1.1.0"
fi

PACKAGE_VERSION="${COSMICSCANNER_PACKAGE_VERSION:-${DEFAULT_VERSION}}"
PACKAGE_NAME="CosmicScanner-${PACKAGE_VERSION}-source"
PACKAGE_ROOT="${COSMICSCANNER_SOURCE_PACKAGE_ROOT:-${ROOT_DIR}/dist}"
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
# RELEASE_NOTES intentionally excluded; release notes are in the GitHub Release

cat > "${PREFIX}/CHANGELOG.md" <<'EOF'
# CosmicScanner Source Package Changelog

## Version history note

- `v1.1.0` is the first public CosmicScanner focused release.
- `1.0.0` and `1.0.5` below are retrospective internal milestones inferred from
  Git history, not earlier public releases or Git tags.

## v1.1.0-source (2026-07-14)

- 首发源码包只保留单性状 GWAS、多性状 GWAS、重复力 GWAS 和 A/G/H 矩阵构建四项功能。
- 源码目录统一为单一扁平 `src/`，不再保留 `libCosmicCore/`、`src/core/` 或
  `src/models/` 子目录。
- 删除首发源码包无关模型和工具源码，包括 RRM、PSR/count/ordinal、ssGWAS
  关联扫描、maternal、survival、GxE、rare-variant、biobank/REGENIE、
  BSLMM/DAP 和 summary-stat 工具。
- 新增 focused CLI 和 focused CMake，仅构建四项功能所需入口、I/O、GRM、
  LMM/MvLMM、repeatability、A/G/H 和输出代码。
- 包内 release validation 覆盖单性状、多性状、重复力 semi/fast 和 A/G/H 构建。
- 重复力验证包含个体内中心化记录级协变量，覆盖 VCE 的零右端项求解回归。
- 默认 repeatability Woodbury 路径跳过无用的完整 GRM 求逆，并在 null summary
  报告实际分析受试者数。
- 新增文章范围、方法、输入输出、验证、复现、源码索引和引用元数据。

## 1.0.5 - retrospective internal milestone (2026-05-22 to 2026-06-24)

- Dedicated repeatability model with exact and fast analysis routes.
- Expanded multivariate compatibility and reference-validation paths.
- Consolidated single-trait, multivariate, and repeatability regression tests.

## 1.0.0 - retrospective internal milestone (2026-04-07)

- Independent GWAS baseline with single-trait LMM/LOCO and an initial
  multivariate LMM path.
- Early repeated-measure analysis used subject-level aggregation rather than
  the later exact record-level model.
EOF

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
  "${PREFIX}/src/h_matrix.cpp" "${PREFIX}/src/h_matrix.h" \
  "${PREFIX}/src/apy.cpp" "${PREFIX}/src/apy.h" \
  "${PREFIX}/src/rare_variant_test.cpp" "${PREFIX}/src/rare_variant_test.h" \
  "${PREFIX}/src/ridge_regenie.cpp" "${PREFIX}/src/ridge_regenie.h" \
  "${PREFIX}/src/models/count_model.cpp" "${PREFIX}/src/models/count_model.h" \
  "${PREFIX}/src/models/gemma_bslmm_model.cpp" "${PREFIX}/src/models/gemma_bslmm_model.h" \
  "${PREFIX}/src/models/gemma_vc_model.cpp" "${PREFIX}/src/models/gemma_vc_model.h" \
  "${PREFIX}/src/models/gxe_model.cpp" "${PREFIX}/src/models/gxe_model.h" \
  "${PREFIX}/src/models/maternal_model.cpp" "${PREFIX}/src/models/maternal_model.h" \
  "${PREFIX}/src/models/ordinal_model.cpp" "${PREFIX}/src/models/ordinal_model.h" \
  "${PREFIX}/src/models/rrm_model.cpp" "${PREFIX}/src/models/rrm_model.h" \
  "${PREFIX}/src/models/survival_model.cpp" "${PREFIX}/src/models/survival_model.h" \
  "${PREFIX}/src/models/threshold_model.cpp" "${PREFIX}/src/models/threshold_model.h" \
  "${PREFIX}/src/models/multivariate_gemma_gxe.cpp" \
  "${PREFIX}/src/models/multivariate_hybrid.cpp"

if [[ -d "${PREFIX}/src/core" ]]; then
  cp -a "${PREFIX}/src/core/." "${PREFIX}/src/"
  rm -rf "${PREFIX}/src/core"
fi
if [[ -d "${PREFIX}/src/models" ]]; then
  cp -a "${PREFIX}/src/models/." "${PREFIX}/src/"
  rm -rf "${PREFIX}/src/models"
fi
rm -f \
  "${PREFIX}/src/h_matrix.cpp" "${PREFIX}/src/h_matrix.h" \
  "${PREFIX}/src/apy.cpp" "${PREFIX}/src/apy.h" \
  "${PREFIX}/src/allele_harmonizer.cpp" "${PREFIX}/src/allele_harmonizer.h" \
  "${PREFIX}/src/genotype_converter.cpp" "${PREFIX}/src/genotype_converter.h" \
  "${PREFIX}/src/grm_factory.cpp" "${PREFIX}/src/grm_factory.h" \
  "${PREFIX}/src/matrix_inversion.cpp" "${PREFIX}/src/matrix_inversion.h" \
  "${PREFIX}/src/qc_engine.cpp" "${PREFIX}/src/qc_engine.h" \
  "${PREFIX}/src/schur.cpp" "${PREFIX}/src/schur.h" \
  "${PREFIX}/src/spa_test.cpp" "${PREFIX}/src/spa_test.h" \
  "${PREFIX}/src/sparse_grm_reader.cpp" "${PREFIX}/src/sparse_grm_reader.h" \
  "${PREFIX}/src/string_utils.h" \
  "${PREFIX}/src/upper_packed_binary.cpp" "${PREFIX}/src/upper_packed_binary.h"

cat > "${PREFIX}/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.15)
project(CosmicScanner VERSION 1.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(GNUInstallDirs)

find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)
find_package(zstd REQUIRED)
if(TARGET zstd::libzstd)
  set(COSMICSCANNER_ZSTD_TARGET zstd::libzstd)
elseif(TARGET zstd::libzstd_shared)
  set(COSMICSCANNER_ZSTD_TARGET zstd::libzstd_shared)
elseif(TARGET zstd::libzstd_static)
  set(COSMICSCANNER_ZSTD_TARGET zstd::libzstd_static)
else()
  message(FATAL_ERROR "A supported zstd CMake target was not found")
endif()
option(COSMICSCANNER_ENABLE_BLAS "Enable BLAS acceleration" OFF)

file(GLOB GWAS_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

if(COSMICSCANNER_ENABLE_BLAS)
  find_package(BLAS)
  if(BLAS_FOUND)
    list(APPEND COSMICSCANNER_EXTRA_LIBS ${BLAS_LIBRARIES})
    list(APPEND COSMICSCANNER_EXTRA_DEFS EIGEN_USE_BLAS)
  endif()
endif()

add_executable(cosmicscanner ${GWAS_SOURCES})
target_include_directories(cosmicscanner PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(cosmicscanner PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX ZLIB::ZLIB ${COSMICSCANNER_ZSTD_TARGET} ${COSMICSCANNER_EXTRA_LIBS})
target_compile_definitions(cosmicscanner PRIVATE ${COSMICSCANNER_EXTRA_DEFS})
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
    std::string scan_tier;
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
    bool multivariate_method_explicit = false;
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
    bool single_step_model = false;
    bool make_agh = false;
    bool write_agh = false;
    std::string ped_file;
    std::string genotyped_list;
    double h_tau = 1.0;
    double h_omega = 1.0;
    double blend_g = 0.05;
    bool write_ebv = false;

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
#include "multivariate_model.h"
#include "repeatability_model.h"
#include "single_step_model.h"
#include "univariate_model.h"

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
    std::cout << "CosmicScanner focused source build\n"
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
              << "  --gwas --scan-tier approx|semi|exact [--no-loco]\n"
              << "  --lm <1|2|3|4>               Linear model Wald/LRT/Score/all\n"
              << "  --lmm-gemma-mode <1|2|3|4>   LMM Wald/LRT/Score/all\n\n"
              << "Multivariate GWAS:\n"
              << "  --multivariate --mpheno-list 1,2[,3...] --scan-tier approx|semi|exact\n\n"
              << "Repeatability GWAS:\n"
              << "  --repeatability-model --tpheno <time_col> --mpheno <trait_col>\n"
              << "  --repeatability-mode exact --scan-tier approx|semi|exact\n\n"
              << "A/G/H matrices:\n"
              << "  --make-agh --ped <file> [--genotyped-list <file>]\n\n"
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
            else if (a == "--scan-tier") {
                cfg.scan_tier = need("--scan-tier");
                if (cfg.scan_tier != "approx" && cfg.scan_tier != "semi" && cfg.scan_tier != "exact")
                    throw std::runtime_error("--scan-tier must be approx, semi, or exact");
            }
            else if (a == "--no-loco") cfg.no_loco = true;
            else if (a == "--fast-lmm") cfg.fast_lmm = true;
            else if (a == "--lm" || a == "-lm" || a == "--lm-mode") { cfg.linear_model = true; cfg.linear_model_mode = std::stoi(need(a.c_str())); }
            else if (a == "--lmm-gemma-mode" || a == "-lmm") cfg.lmm_gemma_mode = std::stoi(need(a.c_str()));
            else if (a == "--lmm-gemma-refit") cfg.lmm_gemma_refit = true;
            else if (a == "--multivariate") cfg.multivariate = true;
            else if (a == "--mv-method") { cfg.multivariate_method = need("--mv-method"); cfg.multivariate_method_explicit = true; }
            else if (a == "--mv-output-vbeta") cfg.multivariate_output_vbeta = true;
            else if (a == "--emmvgwas") { cfg.multivariate = true; cfg.multivariate_method = "emmv"; cfg.multivariate_compatible = "emmvgwas"; cfg.no_loco = true; }
            else if (a == "--repeatability-model") cfg.repeatability_model = true;
            else if (a == "--repeatability-mode") cfg.repeatability_mode = need("--repeatability-mode");
            else if (a == "--repeatability-grm-cache") cfg.repeatability_grm_cache_prefix = need("--repeatability-grm-cache");
            else if (a == "--repeatability-weight") cfg.repeatability_weight_file = need("--repeatability-weight");
            else if (a == "--repeatability-refit-limit") cfg.repeatability_refit_limit = std::stoi(need("--repeatability-refit-limit"));
            else if (a == "--make-agh") { cfg.single_step_model = true; cfg.make_agh = true; cfg.write_agh = true; cfg.run_gwas = true; }
            else if (a == "--ped") cfg.ped_file = need("--ped");
            else if (a == "--genotyped-list") cfg.genotyped_list = need("--genotyped-list");
            else if (a == "--h2-tau") cfg.h_tau = std::stod(need("--h2-tau"));
            else if (a == "--h2-omega") cfg.h_omega = std::stod(need("--h2-omega"));
            else if (a == "--blend-g") cfg.blend_g = std::stod(need("--blend-g"));
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

        if (!cfg.run_gwas) throw std::runtime_error("This focused build requires --gwas or --make-agh.");
        if (cfg.repeatability_mode != "auto" && cfg.repeatability_mode != "exact" && cfg.repeatability_mode != "fast") {
            throw std::runtime_error("Focused source build supports --repeatability-mode auto, exact, or fast only.");
        }
        if (cfg.num_threads > 0) omp_set_num_threads(cfg.num_threads);

        #ifndef COSMICSCANNER_VERSION
        #define COSMICSCANNER_VERSION "unknown"
        #endif
        std::cout << "CosmicScanner focused source v" << COSMICSCANNER_VERSION << "\n";

        PlinkReader plink_reader;
        BgenReader bgen_reader;
        PgenReader pgen_reader;
        auto meta = DataLoader::loadGenotype(cfg.bfile_prefix, cfg.bgen_file, cfg.sample_file,
                                             cfg.pfile_prefix, cfg.pgen_file, cfg.pvar_file, cfg.psam_file,
                                             plink_reader, bgen_reader, pgen_reader);

        auto start = std::chrono::high_resolution_clock::now();
        if (cfg.single_step_model) {
            auto opts = SingleStepModel::Options::fromConfig(cfg);
            SingleStepModel model(opts);
            model.prepare(meta.fids, meta.iids, meta.n_samples, plink_reader, bgen_reader, pgen_reader, meta.is_bgen, meta.is_pgen);
            model.run();
        } else if (cfg.repeatability_model) {
            auto opts = RepeatabilityModel::Options::fromConfig(cfg);
            RepeatabilityModel model(opts);
            model.prepare(meta.fids, meta.iids, meta.n_samples, plink_reader, bgen_reader, pgen_reader, meta.is_bgen, meta.is_pgen);
            model.run();
        } else if (cfg.multivariate) {
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
import re
root = Path(__import__('sys').argv[1])
for path in root.glob("src/*"):
    if path.is_file() and path.suffix in {".cpp", ".h"}:
        path.write_text(path.read_text().replace('"models/', '"'))

uv = root / "src/univariate_model.cpp"
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

lmm_h = root / "src/lmm_engine.h"
text = lmm_h.read_text()
text = text.replace('#include "spa_test.h"  // Saddle Point Approximation\n', '')
text = text.replace('''    // Load Sparse GRM
    void loadSparseGRM(const std::string& prefix);

    Eigen::SparseMatrix<double> sparse_grm;
    bool has_sparse_grm = false;

''', '')
lmm_h.write_text(text)

lmm_cpp = root / "src/lmm_engine.cpp"
text = lmm_cpp.read_text()
text = text.replace('#include "sparse_grm_reader.h"\n', '')
text = text.replace('''void LmmEngine::loadSparseGRM(const std::string& prefix) {
    cosmic::SparseGrmReader::read(prefix, this->sparse_grm);
    this->has_sparse_grm = true;
}

''', '')
text = re.sub(
    r'    SPATest::Options spa_opts;\n.*?    }\n\n',
    '', text, flags=re.DOTALL
)
text = re.sub(
    r'\n\s*if \(options\.use_spa && pval < spa_opts\.pval_threshold\) \{\n.*?\n\s*\}',
    '', text, flags=re.DOTALL
)
lmm_cpp.write_text(text)

mv = root / "src/multivariate_model.cpp"
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
new = '''    if (options_.use_gemma_gxe || options_.use_hybrid) {
        throw std::runtime_error("Focused source build does not include multivariate GxE or hybrid routes.");
    }

    if (options_.use_gemma_refit) {
        runGemmaRefit();
        return;
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

cmake = root / "CMakeLists.txt"
text = cmake.read_text()
sources = "set(GWAS_SOURCES\n" + "".join(
    f"  ${{CMAKE_CURRENT_SOURCE_DIR}}/src/{path.name}\n"
    for path in sorted((root / "src").glob("*.cpp"))
) + ")"
text = text.replace('file(GLOB GWAS_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")', sources)
cmake.write_text(text)
PY

mkdir -p "${PREFIX}/scripts" "${PREFIX}/docs" "${PREFIX}/examples" "${PREFIX}/test"
cat > "${PREFIX}/.gitignore" <<'EOF'
/build/
/test/results/
*.log
EOF
cat > "${PREFIX}/scripts/build.sh" <<'EOF'
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
"${BUILD_DIR}/cosmicscanner" --help >/dev/null

echo "Build passed: ${BUILD_DIR}/cosmicscanner"
EOF
chmod +x "${PREFIX}/scripts/build.sh"

cat > "${PREFIX}/scripts/validate_release.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${COSMICSCANNER_BUILD_DIR:-${ROOT_DIR}/build}"
OUT_DIR="${COSMICSCANNER_VALIDATION_OUT:-${ROOT_DIR}/test/results}"

"${ROOT_DIR}/scripts/build.sh"
BIN="${BUILD_DIR}/cosmicscanner"
rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

run_case() {
  local name="$1"
  shift
  echo "[validate] START ${name}"
  "$@"
  echo "[validate] PASS  ${name}"
}

run_case "single-trait LMM" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/test_data.pheno" \
  --mpheno 1 --gwas --scan-tier semi --no-loco \
  --chr 1 --from-bp 10000 --to-bp 20000 \
  --out "${OUT_DIR}/single_trait"
test -s "${OUT_DIR}/single_trait.assoc.txt"
test "$(wc -l < "${OUT_DIR}/single_trait.assoc.txt")" -eq 12

run_case "multivariate LMM" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_multivariate.pheno" \
  --mpheno-list 1,2 \
  --multivariate --gwas --scan-tier semi --no-loco \
  --out "${OUT_DIR}/multivariate"
test -s "${OUT_DIR}/multivariate.assoc.mv.txt"

run_case "repeatability semi" \
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

run_case "repeatability fast comparison" \
  "${BIN}" \
  --bfile "${ROOT_DIR}/examples/test_data" \
  --pheno "${ROOT_DIR}/examples/demo_repeatability.pheno" \
  --qcovar "${ROOT_DIR}/examples/demo_repeatability.qcovar" \
  --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode fast --no-loco --gwas \
  --out "${OUT_DIR}/repeatability_fast"
test -s "${OUT_DIR}/repeatability_fast.assoc.txt"

run_case "A/G/H matrix construction" \
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

echo "CosmicScanner release validation passed."
EOF
chmod +x "${PREFIX}/scripts/validate_release.sh"

cat > "${PREFIX}/README.md" <<'EOF'
# CosmicScanner Source Package

CosmicScanner v1.1.0 is a focused, reproducible source release for the first
manuscript. It contains four supported functions:

- Single-trait GWAS
- Multivariate GWAS
- Repeatability GWAS
- A-inverse, G, and H-inverse construction

All retained C++ files are placed directly under a single flat `src/` directory.
The command-line interface exposes only the four supported functions. Shared
numerical implementation files required by those workflows remain internal and
do not define additional supported commands.

## Quick Start

```bash
bash scripts/build.sh
bash scripts/validate_release.sh
```

The release defaults to Eigen's native solvers. BLAS/LAPACK are opt-in because
results depend on the selected vendor runtime:

```bash
COSMICSCANNER_ENABLE_BLAS=ON bash scripts/build.sh
```

## Documentation

- `docs/USER_MANUAL_EN.md`: comprehensive English user manual.
- `docs/USER_MANUAL_CN.md`: comprehensive Chinese user manual.
- `docs/TROUBLESHOOTING.md`: solutions for common build, runtime, and output issues.

## Dependencies

- Linux x86_64
- CMake 3.15 or newer
- C++17 compiler
- Eigen3
- OpenMP
- zlib and zstd

If CMake cannot find Eigen3, set `CMAKE_PREFIX_PATH` to the prefix containing
`share/eigen3/cmake/Eigen3Config.cmake`.

## Citation

Citation metadata are provided in `CITATION.cff`.
EOF

cat > "${PREFIX}/CITATION.cff" <<'EOF'
cff-version: 1.2.0
message: "If you use CosmicScanner, please cite this software and the associated manuscript."
title: "CosmicScanner: Gaussian mixed-model association analysis and relationship-matrix construction"
type: software
version: "1.1.0"
date-released: "2026-07-14"
authors:
  - family-names: "Huang"
    given-names: "Jie"
    email: "huangjie2025@zju.edu.cn"
license: MIT
repository-code: "https://github.com/huangjie2025/CosmicScanner"
keywords:
  - GWAS
  - linear mixed model
  - multivariate analysis
  - repeated measurements
  - Woodbury identity
EOF

cat > "${PREFIX}/docs/USER_MANUAL_EN.md" <<'EOF'
# CosmicScanner v1.1.0 User Manual

## Overview

CosmicScanner v1.1.0 is a focused source release supporting four core workflows
for the first manuscript:

- **Single-trait GWAS**: Linear mixed model association analysis for a single
  continuous trait
- **Multivariate GWAS**: Joint association testing of multiple correlated
  continuous traits
- **Repeatability GWAS**: Mixed-model association analysis for repeated-measure
  phenotypes, with exact record-level analysis
- **A/G/H Matrix Construction**: Pedigree-based A⁻¹, genomic G, and H⁻¹ matrix
  building

All workflows share a unified infrastructure for genotype reading, sample
matching, covariate design, variant filtering, and GRM construction.

## Installation

### System Requirements

- Linux x86_64
- CMake 3.15 or newer
- C++17 compiler (GCC 8+, Clang 10+)
- Eigen 3
- OpenMP
- zlib and zstd

### Building from Source

```bash
tar -xzf CosmicScanner-v1.1.0-source.tar.gz
cd CosmicScanner-v1.1.0-source

# Default build (Eigen native solvers)
bash scripts/build.sh

# Optional BLAS/LAPACK acceleration
COSMICSCANNER_ENABLE_BLAS=ON bash scripts/build.sh

# Custom build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target cosmicscanner -j$(nproc)
```

If CMake cannot find Eigen3, set CMAKE_PREFIX_PATH to the directory
containing share/eigen3/cmake/Eigen3Config.cmake.

### Using the Linux Binary Package

```bash
tar -xzf CosmicScanner-v1.1.0-linux-x86_64.tar.gz
cd CosmicScanner-v1.1.0-linux-x86_64
export PATH=$PWD/bin:$PATH
cosmicscanner --help
```

### Build Verification

```bash
bash scripts/validate_release.sh
```

## Input File Formats

### Genotype Data

**PLINK 1 binary (recommended):** `--bfile prefix` for .bed/.bim/.fam.
**PLINK 2 PGEN:** `--pfile prefix` for .pgen/.pvar/.psam.
**BGEN:** `--bgen file.bgen --sample file.sample`.

### Phenotype Files

Single-trait/multivariate: `FID IID TRAIT1 [TRAIT2 ...]`.
Repeatability: `FID IID TIME PHENOTYPE`.
Missing values: NA or -9.

### Covariates

`--qcovar`: quantitative. `--covar`: categorical (dummy-coded).

## CLI Parameters

### Input/Output

| Parameter | Description | Required |
|---|---|---|
| --bfile <prefix> | PLINK 1 binary prefix | One of three |
| --pfile <prefix> | PLINK 2 PGEN prefix | One of three |
| --bgen <file> --sample <f> | BGEN input | Paired |
| --pheno <file> | Phenotype file | Yes |
| --qcovar <file> | Quantitative covariates | No |
| --covar <file> | Categorical covariates | No |
| --out <prefix> | Output prefix | Yes |

### Single-trait GWAS

| Parameter | Description | Default |
|---|---|---|
| --gwas | Enable GWAS scan | Required |
| --scan-tier approx\|semi\|exact | Scan tier | Legacy default |
| --no-loco | Disable LOCO | LOCO enabled |
| --lm <1\|2\|3\|4> | Linear model mode | Off |
| --lmm-gemma-mode <1\|2\|3\|4> | LMM Wald/LRT/Score | Wald |

### Multivariate GWAS

| Parameter | Description | Required |
|---|---|---|
| --multivariate | Enable multivariate | Yes |
| --mpheno-list c1,c2,... | Trait columns | Yes |
| --scan-tier approx\|semi\|exact | Scan tier | Recommended |
| --mv-method <method> | Method (emmv, etc.) | emmv |
| --mv-output-vbeta | Beta covariance | No |

### Repeatability GWAS

| Parameter | Description | Required |
|---|---|---|
| --repeatability-model | Enable model | Yes |
| --tpheno <n> | Time column | Yes |
| --mpheno <n> | Trait column | Yes |
| --repeatability-mode auto\|exact\|fast | Mode | auto |
| --scan-tier approx\|semi\|exact | Scan tier | Recommended |

### A/G/H Matrices

| Parameter | Description | Required |
|---|---|---|
| --make-agh | Build matrices | One of --gwas |
| --ped <file> | Pedigree file | Yes |
| --genotyped-list <file> | Genotyped IDs | No |
| --blend-g <float> | G blend weight | 0.05 |

### Filters

`--maf <x>`, `--hwe <x>`, `--extract <file>`, `--exclude <file>`,
`--chr <chr>`, `--from-bp <pos>`, `--to-bp <pos>`.

### Runtime

`--threads <n>`, `--lmm-null-solver <solver>`, `--reml-maxit <n>`,
`--reml-tol <float>`, `--help`.

## Examples

### Single-trait

```bash
cosmicscanner --bfile data --pheno trait.txt --mpheno 1 \\
  --gwas --scan-tier semi --out results/single
cosmicscanner --bfile data --pheno trait.txt --mpheno 1 \\
  --gwas --scan-tier semi --no-loco --out results/global
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier approx --out app
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier semi --out semi
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier exact --out exact
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier semi \
  --maf 0.01 --hwe 1e-6 --chr 1 --extract keep.txt --threads 8 --out filt
```

### Multivariate

```bash
cosmicscanner --bfile data --pheno traits.txt \\
  --mpheno-list 1,2 --multivariate --gwas --scan-tier semi --no-loco \\
  --out results/mv
cosmicscanner --bfile data --pheno traits.txt \\
  --mpheno-list 1,2,3 --multivariate --gwas --scan-tier semi --no-loco \\
  --out results/mv3
cosmicscanner --bfile data --pheno traits.txt \\
  --mpheno-list 1,2 --multivariate --gwas --scan-tier exact --no-loco \\
  --out results/mv_exact
```

### Repeatability

```bash
cosmicscanner --bfile data --pheno repeated.txt \\
  --tpheno 1 --mpheno 2 \\
  --repeatability-model --repeatability-mode exact \\
  --gwas --scan-tier semi --out results/rep
cosmicscanner --bfile data --pheno repeated.txt \\
  --tpheno 1 --mpheno 2 \\
  --repeatability-model --repeatability-mode fast \\
  --gwas --scan-tier approx --out results/rep_fast
```

### A/G/H

```bash
cosmicscanner --bfile data --make-agh \\
  --ped pedigree.txt --out results/agh
cosmicscanner --bfile data --make-agh \\
  --ped pedigree.txt \\
  --genotyped-list genotyped.txt \\
  --blend-g 0.05 --out results/agh
```

## Scan Tiers

--scan-tier approx|semi|exact:

| Tier | Handling | Use case |
|---|---|---|
| approx | Fixed null parameters | High-throughput screening |
| semi | Fixed ratios, profile scale | Recommended main scan |
| exact | Per-SNP full refit | Small data / fine-mapping |

Limitations: multivariate/repeatability exact require --no-loco.
Repeatability exact requires PLINK input, no record weights.

## Output Formats

### assoc.txt (single-trait, repeatability)

Columns: CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P

### assoc.mv.txt (multivariate)

Columns: CHR SNP BP A1 A2 FREQ WALD_STAT P BETA_1 SE_1 BETA_2 SE_2 ...

### repeatability.null.txt

Key-value: n_records, n_subjects, Vg, Vpe, Ve, Vp, h2, r2, timing, scan mode.

### scan.meta.txt

analysis_model, requested_scan_tier, resolved_scan_tier, snp_vc_strategy.

### Matrix .bin/.id

Lower-triangle row-by-row, 32-bit float, native endian. .id = one ID per line.

### agh.meta.txt

n_pedigree, n_genotyped, grm_algorithm, h_tau, h_omega, blend_g,
id_mapping_verified, hinv_symmetry_max_abs.

## GEMMA Compatibility

--lm <1|2|3|4>, --lmm-gemma-mode <1|2|3|4>, --lmm-gemma-refit,
--lmm-null-solver brent (same as GEMMA default).

## Version

v1.1.0. See CITATION.cff for citation metadata.
EOF

cat > "${PREFIX}/docs/USER_MANUAL_CN.md" <<'EOF'
# CosmicScanner v1.1.0 用户手册

## 概述

CosmicScanner v1.1.0 提供四项核心功能：

- **单性状 GWAS**：针对单个连续性状的线性混合模型关联分析
- **多性状 GWAS**：多个相关连续性状的联合关联检验
- **重复力 GWAS**：重复记录表型的混合模型关联分析
- **A/G/H 矩阵构建**：基于系谱的 A⁻¹、G 和 H⁻¹ 矩阵构建

## 安装

### 系统要求
Linux x86_64, CMake 3.15+, C++17 编译器, Eigen3, OpenMP, zlib, zstd.

### 源码构建
```bash
tar -xzf CosmicScanner-v1.1.0-source.tar.gz
cd CosmicScanner-v1.1.0-source
bash scripts/build.sh
```

### Linux 二进制包
```bash
tar -xzf CosmicScanner-v1.1.0-linux-x86_64.tar.gz
cd CosmicScanner-v1.1.0-linux-x86_64
export PATH=$PWD/bin:$PATH
cosmicscanner --help
```

## 输入格式

### 基因型
PLINK 1 BED: `--bfile <前缀>`。PLINK 2 PGEN: `--pfile <前缀>`。BGEN: `--bgen <文件> --sample <文件>`。

### 表型
单性状/多性状: `FID IID 性状1 [性状2 ...]`。重复力: `FID IID 时间 表型`。缺失值: NA 或 -9。

### 协变量
`--qcovar`: 定量协变量。`--covar`: 分类协变量（虚拟编码）。

## 主要参数

| 参数 | 说明 | 必需 |
|---|---|---|
| --bfile / --pfile / --bgen | 基因型输入 | 三选一 |
| --pheno <文件> | 表型文件 | 是 |
| --out <前缀> | 输出前缀 | 是 |
| --gwas | 启用 GWAS | 或 --make-agh |
| --scan-tier approx\|semi\|exact | 扫描档位 | 推荐 |
| --no-loco | 禁用 LOCO | LOCO 默认 |
| --multivariate | 多性状 GWAS | 配合 --mpheno-list |
| --repeatability-model | 重复力 GWAS | 配合 --tpheno |
| --make-agh | A/G/H 构建 | 配合 --ped |
| --maf / --hwe | 频率/HWE 过滤 | 否 |
| --extract / --exclude | SNP 选择 | 否 |
| --chr / --from-bp / --to-bp | 区间过滤 | 否 |
| --threads <n> | 线程数 | 系统默认 |

## 使用示例

### 单性状
```bash
cosmicscanner --bfile data --pheno trait.txt --mpheno 1 --gwas --scan-tier semi --out out
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier approx --out out.approx
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier semi --out out.semi
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier exact --out out.exact
```

### 多性状
```bash
cosmicscanner --bfile data --pheno traits.txt --mpheno-list 1,2 --multivariate --gwas --scan-tier semi --no-loco --out mv
cosmicscanner --bfile data --pheno traits.txt --mpheno-list 1,2 --multivariate --gwas --scan-tier exact --no-loco --out mv.exact
```

### 重复力
```bash
cosmicscanner --bfile data --pheno repeated.txt --tpheno 1 --mpheno 2 --repeatability-model --repeatability-mode exact --gwas --scan-tier semi --out rep
cosmicscanner --bfile data --pheno repeated.txt --tpheno 1 --mpheno 2 --repeatability-model --repeatability-mode fast --gwas --scan-tier approx --out rep.fast
```

### A/G/H
```bash
cosmicscanner --bfile data --make-agh --ped ped.txt --genotyped-list geno.txt --out agh
```

## 输出文件

| 文件 | 内容 |
|---|---|
| .assoc.txt | CHR, SNP, BP, A1, A2, FREQ, N, BETA, SE, P, LOG10P |
| .assoc.mv.txt | 多性状: WALD_STAT, P, BETA_t, SE_t |
| .repeatability.null.txt | Vg, Vpe, Ve, h2, r2, 时间统计 |
| .scan.meta.txt | analysis_model, scan_tier, vc_strategy |
| .Ainv.bin/.id | A⁻¹ 矩阵下三角 + ID |
| .G.bin/.id | G 矩阵下三角 + ID |
| .Hinv.bin/.id | H⁻¹ 矩阵下三角 + ID |
| .agh.meta.txt | 构建参数和验证信息 |

## 扫描档位

| 档位 | 策略 | 推荐场景 |
|---|---|---|
| approx | 固定 null 方差 | 高通量筛选 |
| semi | 固定比例, profile 残差 | 全基因组主扫描 |
| exact | 每 SNP 重估方差 | 小数据/候选位点 |

多性状/重复力 exact 要求 --no-loco；重复力 exact 额外要求 PLINK 输入和无记录权重。

## 版本

v1.1.0。引用信息见 CITATION.cff。
EOF

cat > "${PREFIX}/docs/TROUBLESHOOTING.md" <<'EOF'
# Troubleshooting Guide

## Build Issues

### CMake cannot find Eigen3
Set CMAKE_PREFIX_PATH to the Eigen3 prefix.

### Compiler does not support C++17
Use GCC 8+ or Clang 10+. Set CXX=g++-10 or equivalent.

### OpenMP not found
Install libomp-dev (Ubuntu) or libomp-devel (CentOS).

## Runtime Errors

### "This focused build requires --gwas or --make-agh"
Add --gwas or --make-agh to the command line.

### "Unsupported option in focused source build"
The parameter is not in v1.1.0 scope. Run cosmicscanner --help for supported
options. RRM, threshold, survival, GxE, SKAT, BSLMM, biobank mode are excluded.

### "Missing value for --<option>"
Provide the required argument.

### Genotype file errors
PLINK BED must be SNP-major (0x01). BGEN requires --sample. PGEN requires
.pgen/.pvar/.psam.

### Phenotype file errors
Minimum 3 columns (FID, IID, trait). Check --mpheno/--tpheno column numbers.

### Repeatability exact tier limitations
Requires --no-loco, PLINK input, and no --repeatability-weight.

### Multivariate exact requires --no-loco
Add --no-loco or use --scan-tier semi.

### GRM fails: "No valid SNPs found"
Loosen MAF/HWE thresholds. Check extract/exclude files.

### VCE does not converge
Increase --reml-maxit or loosen --reml-tol.

## Output Issues

### Empty output
All SNPs filtered out; check MAF/HWE/extract/exclude/chr settings.

### NaN p-values
Small sample size, monomorphic SNP, or collinearity with covariates.

## Performance

### Slow runtime
Increase --threads, use --scan-tier approx, reduce SNPs with --extract/--chr.

### Out of memory
Reduce samples/SNPs, use --chr to limit, avoid --scan-tier exact on large data.

## Platform

### Linux binary GLIBC error
Build from source on the target system.

### Windows/macOS
Not supported in v1.1.0. Requires conditional compilation for sys/resource.h.

## Getting Help

Open an issue at https://github.com/huangjie2025/CosmicScanner/issues with the
full command line, error output, and system information.
EOF



python3 - "${PREFIX}" <<'PY'
import math
import random
import sys
from pathlib import Path

root = Path(sys.argv[1])
examples = root / "examples"
examples.mkdir(parents=True, exist_ok=True)
rng = random.Random(20260714)
n, m = 120, 240
fids = ["SYN"] * n
iids = [f"SYN{i + 1:03d}" for i in range(n)]

genotypes = []
for j in range(m):
    p = 0.08 + 0.37 * ((j % 29) / 28.0)
    row = []
    for _ in range(n):
        g = int(rng.random() < p) + int(rng.random() < p)
        row.append(None if rng.random() < 0.005 else g)
    genotypes.append(row)

with (examples / "test_data.bed").open("wb") as handle:
    handle.write(bytes((0x6C, 0x1B, 0x01)))
    bytes_per_snp = (n + 3) // 4
    for row in genotypes:
        packed = bytearray(bytes_per_snp)
        for i, g in enumerate(row):
            code = 1 if g is None else {2: 0, 1: 2, 0: 3}[g]
            packed[i // 4] |= code << (2 * (i % 4))
        handle.write(packed)

alleles = [("A", "G"), ("C", "T"), ("G", "C"), ("T", "A")]
(examples / "test_data.bim").write_text("".join(
    f"{1 + j // 120}\trs{j + 1}\t0\t{1000 * (1 + j % 120)}\t{alleles[j % 4][0]}\t{alleles[j % 4][1]}\n"
    for j in range(m)
))
(examples / "test_data.fam").write_text("".join(
    f"{fids[i]}\t{iids[i]}\t0\t0\t{1 + i % 2}\t-9\n" for i in range(n)
))

def dosage(j, i):
    value = genotypes[j][i]
    if value is not None:
        return float(value)
    observed = [g for g in genotypes[j] if g is not None]
    return sum(observed) / len(observed)

y1, y2 = [], []
for i in range(n):
    polygenic = 0.35 * (dosage(4, i) - 1.0) - 0.25 * (dosage(37, i) - 1.0)
    shared = rng.gauss(0.0, 0.7)
    y1.append(1.2 + polygenic + shared + rng.gauss(0.0, 0.55))
    y2.append(-0.4 + 0.30 * polygenic + 0.45 * shared + rng.gauss(0.0, 0.75))

(examples / "test_data.pheno").write_text("".join(
    f"{fids[i]}\t{iids[i]}\t{y1[i]:.8f}\n" for i in range(n)
))
(examples / "demo_multivariate.pheno").write_text("".join(
    f"{fids[i]}\t{iids[i]}\t{y1[i]:.8f}\t{y2[i]:.8f}\n" for i in range(n)
))

rep_rows, qcov_rows = [], []
for i in range(60):
    permanent = rng.gauss(0.0, 0.35)
    for time in range(1, 4):
        qcov = (time - 2) + rng.gauss(0.0, 0.15)
        value = y1[i] + permanent + 0.18 * qcov + rng.gauss(0.0, 0.30)
        rep_rows.append(f"{fids[i]}\t{iids[i]}\t{time}\t{value:.8f}\n")
        qcov_rows.append(f"{fids[i]}\t{iids[i]}\t{time}\t{qcov:.8f}\n")
(examples / "demo_repeatability.pheno").write_text("".join(rep_rows))
(examples / "demo_repeatability.qcovar").write_text("".join(qcov_rows))

ped_ids = [f"{fids[i]}:{iids[i]}" for i in range(60)]
ped_rows = []
for i, animal in enumerate(ped_ids):
    if i < 20:
        sire, dam = "0", "0"
    else:
        sire = ped_ids[(i - 20) % 10]
        dam = ped_ids[10 + (i - 20) % 10]
    ped_rows.append(f"{animal}\t{sire}\t{dam}\n")
(examples / "demo.ped").write_text("".join(ped_rows))
(examples / "demo.genotyped.txt").write_text("".join(f"{animal}\n" for animal in ped_ids))
(examples / "README.md").write_text("""# Synthetic Example Data

All files in this directory are deterministic synthetic fixtures generated by
`scripts/package_source.sh` with Python's random seed `20260714`. They contain
no human, animal, cohort, or real-world identifiers or measurements.

- `test_data.*`: PLINK 1 binary data with 120 synthetic samples and 240 SNPs.
- `demo_multivariate.pheno`: two correlated, non-collinear continuous traits.
- `demo_repeatability.*`: three noisy records for each of 60 subjects.
- `demo.ped` and `demo.genotyped.txt`: a synthetic ordered pedigree and subset.

These compact files are intended only for build and smoke validation; they are
not evidence for statistical power or biological conclusions.
""")
PY

cat > "${PREFIX}/test/run_source_smoke.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec "${ROOT_DIR}/scripts/validate_release.sh"
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
