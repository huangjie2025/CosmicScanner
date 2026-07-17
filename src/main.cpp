#include <chrono>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif
#include <omp.h>

#include "config.h"
#include "data_loader.h"
#include "multivariate_model.h"
#include "repeatability_model.h"
#include "single_step_model.h"
#include "univariate_model.h"

#ifndef COSMICSCANNER_VERSION
#define COSMICSCANNER_VERSION "unknown"
#endif

namespace {

const char* kScanTierValues = "approx, semi, exact";
const char* kLmmGemmaModeValues = "1=Wald, 2=LRT, 3=Score, 4=all";
const char* kNullSolverValues = "brent, idul, idul_plus (or idul+), aidul";
const char* kMvMethodValues = "wald, emmv, gemma, gemma-full, gemma-exact-full, eb";
const char* kRepeatabilityModeValues = "auto, exact, fast, stcg";

bool is_int_mode(int mode) {
    return mode >= 1 && mode <= 4;
}

bool is_supported_mv_method(const std::string& method) {
    return method == "wald" || method == "emmv" || method == "gemma"
        || method == "gemma-full" || method == "gemma-exact"
        || method == "gemma-exact-full" || method == "eb"
        || method == "eb-lmm" || method == "lmm-eb";
}

bool is_unavailable_mv_method(const std::string& method) {
    return method == "hybrid" || method == "gemma-gxe" || method == "gxe"
        || method == "gemma-gxe-full" || method == "gxe-full";
}

double peak_memory_mb() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS counters{};
    if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters))) {
        return static_cast<double>(counters.PeakWorkingSetSize) / (1024.0 * 1024.0);
    }
    return 0.0;
#else
    struct rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return 0.0;
#ifdef __APPLE__
    return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
    return static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
#endif
}

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
              << "  --gwas --scan-tier <tier>     LMM scan tier. Values: " << kScanTierValues << "\n"
              << "                                Default: approx when omitted\n"
              << "                                approx=fixed-null; semi=fixed-ratio profile scale;\n"
              << "                                exact=per-SNP variance-component refit\n"
              << "  --lm <mode>                   Linear model test. Values: " << kLmmGemmaModeValues << "\n"
              << "  --lmm-gemma-mode <mode>       LMM primary/output P. Values: " << kLmmGemmaModeValues << "\n"
              << "                                Default: 0 legacy P column; mode 4 emits P_WALD/P_LRT/P_SCORE\n"
              << "  --lmm-null-solver <solver>    Null-model solver only. Values: " << kNullSolverValues << "\n"
              << "                                Default: brent. IDUL-family solvers fall back to brent if needed.\n\n"
              << "Multivariate GWAS:\n"
              << "  --multivariate --mpheno-list 1,2[,3...] --mv-method <method>\n"
              << "  --mv-method <method>          Values: " << kMvMethodValues << "\n"
              << "                                Default: emmv. gemma-exact is accepted as an internal alias.\n"
              << "  --scan-tier <tier>            Shortcut: approx->wald, semi->emmv, exact->gemma-exact\n"
              << "                                Multivariate exact currently requires --no-loco.\n"
              << "                                Hybrid and GxE are unavailable in this focused release.\n\n"
              << "Repeatability GWAS:\n"
              << "  --repeatability-model --tpheno <time_col> --mpheno <trait_col>\n"
              << "  --repeatability-mode <mode>   Values: " << kRepeatabilityModeValues << "; default: auto\n"
              << "  --repeatability-mode exact --scan-tier <tier>  Values: " << kScanTierValues << "\n"
              << "  --repeatability-mode fast|stcg                 Matrix-accelerated repeatability scans\n\n"
              << "A/G/H matrices:\n"
              << "  --make-agh --ped <file> [--genotyped-list <file>]\n\n"
              << "LOCO:\n"
              << "  Default: use LOCO GRMs when multiple chromosomes are present.\n"
              << "  --no-loco                    Use one global GRM for all scanned SNPs.\n\n"
              << "Common filters:\n"
              << "  --maf <x> --hwe <x> --extract <file> --exclude <file>\n"
               << "  --chr <chr> --from-bp <bp> --to-bp <bp> --threads <n>\n";
    std::cout << "  --version                    Print version and exit\n";
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
            else if (a == "--version") { std::cout << "CosmicScanner " << COSMICSCANNER_VERSION << "\n"; return 0; }
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
                    throw std::runtime_error(std::string("--scan-tier must be one of: ") + kScanTierValues);
            }
            else if (a == "--no-loco") cfg.no_loco = true;
            else if (a == "--fast-lmm") cfg.fast_lmm = true;
            else if (a == "--lm" || a == "-lm" || a == "--lm-mode") {
                cfg.linear_model = true;
                cfg.linear_model_mode = std::stoi(need(a.c_str()));
                if (!is_int_mode(cfg.linear_model_mode)) {
                    throw std::runtime_error(std::string("--lm must be one of: ") + kLmmGemmaModeValues);
                }
            }
            else if (a == "--lmm-gemma-mode" || a == "-lmm") {
                cfg.lmm_gemma_mode = std::stoi(need(a.c_str()));
                if (!is_int_mode(cfg.lmm_gemma_mode)) {
                    throw std::runtime_error(std::string("--lmm-gemma-mode must be one of: ") + kLmmGemmaModeValues);
                }
            }
            else if (a == "--lmm-gemma-refit") cfg.lmm_gemma_refit = true;
            else if (a == "--multivariate") cfg.multivariate = true;
            else if (a == "--mv-method") {
                cfg.multivariate_method = need("--mv-method");
                cfg.multivariate_method_explicit = true;
                if (is_unavailable_mv_method(cfg.multivariate_method)) {
                    throw std::runtime_error("--mv-method " + cfg.multivariate_method
                        + " is unavailable in this focused release. Supported values: " + kMvMethodValues);
                }
                if (!is_supported_mv_method(cfg.multivariate_method)) {
                    throw std::runtime_error("Unknown --mv-method: " + cfg.multivariate_method
                        + ". Supported values: " + kMvMethodValues);
                }
            }
            else if (a == "--mv-output-vbeta") cfg.multivariate_output_vbeta = true;
            else if (a == "--mv-null-only") { cfg.multivariate = true; cfg.multivariate_null_only = true; cfg.run_gwas = true; }
            else if (a == "--mv-k-file") cfg.multivariate_k_file = need("--mv-k-file");
            else if (a == "--mv-compare-vg") cfg.multivariate_compare_vg = need("--mv-compare-vg");
            else if (a == "--mv-compare-ve") cfg.multivariate_compare_ve = need("--mv-compare-ve");
            else if (a == "--emmvgwas") { cfg.multivariate = true; cfg.multivariate_method = "emmv"; cfg.multivariate_compatible = "emmvgwas"; cfg.no_loco = true; }
            else if (a == "--repeatability-model") cfg.repeatability_model = true;
            else if (a == "--repeatability-mode") {
                cfg.repeatability_mode = need("--repeatability-mode");
                if (cfg.repeatability_mode != "auto" && cfg.repeatability_mode != "exact"
                    && cfg.repeatability_mode != "fast" && cfg.repeatability_mode != "stcg") {
                    throw std::runtime_error(std::string("--repeatability-mode must be one of: ") + kRepeatabilityModeValues);
                }
            }
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
            else if (a == "--lmm-null-solver") {
                cfg.lmm_null_solver = need("--lmm-null-solver");
                if (cfg.lmm_null_solver != "brent" && cfg.lmm_null_solver != "idul"
                    && cfg.lmm_null_solver != "idul_plus" && cfg.lmm_null_solver != "idul+"
                    && cfg.lmm_null_solver != "aidul") {
                    throw std::runtime_error(std::string("--lmm-null-solver must be one of: ") + kNullSolverValues);
                }
            }
            else if (a == "--reml-maxit") cfg.reml_maxit = std::stoi(need("--reml-maxit"));
            else if (a == "--reml-tol") cfg.reml_tol = std::stod(need("--reml-tol"));
            else throw std::runtime_error("Unsupported option in focused source build: " + a);
        }

        if (!cfg.run_gwas) throw std::runtime_error("This focused build requires --gwas or --make-agh.");
        const int genotype_sources = (!cfg.bfile_prefix.empty() ? 1 : 0)
            + (!cfg.bgen_file.empty() ? 1 : 0)
            + ((!cfg.pfile_prefix.empty() || !cfg.pgen_file.empty()) ? 1 : 0);
        if (genotype_sources != 1) {
            throw std::runtime_error("Specify exactly one genotype input: --bfile, --bgen, or --pfile/--pgen");
        }
        if (cfg.repeatability_mode != "auto" && cfg.repeatability_mode != "exact"
            && cfg.repeatability_mode != "fast" && cfg.repeatability_mode != "stcg") {
            throw std::runtime_error("Focused source build supports --repeatability-mode auto, exact, fast, or stcg only.");
        }
        if (cfg.scan_from_bp < -1 || cfg.scan_to_bp < -1) {
            throw std::runtime_error("--from-bp and --to-bp must be non-negative");
        }
        if (cfg.scan_from_bp >= 0 && cfg.scan_to_bp >= 0 && cfg.scan_from_bp > cfg.scan_to_bp) {
            throw std::runtime_error("--from-bp cannot be greater than --to-bp");
        }
        if (cfg.num_threads > 0) omp_set_num_threads(cfg.num_threads);

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
        std::cout << "Total computation time: " << elapsed.count() << " seconds\n";
        std::cout << "Peak memory usage: " << peak_memory_mb() << " MB\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
