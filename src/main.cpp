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

#ifndef COSMICSCANNER_VERSION
#define COSMICSCANNER_VERSION "unknown"
#endif

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
              << "  --repeatability-mode exact --scan-tier approx|semi|exact\n"
              << "  --repeatability-mode fast|stcg\n\n"
              << "A/G/H matrices:\n"
              << "  --make-agh --ped <file> [--genotyped-list <file>]\n\n"
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
