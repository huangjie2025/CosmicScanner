#pragma once

#include "config.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cosmic {

struct CliInvocation {
    std::vector<std::string> original;
    std::vector<std::string> args;
    std::string task;
    std::string help_topic;
    bool help = false;
    bool pedigree_alias = false;
    bool ped_genotype_alias = false;
};

inline CliInvocation prepareCliInvocation(int argc, char* argv[]) {
    CliInvocation result;
    for (int i = 0; i < argc; ++i) result.original.emplace_back(argv[i]);
    result.args.push_back(argv[0]);
    int start = 1;
    if (argc > 1 && argv[1][0] != '-') {
        result.task = argv[1];
        start = 2;
        if (result.task == "help") {
            result.help = true;
            if (argc > 3) throw std::runtime_error("Usage: cosmicscanner help [topic]");
            result.help_topic = argc == 3 ? argv[2] : "overview";
            return result;
        }
        if (result.task == "version") {
            if (argc != 2) throw std::runtime_error("Usage: cosmicscanner version");
            result.args.push_back("--version"); return result;
        }
        if (result.task == "gwas") result.args.push_back("--gwas");
        else if (result.task == "grm") result.args.push_back("--make-grm");
        else throw std::runtime_error("Unknown command '" + result.task + "'. Use gwas, grm, help, or version.");
    }
    for (int i = start; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help-all") { result.help = true; result.help_topic = "all"; }
        else if (arg == "--help" || arg == "-h") {
            result.help = true;
            result.help_topic = result.task.empty() ? "overview" : result.task;
            if (i + 1 < argc && argv[i + 1][0] != '-') result.help_topic = argv[++i];
        } else if (arg.rfind("--help=", 0) == 0) {
            result.help = true; result.help_topic = arg.substr(7);
        } else {
            if (arg == "--pedigree") { result.pedigree_alias = true; arg = "--ped"; }
            if (arg == "--ped-genotypes") { result.ped_genotype_alias = true; arg = "--ped"; }
            result.args.push_back(arg);
        }
    }
    if (argc < 2) { result.help = true; result.help_topic = "overview"; }
    return result;
}

inline void printCliHelp(const std::string& topic, const std::string& program) {
    std::cout << "CosmicScanner: genetic association analysis\n"
              << "Usage: " << program << " [gwas|grm] [options]\n\n";
    if (topic == "overview" || topic == "gwas") {
        std::cout <<
            "Common options (PLINK/GCTA-style long flags):\n"
            "  --bfile PREFIX         PLINK BED/BIM/FAM input\n"
            "  --pfile PREFIX         PLINK PGEN/PVAR/PSAM input\n"
            "  --bgen FILE --sample FILE  BGEN input (model-specific support)\n"
            "  --pheno FILE           FID IID followed by phenotype columns\n"
            "  --mpheno N             1-based column AFTER FID/IID; default 1\n"
            "  --mpheno-list 1,2      Phenotype columns for multivariate analysis\n"
            "  --covar FILE           Categorical covariates (FID IID ...)\n"
            "  --qcovar FILE          Quantitative covariates (FID IID ...)\n"
            "  --gwas                 Association task; optional after 'gwas' command\n"
            "  --model NAME           lmm, lm, multivariate, repeatability, ...\n"
            "  --scan-tier TIER       approx|semi|exact; see help models for limits\n"
            "  --no-loco              Global scan; otherwise LOCO where supported\n"
            "  --maf F --hwe P        Variant filters\n"
            "  --extract FILE --exclude FILE  Variant-ID inclusion/exclusion\n"
            "  --threads N            Number of threads\n"
            "  --out PREFIX           Output file prefix\n\n"
            "Example:\n"
            "  cosmicscanner --bfile data --pheno trait.txt --gwas --model lmm --out run\n"
            "  cosmicscanner gwas --model multivariate --bfile data --pheno traits.txt --mpheno-list 1,2 --scan-tier approx --reml-maxit 1000 --out mv\n\n";
    } else if (topic == "models" || topic == "lmm" || topic == "lm" || topic == "multivariate") {
        std::cout <<
            "--model NAME selects one model family; old model flags still work.\n"
            "  lmm                 Gaussian LMM (default; scoped reference)\n"
            "  lm                  Linear regression\n"
            "  multivariate        Multivariate LMM; --mpheno-list is required\n"
            "  repeatability       Repeated records; see help repeatability\n"
            "  ssgwas              Method 1 H-subset LMM (experimental)\n"
            "  rrm, animal, maternal, common-env, threshold, count, ordinal, gxe\n"
            "                      Experimental, model-specific contracts\n"
            "  lmm-eb / eb-lmm      Legacy EB modifier; experimental\n\n"
            "--gwas-model additive|additive_dominance|dominance selects SNP/random\n"
            "effect structure WITHIN the single-trait model, not the model family.\n"
            "approx fixes null variance; semi profiles residual scale/covariance;\n"
            "exact refits SNP alternative variance, not finite-sample exact P.\n"
            "Native multivariate semi is unsupported; compatibility semi requires\n"
            "--mv-compatible emmvgwas. Multivariate/repeatability exact needs no-LOCO.\n";
    } else if (topic == "input") {
        std::cout <<
            "Choose ONE genotype source: --bfile, --pfile, --bgen, --vcf, ...\n"
            "Reader availability does not imply support by every model.\n"
            "--pheno FILE: FID IID VALUE...; --mpheno 1 selects the third file column.\n"
            "--mpheno-list 1,2 selects third/fourth file columns.\n"
            "Long records: FID IID TIME VALUE; use --tpheno 1 --mpheno 2.\n"
            "These are phenotype-relative indices, NOT PLINK2 absolute column numbers.\n"
            "--pedigree FILE is ID SIRE DAM, for animal/ssGWAS/AGH models.\n"
            "--ped-genotypes FILE --map FILE is PED/MAP genotype input.\n"
            "Legacy --ped is preserved but context-dependent.\n";
    } else if (topic == "grm") {
        std::cout <<
            "Build a genetic relationship matrix:\n"
            "  cosmicscanner grm --bfile data --threads 4 --out data_grm\n"
            "Equivalent legacy task: --make-grm. --make-grm-alg keeps its existing\n"
            "algorithm numbering and normalization. See --help-all for details.\n";
    } else if (topic == "repeatability") {
        std::cout <<
            "  cosmicscanner gwas --model repeatability --bfile data --pheno long.txt --tpheno 1 --mpheno 2 --repeatability-mode exact --scan-tier approx --out rep\n"
            "--repeatability-mode auto|exact|fast|stcg selects record computation.\n"
            "--scan-tier approx|semi|exact separately selects per-SNP variance handling.\n"
            "Tier exact: BED, no-LOCO, no record weights. Fast/STCG are experimental.\n"
            "auto can change computation with record count; set the mode explicitly\n"
            "for reproducible comparisons. --repeatability-grm-cache PREFIX is optional.\n";
    } else if (topic == "rrm") {
        std::cout <<
            "Experimental RRM: --model rrm --rrm-order N; balanced wide and unbalanced\n"
            "long formats retain their existing definitions. For unbalanced records:\n"
            "  --rrm-unbalance --long-pheno FILE --tpheno 1 --mpheno 2\n"
            "--rrm-gwas-method record|coefficient selects different statistical tests.\n"
            "No automatic format or estimator conversion is performed.\n";
    } else if (topic == "ssgwas") {
        std::cout <<
            "--model ssgwas --pedigree FILE selects Method 1 H-subset fixed-SNP LMM\n"
            "by default (experimental); it does not provide EBV output.\n"
            "--single-step-method backsolve selects Method 2 full-pedigree ssGBLUP.\n"
            "Reference: BED + explicit positive --reml-vg/--reml-ve + LDLT.\n"
            "P-values default off; APY/PCG/estimated-VC routes have separate limits.\n"
            "Use --help-all and the scope contract before selecting extensions.\n";
    } else {
        throw std::runtime_error("Unknown help topic '" + topic + "'. Try models, input, grm, repeatability, rrm, ssgwas, or all.");
    }
    std::cout << "\nHelp: --help [topic], help [topic], --help-all (complete legacy/advanced reference)\n"
              << "Manual: docs/CLI_GUIDE_EN.md / docs/CLI_GUIDE_CN.md\n"
              << "Existing flat commands, aliases and statistical defaults are retained.\n";
}

inline void applyCliModel(GwasConfig& cfg, const std::string& model, const CliInvocation& invocation) {
    const bool legacy_modifier = model == "lmm" || model == "fixed-lmm" || model == "lmm-eb" || model == "eb-lmm";
    const int selected = int(cfg.linear_model) + int(cfg.multivariate) + int(cfg.repeatability_model) +
        int(cfg.rrm_model) + int(cfg.animal_model) + int(cfg.maternal_model) + int(cfg.single_step_model) +
        int(cfg.common_env_model) + int(cfg.threshold_model) + int(cfg.count_model) + int(cfg.ordinal_model) +
        int(cfg.survival_model) + int(cfg.gxe_model) + int(cfg.farmcpu) + int(cfg.blink) +
        int(cfg.firth_logistic) + int(!cfg.logistic_test.empty()) + int(cfg.logistic_glmm) +
        int(cfg.biobank_mode) + int(cfg.saige_style) + int(cfg.run_skat) + int(cfg.run_burden) +
        int(cfg.skat_o) + int(cfg.acat_v);
    if (!model.empty()) {
        if ((!legacy_modifier || !invocation.task.empty()) &&
            (cfg.make_grm || cfg.make_grm_sparse || cfg.make_agh || !cfg.convert_from.empty() ||
             !cfg.batch_models.empty() || cfg.pca || cfg.ld_prune || !cfg.ld_clump_file.empty()))
            throw std::runtime_error("--model is for association; do not combine it with another task.");
        if (legacy_modifier) {
            if (selected > 0 && invocation.task == "gwas")
                throw std::runtime_error("Conflicting model selectors: choose one --model family or legacy model flag.");
            // Historical --model lmm[-eb] was an EB modifier. Preserve that
            // behavior in existing flat scripts with another model flag.
        } else {
            bool same = (model == "lm" && cfg.linear_model) || (model == "multivariate" && cfg.multivariate) ||
                (model == "repeatability" && cfg.repeatability_model) || (model == "rrm" && cfg.rrm_model) ||
                (model == "animal" && cfg.animal_model) || (model == "maternal" && cfg.maternal_model) ||
                (model == "ssgwas" && cfg.single_step_model) || (model == "common-env" && cfg.common_env_model) ||
                (model == "threshold" && cfg.threshold_model) || (model == "count" && cfg.count_model) ||
                (model == "ordinal" && cfg.ordinal_model) || (model == "gxe" && cfg.gxe_model);
            if (selected > int(same)) throw std::runtime_error("Conflicting model selectors for --model " + model);
            if (model == "lm") cfg.linear_model = true;
            else if (model == "multivariate") cfg.multivariate = true;
            else if (model == "repeatability") cfg.repeatability_model = true;
            else if (model == "rrm") cfg.rrm_model = true;
            else if (model == "animal") cfg.animal_model = true;
            else if (model == "maternal") cfg.maternal_model = true;
            else if (model == "ssgwas") cfg.single_step_model = true;
            else if (model == "common-env") cfg.common_env_model = true;
            else if (model == "threshold") cfg.threshold_model = true;
            else if (model == "count") cfg.count_model = true;
            else if (model == "ordinal") cfg.ordinal_model = true;
            else if (model == "gxe") cfg.gxe_model = true;
            else throw std::runtime_error("Unknown model '" + model + "'. See --help models.");
            cfg.run_gwas = true;
        }
        if (model == "multivariate" && cfg.mpheno_list.empty())
            throw std::runtime_error("--model multivariate requires --mpheno-list (e.g. 1,2).");
    }
    if (invocation.task == "grm" && (selected > 0 || cfg.run_gwas))
        throw std::runtime_error("grm is a matrix-only task; remove association/model options.");
    if (invocation.task == "gwas" && (cfg.make_grm || cfg.make_grm_sparse || cfg.make_agh ||
        !cfg.convert_from.empty() || !cfg.batch_models.empty() || cfg.pca || cfg.ld_prune ||
        !cfg.ld_clump_file.empty() || cfg.qc_requested))
        throw std::runtime_error("gwas cannot be combined with a standalone conversion/QC/matrix/batch task.");
    if (invocation.pedigree_alias && invocation.ped_genotype_alias)
        throw std::runtime_error("Choose --pedigree or --ped-genotypes, not both.");
    if ((invocation.pedigree_alias || invocation.ped_genotype_alias) &&
        std::count(invocation.args.begin(), invocation.args.end(), "--ped") > 1)
        throw std::runtime_error("Specify the pedigree/PED input only once.");
    const bool uses_pedigree = cfg.single_step_model || cfg.animal_model || cfg.make_agh;
    if (invocation.pedigree_alias && !uses_pedigree)
        throw std::runtime_error("--pedigree requires an animal, ssGWAS, or AGH model.");
    if (invocation.ped_genotype_alias && uses_pedigree)
        throw std::runtime_error("--ped-genotypes is genotype input, not an animal/ssGWAS pedigree.");
}

} // namespace cosmic
