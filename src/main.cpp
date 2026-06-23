#include <iostream>
#include <string>
#include <sstream>
#include "config.h"
#include "gwas_app.h"

using namespace std;
using namespace cosmic;

void print_usage(const char* prog_name) {
    cout << "Usage: " << prog_name << " [options]\n";
    cout << "Input Options:\n";
    cout << "  --bfile <prefix>    PLINK binary prefix\n";
    cout << "  --bgen <file>       BGEN file (Layout 2, zlib)\n";
    cout << "  --sample <file>     Sample file for BGEN (Oxford format)\n";
    cout << "  --pfile <prefix>    PLINK 2.0 prefix (.pgen/.pvar/.psam)\n";
    cout << "  --pgen <file>       PLINK 2.0 genotype file\n";
    cout << "  --pvar <file>       PLINK 2.0 variant file\n";
    cout << "  --psam <file>       PLINK 2.0 sample file\n";
    cout << "Output Options:\n";
    cout << "  --out <prefix>      Output prefix (default: cosmic_gwas)\n";
    cout << "  --out-ssf           Output in EBI GWAS-SSF format\n";
    cout << "Phenotype/Covariate:\n";
    cout << "  --pheno <file>      Phenotype file (FID IID Pheno)\n";
    cout << "  --mpheno <n>        Phenotype column number (default: 1)\n";
    cout << "  --mpheno-list <list> Comma-separated list of phenotype columns (e.g. 1,2,3)\n";
    cout << "  --multivariate      Run Multivariate LMM (requires --mpheno-list)\n";
    cout << "  --repeatability-model Run Repeatability GWAS on long-format phenotypes\n";
    cout << "  --covar <file>      Categorical Covariate file (FID IID Cat1 ...)\n";
    cout << "  --qcovar <file>     Quantitative Covariate file (FID IID Val1 ...)\n";
    cout << "Analysis Options:\n";
    cout << "  --make-grm          Calculate GRM and exit (PLINK only)\n";
    cout << "  --make-grm-alg <n>  GRM Algorithm (1: VanRaden, 3: Yang, 5: TrueVar/Yang2010). Default: 1\n";
    cout << "  --gwas              Run GWAS scan\n";
    cout << "  --fast-lmm          Use approximate LMM (fixed delta from null model)\n";
    cout << "  --repeatability-mixed-precision Use Float32 kernels for repeatability exact scan\n";
    cout << "  --repeatability-grm-cache <prefix> Save/load repeatability subject GRM cache\n";
    cout << "  --repeatability-refit-validation Run per-SNP variance-component refit validation for repeatability GWAS (small datasets)\n";
    cout << "  --repeatability-refit-limit <n> Limit per-SNP refit validation to first n SNPs (0: all)\n";
    cout << "  --repeatability-mode <mode> Repeatability strategy: auto, exact, fast, stcg (default: auto)\n";
    cout << "  --spa               Use Saddle Point Approximation (for binary traits)\n";
    cout << "  --threshold-model   Binary trait PSR-GWAS (PQL + probit/logit + weighted LMM)\n";
    cout << "  --link <type>       Threshold link: probit or logit (default: probit)\n";
    cout << "  --count-model       Count/binomial PSR-GWAS (PQL + weighted LMM)\n";
    cout << "  --family <type>     Count family: poisson, negbin, or binomial (default: poisson)\n";
    cout << "  --offset <file>     One numeric log-offset/exposure value per sample\n";
    cout << "  --binomial-pheno <file> Binomial FID IID EVENT TRIAL file\n";
    cout << "  --ordinal-model     Ordinal PSR-GWAS MVP (experimental cumulative-binary)\n";
    cout << "  --ordinal-pheno <file> Ordinal FID IID ORDINAL phenotype file\n";
    cout << "  --psr-maxit <n>     Max PSR/PQL iterations (default: 50)\n";
    cout << "  --psr-tol <float>   PSR/PQL convergence tolerance (default: 1e-6)\n";
    cout << "  --psr-mu-eps <float> Mu clipping epsilon for PSR (default: 1e-6)\n";
    cout << "  --psr-min-weight <float> Minimum PSR working weight (default: 1e-8)\n";
    cout << "  --psr-max-weight <float> Maximum PSR working weight (default: 1e8)\n";
    cout << "  --skat              Run SKAT (Sequence Kernel Association Test) for rare variants\n";
    cout << "  --burden            Run Burden Test for rare variants\n";
    cout << "  --set-file <file>   Define variant sets (e.g., genes) for rare variant tests\n";
    cout << "  --biobank-mode      Run REGENIE-style Ridge Regression for ultra-large datasets\n";
    cout << "  --block-size <int>  Block size for REGENIE Step 1 (default: 1000)\n";
    cout << "  --snp-weight <file> SNP weight file for weighted GRM\n";
    cout << "  --pop-class <file>  Population classification file for mixed GRM\n";
    cout << "  --maf <float>       Minimum Allele Frequency filter (default: 0.01)\n";
    cout << "  --hwe <float>       Minimum HWE p-value filter (default: 0.0, disabled)\n";
    cout << "VCE Options:\n";
    cout << "  --reml-maxit <n>    Maximum VCE iterations (default: 50)\n";
    cout << "  --reml-tol <float>  VCE convergence tolerance (default: 1e-6)\n";
    cout << "  --reml-vg <float>   Genetic variance (skip VCE, use this value)\n";
    cout << "  --reml-vpe <float>  Permanent environment variance (skip VCE)\n";
    cout << "  --reml-ve <float>   Residual variance (skip VCE)\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    GwasConfig cfg;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--bfile" && i + 1 < argc) { cfg.bfile_prefix = argv[++i]; }
        else if (arg == "--bgen" && i + 1 < argc) { cfg.bgen_file = argv[++i]; }
        else if (arg == "--sample" && i + 1 < argc) { cfg.sample_file = argv[++i]; }
        else if (arg == "--pfile" && i + 1 < argc) { cfg.pfile_prefix = argv[++i]; }
        else if (arg == "--pgen" && i + 1 < argc) { cfg.pgen_file = argv[++i]; }
        else if (arg == "--pvar" && i + 1 < argc) { cfg.pvar_file = argv[++i]; }
        else if (arg == "--psam" && i + 1 < argc) { cfg.psam_file = argv[++i]; }
        else if (arg == "--out" && i + 1 < argc) { cfg.out_prefix = argv[++i]; }
        else if (arg == "--out-ssf") { cfg.out_ssf = true; }
        else if (arg == "--pheno" && i + 1 < argc) { cfg.pheno_file = argv[++i]; }
        else if (arg == "--long-pheno" && i + 1 < argc) { cfg.long_pheno_file = argv[++i]; cfg.rrm_model = true; }
        else if (arg == "--mpheno" && i + 1 < argc) { cfg.mpheno = stoi(argv[++i]); }
        else if (arg == "--tpheno" && i + 1 < argc) { cfg.tpheno = stoi(argv[++i]); }
        else if (arg == "--rrm-order" && i + 1 < argc) { cfg.rrm_order = stoi(argv[++i]); }
        else if (arg == "--mpheno-list" && i + 1 < argc) {
            string s = argv[++i];
            stringstream ss(s);
            string token;
            while(getline(ss, token, ',')) { cfg.mpheno_list.push_back(stoi(token)); }
        }
        else if (arg == "--multivariate") { cfg.multivariate = true; }
        else if (arg == "--repeatability-model") { cfg.repeatability_model = true; }
        else if (arg == "--threshold-model") { cfg.threshold_model = true; }
        else if (arg == "--link" && i + 1 < argc) { cfg.link_type = argv[++i]; }
        else if (arg == "--count-model") { cfg.count_model = true; }
        else if (arg == "--family" && i + 1 < argc) { cfg.family_type = argv[++i]; }
        else if (arg == "--offset" && i + 1 < argc) { cfg.offset_file = argv[++i]; }
        else if (arg == "--binomial-pheno" && i + 1 < argc) { cfg.binomial_file = argv[++i]; }
        else if (arg == "--ordinal-model") { cfg.ordinal_model = true; }
        else if (arg == "--ordinal-pheno" && i + 1 < argc) { cfg.ordinal_file = argv[++i]; }
        else if (arg == "--psr-maxit" && i + 1 < argc) { cfg.psr_maxit = stoi(argv[++i]); }
        else if (arg == "--psr-tol" && i + 1 < argc) { cfg.psr_tol = stod(argv[++i]); }
        else if (arg == "--psr-mu-eps" && i + 1 < argc) { cfg.psr_mu_eps = stod(argv[++i]); }
        else if (arg == "--psr-min-weight" && i + 1 < argc) { cfg.psr_min_weight = stod(argv[++i]); }
        else if (arg == "--psr-max-weight" && i + 1 < argc) { cfg.psr_max_weight = stod(argv[++i]); }
        else if (arg == "--covar" && i + 1 < argc) { cfg.covar_file = argv[++i]; }
        else if (arg == "--qcovar" && i + 1 < argc) { cfg.qcovar_file = argv[++i]; }
        else if (arg == "--make-grm") { cfg.make_grm = true; }
        else if (arg == "--make-grm-alg" && i + 1 < argc) { cfg.grm_alg = stoi(argv[++i]); }
        else if (arg == "--make-grm-sparse") { cfg.make_grm_sparse = true; }
        else if (arg == "--sparse-cutoff" && i + 1 < argc) { cfg.sparse_cutoff = stod(argv[++i]); }
        else if (arg == "--gwas") { cfg.run_gwas = true; }
        else if (arg == "--fast-lmm") { cfg.fast_lmm = true; }
        else if (arg == "--repeatability-mixed-precision") { cfg.repeatability_mixed_precision = true; }
        else if (arg == "--repeatability-grm-cache" && i + 1 < argc) { cfg.repeatability_grm_cache_prefix = argv[++i]; }
        else if (arg == "--repeatability-refit-validation") { cfg.repeatability_refit_validation = true; }
        else if (arg == "--repeatability-refit-limit" && i + 1 < argc) { cfg.repeatability_refit_limit = stoi(argv[++i]); }
        else if (arg == "--repeatability-mode" && i + 1 < argc) { cfg.repeatability_mode = argv[++i]; }
        else if (arg == "--spa") { cfg.use_spa = true; }
        else if (arg == "--skat") { cfg.run_skat = true; }
        else if (arg == "--burden") { cfg.run_burden = true; }
        else if (arg == "--set-file" && i + 1 < argc) { cfg.set_file = argv[++i]; }
        else if (arg == "--biobank-mode") { cfg.biobank_mode = true; }
        else if (arg == "--block-size" && i + 1 < argc) { cfg.block_size = stoi(argv[++i]); }
        else if (arg == "--no-covar-proj") { cfg.use_covariate_projection = false; }
        else if (arg == "--force-slow") { cfg.force_slow = true; }
        else if (arg == "--no-loco") { cfg.no_loco = true; }
        else if (arg == "--maf" && i + 1 < argc) { cfg.min_maf = stod(argv[++i]); }
        else if (arg == "--hwe" && i + 1 < argc) { cfg.min_hwe = stod(argv[++i]); }
        else if (arg == "--snp-weight" && i + 1 < argc) { cfg.snp_weight_file = argv[++i]; }
        else if (arg == "--pop-class" && i + 1 < argc) { cfg.pop_class_file = argv[++i]; }
        else if (arg == "--reml-maxit" && i + 1 < argc) { cfg.reml_maxit = stoi(argv[++i]); }
        else if (arg == "--reml-tol" && i + 1 < argc) { cfg.reml_tol = stod(argv[++i]); }
        else if (arg == "--reml-vg" && i + 1 < argc) { cfg.reml_vg = stod(argv[++i]); }
        else if (arg == "--reml-vpe" && i + 1 < argc) { cfg.reml_vpe = stod(argv[++i]); }
        else if (arg == "--reml-ve" && i + 1 < argc) { cfg.reml_ve = stod(argv[++i]); }
        else if (arg == "--help" || arg == "-h") { print_usage(argv[0]); return 0; }
    }

    if (cfg.bfile_prefix.empty() && cfg.bgen_file.empty() && cfg.pfile_prefix.empty() && cfg.pgen_file.empty()) {
        print_usage(argv[0]);
        return 1;
    }

    GwasApp app(cfg);
    return app.run();
}
