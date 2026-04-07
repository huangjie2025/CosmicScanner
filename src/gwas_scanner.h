#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "data_filter.h"
#include "lmm_engine.h"
#include "grm_builder.h"
#include "pheno_reader.h"
#include "rare_variant_test.h"
#include "ridge_regenie.h"

namespace cosmic {

class GwasScanner {
public:
    struct Options {
        std::string out_prefix = "cosmic_gwas";
        bool out_ssf = false; // Add switch for EBI GWAS-SSF output
        double min_maf = 0.01;
        double min_hwe = 0.0;
        bool fast_lmm = false;
        int n_threads = 1;
        bool no_loco = false;
        int grm_alg = 1; // 1: VanRaden, 0: Yang et al
        bool biobank_mode = false;
        bool use_covariate_projection = true; // Use exact covariate projection for LMM acceleration
        bool use_spa = false; // Use Saddle Point Approximation for binary traits
        int regenie_block_size = 1000;
        bool run_skat = false;
        bool run_burden = false;
        std::string set_file;
        
        Eigen::VectorXd obs_weights; // Optional weights for observations (e.g., Repeatability Model D^{-1/2})
    };

    GwasScanner(const Options& opts);

    // Run standard LOCO scan (or single chromosome scan if only 1 chrom)
    void runLocoScan(
        const Eigen::VectorXd& y,
        const Eigen::MatrixXf& X,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        bool is_bgen,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );

    // Run Biobank mode scan (Matrix-Free RHE + PCG)
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
    
    // Run Rare Variant Association Test (SKAT / Burden)
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

    // Run Random Regression Model (RRM) scan for longitudinal data
    void runRrmScan(
        const std::vector<PhenoReader::PhenoRecord>& records,
        int rrm_order,
        const Eigen::MatrixXf& X_base,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        bool is_bgen,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );

private:
    Options options;

    // Helper to append results to file
    void appendResults(const std::string& filename, const std::vector<GwasResult>& results);
    
    // Helper to write header
    void writeHeader(const std::string& filename);
};

} // namespace cosmic
