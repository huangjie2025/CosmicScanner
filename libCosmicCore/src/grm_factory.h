#pragma once
/// @file grm_factory.h
/// Unified GRM building helper shared by CosmicSolver and CosmicScanner.
/// Encapsulates the common GrmBuilder::Options + compute + getGrm + extract IDs pattern.
/// Eliminates code duplication across repeatability_model, gwas_scanner, multivariate_model,
/// and matrix_builder.

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "grm_builder.h"
#include "plink_reader.h"
#include "pgen_reader.h"
#include "bgen_reader.h"

namespace cosmic {

/// Configuration for a single GRM build operation.
struct GrmBuildConfig {
    int threads = 1;
    float min_maf = 0.01f;
    float min_hwe = 1e-6f;
    int algorithm = 0;           // 0=VanRaden, 1=Yang, 2=Centered, etc.
    bool compute_homo_hete = false;
    std::string snp_weight_file;
    std::string pop_class_file;
    int block_size = -1;
    std::string chrom_filter;    // empty = all chromosomes (for LOCO)
};

/// Result of a single GRM build operation.
struct GrmBuildResult {
    Eigen::MatrixXd G;          ///< The computed GRM (double)
    std::vector<std::string> ids;   ///< Sample IDs in order
    double total_snps = 0;      ///< Total SNPs used (algorithm-dependent)
    double sum2pq = 0;          ///< Sum(2pq) for Yang normalization
    std::set<std::string> observed_chroms;  ///< Chromosomes observed
};

/// Factory class for building GRMs from various genotype sources.
class GrmFactory {
public:
    /// Build a global GRM from any genotype source.
    /// @param cfg       Build configuration
    /// @param bfile     PLINK prefix (if using PLINK)
    /// @param pfile     PGEN prefix (if using PGEN)
    /// @param plink     Pre-loaded PlinkReader (optional, avoids re-loading)
    /// @param pgen      Pre-loaded PgenReader (optional)
    /// @param bgen      Pre-loaded BgenReader (optional)
    static GrmBuildResult buildGrm(const GrmBuildConfig& cfg,
                                   const std::string& bfile = "",
                                   const std::string& pfile = "",
                                   PlinkReader* plink = nullptr,
                                   PgenReader* pgen = nullptr,
                                   BgenReader* bgen = nullptr);

    /// Build a chromosome-specific GRM for LOCO (Leave-One-Chromosome-Out).
    /// Uses the same options but filters to a specific chromosome.
    static GrmBuildResult buildChrGrm(const GrmBuildConfig& base_cfg,
                                      const std::string& chrom,
                                      PlinkReader* plink = nullptr,
                                      PgenReader* pgen = nullptr,
                                      BgenReader* bgen = nullptr);

    /// Extract sample IDs from a genotype reader into a vector.
    static std::vector<std::string> extractIds(PlinkReader& reader);
    static std::vector<std::string> extractIds(PgenReader& reader);
    static std::vector<std::string> extractIds(BgenReader& reader);

private:
    /// Internal: apply config to GrmBuilder and compute
    template<typename ReaderT>
    static GrmBuildResult doCompute(const GrmBuilder::Options& opts, ReaderT& reader) {
        GrmBuilder builder(opts);
        builder.compute(reader);
        GrmBuildResult result;
        result.G = builder.getGrm();
        result.ids = extractIds(reader);
        result.total_snps = opts.algorithm == 0 ? (double)builder.getTotalSnps() : builder.getSum2pq();
        result.sum2pq = builder.getSum2pq();
        result.observed_chroms = builder.getObservedChroms();
        return result;
    }
};

} // namespace cosmic
