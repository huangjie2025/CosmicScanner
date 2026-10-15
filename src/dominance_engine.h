#pragma once
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>
#include <Eigen/Dense>
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "eigen_vce_d.h"

namespace cosmic {

// Per-SNP result of a dominance-only GWAS scan
// (d+D model; y = Xd*bd + u_D + e).
struct DomResult {
    std::string chrom;
    std::string snp_id;
    long long pos = 0;
    std::string a1;                // Effect (minor) allele
    std::string a2;                // Reference (major) allele
    double freq = 0.0;             // Frequency of A1
    int n_obs = 0;
    double beta_dom = 0.0;
    double se_dom = 0.0;
    double p_dom = std::numeric_limits<double>::quiet_NaN();
    double log10_p_dom = 0.0;
    bool rank_deficient = true;    // true for array slots skipped by filters
};

// Dominance-only LMM scan engine.
//
// Consumes an already-solved EigenVceD null model (sigma_D^2, sigma_e^2, K_D).
// For every biallelic SNP it forms the centered dominance design column
//   g_d = I(heterozygote) - 2p(1-p)  (centered heterozygote indicator)
// rotates into the V eigenbasis (W_proj n x m) and computes the
// P-adjusted 1-column score / information
//   score_d = G_d' P y ,   I_d = G_d' P G_d ,   P = V^{-1} - V^{-1}X(X'V^{-1}X)^{-1}X'V^{-1}
// and the conditional Wald test
//   beta_d = score_d / I_d,  Var(beta_d) = I_d^{-1}
//   P_D  = chi2_1( beta_d^2 / Var(beta_d) )
//
// BED and PGEN use exact hardcall heterozygosity. BGEN uses posterior
// heterozygosity P(AB), centered over observed analysis samples. A one-hot
// BGEN therefore produces the identical marker column.
class LmmEngineD {
public:
    struct Options {
        int thread_num = 1;
        double min_maf = 0.01;
        double min_hwe = 0.0;
        double bgen_dom_max_missing = 0.05;
        double bgen_dom_min_info = 0.01;
        double bgen_dom_entropy_warning = 0.50;
        std::unordered_set<std::string> extract_set;
        std::unordered_set<std::string> exclude_set;
        int64_t scan_from_bp = -1;
        int64_t scan_to_bp = -1;
    };

    LmmEngineD() {}
    LmmEngineD(const Options& opts) : options(opts) {}

    void prepareFromVce(const EigenVceD& vce);

    double getSigma2d() const { return sigma2_d; }
    double getSigma2e() const { return sigma2_e; }

    // Hardcall scan backends. Both index results by the reader's own variant
    // order, so results.size() equals that reader's variant count.
    void scan(const PlinkReader& reader, std::vector<DomResult>& results,
              const std::vector<int>& keep_indices,
              const std::string& chrom_filter = "");
    void scan(PgenReader& reader, std::vector<DomResult>& results,
              const std::vector<int>& keep_indices,
              const std::string& chrom_filter = "");
    void scan(BgenReader& reader, std::vector<DomResult>& results,
              const std::vector<int>& keep_indices,
              const std::string& chrom_filter = "");

private:
    // Per-variant output metadata, filled identically for both backends.
    // a1 is the counted allele (BED A1 / PGEN ALT), a2 is the other allele.
    struct VariantMeta {
        std::string chrom;
        std::string snp_id;
        std::string a1;
        std::string a2;
        uint32_t pos = 0;
    };

    struct ScanWorkspace {
        Eigen::VectorXd p_y;
        Eigen::VectorXd v_inv_diag;
    };

    // Shared null-model projection quantities; backend independent.
    ScanWorkspace prepareWorkspace() const;

    // Shared per-block Wald evaluation; guarantees BED and PGEN produce
    // identical statistics for identical hardcall input.
    void evaluateBlock(const ScanWorkspace& workspace,
                       const Eigen::MatrixXd& gd_block,
                       const std::vector<int>& snp_indices,
                       const std::vector<double>& freqs,
                       const std::vector<VariantMeta>& meta,
                       std::vector<DomResult>& results,
                       long long& tested,
                       long long& rank_deficient) const;

    Options options;

    Eigen::VectorXd D_V;      // Eigenvalues of V
    Eigen::MatrixXd W_proj;   // n x m full rotation (A * U_D * W_local)
    Eigen::VectorXd Uty;
    Eigen::MatrixXd UtX;      // m x n_covars in the projected space

    double sigma2_d = 0.0;
    double sigma2_e = 0.0;
};

} // namespace cosmic
