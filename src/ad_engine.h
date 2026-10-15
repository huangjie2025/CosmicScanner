#pragma once
#include <cmath>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>
#include <Eigen/Dense>
#include "plink_reader.h"
#include "pgen_reader.h"
#include "eigen_vce_ad.h"

namespace cosmic {

// Per-SNP result of an additive + dominance joint scan
// (a + d + A + D model; y = Xa*ba + Xd*bd + u_A + u_D + e).
struct AdResult {
    std::string chrom;
    std::string snp_id;
    long long pos = 0;
    std::string a1;                // Effect (minor) allele
    std::string a2;                // Reference (major) allele
    double freq = 0.0;             // Frequency of A1
    int n_obs = 0;
    double beta_add = 0.0;
    double se_add = 0.0;
    double p_add = std::numeric_limits<double>::quiet_NaN();
    double beta_dom = 0.0;
    double se_dom = 0.0;
    double p_dom = std::numeric_limits<double>::quiet_NaN();
    double log10_p_dom = 0.0;
    double joint_chisq = std::numeric_limits<double>::quiet_NaN();  // 2-df Wald/score statistic
    double p_joint = std::numeric_limits<double>::quiet_NaN();      // joint 2-df P
    double log10_p_joint = 0.0;
    bool rank_deficient = true;    // true for array slots skipped by filters
};

// Additive-Dominance LMM scan engine.
//
// Consumes an already-solved EigenVceAD null model (sigma_A^2, sigma_D^2,
// sigma_e^2, K_A, K_D). For every biallelic SNP it forms the two centered
// design columns
//   g_a = dosage - 2p                        (additive)
//   g_d = I(heterozygote) - mean(I(het))     (dominance)
// rotates both into the V eigenbasis (W_proj n x m) and computes the
// P-adjusted 2-column score / information
//   score_j = G_j' P y ,   I_j = G_j' P G_j ,   P = V^{-1} - V^{-1}X(X'V^{-1}X)^{-1}X'V^{-1}
// and the conditional joint-model Wald tests
//   beta_j = solve(I_j, score_j),  Cov(beta_j) = I_j^{-1}
//   P_ADD  = chi2_1( beta_add^2 / [I_j^{-1}]_11 )
//   P_DOM  = chi2_1( beta_dom^2 / [I_j^{-1}]_22 )
//   joint  = score_j' I_j^{-1} score_j ~ chi2_2
//
// Only hardcall input is supported; the dominance design column is
// undefined for imputed expected dosages. BED and PGEN both supply exact
// hardcall allele counts (PGEN decodes AA/AB/BB/missing per variant), so
// both backends share the identical design columns, centering rules, and
// scoring code below.
class LmmEngineAD {
public:
    struct Options {
        int thread_num = 1;
        double min_maf = 0.01;
        double min_hwe = 0.0;
        std::unordered_set<std::string> extract_set;
        std::unordered_set<std::string> exclude_set;
        int64_t scan_from_bp = -1;
        int64_t scan_to_bp = -1;
    };

    LmmEngineAD() {}
    LmmEngineAD(const Options& opts) : options(opts) {}

    void prepareFromVce(const EigenVceAD& vce);

    double getSigma2a() const { return sigma2_a; }
    double getSigma2d() const { return sigma2_d; }
    double getSigma2e() const { return sigma2_e; }

    // Hardcall scan backends. Both index results by the reader's own variant
    // order, so results.size() equals that reader's variant count.
    void scan(const PlinkReader& reader, std::vector<AdResult>& results,
              const std::vector<int>& keep_indices,
              const std::string& chrom_filter = "");
    void scan(PgenReader& reader, std::vector<AdResult>& results,
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

    // Shared per-block two-column joint Wald evaluation; guarantees BED and
    // PGEN produce identical statistics for identical hardcall input.
    void evaluateBlock(const ScanWorkspace& workspace,
                       const Eigen::MatrixXd& ga_block,
                       const Eigen::MatrixXd& gd_block,
                       const std::vector<int>& snp_indices,
                       const std::vector<double>& freqs,
                       const std::vector<VariantMeta>& meta,
                       std::vector<AdResult>& results,
                       long long& tested,
                       long long& rank_deficient) const;

    Options options;

    Eigen::VectorXd D_V;
    Eigen::MatrixXd W_proj;   // n x m full rotation (A * U_A * W_local)
    Eigen::VectorXd Uty;
    Eigen::MatrixXd UtX;      // 0 columns in the projection-based AD null

    double sigma2_a = 0.0;
    double sigma2_d = 0.0;
    double sigma2_e = 0.0;
};

} // namespace cosmic
