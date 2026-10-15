#pragma once
#include <string>
#include <iosfwd>
#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"
#include "data_filter.h"
#include "pedigree_relationship.h"
#include "apy_ginv.h"
#include "apy_core_pev_operator.h"
#include "apy_inverse_operator.h"
#include "block_streamed_grm.h"
#include "cross_genomic_operator.h"
#include "a22_inverse_operator.h"
#include "design.h"

namespace cosmic {
struct SnpUncertaintyScanResult;

class SingleStepModel {
public:
    enum class Method { Lmm, Backsolve };

    struct Options {
        std::string out_prefix = "cosmic_ssgwas";
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int grm_alg = 1;
        int n_threads = 1;
        bool no_loco = false;
        bool fast_lmm = false;
        bool use_spa = false;
        int reml_maxit = 100;
        double reml_tol = 1e-6;
        int mpheno = 1;
        std::string pheno_file;
        std::string qcovar_file;
        std::string covar_file;
        std::string ped_file;
        std::string genotyped_list;
        double h_tau = 1.0;
        double h_omega = 1.0;
        double blend_g = 0.05;
        bool write_ebv = false;
        bool make_agh = false;
        bool write_agh = false;
        std::string extract_file;
        std::string exclude_file;

        // Method 2 (ssGBLUP-derived back-solving), Phase 6
        Method method = Method::Lmm;
        int ssgwas_window_snps = 10;
        bool ssgwas_tune_g = false;
        bool ssgwas_blend_a22 = false;       // blend G* with A22 instead of I (HIBLUP/preGSf90 reference); backsolve only
        std::string ssgwas_solver = "ldlt";  // "ldlt" | "pcg"
        std::string ssgwas_vce;              // "" = fixed VC path
        double reml_vg = -1.0;               // fixed VC (used when ssgwas_vce empty)
        double reml_ve = -1.0;
        double ssgwas_pcg_tol = 1e-8;
        int ssgwas_pcg_maxit = 2000;
        // Phase C1: frequentist P-values for back-solved SNP effects (Aguilar 2019).
        std::string ssgwas_pvalue_method = "off";  // "off" | "exact" | "apy-exact" | "apy-approx"
        std::string ssgwas_apy_exact_estimator = "legacy-hybrid-v1";
        bool ssgwas_apy_exact_estimator_explicit = false;
        bool allow_large_exact_ssgwas_pvalue = false;
        // Phase E2: block-streamed VanRaden G (no full n_g^2 G).
        // "auto" | "dense" | "streamed".
        std::string ssgwas_grm_mode = "auto";
        // Phase E3-B1: Gcc solver backend.
        std::string ssgwas_gcc_solver = "dense";  // "dense" | "pcg"
        std::string ssgwas_apy_inverse_mode = "dense";  // "dense" | "operator"
        std::string ssgwas_mme_apy_mode = "dense";  // "dense" | "operator"
        std::string ssgwas_pev_apy_mode = "dense";  // "dense" | "operator"
        std::string ssgwas_hinv_mode = "dense";  // "dense" | "operator"
        std::string ssgwas_pev_scan_mode = "batched";  // "scalar" | "batched"
        bool ssgwas_pev_scan_mode_explicit = false;
        int ssgwas_pev_scan_block_size = 256;
        // Phase E3-B2: PCG preconditioner.
        std::string ssgwas_pcg_preconditioner = "jacobi";  // "jacobi" | "lowrank"
        int ssgwas_pcg_preconditioner_rank = 50;
        // Phase C2: APY core partition (Leite et al. 2024, GSE 56:59).
        int ssgwas_apy_core_size = 0;
        std::string ssgwas_apy_core_file;
        int ssgwas_apy_seed = 12345;
        double ssgwas_apy_core_target_ratio = 0.0;  // harness-assigned label (0 = unset)
        double ssgwas_apy_core_variance = 0.99;    // auto core-size target (C98=0.98, C99=0.99)
        std::string ssgwas_apy_core_eigensolver = "auto";  // "auto" | "exact" | "lanczos"

        static Options fromConfig(const GwasConfig& cfg);
    };

    SingleStepModel(const Options& opts);

    void prepare(const std::vector<std::string>& fids,
                 const std::vector<std::string>& iids,
                 int n_samples,
                 PlinkReader& plink_reader,
                 BgenReader& bgen_reader,
                 PgenReader& pgen_reader,
                 VcfReader& vcf_reader,
                 bool is_bgen, bool is_pgen, bool is_vcf);

    void run();

private:
    Options options_;
    int n_samples_ = 0;
    std::vector<std::string> fids_, iids_;
    bool is_bgen_ = false, is_pgen_ = false, is_vcf_ = false;
    PlinkReader* plink_reader_ = nullptr;
    BgenReader* bgen_reader_ = nullptr;
    PgenReader* pgen_reader_ = nullptr;
    VcfReader* vcf_reader_ = nullptr;

    Eigen::VectorXd y_;
    Eigen::MatrixXf X_;
    DataFilter::KeepInfo keep_;

    // Pedigree data (shared PedigreeRelationship builder; also used by the
    // Pedigree Animal Model so that A / A22 construction is common code)
    PedigreeRelationship pedigree_;

    // A⁻¹ in sparse format
    Eigen::SparseMatrix<double> A_inv_;

    // GRM for genotyped individuals
    Eigen::MatrixXd G_;
    Eigen::MatrixXd A22_inv_;
    std::unique_ptr<A22InverseOperator> a22_inverse_operator_;
    int n_genotyped_ = 0;
    std::vector<int> genotyped_indices_;
    std::unordered_set<std::string> requested_genotyped_ids_;
    std::vector<std::string> genotyped_ids_;
    bool has_a22_correction_ = false;
    std::string hinv_preconditioner_diagonal_mode_ = "dense_exact";
    double h_inverse_symmetry_error_ = 0.0;

    // H⁻¹ = A⁻¹ + correction
    Eigen::SparseMatrix<double> H_inv_;
    Eigen::MatrixXd H_sub_;

    double sigma2_g_ = 0.0;
    double sigma2_e_ = 0.0;

    void loadPhenotype();
    void loadCovariates();
    void loadPedigree();
    void buildAInverse();
    void buildGRM();
    void computeA22();
    void computeHInverse();
    Eigen::MatrixXd buildHSubForG(const Eigen::MatrixXd& G) const;
    void runVCE();
    void scanGwas();
    void writeEbv();
    void writeAghMatrices() const;
    void writeNullSummary() const;

    // --- Method 2: ssGBLUP genetic evaluation + back-solving (Phase 6) ---
    // Full-animal MME: y = Xb + Za + e, a ~ N(0, H sigma_a^2), over the complete
    // pedigree animal set. Requires fixed VC (--reml-vg/--reml-ve) or estimated
    // VC via libCosmicCore VCESolver (--ssgwas-vce). Optional Phase C1
    // frequentist SNP P-values (--ssgwas-pvalue-method exact, Aguilar 2019).
    void runBacksolve();
    void buildBacksolveGrm();   // split-path G (+ S = sum 2p(1-p)) with optional tuneG
    void setupGenotypedSamples();  // shared genotyped selection + pedigree mapping (steps 1-2)
    void resolveStreamedGrmMode(); // Phase E2: set grm_streamed_ from --ssgwas-grm-mode / n_g
    void buildBacksolveGrmStreamed();      // Phase E2 pass 1: stats, S, tune scalars, A22 cache
    void buildBacksolveGrmStreamBlocks();  // Phase E2 pass 2: Gcc/Gcn/diag(Gnn) + APY inverses
    void materializeApyMmeReference();
    int runStreamedK98() const;            // Phase E2: Lanczos K98 of the final blended G
    double gwCc(int a, int b) const;       // final blended G core block element (dense or streamed)
    void loadBacksolvePhenotypes();  // phenotype loading keyed by pedigree ID
    void solveGeneticEvaluation();          // MME -> beta_hat, a_hat
    void backsolveSnpEffects();             // q = Ginv*a_g; Z'q; scale conversion
    void computeSnpPvaluesExact();          // Phase C1: selected MME inverse -> Var(u_hat_g) -> SNP PEC
    void computeSnpPvaluesApyExact();       // Phase C2: APY G^-1 MME + core-slice selected inverse (Eq. 7)
    void computeSnpPvaluesApyApprox();      // Phase D: APY G^-1 MME + Eq. 8/9 approximate core PEV (no LHS inverse)
    void acceptExactScanResult(SnpUncertaintyScanResult&& result);
    void writeExactScanAudit(std::ostream& out) const;
    void buildApyCoreSelection();           // Phase C2: core/noncore partition (file or seeded random)
    void writeEbvFile() const;
    void writeSnpFile() const;
    void writePvalueFile() const;
    void writeWindowFile() const;
    void writeMetaFile() const;

    void computeEffectiveRecordDiagonals(const std::vector<GenRecord>& recs,
                                         const FixedDesignG& fd);

    // Method-2 state
    Eigen::MatrixXd G_bs_;        // back-solving GRM (split path, tuneG/blend-applied)
    Eigen::MatrixXd Ginv_bs_;     // ridge inverse of G_bs_
    Eigen::MatrixXd G_bs_tuned_;  // tuned-but-unblended G* (effect back-solve matrix when A22-blend active)
    Eigen::MatrixXd Ginv_bs_tuned_; // ridge inverse of G_bs_tuned_
    double S_2pq_ = 0.0;          // sum_j 2p_j(1-p_j) over the GRM SNP set (VanRaden denominator)
    double tune_a_ = 0.0;         // G = b*G + a alignment factors (0/1 when tuneG off)
    double tune_b_ = 1.0;
    double a22_mean_diag_ = 0.0;  // A22 (pedigree block of genotyped) diagonal/off-diagonal means
    double a22_mean_off_ = 0.0;
    double g_raw_mean_diag_ = 0.0;  // raw G (pre-tune, pre-blend) and final blended G mean summaries
    double g_raw_mean_off_ = 0.0;
    double g_final_mean_diag_ = 0.0;
    double g_final_mean_off_ = 0.0;
    std::string blend_target_name_ = "identity";  // "identity" | "a22"
    std::string effect_ginv_name_ = "blended";    // "blended" | "tuned"
    // Phase E2: block-streamed G state.  When grm_streamed_ is true, G_bs_ /
    // G_bs_tuned_ are NOT allocated; the APY-relevant blocks (Gcc / Gcn /
    // diag(Gnn), tuned-unblended and final blended) are held instead.
    bool grm_streamed_ = false;
    std::unique_ptr<BlockStreamedGrmBuilder> grm_stream_builder_;
    Eigen::MatrixXd gcc_stream_final_;   // final blended core block
    Eigen::VectorXd gnn_diag_stream_final_;
    Eigen::MatrixXd gcc_stream_tuned_;   // tuned-unblended core block
    Eigen::VectorXd gnn_diag_stream_tuned_;
    Eigen::MatrixXd A22_geno_;     // full A22 of genotyped animals (tune/blend-a22 cache)
    std::vector<int> geno_full_idx_;    // genotyped sample indices in sample order (PLINK keep)
    // Phase E3-D2D: matrix-free Gcn cross operators for the streamed operator
    // path (only when the K x Nn Gcn block is NOT materialized).  gcn_cross_final_
    // serves the MME / PEV operators; gcn_cross_tuned_ serves the effect
    // back-solve operator under the A22-blend reference path.
    std::unique_ptr<CrossGenomicOperator> gcn_cross_final_;
    std::unique_ptr<CrossGenomicOperator> gcn_cross_tuned_;
    long long gcn_regen_after_construct_ = 0;   // E3-D2D no-regen diagnostic
    std::vector<int> pheno_ped_rows_;   // pedigree rows with a phenotype (Method 2)
    std::vector<double> pheno_y_;       // phenotypes aligned to pheno_ped_rows_
    std::vector<int> snp_keep_indices_; // SNP indices (into bim/zt outputs) kept for reporting
    Eigen::VectorXd beta_hat_;    // fixed effects (p x 1)
    Eigen::VectorXd a_hat_;       // breeding values for ALL pedigree animals
    Eigen::VectorXd a_hat_g_;     // breeding values for genotyped animals (genotyped order)
    Eigen::VectorXd snp_effect_allele_;   // per-allele-substitution scale (M = x - 2p convention)
    Eigen::VectorXd snp_effect_std_;      // per standardized-genotype SD: EA_j * sd_j
                                          // (Z = (x-mean)/sd); Z*std = M*EA = GEBV(typed)
    Eigen::VectorXd snp_variance_;        // 2p(1-p) * effect_allele^2
    double total_snp_variance_ = 0.0;     // denominator for PCT_GENETIC_VARIANCE
    // timing (seconds)
    double t_parse_ = 0.0, t_ainv_ = 0.0, t_g_ = 0.0, t_ginv_ = 0.0;
    long long point_multiply_zt_v_calls_ = 0;
    std::string exact_scan_actual_mode_ = "SCALAR_ONEHOT_REFERENCE";
    double t_hnv_ = 0.0, t_mme_ = 0.0, t_solve_ = 0.0, t_bs_ = 0.0;
    double t_snp_ = 0.0, t_win_ = 0.0;
    int n_ped_ = 0, n_phenotyped_ = 0;
    int n_genotyped_with_pheno_ = 0, n_ungenotyped_with_pheno_ = 0;
    int solver_iterations_ = 0;
    double solver_residual_ = 0.0;
    bool solver_converged_ = true;
    std::string vce_mode_used_;
    double sigma2_e_backsolve_ = 0.0;
    double sigma2_a_backsolve_ = 0.0;

    // --- Phase C1: exact SNP P-value state (--ssgwas-pvalue-method exact) ---
    // Convention (docs/SSGWAS_SNP_PVALUE_METHOD.md, C0 gate):
    //   MME (lambda-form, no sigma_e^-2): C_lambda theta = rhs
    //   PEV(a)          = sigma_e2 * C_lambda^-1            (breeding-value PEV)
    //   Var(u_hat_g)    = G_eval * sigma_a2 - sigma_e2 * Cinv_gg
    //                     (variance of the predictor; G_eval = final blended Gw,
    //                      Cinv_gg = genotyped block of C_lambda^-1)
    //   Var(EA_j)       = B_j * Var(u_hat_g) * B_j', B_j = M_j' G_backsolve^-1 / S
    //                     (sampling variance of the SNP-effect estimate; G_backsolve
    //                      = tuned-unblended G* when blend-a22, else G_eval)
    Eigen::SparseMatrix<double> mme_lhs_;   // C_lambda (lambda-form), kept for PEV factorization
    Eigen::MatrixXd cinv_gg_;               // genotyped block of C_lambda^-1 (selected solves)
    Eigen::MatrixXd vhat_g_;                // Var(u_hat_g) = G_eval*sigma_a2 - sigma_e2*Cinv_gg

    // --- Phase C2: APY exact SNP P-value state (--ssgwas-pvalue-method apy-exact) ---
    // Convention (docs/SSGWAS_APY_PVALUE_METHOD.md; Leite et al. 2024 Eqs. 5-8):
    //   MME genotyped block   = tau*G_APY^-1 - omega*A22^-1     (APY inverse of Gw)
    //   Var(u_hat_core)       = (Gw)_cc*sigma_a2 - sigma_e2*(C_lambda^-1)_cc
    //   Var(EA_j)             = w_c' Var(u_hat_core) w_c / S^2,
    //                           w_c = (G_backsolve^-1)_cc M_cj   (core-only, Eq. 7)
    std::vector<int> core_indices_;         // APY core (genotyped order, ascending)
    std::vector<int> noncore_indices_;      // APY noncore (genotyped order, ascending)
    std::string apy_core_selection_mode_;   // "RANDOM_SEEDED" | "FILE"
    std::unique_ptr<ApyGinv> apy_ginv_mme_;   // APY inverse of the final Gw (MME block)
    std::unique_ptr<ApyGinv> apy_ginv_eff_;   // APY inverse of the effect matrix (G* when blend-a22, else Gw)
    std::unique_ptr<ApyInverseOperator> apy_inverse_eff_operator_;
    std::unique_ptr<ApyInverseOperator> apy_inverse_mme_operator_;
    Eigen::MatrixXd cinv_cc_;               // core slice of C_lambda^-1 ((C_lambda^-1)_cc, selected solves)
    Eigen::MatrixXd vhat_c_;                // Var(u_hat_core) = (Gw)_cc*sigma_a2 - sigma_e2*Cinv_cc
    std::string cinv_cc_max_asym_;          // max |Cinv_cc - Cinv_cc'| before symmetrization
    std::string apy_gcc_min_diag_d_, apy_gcc_max_diag_d_;  // Gcc LDLT pivot diagnostics
    double apy_mnn_min_ = 0.0, apy_mnn_max_ = 0.0;
    int apy_n_core_ = 0, apy_n_noncore_ = 0;
    bool apy_core_auto_ = false;                // true if --ssgwas-apy-core-size auto was used
    int apy_core_auto_K98_ = -1;                // K98 value when auto mode
    int apy_core_auto_K99_ = -1;                // K99 value when auto mode
    double apy_core_auto_variance_ = 0.0;       // target variance when auto mode
    std::string apy_core_eigensolver_used_;     // "dense" | "lanczos"
    int apy_core_eig_iterations_ = 0;           // Lanczos iterations
    double apy_core_trace_g_ = 0.0;             // trace(G) computed for auto core

    // --- Phase D: approximate SNP P-value state (--ssgwas-pvalue-method
    // apy-approx; Leite et al. 2024 Eq. 8/9, single Schur/Woodbury form) ---
    //   C_core_approx^-1 = D_c + lambda*Gcc^-1 + Q*S*Q'  (S_j = lambda d_j/(lambda + d_j m_jj))
    //   D (Eq. 9)        = diag(W'(I - X(X'X)^-1 X')W)  (fixed-effect absorption)
    //   Var(EA_j)        = w_c' (Gw_cc*sigma_a2 - sigma_e2*C_core_approx) w_c / S^2
    // No dense n_g^2 uncertainty block is ever allocated in this mode.
    std::vector<double> d_genotyped_;       // Eq. 9 diagonal effective records, genotyped order (n_g)
    int k98_g_ = -1;                        // eigen-spectrum K98 of final G_bs_ (cumulative 98% variance)
    double k98_cumfraction_ = 0.98;
    Eigen::MatrixXd pec_c_appx_;            // C_core_approx (n_c^2), Eq. 8 Schur form
    Eigen::MatrixXd vhat_c_appx_;           // Var(u_hat_core) approx = Gw_cc*sigma_a2 - sigma_e2*Pec_c_appx
    int dense_uncertainty_allocs_ = 0;      // debug counter; only exact/apy-exact paths increment
    double estimated_peak_memory_bytes_ = 0.0;
    std::string apy_pev_backend_ = "NA";
    double apy_pev_max_relative_residual_ = 0.0;

    std::vector<double> snp_se_allele_;
    std::vector<double> snp_z_;
    std::vector<double> snp_p_;
    std::vector<double> snp_log10p_;
    std::vector<std::string> snp_pvalue_status_;
    int n_pec_clamped_ = 0;                 // PEC in [-tol, 0] -> 0 (NEGATIVE_VARIANCE_CLAMPED)
    int n_pec_failure_ = 0;                 // PEC < -tol          -> NUMERICAL_FAILURE
    int n_zero_se_ = 0;                     // SE <= threshold     -> ZERO_SE
    std::string cinv_gg_max_asym_;          // max |Cinv - Cinv'| before symmetrization
    std::string var_ghat_min_eig_, var_ghat_max_eig_, var_ghat_eig_mode_;
    double t_pev_inverse_ = 0.0, t_pev_vhat_ = 0.0, t_pev_scan_ = 0.0, t_pev_out_ = 0.0;
    double t_pev_scan_setup_ = 0.0, t_pev_scan_onehot_ = 0.0;
    double t_pev_scan_genotype_ = 0.0, t_pev_scan_core_gather_ = 0.0;
    double t_pev_scan_gcc_solve_ = 0.0, t_pev_scan_gw_matvec_ = 0.0;
    double t_pev_scan_pec_matvec_ = 0.0, t_pev_scan_quadratic_ = 0.0;
    double t_pev_scan_statistics_ = 0.0, t_pev_scan_other_ = 0.0;
    long long pev_scan_multiply_z_v_calls_ = 0, pev_scan_genotype_passes_ = 0;
    long long pev_scan_bed_bytes_read_estimate_ = 0;
    long long pev_scan_variants_tested_ = 0, pev_scan_allocations_ = 0;
    long long pev_scan_blocks_ = 0, pev_scan_variants_decoded_ = 0;
    long long pev_scan_block_reallocations_ = 0, pev_scan_peak_transient_bytes_ = 0;
    long long pev_scan_gcc_solve_calls_ = 0, pev_scan_gcc_rhs_columns_ = 0;
    long long pev_scan_gw_gemm_calls_ = 0, pev_scan_pec_gemm_calls_ = 0;
    long long pev_scan_missing_to_zero_ = 0;
    int pev_scan_effective_block_size_ = 1;
    double estimated_cinv_gg_bytes_ = 0.0;
};

} // namespace cosmic
