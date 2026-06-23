#pragma once

#include <map>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"

namespace cosmic {

class OrdinalModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_ordinal";
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int grm_alg = 1;
        int n_threads = 1;
        bool no_loco = false;
        bool fast_lmm = false;
        std::string link_type = "logit";
        int pql_maxit = 50;
        double pql_tol = 1e-6;
        double mu_eps = 1e-6;
        double min_weight = 1e-8;
        double max_weight = 1e8;
        int reml_maxit = 50;
        double reml_tol = 1e-6;
        int mpheno = 1;
        std::string ordinal_mode = "experimental-cumulative-binary";
        std::string ordinal_stacked_engine = "fixed-gls";
        std::string ordinal_vc_method = "fixed";
        double ordinal_sigma2_g = 0.5;
        double ordinal_vc_grid_min = 0.0;
        double ordinal_vc_grid_max = 2.0;
        int ordinal_vc_grid_points = 21;
        bool ordinal_allow_fallback = false;
        std::string ordinal_test = "wald";
        std::string pheno_file;
        std::string ordinal_file;
        std::string qcovar_file;
        std::string covar_file;

        static Options fromConfig(const GwasConfig& cfg);
    };

    explicit OrdinalModel(const Options& opts);

    void prepare(const std::vector<std::string>& fids,
                 const std::vector<std::string>& iids,
                 int n_samples,
                 PlinkReader& plink_reader,
                 BgenReader& bgen_reader,
                 PgenReader& pgen_reader,
                 bool is_bgen,
                 bool is_pgen);

    void run();

private:
    Options options_;
    std::vector<std::string> fids_, iids_;
    int n_samples_ = 0;
    bool is_bgen_ = false, is_pgen_ = false;
    PlinkReader* plink_reader_ = nullptr;
    BgenReader* bgen_reader_ = nullptr;
    PgenReader* pgen_reader_ = nullptr;
    Eigen::VectorXd y_;
    std::vector<double> levels_;
    std::map<double, int> level_counts_;

    void loadOrdinalPhenotype();
    void loadCovariates();
    void validateInputs() const;
    std::string writeThresholdPhenotype(double threshold, int threshold_index) const;
    void writeThresholdSummary() const;
    void writeCombinedConvergence(const std::vector<std::string>& threshold_prefixes) const;
    void writeCombinedGwas(const std::vector<std::string>& threshold_prefixes, double runtime_sec) const;

    void runStackedPsr();
    void buildStackedData(int n_thresh, int n_stacked,
                          Eigen::VectorXd& z_stack,
                          Eigen::MatrixXd& X_stack,
                          Eigen::VectorXi& row2indiv) const;
    Eigen::MatrixXd buildGrmMatrix() const;
    void writeStackedThresholds(int n_thresh) const;
    void writeStackedConvergence(const std::vector<std::string>& rows) const;
    void writeStackedBlup(const Eigen::VectorXd& u_hat,
                          const Eigen::VectorXd& sigma2_g) const;

    struct GridRemlResult {
        double sigma2_g = 0.5;
        double reml_loglik = 0.0;
        bool boundary = false;
        int n_valid = 0;
        bool jitter_used = false;
        double jitter_value = 0.0;
        std::string status = "OK";
    };
    GridRemlResult estimateSigma2GridReml(
        const Eigen::VectorXd& y_tilde,
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& K_stack,
        const Eigen::VectorXd& w,
        double current_sigma2) const;
    void writeStackedVc(const std::vector<std::string>& vc_rows) const;

    Eigen::MatrixXd X_covariates_;
    bool has_covariates_ = false;

    // Stored for score diagnostics
    Eigen::VectorXd z_stack_;
    Eigen::VectorXd mu_;
    Eigen::VectorXd w_;
    Eigen::VectorXi row2indiv_;
};

} // namespace cosmic
