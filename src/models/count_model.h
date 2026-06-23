#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "data_filter.h"
#include "eigen_vce.h"

namespace cosmic {

enum class CountFamily {
    Poisson,
    NegativeBinomial,
    Binomial
};

class CountModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_count";
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int grm_alg = 1;
        int n_threads = 1;
        bool no_loco = false;
        bool fast_lmm = false;
        bool use_spa = false;
        CountFamily family = CountFamily::Poisson;
        bool estimate_overdispersion = false;
        int pql_maxit = 15;
        double pql_tol = 1e-6;
        double mu_eps = 1e-6;
        double min_weight = 1e-8;
        double max_weight = 1e8;
        int reml_maxit = 100;
        double reml_tol = 1e-6;
        int mpheno = 1;
        std::string pheno_file;
        std::string qcovar_file;
        std::string covar_file;
        std::string offset_file;
        std::string binomial_file;

        static Options fromConfig(const GwasConfig& cfg);
    };

    CountModel(const Options& opts);

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
    int n_samples_ = 0;
    std::vector<std::string> fids_, iids_;
    bool is_bgen_ = false, is_pgen_ = false;
    PlinkReader* plink_reader_ = nullptr;
    BgenReader* bgen_reader_ = nullptr;
    PgenReader* pgen_reader_ = nullptr;

    Eigen::VectorXd y_;
    Eigen::VectorXd offset_;
    Eigen::VectorXd trials_;
    Eigen::MatrixXf X_;
    DataFilter::KeepInfo keep_;

    double sigma2_g_ = 0.0;
    double sigma2_e_ = 1.0;
    double overdispersion_ = 0.1;

    Eigen::VectorXd mu_;
    Eigen::VectorXd working_y_;
    Eigen::VectorXd weights_;
    Eigen::VectorXd eta_;

    void loadPhenotype();
    void loadCovariates();
    void loadOffset();
    void validateInputs() const;
    void runPQL();
    double linkFun(double x) const;
    double linkDeriv(double x) const;
    double varianceFun(double mu) const;
    void scanGwas();
    void writeNullSummary() const;
    void writeConvergenceLog(const std::vector<std::string>& rows) const;
};

} // namespace cosmic
