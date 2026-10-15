#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"
#include "data_filter.h"
#include "eigen_vce.h"

namespace cosmic {

enum class ThresholdLink {
    Probit,
    Logit
};

class ThresholdModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_threshold";
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int grm_alg = 1;
        int n_threads = 1;
        bool no_loco = false;
        bool fast_lmm = false;
        ThresholdLink link = ThresholdLink::Probit;
        double prevalence = -1.0;
        bool liability_scale = false;
        int pql_maxit = 15;
        double pql_tol = 1e-6;
        double mu_eps = 1e-6;
        double min_weight = 1e-8;
        double max_weight = 1e8;
        int reml_maxit = 100;
        double reml_tol = 1e-6;
        bool psr_null_only = false;
        int mpheno = 1;
        std::string pheno_file;
        std::string qcovar_file;
        std::string covar_file;
        // Explicit opt-in to dichotomise a non-binary response at 0.5 (see --threshold-dichotomize).
        bool allow_dichotomize = false;

        static Options fromConfig(const GwasConfig& cfg);
    };

    ThresholdModel(const Options& opts);

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

    double sigma2_g_ = 0.0;
    double sigma2_e_ = 1.0;

    Eigen::VectorXd mu_;
    Eigen::VectorXd working_y_;
    Eigen::VectorXd weights_;
    Eigen::VectorXd eta_;

    void loadPhenotype();
    void loadCovariates();
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
