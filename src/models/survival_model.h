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

namespace cosmic {

class SurvivalModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_survival";
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int grm_alg = 1;
        int n_threads = 1;
        bool no_loco = false;
        bool use_spa = false;
        int time_col = 3;
        int event_col = 4;
        int reml_maxit = 100;
        double reml_tol = 1e-6;
        int mpheno = 1;
        std::string pheno_file;
        std::string qcovar_file;
        std::string covar_file;
        std::string extract_file;
        std::string exclude_file;

        static Options fromConfig(const GwasConfig& cfg);
    };

    SurvivalModel(const Options& opts);

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
    Eigen::VectorXd survival_time_;
    Eigen::VectorXd event_flag_;
    DataFilter::KeepInfo keep_;

    double sigma2_g_ = 0.0;
    double sigma2_e_ = 1.0;

    void loadPhenotype();
    void loadCovariates();
    void runCoxNull();
    void martingaleResiduals(Eigen::VectorXd& M, Eigen::VectorXd& eta);
    void scanGwas();
    void writeNullSummary() const;
};

} // namespace cosmic
