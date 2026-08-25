#pragma once

#include <string>
#include <vector>
#include <Eigen/Dense>

#include "bgen_reader.h"
#include "config.h"
#include "data_filter.h"
#include "pgen_reader.h"
#include "plink_reader.h"

namespace cosmic {

class GemmaVcModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_gwas";
        std::string pheno_file;
        std::string qcovar_file;
        std::string covar_file;
        std::string extract_file;
        std::string exclude_file;
        int mpheno = 1;
        int mode = 2;
        int grm_alg = 1;
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int n_threads = 1;
        int reml_maxit = 100;
        double reml_tol = 1e-6;
        bool use_covariate_projection = true;
        std::string lmm_null_solver = "brent";

        static Options fromConfig(const GwasConfig& cfg);
    };

    explicit GemmaVcModel(const Options& opts);

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
    struct VcResult {
        std::string method;
        double sigma2_g = 0.0;
        double sigma2_e = 0.0;
        double pve = 0.0;
        double delta = 0.0;
        int n = 0;
        int n_covars = 0;
        int n_snps = 0;
    };

    VcResult runHe(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const Eigen::MatrixXd& G) const;
    VcResult runReml(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const Eigen::MatrixXd& G) const;
    Eigen::MatrixXd buildGrm() const;
    void writeResult(const VcResult& res) const;

    Options options_;
    int n_samples_ = 0;
    std::vector<std::string> fids_, iids_;
    bool is_bgen_ = false, is_pgen_ = false;
    PlinkReader* plink_reader_ = nullptr;
    BgenReader* bgen_reader_ = nullptr;
    PgenReader* pgen_reader_ = nullptr;
    Eigen::VectorXd y_;
    Eigen::MatrixXf X_;
    DataFilter::KeepInfo keep_;
};

} // namespace cosmic
