#pragma once
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <Eigen/Dense>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "data_filter.h"
#include "eigen_vce.h"

namespace cosmic {

class MaternalModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_maternal";
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int grm_alg = 1;
        int n_threads = 1;
        bool no_loco = false;
        bool fast_lmm = false;
        bool use_spa = false;
        bool maternal_pe = true;
        bool direct_maternal_cov = true;
        int reml_maxit = 100;
        double reml_tol = 1e-6;
        int mpheno = 1;
        int maternal_id_col = 3;
        std::string maternal_snp_test = "direct"; // direct | maternal | joint
        std::string pheno_file;
        std::string qcovar_file;
        std::string covar_file;
        std::string maternal_grm_file;
        std::string direct_grm_file;
        std::string extract_file;
        std::string exclude_file;

        static Options fromConfig(const GwasConfig& cfg);
    };

    MaternalModel(const Options& opts);

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
    Eigen::MatrixXf X_;
    DataFilter::KeepInfo keep_;

    std::vector<int> record_to_animal_;
    std::vector<int> record_to_dam_;
    std::unordered_map<std::string, int> iid_to_subject_idx_;

    Eigen::MatrixXd G_direct_;
    Eigen::MatrixXd G_maternal_;
    Eigen::MatrixXd G_dm_;

    int n_dams_ = 0;
    int n_animals_ = 0;

    double sigma2_ad_ = 0.0;
    double sigma2_am_ = 0.0;
    double sigma_adm_ = 0.0;
    double sigma2_c_ = 0.0;
    double sigma2_e_ = 0.0;

    void loadPhenotype();
    void loadCovariates();
    void buildGRMs();
    void loadDamMap();
    void runVCE();
    void scanGwas();
    // Model-specific maternal SNP test scan (maternal or joint design).
    // Builds the eligible analysis subset (animals whose dam has a genotype row),
    // fits the shared null covariance once on that subset, then scans with a
    // 1-column (maternal) or 2-column (direct+maternal, joint 2-df Wald) design.
    void scanMaternalSnp();
    void runLocoScan(const Eigen::VectorXd& y_agg,
                     const Eigen::MatrixXf& X_agg,
                     double vg_eff, double ve_eff,
                     int n_animal,
                     const std::string& out_file);
    void writeNullSummary() const;
};

} // namespace cosmic
