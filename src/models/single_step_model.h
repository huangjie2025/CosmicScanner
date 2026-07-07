#pragma once
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "data_filter.h"

namespace cosmic {

class SingleStepModel {
public:
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
        std::string extract_file;
        std::string exclude_file;

        static Options fromConfig(const GwasConfig& cfg);
    };

    SingleStepModel(const Options& opts);

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

    // Pedigree data
    struct PedRecord {
        std::string id;
        std::string sire;
        std::string dam;
    };
    std::vector<PedRecord> pedigree_;
    std::unordered_map<std::string, int> ped_id_to_idx_;

    // A⁻¹ in sparse format
    Eigen::SparseMatrix<double> A_inv_;

    // GRM for genotyped individuals
    Eigen::MatrixXd G_;
    Eigen::MatrixXd A22_inv_;
    int n_genotyped_ = 0;
    std::vector<int> genotyped_indices_;
    bool has_a22_correction_ = false;
    double h_inverse_symmetry_error_ = 0.0;

    // H⁻¹ = A⁻¹ + correction
    Eigen::SparseMatrix<double> H_inv_;

    double sigma2_g_ = 0.0;
    double sigma2_e_ = 0.0;

    void loadPhenotype();
    void loadCovariates();
    void loadPedigree();
    void buildAInverse();
    void buildGRM();
    void computeA22();
    void computeHInverse();
    void runVCE();
    void scanGwas();
    void writeEbv();
    void writeNullSummary() const;
};

} // namespace cosmic
