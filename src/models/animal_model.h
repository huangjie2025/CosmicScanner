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
#include "pedigree_relationship.h"

namespace cosmic {

/// Pedigree-based Simple Animal Model GWAS:
///
///   y = Xb + x_j beta_j + Za + e
///   a ~ N(0, A sigma_a^2),  e ~ N(0, I sigma_e^2)
///
/// A is the full pedigree numerator relationship matrix; the analysis
/// (genotyped with phenotype) individuals use the A22 block. VCE via EigenVce,
/// per-SNP scan via LmmEngine with the same output format as standard GWAS.
class AnimalModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_animal_model";
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int reml_maxit = 100;
        double reml_tol = 1e-6;
        int mpheno = 1;
        std::string pheno_file;
        std::string qcovar_file;
        std::string covar_file;
        std::string ped_file;
        std::string extract_file;
        std::string exclude_file;
        bool write_ped_kernel = false;
        int n_threads = 1;
        bool no_loco = false;

        static Options fromConfig(const GwasConfig& cfg);
    };

    AnimalModel(const Options& opts);

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
    void loadPhenotype();
    void loadCovariates();
    void buildPedigreeRelationship();
    void runVCE();
    void scanGwas();
    void writeKernel() const;
    void writeNullSummary() const;

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

    PedigreeRelationship pedigree_;
    std::vector<int> analysis_pedigree_indices_;
    Eigen::MatrixXd A22_;

    double sigma2_a_ = 0.0;
    double sigma2_e_ = 0.0;
    bool vce_converged_ = false;
    int vce_iterations_ = 0;

    double t_ped_parse_ = 0.0;
    double t_ainv_ = 0.0;
    double t_a22_ = 0.0;
    double t_vce_ = 0.0;
    double t_scan_ = 0.0;
};

} // namespace cosmic
