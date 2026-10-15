#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"
#include "pheno_reader.h"
#include "data_filter.h"

namespace cosmic {

class RrmModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_rrm";
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int grm_alg = 1;
        int n_threads = 1;
        bool no_loco = false;
        bool fast_lmm = false;
        bool use_spa = false;
        int rrm_order = 1;
        int rrm_curve_grid = 21;
        int rrm_top = 100;
        int rrm_ai_iter = -1;
        double rrm_sig_threshold = 0.05;
        bool rrm_unbalance = false;
        bool rrm_loco_global_null = false;
        std::string rrm_basis = "legendre";
        std::string rrm_gwas_method = "coefficient";
        int tpheno = 1;
        int mpheno = 1;
        std::string long_pheno_file;
        std::string pheno_file;
        std::string bfile_prefix;
        std::string rrm_compatible;
        std::string rrm_compatible_varcom;
        std::vector<int> mpheno_list;
        std::string qcovar_file;
        std::string covar_file;
        std::string tpoint_col_name;
        std::string tclass_names;
        std::string tcovar_names;

        static Options fromConfig(const GwasConfig& cfg);
    };

    RrmModel(const Options& opts);

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
    Eigen::MatrixXf X_base_;
    Eigen::MatrixXd Y_wide_;
    std::vector<PhenoReader::PhenoRecord> records_;
    DataFilter::KeepInfo keep_;

    void runGmaFixregCompatible();
    void runUnbalanced();
};

} // namespace cosmic
