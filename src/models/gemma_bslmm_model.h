#pragma once

#include <string>
#include <vector>
#include <Eigen/Dense>

#include "config.h"
#include "data_filter.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"

namespace cosmic {

class GemmaBslmmModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_gwas";
        std::string pheno_file;
        std::string qcovar_file;
        std::string covar_file;
        int mpheno = 1;
        int mode = 2;
        int top_sparse = 100;
        double h2 = 0.5;
        double min_maf = 0.01;
        std::string extract_file;
        std::string exclude_file;
        std::string bf_file;
        std::string hyp_file;
        int h_grid = 5;
        int rho_grid = 5;
        double h_min = 0.05;
        double h_max = 0.95;
        double rho_min = 0.05;
        double rho_max = 0.95;
        double p_min = -4.0; // log10(pi)
        double p_max = 0.0;  // log10(pi)
        bool p_min_set = false;
        bool p_max_set = false;
        double gamma_mean = 0.0;
        double h_scale = 1.0;
        double rho_scale = 1.0;
        double p_scale = 1.0;
        int mh_steps = 0;
        int s_min = 0;
        int rpace = 10;
        int wpace = 1000;
        int burnin = 100000;
        int sampling = 1000000;
        long long seed = 0;

        static Options fromConfig(const GwasConfig& cfg);
    };

    explicit GemmaBslmmModel(const Options& options);

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
    struct VariantMeta {
        std::string chrom;
        std::string id;
        long long pos = 0;
        std::string a1;
        std::string a2;
        double freq = 0.0;
        double mean = 0.0;
        double sd = 1.0;
    };

    Options options_;
    int n_samples_ = 0;
    std::vector<std::string> fids_;
    std::vector<std::string> iids_;
    bool is_bgen_ = false;
    bool is_pgen_ = false;
    PlinkReader* plink_reader_ = nullptr;
    BgenReader* bgen_reader_ = nullptr;
    PgenReader* pgen_reader_ = nullptr;
    Eigen::VectorXd y_;
    Eigen::MatrixXf X_;
    DataFilter::KeepInfo keep_;

    void loadDesign(Eigen::MatrixXd& Z, std::vector<VariantMeta>& variants) const;
    void runDapFineMap() const;
};

} // namespace cosmic
