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
    void validateInputs() const;
    std::string writeThresholdPhenotype(double threshold, int threshold_index) const;
    void writeThresholdSummary() const;
    void writeCombinedConvergence(const std::vector<std::string>& threshold_prefixes) const;
    void writeCombinedGwas(const std::vector<std::string>& threshold_prefixes, double runtime_sec) const;
};

} // namespace cosmic
