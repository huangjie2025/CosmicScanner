#pragma once

#include <string>
#include <vector>
#include <Eigen/Dense>

#include "bgen_reader.h"
#include "config.h"
#include "pgen_reader.h"
#include "plink_reader.h"
#include "vcf/vcf_reader.h"

namespace cosmic {

class GemmaCalccor {
public:
    explicit GemmaCalccor(const GwasConfig& cfg);

    void run(PlinkReader& plink_reader,
             BgenReader& bgen_reader,
             PgenReader& pgen_reader,
             VcfReader& vcf_reader,
             bool is_bgen,
             bool is_pgen,
             bool is_vcf);

private:
    struct VariantRecord {
        std::string chr;
        std::string rs;
        long long bp = -9;
        double cm = -9.0;
        std::string allele1;
        std::string allele0;
        int n_miss = 0;
        int n_obs = 0;
        double af = 0.0;
        double var = 0.0;
        Eigen::VectorXd centered;
    };

    const GwasConfig& cfg_;

    bool passesIdFilters(const std::string& rs) const;
    bool passesRegionFilters(const std::string& chr, long long bp) const;
    bool normalizeVariant(const std::string& chr,
                          const std::string& rs,
                          long long bp,
                          double cm,
                          const std::string& allele1,
                          const std::string& allele0,
                          const Eigen::VectorXd& dosage,
                          VariantRecord& rec) const;
std::vector<VariantRecord> loadVariants(PlinkReader& plink_reader,
                                             BgenReader& bgen_reader,
                                             PgenReader& pgen_reader,
                                             VcfReader& vcf_reader,
                                             bool is_bgen,
                                             bool is_pgen,
                                             bool is_vcf) const;
    void writeCor(const std::vector<VariantRecord>& variants) const;
};

} // namespace cosmic
