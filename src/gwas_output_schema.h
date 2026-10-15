#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cosmic::gwas_output {

inline constexpr const char* kUnifiedHeader =
    "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN"
    "\tBETA_A\tSE_A\tP_A\tLOG10P_A"
    "\tBETA_D\tSE_D\tP_D\tLOG10P_D"
    "\tP_JOINT\tLOG10P_JOINT\n";

inline double negativeLog10(double p) {
    return -std::log10(std::max(p, std::numeric_limits<double>::min()));
}

template <typename Result>
inline void writePrefix(std::ostream& output, const Result& result) {
    output << result.chrom << '\t' << result.snp_id << '\t' << result.pos << '\t'
           << result.a1 << '\t' << result.a2 << '\t' << result.freq << '\t'
           << result.n_obs;
}

template <typename Result>
inline void writeAdditiveRow(std::ostream& output, const Result& result) {
    writePrefix(output, result);
    output << '\t' << result.beta << '\t' << result.se << '\t' << result.p_value
           << '\t' << result.log10_p
           << "\tNA\tNA\tNA\tNA\tNA\tNA\n";
}

template <typename Result>
inline void writeDominanceRow(std::ostream& output, const Result& result) {
    writePrefix(output, result);
    output << "\tNA\tNA\tNA\tNA"
           << '\t' << result.beta_dom << '\t' << result.se_dom << '\t' << result.p_dom
           << '\t' << result.log10_p_dom << "\tNA\tNA\n";
}

template <typename Result>
inline void writeAdRow(std::ostream& output, const Result& result) {
    writePrefix(output, result);
    output << '\t' << result.beta_add << '\t' << result.se_add << '\t' << result.p_add
           << '\t' << negativeLog10(result.p_add)
           << '\t' << result.beta_dom << '\t' << result.se_dom << '\t' << result.p_dom
           << '\t' << result.log10_p_dom
           << '\t' << result.p_joint << '\t' << result.log10_p_joint << '\n';
}

inline void appendCommonMetadata(const std::string& filename,
                                 const std::string& model,
                                 int n,
                                 bool loco,
                                 bool covariate_projection,
                                 bool vce_converged) {
    std::ofstream metadata(filename, std::ios::app);
    if (!metadata) throw std::runtime_error("Cannot append GWAS scan metadata: " + filename);
    metadata << "MODEL\t" << model << "\n"
             << "N\t" << n << "\n"
             << "N_SEMANTICS\tanalysis_samples_after_phenotype_and_covariate_filtering;"
                "per-SNP_missing_genotypes_mean-imputed\n"
             << "FREQ_SEMANTICS\tA1_allele_frequency_(not_MAF)\n"
             << "LOCO\t" << (loco ? "true" : "false") << "\n"
             << "COVARIATE_PROJECTION\t" << (covariate_projection ? "true" : "false") << "\n"
             << "VCE_CONVERGED\t" << (vce_converged ? "true" : "false") << "\n";
}

inline void appendVarianceMetadata(const std::string& filename,
                                   const std::string& name,
                                   const std::vector<std::string>& labels,
                                   const std::vector<double>& values) {
    if (labels.size() != values.size() || values.empty()) return;
    std::ofstream metadata(filename, std::ios::app);
    if (!metadata) throw std::runtime_error("Cannot append variance metadata: " + filename);
    if (values.size() == 1) {
        metadata << name << '\t' << values.front() << '\n';
    } else {
        metadata << name << "\tPER_NULL_MODEL\n";
        for (size_t i = 0; i < values.size(); ++i) {
            metadata << name << '_' << labels[i] << '\t' << values[i] << '\n';
        }
    }
}

} // namespace cosmic::gwas_output
