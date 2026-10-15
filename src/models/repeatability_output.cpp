#include "repeatability_helpers.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace cosmic {

void writeAssocHeader(const std::string& filename, bool out_ssf, bool eb_lmm) {
    std::ofstream ofs(filename);
    if (eb_lmm && !out_ssf) {
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tEB_BETA\tEB_SE\tEB_WALD\tP\tLOG10P"
            << "\tMARKER_LAMBDA\tPRIOR_VAR\tPOST_VAR\tEB_CONFIDENCE"
            << "\tEB_LRT\tEB_LRT_P\tEB_SIGMA2\tEB_CONV\n";
    } else if (out_ssf) {
        ofs << "chromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tsample_size\tbeta\tstandard_error\tp_value\tvariant_id\n";
    } else {
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tLOG10P\n";
    }
}

void appendAssocResult(const std::string& filename, const GwasResult& res, bool out_ssf, bool eb_lmm) {
    std::ofstream ofs(filename, std::ios::app);
    if (std::isnan(res.p_value)) return;
    if (eb_lmm && !out_ssf) {
        ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
            << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
            << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
            << res.eb_wald << "\t" << res.p_value << "\t" << res.log10_p << "\t"
            << res.eb_marker_lambda << "\t" << res.eb_prior_var << "\t"
            << res.eb_post_var << "\t" << res.eb_confidence << "\t"
            << res.eb_lrt << "\t" << res.eb_lrt_p << "\t"
            << res.eb_sigma2 << "\t" << res.eb_conv << "\n";
    } else if (out_ssf) {
        ofs << res.chrom << "\t" << res.pos << "\t"
            << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
            << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
            << res.p_value << "\t" << res.snp_id << "\n";
    } else {
        ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
            << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
            << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
            << res.p_value << "\t" << res.log10_p << "\n";
    }
}

void appendAssocResults(const std::string& filename, const std::vector<GwasResult>& results, bool out_ssf, bool eb_lmm) {
    for (const auto& res : results) appendAssocResult(filename, res, out_ssf, eb_lmm);
}

void accumulateEbResults(const std::vector<GwasResult>& results, double& effective_tests, long long& tested_markers) {
    for (const auto& res : results) {
        if (std::isnan(res.p_value)) continue;
        effective_tests += std::max(0.0, std::min(1.0, res.eb_confidence));
        tested_markers++;
    }
}

void writeRepeatabilityEbSummary(const std::string& filename, double effective_tests, long long tested_markers) {
    std::ofstream ofs(filename);
    if (!ofs) throw std::runtime_error("Cannot write repeatability EB summary: " + filename);
    ofs << std::setprecision(12);
    double m_eff = std::max(effective_tests, 1.0);
    ofs << "Metric\tValue\n";
    ofs << "method\tRepeatability EB-LMM random-SNP scan (fast aggregated subject-level)\n";
    ofs << "tested_markers\t" << tested_markers << "\n";
    ofs << "effective_tests\t" << effective_tests << "\n";
    ofs << "effective_tests_per_marker\t" << (tested_markers > 0 ? effective_tests / (double)tested_markers : 0.0) << "\n";
    ofs << "bonferroni_0.05\t" << (tested_markers > 0 ? 0.05 / (double)tested_markers : 0.0) << "\n";
    ofs << "eb_mbonf_0.05\t" << (0.05 / m_eff) << "\n";
    ofs << "note\tExperimental repeatability EB uses aggregated weighted subject means; exact record-level EB is not used in this mode.\n";
}

} // namespace cosmic
