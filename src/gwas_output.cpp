#include "gwas_scanner.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace cosmic {

void GwasScanner::writeHeader(const std::string& filename) {
    std::ofstream ofs(filename);
    eb_effective_tests_ = 0.0;
    eb_tested_markers_ = 0;
    if (options.run_skat || options.run_burden) {
        ofs << "SET_ID\tN_SNPS\tSTAT\tP\n";
    } else if (options.eb_lmm) {
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tEB_BETA\tEB_SE\tEB_WALD\tP\tLOG10P"
            << "\tMARKER_LAMBDA\tPRIOR_VAR\tPOST_VAR\tEB_CONFIDENCE"
            << "\tEB_LRT\tEB_LRT_P\tEB_SIGMA2\tEB_CONV\n";
    } else if (options.lmm_gemma_mode > 0) {
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP";
        if (options.lmm_gemma_mode == 4) ofs << "\tP_WALD\tP_LRT\tP_SCORE";
        ofs << "\tWALD\tLRT\tSCORE\n";
    } else if (options.out_ssf) {
        ofs << "chromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tsample_size\tbeta\tstandard_error\tp_value\tvariant_id\n";
    } else {
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tLOG10P\n";
    }
    ofs.close();
}

void GwasScanner::writeEbSummary(const std::string& filename) {
    if (!options.eb_lmm) return;
    std::ofstream ofs(filename);
    if (!ofs) throw std::runtime_error("Cannot write EB summary: " + filename);
    ofs << std::setprecision(12);
    const double m_eff = std::max(eb_effective_tests_, 1.0);
    ofs << "Metric\tValue\n";
    ofs << "method\tEB-LMM random-SNP scan\n";
    ofs << "tested_markers\t" << eb_tested_markers_ << "\n";
    ofs << "effective_tests\t" << eb_effective_tests_ << "\n";
    ofs << "effective_tests_per_marker\t" << (eb_tested_markers_ > 0 ? eb_effective_tests_ / (double)eb_tested_markers_ : 0.0) << "\n";
    ofs << "bonferroni_0.05\t" << (eb_tested_markers_ > 0 ? 0.05 / (double)eb_tested_markers_ : 0.0) << "\n";
    ofs << "eb_mbonf_0.05\t" << (0.05 / m_eff) << "\n";
    ofs << "note\tValidated EB-LMM auxiliary scan; report standard Bonferroni/FDR alongside EB-mBonf and keep standard LMM results in formal analyses.\n";
}

void GwasScanner::appendResults(const std::string& filename, const std::vector<GwasResult>& results) {
    std::ofstream ofs(filename, std::ios::app);
    for (const auto& res : results) {
        if (std::isnan(res.p_value)) continue;
        if (options.run_skat || options.run_burden) {
            ofs << res.snp_id << "\t" << res.n_obs << "\t" << res.beta << "\t" << res.p_value << "\n";
        } else if (options.out_ssf) {
            ofs << res.chrom << "\t" << res.pos << "\t"
                << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
                << res.p_value << "\t" << res.snp_id << "\n";
        } else if (options.eb_lmm) {
            eb_effective_tests_ += std::max(0.0, std::min(1.0, res.eb_confidence));
            eb_tested_markers_++;
            ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
                << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
                << res.eb_wald << "\t" << res.p_value << "\t" << res.log10_p << "\t"
                << res.eb_marker_lambda << "\t" << res.eb_prior_var << "\t"
                << res.eb_post_var << "\t" << res.eb_confidence << "\t"
                << res.eb_lrt << "\t" << res.eb_lrt_p << "\t"
                << res.eb_sigma2 << "\t" << res.eb_conv << "\n";
        } else if (options.lmm_gemma_mode > 0) {
            double p = res.p_wald;
            if (options.lmm_gemma_mode == 2) p = res.p_lrt;
            else if (options.lmm_gemma_mode == 3) p = res.p_score;
            ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
                << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t" << p;
            if (options.lmm_gemma_mode == 4) {
                ofs << "\t" << res.p_wald << "\t" << res.p_lrt << "\t" << res.p_score;
            }
            ofs << "\t" << res.wald_stat << "\t" << res.lrt_stat << "\t" << res.score_stat << "\n";
        } else {
            ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
                << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
                << res.p_value << "\t" << res.log10_p << "\n";
        }
    }
}

} // namespace cosmic
