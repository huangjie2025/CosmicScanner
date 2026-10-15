#include "gwas_scanner.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <chrono>
#include <sstream>
#include "validation.h"
#include "association/sumstats/sumstats.h"
#include "gzip_writer.h"
#include "gwas_output_schema.h"

namespace cosmic {

namespace {

// Helper to create output stream (plain or gzip)
std::unique_ptr<std::ostream> makeOutputStream(const std::string& filename, bool use_gzip, bool append = false) {
    if (use_gzip) {
        return std::make_unique<GzipOutputStream>(filename, Z_DEFAULT_COMPRESSION, append);
    } else {
        auto stream = std::make_unique<std::ofstream>(filename, std::ios::out |
                                                     (append ? std::ios::app : std::ios::trunc));
        if (!*stream) throw std::runtime_error("Cannot open association output: " + filename);
        return stream;
    }
}

void finishOutputStream(std::ostream& stream, const std::string& filename) {
    stream.flush();
    if (!stream) throw std::runtime_error("Failed while writing association output: " + filename);
    if (auto* gzip = dynamic_cast<GzipOutputStream*>(&stream)) gzip->close();
    else if (auto* file = dynamic_cast<std::ofstream*>(&stream)) file->close();
    if (!stream) throw std::runtime_error("Failed closing association output: " + filename);
}

ValidationReport validate_association_results(const std::vector<GwasResult>& results,
                                              const GwasScanner::Options& options) {
    ValidationReport report;
    if (options.run_skat || options.run_burden) return report;
    std::size_t tested = 0;
    for (const auto& res : results) {
        if (std::isnan(res.p_value)) continue;
        ++tested;
        const double p = options.lmm_gemma_mode == 2 ? res.p_lrt
                       : options.lmm_gemma_mode == 3 ? res.p_score
                                                     : res.p_value;
        validate_probability_value(report, p, "ASSOC_INVALID_PVALUE", "CosmicScanner association output");
        validate_finite_scalar(report, res.beta, "ASSOC_NONFINITE_BETA", "CosmicScanner association output");
        validate_nonnegative_finite_scalar(report, res.se, "ASSOC_INVALID_SE", "CosmicScanner association output");
        if (report.has_fail()) break;
    }
    if (!report.has_fail()) {
        report.pass(ValidationStage::PostRun, "ASSOC_RESULTS_VALID",
                    "CosmicScanner association output",
                    std::to_string(tested) + " tested marker results passed numerical sanity checks");
    }
    std::cout << report.summary(false);
    report.throw_if_fail();
    return report;
}

} // namespace

void GwasScanner::writeHeader(const std::string& filename) {
    auto ofs = makeOutputStream(filename, options.out_gz);
    eb_effective_tests_ = 0.0;
    eb_tested_markers_ = 0;
    if (options.run_skat || options.run_burden) {
        *ofs << "SET_ID\tN_SNPS\tSTAT\tP\n";
    } else if (options.eb_lmm) {
        *ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tEB_BETA\tEB_SE\tEB_WALD\tP\tLOG10P"
            << "\tMARKER_LAMBDA\tPRIOR_VAR\tPOST_VAR\tEB_CONFIDENCE"
            << "\tEB_LRT\tEB_LRT_P\tEB_SIGMA2\tEB_CONV\n";
    } else if (options.lmm_gemma_mode > 0) {
        *ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP";
        if (options.lmm_gemma_mode == 4) *ofs << "\tP_WALD\tP_LRT\tP_SCORE";
        *ofs << "\tWALD\tLRT\tSCORE\n";
    } else if (options.out_ssf) {
        *ofs << "chromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tsample_size\tbeta\tstandard_error\tp_value\tvariant_id\n";
    } else {
        *ofs << gwas_output::kUnifiedHeader;
    }
    finishOutputStream(*ofs, filename);
}

void GwasScanner::writeEbSummary(const std::string& filename) {
    if (!options.eb_lmm) return;
    auto ofs = makeOutputStream(filename, options.out_gz);
    *ofs << std::setprecision(12);
    const double m_eff = std::max(eb_effective_tests_, 1.0);
    *ofs << "Metric\tValue\n";
    *ofs << "method\tEB-LMM random-SNP scan\n";
    *ofs << "tested_markers\t" << eb_tested_markers_ << "\n";
    *ofs << "effective_tests\t" << eb_effective_tests_ << "\n";
    *ofs << "effective_tests_per_marker\t" << (eb_tested_markers_ > 0 ? eb_effective_tests_ / (double)eb_tested_markers_ : 0.0) << "\n";
    *ofs << "bonferroni_0.05\t" << (eb_tested_markers_ > 0 ? 0.05 / (double)eb_tested_markers_ : 0.0) << "\n";
    *ofs << "eb_mbonf_0.05\t" << (0.05 / m_eff) << "\n";
    *ofs << "note\tValidated EB-LMM auxiliary scan; report standard Bonferroni/FDR alongside EB-mBonf and keep standard LMM results in formal analyses.\n";
    finishOutputStream(*ofs, filename);
}

void GwasScanner::appendResults(const std::string& filename, const std::vector<GwasResult>& results) {
    validation_report_.merge(validate_association_results(results, options));
    // Accumulate EB counters during scan
    if (options.eb_lmm) {
        for (const auto& res : results) {
            if (std::isnan(res.p_value)) continue;
            eb_effective_tests_ += std::max(0.0, std::min(1.0, res.eb_confidence));
            eb_tested_markers_++;
        }
    }
    // Accumulate for single-pass output
    assoc_results_.insert(assoc_results_.end(), results.begin(), results.end());
    // Accumulate for Manhattan/QQ plot TSVs
    if (options.manhattan_tsv || options.qq_tsv) {
        plot_results_.insert(plot_results_.end(), results.begin(), results.end());
    }
}

void GwasScanner::flushAssocResults(const std::string& filename) {
    if (assoc_results_.empty()) return;
    auto ofs = makeOutputStream(filename, options.out_gz);

    // Write header
    if (options.run_skat || options.run_burden) {
        *ofs << "SET_ID\tN_SNPS\tSTAT\tP\n";
    } else if (options.eb_lmm) {
        *ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tEB_BETA\tEB_SE\tEB_WALD\tP\tLOG10P"
            << "\tMARKER_LAMBDA\tPRIOR_VAR\tPOST_VAR\tEB_CONFIDENCE"
            << "\tEB_LRT\tEB_LRT_P\tEB_SIGMA2\tEB_CONV\n";
    } else if (options.lmm_gemma_mode > 0) {
        *ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP";
        if (options.lmm_gemma_mode == 4) *ofs << "\tP_WALD\tP_LRT\tP_SCORE";
        *ofs << "\tWALD\tLRT\tSCORE\n";
    } else if (options.out_ssf) {
        *ofs << "chromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tsample_size\tbeta\tstandard_error\tp_value\tvariant_id\n";
    } else {
        *ofs << gwas_output::kUnifiedHeader;
    }

    // Write data
    for (const auto& res : assoc_results_) {
        if (std::isnan(res.p_value)) continue;
        if (options.run_skat || options.run_burden) {
            *ofs << res.snp_id << "\t" << res.n_obs << "\t" << res.beta << "\t" << res.p_value << "\n";
        } else if (options.out_ssf) {
            *ofs << res.chrom << "\t" << res.pos << "\t"
                << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
                << res.p_value << "\t" << res.snp_id << "\n";
        } else if (options.eb_lmm) {
            *ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
                << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
                << res.eb_wald << "\t" << res.p_value << "\t" << res.log10_p << "\t"
                << std::setprecision(6) << res.eb_marker_lambda << "\t"
                << std::setprecision(6) << res.eb_prior_var << "\t"
                << std::setprecision(6) << res.eb_post_var << "\t"
                << std::setprecision(12) << res.eb_confidence << "\t"
                << std::setprecision(6) << res.eb_lrt << "\t"
                << std::setprecision(6) << res.eb_lrt_p << "\t"
                << std::setprecision(6) << res.eb_sigma2 << "\t"
                << res.eb_conv << "\n";
        } else if (options.lmm_gemma_mode > 0) {
            double p = res.p_wald;
            if (options.lmm_gemma_mode == 2) p = res.p_lrt;
            else if (options.lmm_gemma_mode == 3) p = res.p_score;
            *ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
                << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t" << p;
            if (options.lmm_gemma_mode == 4) {
                *ofs << "\t" << res.p_wald << "\t" << res.p_lrt << "\t" << res.p_score;
            }
            *ofs << "\t" << res.wald_stat << "\t" << res.lrt_stat << "\t" << res.score_stat << "\n";
        } else {
            gwas_output::writeAdditiveRow(*ofs, res);
        }
    }
    finishOutputStream(*ofs, filename);
    assoc_results_.clear();
}

// Append canonical CosmicScanner summary statistics (.cosmic.sumstats)
void GwasScanner::appendCosmicSumstats(const std::string& filename,
                                        const std::vector<GwasResult>& results,
                                        const std::string& model_name) {
    using namespace cosmic::association::sumstats;

    // First chromosome replaces a previous run; subsequent chromosomes append.
    // gzip append writes concatenated members, supported by gzread/gzip readers.
    auto ofs = makeOutputStream(filename, options.out_gz, sumstats_header_written_);

    if (!sumstats_header_written_) {
        SumstatsMetadata meta;
        meta.model = model_name;
        meta.n_samples = 0;
        meta.n_variants = 0;

        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::gmtime(&in_time_t), "%Y-%m-%d %H:%M:%S UTC");
        meta.date = ss.str();

        *ofs << "# CosmicScanner Summary Statistics v" << meta.version << "\n";
        *ofs << "# Created: " << meta.date << "\n";
        if (!meta.model.empty()) *ofs << "# Model: " << meta.model << "\n";
        *ofs << "#\n";

        *ofs << "CHR\tSNP\tPOS\tA1\tA2\tBETA\tSE\tP\tN\tAF\tLOG10P\tMODEL";
        *ofs << "\tZ\tWEIGHT\tINFO\tWARNING\n";
        sumstats_header_written_ = true;
    }

    for (const auto& res : results) {
        if (std::isnan(res.p_value)) continue;

        *ofs << res.chrom << '\t'
            << res.snp_id << '\t'
            << res.pos << '\t'
            << res.a1 << '\t'
            << res.a2 << '\t'
            << std::setprecision(10) << res.beta << '\t'
            << std::setprecision(10) << res.se << '\t'
            << std::setprecision(10) << res.p_value << '\t'
            << res.n_obs << '\t'
            << std::setprecision(6) << res.freq << '\t'
            << std::setprecision(6) << res.log10_p << '\t'
            << model_name << '\t'
            << std::setprecision(6) << ((res.se > 0) ? (res.beta / res.se) : 0.0) << '\t'
            << std::setprecision(6) << ((res.se > 0) ? (1.0 / (res.se * res.se)) : 0.0) << '\t'
            << std::setprecision(4) << 1.0 << '\t'
            << ".\n";
    }
    finishOutputStream(*ofs, filename);
}

void GwasScanner::writePlotTsv(const std::string& prefix,
                                const std::vector<GwasResult>& results) {
    struct RankItem { std::string chrom; std::string snp; long long bp; double p; double log10p; };
    std::vector<RankItem> items;
    items.reserve(results.size());
    for (const auto& r : results) {
        if (std::isnan(r.p_value) || r.p_value <= 0.0) continue;
        RankItem ri;
        ri.chrom = r.chrom;
        ri.snp = r.snp_id;
        ri.bp = r.pos;
        ri.p = r.p_value;
        ri.log10p = r.log10_p;
        items.push_back(ri);
    }
    if (items.empty()) return;

    // Sort by p-value ascending for QQ
    std::vector<RankItem> sorted_by_p = items;
    std::sort(sorted_by_p.begin(), sorted_by_p.end(),
              [](const RankItem& a, const RankItem& b) { return a.p < b.p; });

    const int m = static_cast<int>(sorted_by_p.size());

    // Manhattan: write cumulative-BP TSV (CHR, SNP, BP, P, LOG10P, CUMULATIVE_BP)
    if (options.manhattan_tsv) {
        // Compute cumulative offsets per chromosome
        std::map<std::string, long long> chrom_max_bp;
        for (const auto& it : items) {
            auto it2 = chrom_max_bp.find(it.chrom);
            if (it2 == chrom_max_bp.end() || it.bp > it2->second)
                chrom_max_bp[it.chrom] = it.bp;
        }
        // Sort chromosomes
        std::vector<std::string> chroms;
        for (auto& kv : chrom_max_bp) chroms.push_back(kv.first);
        std::sort(chroms.begin(), chroms.end());
        std::map<std::string, long long> chrom_offset;
        long long cum = 0;
        for (const auto& c : chroms) {
            chrom_offset[c] = cum;
            cum += chrom_max_bp[c] + 1;
        }

        std::string manhattan_file = prefix + ".manhattan.tsv";
        auto ofs = makeOutputStream(manhattan_file, options.out_gz);
        *ofs << "CHR\tSNP\tBP\tP\tLOG10P\tCUMULATIVE_BP\n";
        *ofs << std::setprecision(15);
        for (const auto& it : items) {
            long long cum_bp = chrom_offset[it.chrom] + it.bp;
            *ofs << it.chrom << "\t" << it.snp << "\t" << it.bp << "\t"
                 << it.p << "\t" << it.log10p << "\t" << cum_bp << "\n";
        }
        finishOutputStream(*ofs, manhattan_file);
    }

    // QQ: write expected vs observed -log10(p)
    if (options.qq_tsv) {
        std::string qq_file = prefix + ".qq.tsv";
        auto ofs = makeOutputStream(qq_file, options.out_gz);
        *ofs << "EXPECTED_LOG10P\tOBSERVED_LOG10P\tP\tSNP\tCHR\tBP\n";
        *ofs << std::setprecision(15);
        for (int i = 0; i < m; ++i) {
            // Expected: uniform order statistic median approximation
            double expected_p = (i + 0.5) / m;
            double expected_log10p = -std::log10(expected_p);
            *ofs << expected_log10p << "\t" << sorted_by_p[i].log10p << "\t"
                 << sorted_by_p[i].p << "\t" << sorted_by_p[i].snp << "\t"
                 << sorted_by_p[i].chrom << "\t" << sorted_by_p[i].bp << "\n";
        }
        finishOutputStream(*ofs, qq_file);
    }
}

} // namespace cosmic
