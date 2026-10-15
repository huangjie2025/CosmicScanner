#include "association/fine_mapping/fine_mapping.h"

#include <Eigen/Cholesky>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace cosmic::association::fine_mapping {

namespace {

double twoSidedNormalP(double z) {
    return std::erfc(std::abs(z) / std::sqrt(2.0));
}

sumstats::SumstatsRecord toSumstats(const cojo::SummaryRecord& source) {
    sumstats::SumstatsRecord result;
    result.chrom = source.chrom;
    result.snp_id = source.snp_id;
    result.pos = source.position;
    result.a1 = source.a1;
    result.a2 = source.a2;
    result.beta = source.beta;
    result.se = source.se;
    result.p_value = source.p_value;
    result.n = source.n;
    result.freq = source.freq;
    result.z_score = source.beta / source.se;
    result.weight = 1.0 / (source.se * source.se);
    return result;
}

} // namespace

FineMappingEngine::FineMappingEngine(FineMappingOptions options) : options_(options) {
    if (!(options_.prior_variance > 0.0) || !std::isfinite(options_.prior_variance))
        throw std::invalid_argument("Fine-mapping prior variance must be finite and positive");
    if (!(options_.credible_level > 0.0 && options_.credible_level < 1.0))
        throw std::invalid_argument("Fine-mapping credible level must be in (0, 1)");
    if (options_.max_signals < 1)
        throw std::invalid_argument("Fine-mapping maximum signals must be positive");
    if (!(options_.signal_p_threshold > 0.0 && options_.signal_p_threshold < 1.0))
        throw std::invalid_argument("Fine-mapping signal threshold must be in (0, 1)");
}

double FineMappingEngine::logApproximateBayesFactor(double beta, double se, double prior_variance) {
    if (!std::isfinite(beta) || !std::isfinite(se) || !(se > 0.0) ||
        !std::isfinite(prior_variance) || !(prior_variance > 0.0))
        throw std::invalid_argument("Invalid Wakefield ABF input");
    const long double variance = static_cast<long double>(se) * se;
    const long double prior = prior_variance;
    const long double shrinkage = prior / (variance + prior);
    const long double z = static_cast<long double>(beta) / se;
    const long double value = 0.5L * (std::log(variance) - std::log(variance + prior)) +
                              0.5L * z * z * shrinkage;
    return static_cast<double>(std::max(-700.0L, std::min(700.0L, value)));
}

FineMappingResult FineMappingEngine::run(
    const std::vector<cojo::SummaryRecord>& summary_stats,
    const cojo::LDReference& ld_reference) const {
    const auto started = std::chrono::steady_clock::now();
    FineMappingResult result;
    if (ld_reference.ld_matrix.rows() != ld_reference.ld_matrix.cols() ||
        ld_reference.ld_matrix.rows() != static_cast<Eigen::Index>(ld_reference.snp_ids.size()))
        throw std::invalid_argument("Fine-mapping LD reference dimensions are inconsistent");

    std::vector<sumstats::SumstatsRecord> variants;
    std::vector<int> ld_indices;
    std::unordered_set<std::string> seen;
    variants.reserve(summary_stats.size());
    ld_indices.reserve(summary_stats.size());

    for (const auto& source : summary_stats) {
        const auto found = ld_reference.snp_to_idx.find(source.snp_id);
        if (found == ld_reference.snp_to_idx.end() || !seen.insert(source.snp_id).second) {
            ++result.rejected_variants;
            continue;
        }
        const int ld_idx = found->second;
        if (ld_idx < 0 || ld_idx >= static_cast<int>(ld_reference.snp_ids.size()) ||
            source.chrom != ld_reference.chroms[ld_idx] || source.position != ld_reference.positions[ld_idx]) {
            ++result.rejected_variants;
            continue;
        }
        auto harmonized = sumstats::harmonizeAlleles(
            toSumstats(source), ld_reference.a1s[ld_idx], ld_reference.a2s[ld_idx]);
        if (!harmonized.harmonized || !std::isfinite(harmonized.record.beta) ||
            !std::isfinite(harmonized.record.se) || !(harmonized.record.se > 0.0)) {
            ++result.rejected_variants;
            continue;
        }
        harmonized.record.z_score = harmonized.record.beta / harmonized.record.se;
        harmonized.record.weight = 1.0 / (harmonized.record.se * harmonized.record.se);
        variants.push_back(std::move(harmonized.record));
        ld_indices.push_back(ld_idx);
    }
    if (variants.empty()) throw std::runtime_error("Fine-mapping found no harmonized polymorphic variants");

    const int m = static_cast<int>(variants.size());
    Eigen::MatrixXd ld(m, m);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j)
            ld(i, j) = ld_reference.ld_matrix(ld_indices[i], ld_indices[j]);
    if (!ld.allFinite()) throw std::runtime_error("Fine-mapping LD matrix contains non-finite values");
    ld = 0.5 * (ld + ld.transpose());
    ld.diagonal().setOnes();

    Eigen::VectorXd marginal_z(m);
    for (int i = 0; i < m; ++i) marginal_z[i] = variants[i].beta / variants[i].se;
    std::vector<int> selected;

    auto append_signal = [&](int signal, const Eigen::VectorXd& beta, const Eigen::VectorXd& se,
                             const Eigen::VectorXd& log_abf, const std::vector<char>& eligible,
                             int lead) {
        double maximum = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < m; ++i) if (eligible[i]) maximum = std::max(maximum, log_abf[i]);
        Eigen::VectorXd pip = Eigen::VectorXd::Zero(m);
        double denominator = 0.0;
        for (int i = 0; i < m; ++i) {
            if (!eligible[i]) continue;
            pip[i] = std::exp(log_abf[i] - maximum);
            denominator += pip[i];
        }
        if (!(denominator > 0.0) || !std::isfinite(denominator))
            throw std::runtime_error("Fine-mapping posterior normalization failed");
        pip /= denominator;

        std::vector<int> order;
        for (int i = 0; i < m; ++i) if (eligible[i]) order.push_back(i);
        std::stable_sort(order.begin(), order.end(), [&](int left, int right) {
            if (pip[left] != pip[right]) return pip[left] > pip[right];
            return variants[left].snp_id < variants[right].snp_id;
        });
        Eigen::VectorXd cumulative = Eigen::VectorXd::Zero(m);
        std::vector<char> member(m, false);
        double running = 0.0;
        for (int idx : order) {
            running = std::min(1.0, running + pip[idx]);
            cumulative[idx] = running;
            if (signal > 0 && running - pip[idx] < options_.credible_level) member[idx] = true;
        }

        for (int i = 0; i < m; ++i) {
            FineMappingRow row;
            row.variant = variants[i];
            row.variant.beta = beta[i];
            row.variant.se = se[i];
            row.abf = eligible[i] ? std::exp(std::min(700.0, log_abf[i])) : 0.0;
            row.pip = pip[i];
            row.signal = signal;
            row.credible_set = member[i] ? signal : 0;
            row.cumulative_posterior = cumulative[i];
            result.rows.push_back(std::move(row));
        }
        for (int idx : order) {
            if (!member[idx]) continue;
            result.credible_sets.push_back({signal, signal, variants[idx].snp_id, pip[idx],
                                            cumulative[idx], variants[lead].snp_id, ld(idx, lead)});
        }
    };

    for (int signal = 1; signal <= options_.max_signals; ++signal) {
        Eigen::VectorXd conditional_z = marginal_z;
        Eigen::VectorXd conditional_beta(m);
        Eigen::VectorXd conditional_se(m);
        for (int j = 0; j < m; ++j) {
            conditional_beta[j] = variants[j].beta;
            conditional_se[j] = variants[j].se;
        }
        std::vector<char> eligible(m, true);
        for (int idx : selected) eligible[idx] = false;

        if (!selected.empty()) {
            const int k = static_cast<int>(selected.size());
            Eigen::MatrixXd selected_ld(k, k);
            Eigen::VectorXd selected_z(k);
            for (int i = 0; i < k; ++i) {
                selected_z[i] = marginal_z[selected[i]];
                for (int j = 0; j < k; ++j) selected_ld(i, j) = ld(selected[i], selected[j]);
            }
            Eigen::LDLT<Eigen::MatrixXd> solver;
            double ridge = 0.0;
            for (int attempt = 0; attempt < 5; ++attempt) {
                Eigen::MatrixXd regularized = selected_ld;
                regularized.diagonal().array() += ridge;
                solver.compute(regularized);
                if (solver.info() == Eigen::Success && solver.vectorD().array().abs().minCoeff() > 1e-10) break;
                ridge = ridge == 0.0 ? 1e-10 : ridge * 10.0;
            }
            if (solver.info() != Eigen::Success || solver.vectorD().array().abs().minCoeff() <= 1e-10) break;
            const Eigen::VectorXd solved_z = solver.solve(selected_z);
            for (int j = 0; j < m; ++j) {
                if (!eligible[j]) continue;
                Eigen::VectorXd correlations(k);
                for (int i = 0; i < k; ++i) correlations[i] = ld(j, selected[i]);
                const double residual_variance = 1.0 - correlations.dot(solver.solve(correlations));
                if (!(residual_variance > 1e-8) || !std::isfinite(residual_variance)) {
                    eligible[j] = false;
                    continue;
                }
                conditional_z[j] = (marginal_z[j] - correlations.dot(solved_z)) /
                                   std::sqrt(residual_variance);
                conditional_se[j] = variants[j].se / std::sqrt(residual_variance);
                conditional_beta[j] = conditional_z[j] * conditional_se[j];
            }
        }

        int lead = -1;
        double lead_abs_z = -1.0;
        Eigen::VectorXd log_abf = Eigen::VectorXd::Constant(m, -std::numeric_limits<double>::infinity());
        for (int j = 0; j < m; ++j) {
            if (!eligible[j] || !std::isfinite(conditional_z[j]) || !std::isfinite(conditional_se[j])) continue;
            log_abf[j] = logApproximateBayesFactor(conditional_beta[j], conditional_se[j], options_.prior_variance);
            if (std::abs(conditional_z[j]) > lead_abs_z) {
                lead_abs_z = std::abs(conditional_z[j]);
                lead = j;
            }
        }
        if (lead < 0) break;
        if (twoSidedNormalP(conditional_z[lead]) > options_.signal_p_threshold) {
            if (signal == 1) append_signal(0, conditional_beta, conditional_se, log_abf, eligible, lead);
            break;
        }
        append_signal(signal, conditional_beta, conditional_se, log_abf, eligible, lead);
        selected.push_back(lead);
        result.number_signals = signal;
    }

    result.number_variants = m;
    result.ld_dimension = m;
    result.runtime_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    return result;
}

FineMappingRegion parseRegion(const std::string& text) {
    const size_t colon = text.find(':');
    const size_t dash = text.find('-', colon == std::string::npos ? 0 : colon + 1);
    if (colon == std::string::npos || dash == std::string::npos || colon == 0 || dash <= colon + 1)
        throw std::invalid_argument("--finemap-region must have form chr:start-end");
    FineMappingRegion region;
    region.chrom = text.substr(0, colon);
    try {
        region.start = std::stoll(text.substr(colon + 1, dash - colon - 1));
        region.end = std::stoll(text.substr(dash + 1));
    } catch (...) {
        throw std::invalid_argument("--finemap-region contains invalid coordinates");
    }
    if (region.start < 0 || region.end < region.start)
        throw std::invalid_argument("--finemap-region coordinates are invalid");
    return region;
}

void writeFineMappingResults(const FineMappingResult& result,
                             const FineMappingOptions& options,
                             const std::string& out_prefix,
                             const std::string& ld_source) {
    std::ofstream variants(out_prefix + ".finemap.tsv");
    if (!variants) throw std::runtime_error("Cannot write fine-mapping output");
    variants << "CHR\tPOS\tSNP\tA1\tA2\tBETA\tSE\tABF\tPIP\tSIGNAL\tCREDIBLE_SET\tCUM_POSTERIOR\n";
    variants << std::setprecision(12);
    for (const auto& row : result.rows) {
        variants << row.variant.chrom << '\t' << row.variant.pos << '\t' << row.variant.snp_id << '\t'
                 << row.variant.a1 << '\t' << row.variant.a2 << '\t' << row.variant.beta << '\t'
                 << row.variant.se << '\t' << row.abf << '\t' << row.pip << '\t' << row.signal << '\t'
                 << row.credible_set << '\t' << row.cumulative_posterior << '\n';
    }

    std::ofstream sets(out_prefix + ".credible_sets.tsv");
    if (!sets) throw std::runtime_error("Cannot write credible-set output");
    sets << "SIGNAL\tCREDIBLE_SET\tSNP\tPIP\tCUM_POSTERIOR\tLEAD_SNP\tLD_WITH_LEAD\n";
    sets << std::setprecision(12);
    for (const auto& row : result.credible_sets) {
        sets << row.signal << '\t' << row.credible_set << '\t' << row.snp_id << '\t' << row.pip << '\t'
             << row.cumulative_posterior << '\t' << row.lead_snp << '\t' << row.ld_with_lead << '\n';
    }

    std::ofstream metadata(out_prefix + ".finemap.meta.tsv");
    if (!metadata) throw std::runtime_error("Cannot write fine-mapping metadata");
    metadata << "method\tCosmicScanner conditional multi-signal Wakefield ABF\n"
             << "prior_variance\t" << options.prior_variance << '\n'
             << "credible_level\t" << options.credible_level << '\n'
             << "max_signals\t" << options.max_signals << '\n'
             << "signal_p_threshold\t" << options.signal_p_threshold << '\n'
             << "ld_source\t" << ld_source << '\n'
             << "number_variants\t" << result.number_variants << '\n'
             << "ld_dimension\t" << result.ld_dimension << '\n'
             << "number_signals\t" << result.number_signals << '\n'
             << "rejected_variants\t" << result.rejected_variants << '\n'
             << "runtime_seconds\t" << std::setprecision(8) << result.runtime_seconds << '\n';
}

} // namespace cosmic::association::fine_mapping
