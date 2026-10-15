#include "association/meta/meta_engine.h"
#include "association/sumstats/sumstats.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace cosmic::association::meta {

namespace {

#ifndef M_1_SQRT2PI
#define M_1_SQRT2PI 0.39894228040143267793994605993438
#endif

double regularizedGammaQ(double shape, double x) {
    if (!(shape > 0.0) || x < 0.0) return 1.0;
    if (x == 0.0) return 1.0;
    if (!std::isfinite(x)) return 0.0;
    constexpr int max_iterations = 10000;
    constexpr double epsilon = 1e-14;
    constexpr double floor = 1e-300;
    const double log_scale = -x + shape * std::log(x) - std::lgamma(shape);
    if (x < shape + 1.0) {
        double term = 1.0 / shape;
        double sum = term;
        double denominator = shape;
        for (int i = 1; i <= max_iterations; ++i) {
            denominator += 1.0;
            term *= x / denominator;
            sum += term;
            if (std::abs(term) <= std::abs(sum) * epsilon) break;
        }
        return std::max(0.0, std::min(1.0, 1.0 - sum * std::exp(log_scale)));
    }
    double b = x + 1.0 - shape;
    double c = 1.0 / floor;
    double d = 1.0 / std::max(b, floor);
    double fraction = d;
    for (int i = 1; i <= max_iterations; ++i) {
        const double a = -static_cast<double>(i) * (static_cast<double>(i) - shape);
        b += 2.0;
        d = a * d + b;
        if (std::abs(d) < floor) d = floor;
        c = b + a / c;
        if (std::abs(c) < floor) c = floor;
        d = 1.0 / d;
        const double delta = d * c;
        fraction *= delta;
        if (std::abs(delta - 1.0) <= epsilon) break;
    }
    return std::max(0.0, std::min(1.0, fraction * std::exp(log_scale)));
}

double chi2Survival(double x, int df) {
    if (x <= 0.0 || df <= 0) return 1.0;
    return regularizedGammaQ(0.5 * df, 0.5 * x);
}

double normal_cdf(double x) {
    return 0.5 * std::erfc(-x * M_SQRT1_2);
}

// Parse canonical sumstats format
StudySummary parseCanonicalSumstatsLine(
    const std::string& line,
    const std::unordered_map<std::string, size_t>& columns) {
    std::istringstream iss(line);
    std::vector<std::string> fields;
    std::string field;
    while (iss >> field) fields.push_back(field);
    
    if (fields.size() < 9) return StudySummary{};
    
    StudySummary s;
    auto get = [&](const std::string& name) -> std::string {
        const auto it = columns.find(name);
        return it != columns.end() && it->second < fields.size() ? fields[it->second] : "";
    };
    auto get_numeric = [&](const std::string& primary, const std::string& fallback) -> std::string {
        const std::string value = get(primary);
        return value.empty() || value == "NA" || value == "." ? get(fallback) : value;
    };
    if (columns.count("CHR") && columns.count("SNP")) {
        s.chrom = get("CHR");
        s.snp_id = get("SNP");
        const std::string position = get_numeric("BP", "POS");
        const std::string beta = get_numeric("BETA_A", "BETA_D");
        const std::string se = get_numeric("SE_A", "SE_D");
        const std::string p = get_numeric("P_A", "P_D");
        s.position = std::stoll(position);
        s.a1 = get("A1");
        s.a2 = get("A2");
        s.beta = std::stod(beta.empty() ? get("BETA") : beta);
        s.se = std::stod(se.empty() ? get("SE") : se);
        s.p_value = std::stod(p.empty() ? get("P") : p);
        s.n = std::stoi(get("N"));
        const std::string freq = get_numeric("FREQ", "AF");
        if (!freq.empty()) s.freq = std::stod(freq);
        if (s.se > 0.0) s.z_score = s.beta / s.se;
        if (s.se > 0.0) s.weight = 1.0 / (s.se * s.se);
        return s;
    }

    // Detect format: canonical sumstats has CHR (numeric) first
    bool is_canonical = false;
    if (!fields[0].empty() && std::isdigit(fields[0][0])) {
        is_canonical = true;
    }
    
    if (is_canonical) {
        // Canonical: CHR, SNP, POS, A1, A2, BETA, SE, P, N, AF, ...
        s.chrom = fields[0];
        s.snp_id = fields[1];
        s.position = std::stoll(fields[2]);
        s.a1 = fields[3];
        s.a2 = fields[4];
        s.beta = std::stod(fields[5]);
        s.se = std::stod(fields[6]);
        s.p_value = std::stod(fields[7]);
        s.n = std::stoi(fields[8]);
        if (fields.size() > 9) s.freq = std::stod(fields[9]);
    } else {
        // Old format: SNP, CHR, POS, A1, A2, BETA, SE, P, N, FREQ
        s.snp_id = fields[0];
        s.chrom = fields[1];
        s.position = std::stoll(fields[2]);
        s.a1 = fields[3];
        s.a2 = fields[4];
        s.beta = std::stod(fields[5]);
        s.se = std::stod(fields[6]);
        s.p_value = std::stod(fields[7]);
        s.n = std::stoi(fields[8]);
        if (fields.size() > 9) s.freq = std::stod(fields[9]);
    }
    if (s.se > 0.0) s.z_score = s.beta / s.se;
    if (s.se > 0.0) s.weight = 1.0 / (s.se * s.se);
    return s;
}

// Keep old parser for backward compatibility
StudySummary parseSummaryLine(const std::string& line) {
    std::istringstream iss(line);
    StudySummary s;
    iss >> s.snp_id >> s.chrom >> s.position >> s.a1 >> s.a2 >> s.beta >> s.se >> s.p_value >> s.n >> s.freq;
    if (s.se > 0.0) s.z_score = s.beta / s.se;
    if (s.se > 0.0) s.weight = 1.0 / (s.se * s.se);
    return s;
}

} // namespace

MetaAnalysisEngine::MetaAnalysisEngine(const MetaAnalysisOptions& opts)
    : options_(opts) {}

void MetaAnalysisEngine::addStudy(const std::string& study_id,
                                   const std::vector<StudySummary>& summaries) {
    if (!study_indices_.count(study_id)) {
        study_indices_[study_id] = static_cast<int>(study_order_.size());
        study_order_.push_back(study_id);
    }
    studies_[study_id] = summaries;
    std::unordered_set<std::string> seen;
    for (auto& s : studies_[study_id]) {
        if (s.snp_id.empty() || !seen.insert(s.snp_id).second)
            throw std::runtime_error("Duplicate SNP in study " + study_id + ": " + s.snp_id);
        s.study_id = study_id;
        if (s.se > 0.0) s.weight = 1.0 / (s.se * s.se);
    }
}

void MetaAnalysisEngine::addStudyFromFile(const std::string& study_id,
                                            const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file) throw std::runtime_error("Cannot open summary file: " + filepath);

    std::string line;
    // Skip all comment lines
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        // Found first non-comment line - this should be the header
        break;
    }
    
    // Skip the header line
    if (line.empty() || line[0] == '#') {
        // No data lines
        addStudy(study_id, {});
        return;
    }
    std::unordered_map<std::string, size_t> columns;
    {
        std::istringstream header(line);
        std::string name;
        size_t index = 0;
        while (header >> name) columns[name] = index++;
    }

    std::vector<StudySummary> summaries;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        StudySummary s = parseCanonicalSumstatsLine(line, columns);
        if (s.se <= 0.0 || !std::isfinite(s.beta) || !std::isfinite(s.se) ||
            !std::isfinite(s.p_value) || s.n < 1) continue;
        summaries.push_back(std::move(s));
    }
    addStudy(study_id, std::move(summaries));
}

MetaResult MetaAnalysisEngine::metaAnalyzeSnp(
    const std::vector<StudySummary>& snp_data) const {

    MetaResult result;
    if (snp_data.empty()) {
        result.valid = false;
        result.warning = "No study data";
        return result;
    }

    const auto& first = snp_data[0];
    result.snp_id = first.snp_id;
    result.chrom = first.chrom;
    result.position = first.position;
    result.a1 = first.a1;
    result.a2 = first.a2;

    std::vector<double> betas, weights;
    std::vector<int> ns;
    double total_weight = 0.0;
    double weighted_freq = 0.0;
    int total_n = 0;

    for (const auto& s : snp_data) {
        if (s.se <= 0.0 || !std::isfinite(s.beta) || !std::isfinite(s.se)) continue;
        double w = s.weight;
        if (options_.use_gc_correction && options_.gc_lambda > 1.0) {
            w /= options_.gc_lambda;
        }
        if (!(w > 0.0) || !std::isfinite(w)) continue;
        betas.push_back(s.beta);
        weights.push_back(w);
        ns.push_back(s.n);
        total_weight += w;
        weighted_freq += s.freq * s.n;
        total_n += s.n;
    }

    int k = static_cast<int>(betas.size());
    if (k < options_.min_studies) {
        result.valid = false;
        result.warning = "Insufficient valid studies: " + std::to_string(k);
        return result;
    }

    result.n_studies = k;
    result.total_n = total_n;
    result.freq = total_n > 0 ? weighted_freq / total_n : 0.0;

    if (!(total_weight > 0.0) || !std::isfinite(total_weight)) {
        result.valid = false;
        result.warning = "Invalid inverse-variance weights";
        return result;
    }

    double beta_sum = 0.0;
    for (int i = 0; i < k; ++i) beta_sum += weights[i] * betas[i];
    double beta_fe = beta_sum / total_weight;
    double se_fe = std::sqrt(1.0 / total_weight);

    result.heterogeneity = computeHeterogeneity(betas, weights, ns);

    double beta_final = beta_fe;
    double se_final = se_fe;

    if (options_.method == MetaAnalysisMethod::RANDOM_EFFECT_DL ||
        options_.method == MetaAnalysisMethod::RANDOM_EFFECT_REML) {

        double tau2 = estimateTau2DL(betas, weights, result.heterogeneity);
        if (options_.method == MetaAnalysisMethod::RANDOM_EFFECT_REML) {
            tau2 = estimateTau2REML(betas, weights, tau2);
        }
        result.heterogeneity.tau2 = tau2;

        if (tau2 > 0.0) {
            std::vector<double> weights_re(k);
            double total_weight_re = 0.0;
            for (int i = 0; i < k; ++i) {
                double v = 1.0 / weights[i] + tau2;
                weights_re[i] = v > 0.0 ? 1.0 / v : 0.0;
                total_weight_re += weights_re[i];
            }
            if (total_weight_re > 0.0) {
                double beta_sum_re = 0.0;
                for (int i = 0; i < k; ++i) beta_sum_re += weights_re[i] * betas[i];
                beta_final = beta_sum_re / total_weight_re;
                se_final = std::sqrt(1.0 / total_weight_re);
            }
        }
    }

    result.beta = beta_final;
    result.se = se_final;
    result.z_score = se_final > 0.0 ? beta_final / se_final : 0.0;
    result.p_value = 2.0 * normal_cdf(-std::abs(result.z_score));
    result.valid = true;

    return result;
}

HeterogeneityStats MetaAnalysisEngine::computeHeterogeneity(
    const std::vector<double>& betas,
    const std::vector<double>& weights,
    const std::vector<int>& ns) const {

    HeterogeneityStats het;
    int k = static_cast<int>(betas.size());
    if (k < 2) return het;

    double beta_sum = 0.0, weight_sum = 0.0;
    for (int i = 0; i < k; ++i) {
        beta_sum += weights[i] * betas[i];
        weight_sum += weights[i];
    }
    double beta_fe = weight_sum > 0.0 ? beta_sum / weight_sum : 0.0;

    double q = 0.0;
    for (int i = 0; i < k; ++i) {
        double diff = betas[i] - beta_fe;
        q += weights[i] * diff * diff;
    }

    het.q_statistic = q;
    het.q_df = k - 1;
    het.q_p_value = chi2Survival(q, k - 1);
    het.significant_heterogeneity = het.q_p_value < options_.heterogeneity_threshold;

    if (q > het.q_df && q > 0.0) {
        het.i2 = 100.0 * (q - het.q_df) / q;
        het.i2 = std::max(0.0, std::min(100.0, het.i2));
    } else {
        het.i2 = 0.0;
    }

    double sum_w2 = 0.0;
    for (int i = 0; i < k; ++i) sum_w2 += weights[i] * weights[i];
    const double c = weight_sum > 0.0 ? weight_sum - sum_w2 / weight_sum : 0.0;
    het.tau2 = c > 0.0 ? std::max(0.0, (q - het.q_df) / c) : 0.0;
    het.h2 = het.q_df > 0.0 ? std::max(1.0, q / het.q_df) : 1.0;

    return het;
}

double MetaAnalysisEngine::estimateTau2DL(
    const std::vector<double>& betas,
    const std::vector<double>& weights,
    const HeterogeneityStats& het) const {

    if (het.q_df <= 0.0) return 0.0;

    double sum_w = 0.0;
    for (double w : weights) sum_w += w;
    double sum_w2 = 0.0;
    for (double w : weights) sum_w2 += w * w;
    const double c = sum_w > 0.0 ? sum_w - sum_w2 / sum_w : 0.0;
    if (!(c > 0.0)) return 0.0;
    return std::max(0.0, (het.q_statistic - het.q_df) / c);
}

double MetaAnalysisEngine::estimateTau2REML(
    const std::vector<double>& betas,
    const std::vector<double>& weights,
    double initial_tau2) const {

    const int k = static_cast<int>(betas.size());
    if (k < 2) return 0.0;
    auto score = [&](double tau2) {
        std::vector<double> random_weights(k);
        double sum_w = 0.0;
        double sum_w2 = 0.0;
        double weighted_mean = 0.0;
        for (int i = 0; i < k; ++i) {
            random_weights[i] = 1.0 / (1.0 / weights[i] + tau2);
            sum_w += random_weights[i];
            sum_w2 += random_weights[i] * random_weights[i];
            weighted_mean += random_weights[i] * betas[i];
        }
        weighted_mean /= sum_w;
        double residual = 0.0;
        for (int i = 0; i < k; ++i) {
            const double difference = betas[i] - weighted_mean;
            residual += random_weights[i] * random_weights[i] * difference * difference;
        }
        return residual - (sum_w - sum_w2 / sum_w);
    };

    if (score(0.0) <= 0.0) return 0.0;
    double lower = 0.0;
    double upper = std::max(1.0, initial_tau2);
    for (int i = 0; i < 100 && score(upper) > 0.0; ++i) upper *= 2.0;
    if (score(upper) > 0.0) return upper;
    for (int i = 0; i < 200; ++i) {
        const double midpoint = 0.5 * (lower + upper);
        if (score(midpoint) > 0.0) lower = midpoint;
        else upper = midpoint;
        if (upper - lower <= 1e-10 * (1.0 + midpoint)) break;
    }
    return 0.5 * (lower + upper);
}

double MetaAnalysisEngine::computeQStatistic(
    const std::vector<double>& betas,
    const std::vector<double>& weights) const {

    double beta_sum = 0.0, weight_sum = 0.0;
    for (size_t i = 0; i < betas.size(); ++i) {
        beta_sum += weights[i] * betas[i];
        weight_sum += weights[i];
    }
    if (weight_sum <= 0.0) return 0.0;

    double beta_fe = beta_sum / weight_sum;
    double q = 0.0;
    for (size_t i = 0; i < betas.size(); ++i) {
        double diff = betas[i] - beta_fe;
        q += weights[i] * diff * diff;
    }
    return q;
}

const std::unordered_map<std::string, std::vector<StudySummary>>&
MetaAnalysisEngine::getStudies() const {
    return studies_;
}

std::vector<MetaResult> MetaAnalysisEngine::run() {
    if (studies_.empty()) return {};

    std::unordered_map<std::string, std::pair<std::string, std::string>> ref_alleles;
    for (const auto& study_id : study_order_) {
        for (const auto& summary : studies_.at(study_id)) {
            ref_alleles.emplace(summary.snp_id, std::make_pair(summary.a1, summary.a2));
        }
    }

    for (const auto& study_id : study_order_) {
        auto& summaries = studies_.at(study_id);
        using namespace cosmic::association::sumstats;
        
        std::vector<SumstatsRecord> sumstats_records;
        sumstats_records.reserve(summaries.size());
        for (const auto& s : summaries) {
            SumstatsRecord r;
            r.snp_id = s.snp_id;
            r.chrom = s.chrom;
            r.pos = s.position;
            r.a1 = s.a1;
            r.a2 = s.a2;
            r.beta = s.beta;
            r.se = s.se;
            r.p_value = s.p_value;
            r.n = s.n;
            r.freq = s.freq;
            r.weight = s.weight;
            sumstats_records.push_back(r);
        }
        
        auto harmonized = harmonizeSumstats(sumstats_records, ref_alleles);
        
        std::vector<StudySummary> harmonized_summaries;
        harmonized_summaries.reserve(summaries.size());
        
        for (const auto& h : harmonized) {
            if (h.harmonized) {
                StudySummary s;
                s.snp_id = h.record.snp_id;
                s.chrom = h.record.chrom;
                s.position = h.record.pos;
                s.a1 = h.record.a1;
                s.a2 = h.record.a2;
                s.beta = h.record.beta;
                s.se = h.record.se;
                s.p_value = h.record.p_value;
                s.n = h.record.n;
                s.freq = h.record.freq;
                s.z_score = h.record.z_score;
                s.weight = h.record.weight;
                harmonized_summaries.push_back(s);
            } else {
            }
        }
        
        summaries = std::move(harmonized_summaries);
    }

    std::unordered_map<std::string, std::vector<StudySummary>> snp_to_data;
    std::vector<std::string> snp_order;
    for (const auto& study_id : study_order_) {
        for (const auto& summary : studies_.at(study_id)) {
            auto [found, inserted] = snp_to_data.emplace(summary.snp_id, std::vector<StudySummary>{});
            if (inserted) snp_order.push_back(summary.snp_id);
            found->second.push_back(summary);
        }
    }

    std::vector<MetaResult> results;
    results.reserve(snp_to_data.size());

    for (const auto& snp_id : snp_order) {
        MetaResult result = metaAnalyzeSnp(snp_to_data.at(snp_id));
        if (result.valid) results.push_back(std::move(result));
    }

    return results;
}

void writeMetaResults(
    const std::vector<MetaResult>& results,
    const std::string& out_prefix,
    const MetaWriterOptions& opts) {

    if (opts.write_metal) {
        std::ofstream out(out_prefix + ".meta");
        out << "SNP\tCHR\tBP\tA1\tA2\tFREQ\tBETA\tSE\tP\tN_STUDIES\tTOTAL_N\tQ\tQ_DF\tQ_P\tI2\tTAU2\n";
        out << std::setprecision(10);
        for (const auto& r : results) {
            out << r.snp_id << '\t' << r.chrom << '\t' << r.position << '\t'
                << r.a1 << '\t' << r.a2 << '\t' << r.freq << '\t'
                << r.beta << '\t' << r.se << '\t' << r.p_value << '\t'
                << r.n_studies << '\t' << r.total_n << '\t'
                << r.heterogeneity.q_statistic << '\t'
                << r.heterogeneity.q_df << '\t'
                << r.heterogeneity.q_p_value << '\t'
                << r.heterogeneity.i2 << '\t'
                << r.heterogeneity.tau2 << '\n';
        }
    }

    if (opts.write_simple) {
        std::ofstream out(out_prefix + ".meta.assoc");
        out << "SNP\tCHR\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\n";
        out << std::setprecision(10);
        for (const auto& r : results) {
            out << r.snp_id << '\t' << r.chrom << '\t' << r.position << '\t'
                << r.a1 << '\t' << r.a2 << '\t' << r.freq << '\t'
                << r.total_n << '\t' << r.beta << '\t' << r.se << '\t'
                << r.p_value << '\n';
        }
    }
}

} // namespace cosmic::association::meta
