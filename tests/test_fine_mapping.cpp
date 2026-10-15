#include "association/fine_mapping/fine_mapping.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using cosmic::association::cojo::LDReference;
using cosmic::association::cojo::SummaryRecord;
using cosmic::association::fine_mapping::FineMappingEngine;
using cosmic::association::fine_mapping::FineMappingOptions;

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

LDReference makeReference(const Eigen::MatrixXd& ld) {
    LDReference ref;
    ref.ld_matrix = ld;
    for (int i = 0; i < ld.rows(); ++i) {
        ref.snp_ids.push_back("rs" + std::to_string(i + 1));
        ref.chroms.push_back("1");
        ref.positions.push_back(100 + i * 100);
        ref.a1s.push_back("A");
        ref.a2s.push_back("C");
        ref.snp_to_idx[ref.snp_ids.back()] = i;
    }
    return ref;
}

std::vector<SummaryRecord> makeSummary(const std::vector<double>& z, double se = 0.1) {
    std::vector<SummaryRecord> result;
    for (size_t i = 0; i < z.size(); ++i) {
        SummaryRecord record;
        record.snp_id = "rs" + std::to_string(i + 1);
        record.chrom = "1";
        record.position = 100 + static_cast<long long>(i) * 100;
        record.a1 = "A";
        record.a2 = "C";
        record.beta = z[i] * se;
        record.se = se;
        record.p_value = std::erfc(std::abs(z[i]) / std::sqrt(2.0));
        record.freq = 0.25;
        record.n = 1000;
        result.push_back(record);
    }
    return result;
}

double signalPipSum(const cosmic::association::fine_mapping::FineMappingResult& result, int signal) {
    double total = 0.0;
    for (const auto& row : result.rows) if (row.signal == signal) total += row.pip;
    return total;
}

} // namespace

int main() {
    const double expected_log_abf = 0.7952810437829499;
    check(std::abs(FineMappingEngine::logApproximateBayesFactor(0.2, 0.1, 0.04) - expected_log_abf) < 1e-12,
          "Wakefield ABF oracle mismatch");
    check(std::isfinite(FineMappingEngine::logApproximateBayesFactor(100.0, 1e-5, 0.04)),
          "extreme-significance ABF is not finite");

    FineMappingOptions options;
    options.signal_p_threshold = 1e-4;
    options.max_signals = 3;

    Eigen::MatrixXd one_ld = Eigen::MatrixXd::Identity(3, 3);
    one_ld(0, 1) = one_ld(1, 0) = 0.7;
    auto one = FineMappingEngine(options).run(makeSummary({8.0, 5.2, 0.2}), makeReference(one_ld));
    check(one.number_signals == 1, "single-causal locus did not yield one signal");
    check(std::abs(signalPipSum(one, 1) - 1.0) < 1e-12, "single-signal PIPs do not normalize");
    check(!one.credible_sets.empty() && one.credible_sets.front().lead_snp == "rs1", "causal SNP is not credible-set lead");
    bool causal_member = false;
    for (const auto& row : one.credible_sets) causal_member |= row.snp_id == "rs1";
    check(causal_member, "causal SNP missing from credible set");

    auto ambiguous_reference = makeReference(one_ld);
    ambiguous_reference.a2s[2] = "T";
    auto ambiguous_summary = makeSummary({8.0, 5.2, 0.2});
    ambiguous_summary[2].a2 = "T";
    auto ambiguity_filtered = FineMappingEngine(options).run(ambiguous_summary, ambiguous_reference);
    check(ambiguity_filtered.rejected_variants == 1 && ambiguity_filtered.number_variants == 2,
          "unsafe ambiguous SNP was not rejected");

    Eigen::MatrixXd two_ld = Eigen::MatrixXd::Identity(4, 4);
    two_ld(0, 1) = two_ld(1, 0) = 0.7;
    two_ld(2, 3) = two_ld(3, 2) = 0.65;
    auto two = FineMappingEngine(options).run(makeSummary({8.0, 5.6, 7.0, 4.55}), makeReference(two_ld));
    check(two.number_signals == 2, "two independent causal signals were not separated");
    check(std::abs(signalPipSum(two, 1) - 1.0) < 1e-12 && std::abs(signalPipSum(two, 2) - 1.0) < 1e-12,
          "multi-signal PIPs do not normalize per signal");
    check(two.credible_sets.front().lead_snp == "rs1", "first independent lead is incorrect");
    bool second_lead = false;
    for (const auto& row : two.credible_sets) second_lead |= row.signal == 2 && row.lead_snp == "rs3";
    check(second_lead, "second independent lead is incorrect");

    Eigen::MatrixXd high_ld = Eigen::MatrixXd::Identity(2, 2);
    high_ld(0, 1) = high_ld(1, 0) = 0.999;
    options.max_signals = 1;
    auto high_ld_result = FineMappingEngine(options).run(makeSummary({7.0, 7.0}), makeReference(high_ld));
    check(high_ld_result.credible_sets.size() == 2, "high-LD credible set should include both variants");
    check(std::abs(high_ld_result.rows[0].pip - 0.5) < 1e-12 && std::abs(high_ld_result.rows[1].pip - 0.5) < 1e-12,
          "high-LD posterior should be shared without artificial certainty");

    options.max_signals = 3;
    options.signal_p_threshold = 5e-8;
    auto null_result = FineMappingEngine(options).run(makeSummary({0.2, -0.5, 1.0}),
                                                      makeReference(Eigen::MatrixXd::Identity(3, 3)));
    check(null_result.number_signals == 0 && null_result.credible_sets.empty(),
          "null locus created unsupported signals: " + std::to_string(null_result.number_signals) +
          "/" + std::to_string(null_result.credible_sets.size()));
    check(null_result.rows.size() == 3 && std::abs(signalPipSum(null_result, 0) - 1.0) < 1e-12,
          "null-locus posterior is not finite and normalized");

    options.signal_p_threshold = 1e-4;
    options.max_signals = 1;
    auto reference = makeReference(one_ld);
    auto original_summary = makeSummary({8.0, 5.2, 0.2});
    auto swapped_summary = original_summary;
    for (auto& record : swapped_summary) {
        std::swap(record.a1, record.a2);
        record.beta = -record.beta;
    }
    auto original = FineMappingEngine(options).run(original_summary, reference);
    auto swapped = FineMappingEngine(options).run(swapped_summary, reference);
    check(original.rows.size() == swapped.rows.size(), "allele swap changed result size");
    for (size_t i = 0; i < original.rows.size(); ++i)
        check(std::abs(original.rows[i].pip - swapped.rows[i].pip) < 1e-12,
              "allele harmonization changed fine-mapping PIP");

    return 0;
}
