#pragma once

#include "association/cojo/cojo_engine.h"
#include "association/sumstats/sumstats.h"

#include <string>
#include <vector>

namespace cosmic::association::fine_mapping {

struct FineMappingOptions {
    double prior_variance = 0.04;
    double credible_level = 0.95;
    int max_signals = 5;
    double signal_p_threshold = 5e-8;
};

struct FineMappingRow {
    sumstats::SumstatsRecord variant;
    double abf = 0.0;
    double pip = 0.0;
    int signal = 0;
    int credible_set = 0;
    double cumulative_posterior = 0.0;
};

struct CredibleSetRow {
    int signal = 0;
    int credible_set = 0;
    std::string snp_id;
    double pip = 0.0;
    double cumulative_posterior = 0.0;
    std::string lead_snp;
    double ld_with_lead = 0.0;
};

struct FineMappingResult {
    std::vector<FineMappingRow> rows;
    std::vector<CredibleSetRow> credible_sets;
    int number_variants = 0;
    int number_signals = 0;
    int rejected_variants = 0;
    int ld_dimension = 0;
    double runtime_seconds = 0.0;
};

struct FineMappingRegion {
    std::string chrom;
    long long start = 0;
    long long end = 0;
};

class FineMappingEngine {
public:
    explicit FineMappingEngine(FineMappingOptions options = {});

    FineMappingResult run(
        const std::vector<cojo::SummaryRecord>& summary_stats,
        const cojo::LDReference& ld_reference) const;

    static double logApproximateBayesFactor(double beta, double se, double prior_variance);

private:
    FineMappingOptions options_;
};

FineMappingRegion parseRegion(const std::string& text);

void writeFineMappingResults(const FineMappingResult& result,
                             const FineMappingOptions& options,
                             const std::string& out_prefix,
                             const std::string& ld_source);

} // namespace cosmic::association::fine_mapping
