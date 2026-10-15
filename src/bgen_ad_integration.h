#pragma once
#include "bgen_ad_null.h"
#include <optional>

namespace cosmic {
struct BgenAdResearchOptions {
    std::string population_file;
    std::string snp_weight_file;
};

struct BgenAdIntegrationNull {
    std::string label;
    BgenAdNullResult fit;
    long long background_markers = 0;
    double denominator_a = 0.0;
    double denominator_d = 0.0;
    long long tested = 0;
    long long rank_deficient = 0;
};

struct BgenAdIntegrationDiagnostics {
    std::optional<BgenAdKernelResult> kernels;
    std::vector<BgenAdIntegrationNull> nulls;
    long long per_snp_null_refits = 0;
};
}
