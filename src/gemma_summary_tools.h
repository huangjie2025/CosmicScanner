#pragma once

#include "bgen_reader.h"
#include "config.h"
#include "pgen_reader.h"
#include "plink_reader.h"

namespace cosmic {

class GemmaSummaryTools {
public:
    explicit GemmaSummaryTools(const GwasConfig& cfg);

    void runGs(PlinkReader& plink_reader,
               BgenReader& bgen_reader,
               PgenReader& pgen_reader,
               bool is_bgen,
               bool is_pgen);
    void runGq();
    void runSummaryVc();
    void runBetaCorVc();
    void runSummaryCi();

private:
    const GwasConfig& cfg_;
};

} // namespace cosmic
