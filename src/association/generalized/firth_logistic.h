#pragma once

#include "association/association_test.h"

#include <vector>

namespace cosmic::association {

struct FirthOptions {
    int max_iterations = 100;
    double tolerance = 1e-8;
    double maximum_step = 5.0;
};

struct FirthResult {
    double beta = 0.0;
    double se = 0.0;
    double p_value = 1.0;
    int iterations = 0;
    bool converged = false;
};

class FirthLogisticTest final : public AssociationTest {
public:
    explicit FirthLogisticTest(FirthOptions options = {});

    AssociationModel getModelType() const override { return AssociationModel::FIRTH_LOGISTIC; }
    AssociationCapabilities capabilities() const override { return {true, false, false}; }
    void prepare(const ScanContext& context) override;
    std::vector<AssociationResult> scanBlock(const GenotypeBlock& block) override;

    const std::vector<FirthResult>& results() const { return results_; }

private:
    FirthOptions options_;
    const ScanContext* context_ = nullptr;
    std::vector<FirthResult> results_;
};

} // namespace cosmic::association
