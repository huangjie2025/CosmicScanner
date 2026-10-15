#pragma once

#include "association/association_test.h"

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace cosmic::association {

enum class MultiLocusMethod { FARMCPU, BLINK, CONDITIONAL, CONDITIONAL_JOINT };

struct MultiLocusOptions {
    MultiLocusMethod method = MultiLocusMethod::FARMCPU;
    int max_iterations = 5;
    int max_pseudo_qtn = 10;
    double ld_r2 = 0.7;
    std::vector<std::string> condition_ids;
    std::vector<double> initial_p_values;
};

struct MultiLocusResult {
    std::string variant_id;
    double beta = 0.0;
    double se = 0.0;
    double p_value = 1.0;
    int n_obs = 0;
    int iteration = 0;
    bool selected_pseudo_qtn = false;
};

MultiLocusOptions farmCpuOptions(int max_iterations, int max_pseudo_qtn);
MultiLocusOptions blinkOptions(int max_pseudo_qtn = 10);
MultiLocusOptions conditionalOptions(std::vector<std::string> condition_ids);
MultiLocusOptions conditionalJointOptions(std::vector<std::string> condition_ids);

class MultiLocusEngine final : public AssociationTest {
public:
    explicit MultiLocusEngine(MultiLocusOptions options = {});

    AssociationModel getModelType() const override;
    AssociationCapabilities capabilities() const override;
    void prepare(const ScanContext& context) override;
    std::vector<AssociationResult> scanBlock(const GenotypeBlock& block) override;

    const std::vector<MultiLocusResult>& results() const { return results_; }
    const std::vector<int>& selectedIndices() const { return selected_; }

private:
    MultiLocusOptions options_;
    const ScanContext* context_ = nullptr;
    std::vector<MultiLocusResult> results_;
    std::vector<int> selected_;
};

} // namespace cosmic::association
