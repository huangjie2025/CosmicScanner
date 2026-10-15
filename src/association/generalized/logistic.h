#pragma once

#include "association/association_test.h"

#include <string>
#include <vector>

namespace cosmic::association {

enum class LogisticTestType { SCORE, WALD };

struct LogisticOptions {
    LogisticTestType test = LogisticTestType::SCORE;
    int max_iterations = 50;
    double tolerance = 1e-8;
    double maximum_step = 10.0;
};

struct LogisticResult {
    double beta = 0.0;
    double se = 0.0;
    double statistic = 0.0;
    double p_value = 1.0;
    int iterations = 0;
    bool converged = false;
    bool separated = false;
};

class LogisticTest final : public AssociationTest {
public:
    explicit LogisticTest(LogisticOptions options = {});

    AssociationModel getModelType() const override { return AssociationModel::LOGISTIC; }
    AssociationCapabilities capabilities() const override { return {true, false, false}; }
    void prepare(const ScanContext& context) override;
    std::vector<AssociationResult> scanBlock(const GenotypeBlock& block) override;

    const std::vector<LogisticResult>& results() const { return results_; }
    const Eigen::VectorXd& nullProbabilities() const { return null_probability_; }

private:
    LogisticResult fitWald(const Eigen::MatrixXd& design) const;

    LogisticOptions options_;
    Eigen::VectorXd phenotype_, null_beta_, null_probability_, null_weight_, null_residual_;
    Eigen::MatrixXd covariates_, null_information_inverse_;
    std::vector<Eigen::Index> observation_to_sample_;
    std::vector<LogisticResult> results_;
};

} // namespace cosmic::association
