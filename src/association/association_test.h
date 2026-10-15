#pragma once

#include <Eigen/Core>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace cosmic {

enum class AssociationModel {
    LMM, MVLMM, RRM, REPEATABILITY, MATERNAL, SINGLE_STEP, ANIMAL,
    THRESHOLD, COUNT, ORDINAL, SURVIVAL, GXE, COMMON_ENV,
    SKAT, BURDEN, FARMCPU, BLINK, CONDITIONAL, CONDITIONAL_JOINT, FIRTH_LOGISTIC, SAIGE_STYLE,
    LOGISTIC, LOGISTIC_GLMM
};

// Borrowed inputs must outlive the prepared test and must not be mutated during use.
struct ScanContext {
    const Eigen::MatrixXd& phenotypes; // observations x traits
    const Eigen::MatrixXd& covariates; // observations x covariates (include intercept explicitly)
    // One zero-based genotype row per observation; repeats and reordering are allowed.
    // Required even for identity mapping; all entries must be valid for each block.
    const std::vector<Eigen::Index>& observation_to_sample;
};

struct GenotypeBlock {
    Eigen::MatrixXd dosages; // samples x variants, indexed by observation_to_sample
    std::vector<std::string> variant_ids; // one per dosage column
};

enum class AssociationStatus {
    NOT_TESTED,
    SUCCESS,
    FILTERED,
    INVALID_INPUT,
    NUMERICAL_FAILURE,
    UNSUPPORTED
};

// One hypothesis, potentially spanning multiple effects, traits, or variants.
// Numeric inference fields are meaningful only on SUCCESS; missing values stay NaN.
struct AssociationResult {
    std::string test_id;
    AssociationStatus status = AssociationStatus::NOT_TESTED;
    std::string message;
    Eigen::Index n_obs = 0;
    Eigen::VectorXd effects;
    Eigen::MatrixXd effect_covariance;
    double statistic = std::numeric_limits<double>::quiet_NaN();
    double degrees_of_freedom = std::numeric_limits<double>::quiet_NaN();
    double p_value = std::numeric_limits<double>::quiet_NaN();
};

struct AssociationCapabilities {
    bool block_scan = false;
    bool iterative_scan = false;
    bool set_test = false;
};

// Shared extension interface; model factories register only explicitly implemented adapters.
// Iterative and set workflows are distinct capabilities, not independent block scans;
// their execution contracts are intentionally deferred until adapters are implemented.
class AssociationTest {
public:
    virtual ~AssociationTest() = default;
    virtual AssociationModel getModelType() const = 0;
    virtual AssociationCapabilities capabilities() const = 0;
    virtual void prepare(const ScanContext& context) = 0;
    // A block-capable adapter must override this; the default explicitly rejects use.
    virtual std::vector<AssociationResult> scanBlock(const GenotypeBlock& block);
};

using AssociationTestPtr = std::unique_ptr<AssociationTest>;

class AssociationTestFactory {
public:
    // FarmCPU, BLINK, Firth logistic, and SAIGE-style have native adapters. Existing
    // production engines remain separate and unavailable adapters throw.
    static AssociationTestPtr create(AssociationModel model);
};

} // namespace cosmic
