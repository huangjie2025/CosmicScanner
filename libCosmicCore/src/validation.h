#pragma once

#include "analysis_dataset.h"

#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cosmic {

enum class ValidationSeverity {
    Pass,
    Warning,
    Fail
};

enum class ValidationStage {
    PreRun,
    Runtime,
    PostRun
};

struct ValidationIssue {
    std::string code;
    ValidationStage stage = ValidationStage::PreRun;
    ValidationSeverity severity = ValidationSeverity::Pass;
    std::string context;
    std::string message;
    std::string field;
    std::string observed_value;
    std::string expected_value;
};

struct ValidationCounts {
    std::size_t passed = 0;
    std::size_t warnings = 0;
    std::size_t failures = 0;
};

class ValidationReport {
public:
    void add(ValidationSeverity severity,
             ValidationStage stage,
             std::string code,
             std::string context,
             std::string message,
             std::string field = {},
             std::string observed_value = {},
             std::string expected_value = {}) {
        issues_.push_back(ValidationIssue{
            std::move(code), stage, severity, std::move(context), std::move(message),
            std::move(field), std::move(observed_value), std::move(expected_value)});
    }

    void pass(ValidationStage stage, std::string code, std::string context, std::string message) {
        add(ValidationSeverity::Pass, stage, std::move(code), std::move(context), std::move(message));
    }

    void warning(ValidationStage stage, std::string code, std::string context, std::string message) {
        add(ValidationSeverity::Warning, stage, std::move(code), std::move(context), std::move(message));
    }

    void fail(ValidationStage stage, std::string code, std::string context, std::string message) {
        add(ValidationSeverity::Fail, stage, std::move(code), std::move(context), std::move(message));
    }

    void merge(const ValidationReport& other) {
        issues_.insert(issues_.end(), other.issues_.begin(), other.issues_.end());
    }

    const std::vector<ValidationIssue>& issues() const { return issues_; }
    bool empty() const { return issues_.empty(); }

    ValidationCounts counts() const {
        ValidationCounts out;
        for (const auto& issue : issues_) {
            switch (issue.severity) {
                case ValidationSeverity::Pass: ++out.passed; break;
                case ValidationSeverity::Warning: ++out.warnings; break;
                case ValidationSeverity::Fail: ++out.failures; break;
            }
        }
        return out;
    }

    bool has_fail() const { return counts().failures > 0; }
    bool has_warning() const { return counts().warnings > 0; }

    ValidationSeverity overall_status() const {
        const auto c = counts();
        if (c.failures > 0) return ValidationSeverity::Fail;
        if (c.warnings > 0) return ValidationSeverity::Warning;
        return ValidationSeverity::Pass;
    }

    std::string summary(bool include_passes = false) const {
        const auto c = counts();
        std::ostringstream os;
        os << "Validation : " << to_string(overall_status()) << "\n";
        os << "Checks     : " << c.passed << " passed, "
           << c.warnings << " warnings, " << c.failures << " failures\n";
        for (const auto& issue : issues_) {
            if (issue.severity == ValidationSeverity::Pass && !include_passes) continue;
            os << "  [" << issue.code << "] " << issue.message;
            if (!issue.context.empty()) os << " (" << issue.context << ")";
            os << "\n";
        }
        return os.str();
    }

    void throw_if_fail() const {
        if (has_fail()) {
            throw std::runtime_error(summary(false));
        }
    }

private:
    static const char* to_string(ValidationSeverity severity) {
        switch (severity) {
            case ValidationSeverity::Pass: return "PASS";
            case ValidationSeverity::Warning: return "WARNING";
            case ValidationSeverity::Fail: return "FAIL";
        }
        return "FAIL";
    }

    std::vector<ValidationIssue> issues_;
};

inline const char* to_string(ValidationSeverity severity) {
    switch (severity) {
        case ValidationSeverity::Pass: return "PASS";
        case ValidationSeverity::Warning: return "WARNING";
        case ValidationSeverity::Fail: return "FAIL";
    }
    return "FAIL";
}

inline const char* to_string(ValidationStage stage) {
    switch (stage) {
        case ValidationStage::PreRun: return "PRE_RUN";
        case ValidationStage::Runtime: return "RUNTIME";
        case ValidationStage::PostRun: return "POST_RUN";
    }
    return "PRE_RUN";
}

inline std::string validation_count_string(std::size_t value) {
    return std::to_string(value);
}

inline void validate_sample_registry(ValidationReport& report,
                                     const SampleRegistry& registry,
                                     const std::string& context,
                                     bool require_phenotype = false,
                                     bool require_covariates = false) {
    const auto& d = registry.diagnostics();
    if (registry.empty()) {
        report.fail(ValidationStage::PreRun, "IO_EMPTY_ANALYSIS_SET", context,
                    "analysis sample set is empty");
    } else {
        report.pass(ValidationStage::PreRun, "IO_SAMPLE_REGISTRY_PRESENT", context,
                    "analysis sample registry is present");
    }
    if (d.duplicate_genotype_ids > 0) {
        report.fail(ValidationStage::PreRun, "IO_DUPLICATE_SAMPLE_ID", context,
                    "duplicate genotype/sample IDs make sample identity ambiguous");
    }
    if (d.duplicate_phenotype_ids > 0) {
        report.fail(ValidationStage::PreRun, "IO_DUPLICATE_PHENOTYPE_ID", context,
                    "duplicate phenotype IDs make phenotype ownership ambiguous");
    }
    if (d.duplicate_covariate_ids > 0) {
        report.fail(ValidationStage::PreRun, "IO_DUPLICATE_COVARIATE_ID", context,
                    "duplicate covariate IDs make covariate ownership ambiguous");
    }
    if (require_phenotype && d.matched_phenotypes == 0) {
        report.fail(ValidationStage::PreRun, "IO_SAMPLE_MISMATCH", context,
                    "no phenotype rows match analysis samples");
    }
    if (require_phenotype && d.genotype_without_phenotype > 0) {
        report.warning(ValidationStage::PreRun, "IO_MISSING_PHENOTYPE", context,
                       validation_count_string(d.genotype_without_phenotype) +
                       " analysis samples have no matched phenotype row");
    }
    if (d.phenotype_without_genotype > 0) {
        report.warning(ValidationStage::PreRun, "IO_UNMATCHED_PHENOTYPE_ID", context,
                       validation_count_string(d.phenotype_without_genotype) +
                       " phenotype rows do not match genotype/matrix samples");
    }
    if (require_covariates && d.matched_covariates == 0) {
        report.fail(ValidationStage::PreRun, "IO_COVARIATE_MISMATCH", context,
                    "no covariate rows match analysis samples");
    }
    if (require_covariates && d.genotype_without_covariate > 0) {
        report.warning(ValidationStage::PreRun, "IO_MISSING_COVARIATE", context,
                       validation_count_string(d.genotype_without_covariate) +
                       " analysis samples have no matched covariate row");
    }
    if (d.covariate_without_genotype > 0) {
        report.warning(ValidationStage::PreRun, "IO_UNMATCHED_COVARIATE_ID", context,
                       validation_count_string(d.covariate_without_genotype) +
                       " covariate rows do not match genotype/matrix samples");
    }
}

inline void validate_analysis_dataset_prerun(ValidationReport& report,
                                             const AnalysisDataset& dataset,
                                             const std::string& context,
                                             bool require_phenotype,
                                             bool require_markers,
                                             bool require_covariates = false,
                                             bool allow_missing_phenotype = true,
                                             bool allow_missing_covariates = true) {
    validate_sample_registry(report, dataset.samples(), context, require_phenotype, require_covariates);

    const auto genotype = dataset.genotype();
    if (genotype) {
        const auto& meta = genotype->metadata();
        if (meta.n_samples != dataset.samples().size()) {
            report.fail(ValidationStage::PreRun, "GENOTYPE_SAMPLE_COUNT_MISMATCH", context,
                        "genotype metadata sample count differs from SampleRegistry size");
        } else {
            report.pass(ValidationStage::PreRun, "GENOTYPE_SAMPLE_COUNT_MATCH", context,
                        "genotype metadata sample count matches SampleRegistry");
        }
        if (require_markers && meta.n_markers == 0) {
            report.fail(ValidationStage::PreRun, "GENOTYPE_EMPTY_MARKER_SET", context,
                        "genotype marker count is zero");
        }
    }

    const auto* phenotype = dataset.phenotype();
    if (require_phenotype && !phenotype) {
        report.fail(ValidationStage::PreRun, "PHENOTYPE_REQUIRED", context,
                    "required phenotype table is absent");
    } else if (phenotype) {
        if (phenotype->traits() == 0 || phenotype->values.size() != phenotype->rows() * phenotype->traits()) {
            report.fail(ValidationStage::PreRun, "PHENOTYPE_DIMENSION_MISMATCH", context,
                        "phenotype table dimensions do not match stored values");
        } else {
            std::size_t missing = 0;
            for (std::size_t i = 0; i < phenotype->values.size(); ++i) {
                const bool marked_missing = i < phenotype->missing_mask.size() && phenotype->missing_mask[i];
                const double v = phenotype->values[i];
                if (marked_missing || std::isnan(v)) {
                    ++missing;
                    continue;
                }
                if (!std::isfinite(v)) {
                    report.fail(ValidationStage::PreRun, "PHENOTYPE_NONFINITE", context,
                                "phenotype contains non-finite value that is not a missing marker");
                    break;
                }
            }
            if (missing > 0) {
                if (allow_missing_phenotype) {
                    report.warning(ValidationStage::PreRun, "IO_MISSING_PHENOTYPE", context,
                                   validation_count_string(missing) +
                                   " phenotype values are missing and will follow existing model filtering");
                } else {
                    report.fail(ValidationStage::PreRun, "IO_MISSING_PHENOTYPE", context,
                                "missing phenotype values are not allowed for this path");
                }
            }
            report.pass(ValidationStage::PreRun, "PHENOTYPE_DIMENSIONS_VALID", context,
                        "phenotype table dimensions are valid");
        }
    }

    const auto* covariates = dataset.covariates();
    if (require_covariates && !covariates) {
        report.fail(ValidationStage::PreRun, "COVARIATE_REQUIRED", context,
                    "required covariate table is absent");
    } else if (covariates) {
        if (covariates->numeric_values.size() != covariates->rows() * covariates->numeric_cols()) {
            report.fail(ValidationStage::PreRun, "COVARIATE_DIMENSION_MISMATCH", context,
                        "covariate table dimensions do not match stored values");
        } else {
            std::size_t missing = 0;
            for (std::size_t i = 0; i < covariates->numeric_values.size(); ++i) {
                const bool marked_missing = i < covariates->numeric_missing_mask.size() &&
                                            covariates->numeric_missing_mask[i];
                const double v = covariates->numeric_values[i];
                if (marked_missing || std::isnan(v)) {
                    ++missing;
                    continue;
                }
                if (!std::isfinite(v)) {
                    report.fail(ValidationStage::PreRun, "COVARIATE_NONFINITE", context,
                                "covariate table contains non-finite value that is not a missing marker");
                    break;
                }
            }
            if (missing > 0) {
                if (allow_missing_covariates) {
                    report.warning(ValidationStage::PreRun, "IO_MISSING_COVARIATE", context,
                                   validation_count_string(missing) +
                                   " covariate values are missing and will follow existing model filtering");
                } else {
                    report.fail(ValidationStage::PreRun, "IO_MISSING_COVARIATE", context,
                                "missing covariate values are not allowed for this path");
                }
            }
            report.pass(ValidationStage::PreRun, "COVARIATE_DIMENSIONS_VALID", context,
                        "covariate table dimensions are valid");
        }
    }
}

inline void validate_matrix_dimensions(ValidationReport& report,
                                       std::size_t rows,
                                       std::size_t cols,
                                       std::size_t expected_rows,
                                       const std::string& context,
                                       bool require_square = true) {
    if (require_square && rows != cols) {
        report.fail(ValidationStage::PreRun, "MATRIX_NOT_SQUARE", context,
                    "matrix must be square for this model");
    }
    if (rows != expected_rows || cols != expected_rows) {
        report.fail(ValidationStage::PreRun, "MATRIX_DIMENSION_MISMATCH", context,
                    "matrix dimensions do not match SampleRegistry size");
    } else {
        report.pass(ValidationStage::PreRun, "MATRIX_DIMENSION_MATCH", context,
                    "matrix dimensions match SampleRegistry size");
    }
}

inline void validate_eigen_vector_finite(ValidationReport& report,
                                         const Eigen::Ref<const Eigen::VectorXd>& values,
                                         const std::string& code,
                                         const std::string& context,
                                         ValidationStage stage = ValidationStage::PostRun) {
    for (Eigen::Index i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values(i))) {
            report.fail(stage, code, context, "vector contains NaN or Inf");
            return;
        }
    }
    report.pass(stage, code + "_FINITE", context, "vector values are finite");
}

inline void validate_variance_components(ValidationReport& report,
                                         const std::vector<double>& variances,
                                         const std::string& context,
                                         ValidationStage stage = ValidationStage::PostRun) {
    for (double v : variances) {
        if (!std::isfinite(v)) {
            report.fail(stage, "VARIANCE_NONFINITE", context, "variance component is NaN or Inf");
            return;
        }
        if (v < 0.0) {
            report.fail(stage, "VARIANCE_NEGATIVE", context, "variance component is negative");
            return;
        }
    }
    report.pass(stage, "VARIANCE_COMPONENTS_VALID", context, "variance components are finite and non-negative");
}

inline void validate_probability_value(ValidationReport& report,
                                       double value,
                                       const std::string& code,
                                       const std::string& context,
                                       ValidationStage stage = ValidationStage::PostRun) {
    if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
        report.fail(stage, code, context, "probability value is non-finite or outside [0,1]");
    }
}

inline void validate_finite_scalar(ValidationReport& report,
                                   double value,
                                   const std::string& code,
                                   const std::string& context,
                                   ValidationStage stage = ValidationStage::PostRun) {
    if (!std::isfinite(value)) {
        report.fail(stage, code, context, "scalar value is NaN or Inf");
    }
}

inline void validate_nonnegative_finite_scalar(ValidationReport& report,
                                               double value,
                                               const std::string& code,
                                               const std::string& context,
                                               ValidationStage stage = ValidationStage::PostRun) {
    if (!std::isfinite(value) || value < 0.0) {
        report.fail(stage, code, context, "scalar value is non-finite or negative");
    }
}

} // namespace cosmic
