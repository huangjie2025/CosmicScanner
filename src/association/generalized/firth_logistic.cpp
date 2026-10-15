#include "association/generalized/firth_logistic.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cosmic::association {
namespace {

double logistic(double value) {
    if (value >= 0.0) {
        const double z = std::exp(-value);
        return 1.0 / (1.0 + z);
    }
    const double z = std::exp(value);
    return z / (1.0 + z);
}

double logOnePlusExp(double value) {
    return value > 0.0 ? value + std::log1p(std::exp(-value))
                       : std::log1p(std::exp(value));
}

bool information(const Eigen::MatrixXd& design, const Eigen::VectorXd& beta,
                 Eigen::VectorXd& probability, Eigen::VectorXd& weight,
                 Eigen::MatrixXd& inverse, double& log_determinant) {
    const Eigen::VectorXd eta = design * beta;
    for (int i = 0; i < eta.size(); ++i) probability[i] = logistic(eta[i]);
    weight = (probability.array() * (1.0 - probability.array())).max(1e-12);
    const Eigen::MatrixXd matrix = design.transpose() * weight.asDiagonal() * design;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(matrix);
    if (eigen.info() != Eigen::Success || eigen.eigenvalues().minCoeff() <= 1e-12) return false;
    inverse = eigen.operatorInverseSqrt() * eigen.operatorInverseSqrt();
    log_determinant = eigen.eigenvalues().array().log().sum();
    return inverse.allFinite() && std::isfinite(log_determinant);
}

double penalizedLikelihood(const Eigen::VectorXd& y, const Eigen::MatrixXd& design,
                           const Eigen::VectorXd& beta) {
    Eigen::VectorXd probability(y.size()), weight(y.size());
    Eigen::MatrixXd inverse;
    double log_determinant = 0.0;
    if (!information(design, beta, probability, weight, inverse, log_determinant))
        return -std::numeric_limits<double>::infinity();
    const Eigen::VectorXd eta = design * beta;
    double likelihood = 0.5 * log_determinant;
    for (int i = 0; i < y.size(); ++i) likelihood += y[i] * eta[i] - logOnePlusExp(eta[i]);
    return likelihood;
}

FirthResult fit(const Eigen::VectorXd& y, const Eigen::MatrixXd& design,
                const FirthOptions& options) {
    FirthResult result;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> rank_check(design);
    if (rank_check.rank() != design.cols() || design.rows() <= design.cols()) return result;
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(design.cols());
    const double cases = y.sum();
    if ((design.col(0).array() - 1.0).abs().maxCoeff() < 1e-12)
        beta[0] = std::log((cases + 0.5) / (y.size() - cases + 0.5));
    double objective = penalizedLikelihood(y, design, beta);
    Eigen::VectorXd probability(y.size()), weight(y.size());
    Eigen::MatrixXd inverse;
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        double log_determinant = 0.0;
        if (!information(design, beta, probability, weight, inverse, log_determinant)) return result;
        Eigen::VectorXd leverage(y.size());
        for (int i = 0; i < y.size(); ++i)
            leverage[i] = weight[i] * (design.row(i) * inverse * design.row(i).transpose())(0, 0);
        const Eigen::VectorXd adjusted = y - probability +
            (leverage.array() * (0.5 - probability.array())).matrix();
        Eigen::VectorXd step = inverse * design.transpose() * adjusted;
        const double largest = step.cwiseAbs().maxCoeff();
        if (largest > options.maximum_step) step *= options.maximum_step / largest;
        double scale = 1.0;
        double candidate_objective = -std::numeric_limits<double>::infinity();
        Eigen::VectorXd candidate;
        while (scale >= std::ldexp(1.0, -20)) {
            candidate = beta + scale * step;
            candidate_objective = penalizedLikelihood(y, design, candidate);
            if (candidate_objective >= objective - 1e-12) break;
            scale *= 0.5;
        }
        if (!std::isfinite(candidate_objective)) return result;
        beta = candidate;
        objective = candidate_objective;
        result.iterations = iteration;
        if ((scale * step).cwiseAbs().maxCoeff() < options.tolerance) {
            result.converged = true;
            break;
        }
    }
    double log_determinant = 0.0;
    if (!result.converged || !information(design, beta, probability, weight, inverse, log_determinant))
        return result;
    const int marker = design.cols() - 1;
    result.beta = beta[marker];
    result.se = std::sqrt(inverse(marker, marker));
    const double z = result.beta / result.se;
    result.p_value = std::erfc(std::abs(z) / std::sqrt(2.0));
    if (!std::isfinite(result.beta) || !std::isfinite(result.se) || result.se <= 0.0 ||
        !std::isfinite(result.p_value)) result.converged = false;
    return result;
}

} // namespace

FirthLogisticTest::FirthLogisticTest(FirthOptions options) : options_(options) {
    if (options_.max_iterations < 1 || !std::isfinite(options_.tolerance) || options_.tolerance <= 0.0 ||
        !std::isfinite(options_.maximum_step) || options_.maximum_step <= 0.0)
        throw std::invalid_argument("Invalid Firth logistic options");
}

void FirthLogisticTest::prepare(const ScanContext& context) {
    if (context.phenotypes.cols() != 1 || context.covariates.rows() != context.phenotypes.rows() ||
        !context.phenotypes.allFinite() || !context.covariates.allFinite())
        throw std::invalid_argument("Firth logistic requires one finite aligned phenotype and covariates");
    int cases = 0;
    for (int i = 0; i < context.phenotypes.rows(); ++i) {
        const double value = context.phenotypes(i, 0);
        if (value != 0.0 && value != 1.0) throw std::invalid_argument("Firth logistic phenotype must be coded 0/1");
        cases += value == 1.0;
    }
    if (cases == 0 || cases == context.phenotypes.rows())
        throw std::invalid_argument("Firth logistic requires both cases and controls");
    context_ = &context;
}

std::vector<AssociationResult> FirthLogisticTest::scanBlock(const GenotypeBlock& block) {
    if (!context_) throw std::logic_error("Firth logistic engine was not prepared");
    if (block.dosages.cols() != static_cast<int>(block.variant_ids.size()) || !block.dosages.allFinite() ||
        static_cast<int>(context_->observation_to_sample.size()) != context_->phenotypes.rows())
        throw std::invalid_argument("Invalid Firth logistic genotype block");
    const int n = context_->phenotypes.rows();
    results_.clear();
    std::vector<AssociationResult> output;
    results_.reserve(block.dosages.cols());
    output.reserve(block.dosages.cols());
    for (int marker = 0; marker < block.dosages.cols(); ++marker) {
        Eigen::MatrixXd design(n, context_->covariates.cols() + 1);
        design.leftCols(context_->covariates.cols()) = context_->covariates;
        for (int i = 0; i < n; ++i) {
            const Eigen::Index row = context_->observation_to_sample[i];
            if (row < 0 || row >= block.dosages.rows())
                throw std::invalid_argument("Firth genotype sample index is out of range");
            design(i, design.cols() - 1) = block.dosages(row, marker);
        }
        FirthResult fit_result = fit(context_->phenotypes.col(0), design, options_);
        AssociationResult association;
        association.test_id = block.variant_ids[marker];
        association.n_obs = n;
        if (fit_result.converged) {
            association.status = AssociationStatus::SUCCESS;
            association.effects = Eigen::VectorXd::Constant(1, fit_result.beta);
            association.effect_covariance = Eigen::MatrixXd::Constant(1, 1, fit_result.se * fit_result.se);
            association.statistic = fit_result.beta * fit_result.beta / (fit_result.se * fit_result.se);
            association.degrees_of_freedom = 1.0;
            association.p_value = fit_result.p_value;
        } else {
            association.status = AssociationStatus::NUMERICAL_FAILURE;
            association.message = "Firth iteration did not converge or design was singular";
        }
        results_.push_back(fit_result);
        output.push_back(std::move(association));
    }
    return output;
}

} // namespace cosmic::association
