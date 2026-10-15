#include "association/generalized/logistic.h"

#include <Eigen/Cholesky>
#include <Eigen/QR>
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

bool solveInformation(const Eigen::MatrixXd& design, const Eigen::VectorXd& weight,
                      const Eigen::VectorXd& score, Eigen::VectorXd& step,
                      Eigen::MatrixXd* inverse = nullptr) {
    const Eigen::MatrixXd information = design.transpose() * weight.asDiagonal() * design;
    Eigen::LDLT<Eigen::MatrixXd> solver(information);
    if (solver.info() != Eigen::Success || !solver.isPositive()) return false;
    step = solver.solve(score);
    if (!step.allFinite()) return false;
    if (inverse) {
        *inverse = solver.solve(Eigen::MatrixXd::Identity(information.rows(), information.cols()));
        if (!inverse->allFinite()) return false;
    }
    return true;
}

double logLikelihood(const Eigen::VectorXd& y, const Eigen::MatrixXd& design,
                     const Eigen::VectorXd& beta) {
    const Eigen::VectorXd eta = design * beta;
    double value = 0.0;
    for (int i = 0; i < eta.size(); ++i) {
        const double softplus = eta[i] > 0.0 ? eta[i] + std::log1p(std::exp(-eta[i]))
                                             : std::log1p(std::exp(eta[i]));
        value += y[i] * eta[i] - softplus;
    }
    return value;
}

bool fitLogistic(const Eigen::VectorXd& y, const Eigen::MatrixXd& design,
                 const LogisticOptions& options, Eigen::VectorXd& beta,
                 Eigen::VectorXd& probability, Eigen::VectorXd& weight,
                 Eigen::MatrixXd& information_inverse, int& iterations,
                 bool& separated) {
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> rank(design);
    if (design.rows() <= design.cols() || rank.rank() != design.cols()) return false;
    beta = Eigen::VectorXd::Zero(design.cols());
    if ((design.col(0).array() - 1.0).abs().maxCoeff() < 1e-12) {
        const double cases = y.sum();
        beta[0] = std::log(cases / (y.size() - cases));
    }
    double objective = logLikelihood(y, design, beta);
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        const Eigen::VectorXd eta = design * beta;
        probability.resize(y.size());
        for (int i = 0; i < y.size(); ++i) probability[i] = logistic(eta[i]);
        weight = (probability.array() * (1.0 - probability.array())).max(1e-15).matrix();
        Eigen::VectorXd step;
        if (!solveInformation(design, weight, design.transpose() * (y - probability), step)) {
            separated = eta.cwiseAbs().maxCoeff() > 20.0;
            return false;
        }
        const double largest = step.cwiseAbs().maxCoeff();
        if (largest > options.maximum_step) step *= options.maximum_step / largest;
        double scale = 1.0;
        Eigen::VectorXd candidate;
        double candidate_objective = -std::numeric_limits<double>::infinity();
        while (scale >= std::ldexp(1.0, -20)) {
            candidate = beta + scale * step;
            candidate_objective = logLikelihood(y, design, candidate);
            if (candidate_objective >= objective - 1e-12) break;
            scale *= 0.5;
        }
        if (!std::isfinite(candidate_objective)) return false;
        beta = candidate;
        objective = candidate_objective;
        iterations = iteration;
        if ((scale * step).cwiseAbs().maxCoeff() < options.tolerance) {
            const Eigen::VectorXd final_eta = design * beta;
            for (int i = 0; i < y.size(); ++i) probability[i] = logistic(final_eta[i]);
            weight = (probability.array() * (1.0 - probability.array())).max(1e-15).matrix();
            Eigen::VectorXd unused;
            return solveInformation(design, weight, Eigen::VectorXd::Zero(design.cols()),
                                    unused, &information_inverse);
        }
    }
    separated = (design * beta).cwiseAbs().maxCoeff() > 30.0;
    return false;
}

} // namespace

LogisticTest::LogisticTest(LogisticOptions options) : options_(options) {
    if (options_.max_iterations < 1 || !(options_.tolerance > 0.0) ||
        !(options_.maximum_step > 0.0) || !std::isfinite(options_.tolerance) ||
        !std::isfinite(options_.maximum_step))
        throw std::invalid_argument("Invalid logistic options");
}

void LogisticTest::prepare(const ScanContext& context) {
    if (context.phenotypes.cols() != 1 || context.covariates.rows() != context.phenotypes.rows() ||
        context.covariates.cols() < 1 || !context.phenotypes.allFinite() ||
        !context.covariates.allFinite() ||
        context.observation_to_sample.size() != static_cast<size_t>(context.phenotypes.rows()))
        throw std::invalid_argument("Logistic regression requires one finite aligned phenotype and covariates");
    int cases = 0;
    for (int i = 0; i < context.phenotypes.rows(); ++i) {
        const double value = context.phenotypes(i, 0);
        if (value != 0.0 && value != 1.0)
            throw std::invalid_argument("Logistic phenotype must be coded exactly 0/1");
        cases += value == 1.0;
    }
    if (cases == 0 || cases == context.phenotypes.rows())
        throw std::invalid_argument("Logistic regression requires both cases and controls");
    phenotype_ = context.phenotypes.col(0);
    covariates_ = context.covariates;
    observation_to_sample_ = context.observation_to_sample;
    int iterations = 0;
    bool separated = false;
    if (!fitLogistic(phenotype_, covariates_, options_, null_beta_, null_probability_,
                     null_weight_, null_information_inverse_, iterations, separated))
        throw std::runtime_error("Logistic null model did not converge or covariates separate the phenotype");
    null_residual_ = phenotype_ - null_probability_;
}

LogisticResult LogisticTest::fitWald(const Eigen::MatrixXd& design) const {
    LogisticResult result;
    Eigen::VectorXd beta, probability, weight;
    Eigen::MatrixXd inverse;
    if (!fitLogistic(phenotype_, design, options_, beta, probability, weight, inverse,
                     result.iterations, result.separated)) return result;
    const int marker = design.cols() - 1;
    result.beta = beta[marker];
    result.se = std::sqrt(inverse(marker, marker));
    result.statistic = result.beta * result.beta / (result.se * result.se);
    result.p_value = std::erfc(std::sqrt(result.statistic / 2.0));
    result.converged = std::isfinite(result.beta) && std::isfinite(result.se) && result.se > 0.0 &&
                       std::isfinite(result.p_value);
    return result;
}

std::vector<AssociationResult> LogisticTest::scanBlock(const GenotypeBlock& block) {
    if (null_probability_.size() == 0) throw std::logic_error("Logistic engine was not prepared");
    if (block.dosages.cols() != static_cast<int>(block.variant_ids.size()) || !block.dosages.allFinite())
        throw std::invalid_argument("Invalid logistic genotype block");
    results_.clear();
    std::vector<AssociationResult> output;
    for (int marker = 0; marker < block.dosages.cols(); ++marker) {
        Eigen::VectorXd genotype(phenotype_.size());
        for (int i = 0; i < genotype.size(); ++i) {
            const Eigen::Index row = observation_to_sample_[i];
            if (row < 0 || row >= block.dosages.rows())
                throw std::invalid_argument("Logistic genotype sample index is out of range");
            genotype[i] = block.dosages(row, marker);
        }
        LogisticResult tested;
        if (options_.test == LogisticTestType::WALD) {
            Eigen::MatrixXd design(phenotype_.size(), covariates_.cols() + 1);
            design.leftCols(covariates_.cols()) = covariates_;
            design.col(design.cols() - 1) = genotype;
            tested = fitWald(design);
        } else {
            const Eigen::VectorXd coefficient = null_information_inverse_ * covariates_.transpose() *
                (null_weight_.array() * genotype.array()).matrix();
            const Eigen::VectorXd adjusted = genotype - covariates_ * coefficient;
            const double variance = (null_weight_.array() * adjusted.array().square()).sum();
            if (variance > 1e-12 && std::isfinite(variance)) {
                const double score = genotype.dot(null_residual_);
                tested.beta = score / variance;
                tested.se = 1.0 / std::sqrt(variance);
                tested.statistic = score * score / variance;
                tested.p_value = std::erfc(std::sqrt(tested.statistic / 2.0));
                tested.converged = std::isfinite(tested.p_value);
            }
        }
        AssociationResult generic;
        generic.test_id = block.variant_ids[marker];
        generic.n_obs = phenotype_.size();
        if (tested.converged) {
            generic.status = AssociationStatus::SUCCESS;
            generic.effects = Eigen::VectorXd::Constant(1, tested.beta);
            generic.effect_covariance = Eigen::MatrixXd::Constant(1, 1, tested.se * tested.se);
            generic.statistic = tested.statistic;
            generic.degrees_of_freedom = 1.0;
            generic.p_value = tested.p_value;
        } else {
            generic.status = AssociationStatus::NUMERICAL_FAILURE;
            generic.message = tested.separated ? "Complete or quasi-complete separation" :
                                                "Logistic fit did not converge or design was singular";
        }
        results_.push_back(tested);
        output.push_back(std::move(generic));
    }
    return output;
}

} // namespace cosmic::association
