#include "bgen_ad_reml.h"
#include <array>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace cosmic::detail {
namespace {

using Vector = Eigen::Vector3d;
using Matrix = Eigen::Matrix3d;
using Face = std::array<bool, 3>;

struct Evaluation {
    double nll = std::numeric_limits<double>::infinity();
    Vector gradient = Vector::Zero();
    Matrix fisher = Matrix::Zero();
    Matrix hessian = Matrix::Zero();
    bool valid = false;
};

struct Problem {
    std::array<Eigen::MatrixXd, 3> kernel;
    Eigen::VectorXd y;
    Vector conversion;
    double gradient_tolerance;
    double likelihood_tolerance;

    Evaluation evaluate(const Vector& q) const {
        Evaluation result;
        Eigen::MatrixXd covariance = q[0] * kernel[0] + q[1] * kernel[1] + q[2] * kernel[2];
        Eigen::LLT<Eigen::MatrixXd> factor(covariance);
        if (factor.info() != Eigen::Success) return result;
        const Eigen::VectorXd inverse_y = factor.solve(y);
        const Eigen::MatrixXd inverse = factor.solve(Eigen::MatrixXd::Identity(y.size(), y.size()));
        const Eigen::VectorXd diagonal = factor.matrixL().toDenseMatrix().diagonal();
        result.nll = diagonal.array().log().sum() + 0.5 * y.dot(inverse_y);
        std::array<Eigen::MatrixXd, 3> inverse_k;
        std::array<Eigen::VectorXd, 3> k_y;
        for (int i = 0; i < 3; ++i) {
            inverse_k[i].noalias() = inverse * kernel[i];
            k_y[i].noalias() = kernel[i] * inverse_y;
            result.gradient[i] = 0.5 * (inverse_k[i].trace() - inverse_y.dot(k_y[i]));
        }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j <= i; ++j) {
                const double trace = (inverse_k[i].array() * inverse_k[j].transpose().array()).sum();
                result.fisher(i, j) = result.fisher(j, i) = 0.5 * trace;
                result.hessian(i, j) = result.hessian(j, i) = k_y[i].dot(inverse * k_y[j]) - 0.5 * trace;
            }
        }
        result.valid = std::isfinite(result.nll) && result.gradient.allFinite() &&
                       result.fisher.allFinite() && result.hessian.allFinite();
        return result;
    }

    double projectedNorm(const Vector& q, const Vector& g, const Face& face) const {
        double norm = 0.0;
        for (int i = 0; i < 3; ++i) {
            if (!face[i]) continue;
            // Residual variance is strictly positive: its score must vanish.
            const double value = (i < 2 && q[i] == 0.0) ? std::min(0.0, g[i]) : g[i];
            norm = std::max(norm, std::abs(value));
        }
        return norm / static_cast<double>(y.size());
    }
};

std::string boundaryStatus(const Vector& q) {
    if (q[0] == 0.0 && q[1] == 0.0) return "AD_ZERO";
    if (q[0] == 0.0) return "A_ZERO";
    if (q[1] == 0.0) return "D_ZERO";
    return "INTERIOR";
}

BgenAdRemlCandidate optimize(const Problem& problem, const Face& face,
                            Vector q, const std::string& model, const std::string& start,
                            int max_iterations) {
    constexpr double residual_floor = 1e-10;
    for (int i = 0; i < 3; ++i) q[i] = face[i] ? std::max(q[i], i == 2 ? residual_floor : 0.0) : 0.0;
    BgenAdRemlCandidate result;
    result.model = model;
    result.start = start;
    double delta = std::numeric_limits<double>::infinity();
    Evaluation current = problem.evaluate(q);
    if (model == "AD_ZERO") {
        // Closed-form residual-only REML: q_e=1 because the numerical scale
        // equals y_p'y_p/r. It still needs full-domain KKT validation below.
        q << 0.0, 0.0, 1.0;
        current = problem.evaluate(q);
        delta = 0.0;
    } else {
        for (int iteration = 0; iteration < max_iterations && current.valid; ++iteration) {
            if (problem.projectedNorm(q, current.gradient, face) <= problem.gradient_tolerance &&
                delta <= problem.likelihood_tolerance * std::max(1.0, std::abs(current.nll))) break;
            std::vector<int> active;
            for (int i = 0; i < 3; ++i) {
                if (face[i] && (i == 2 || q[i] > 0.0 || current.gradient[i] < 0.0)) active.push_back(i);
            }
            if (active.empty()) break;
            const auto direction = [&](const Matrix& information) {
                Eigen::MatrixXd h(active.size(), active.size());
                Eigen::VectorXd g(active.size());
                for (size_t i = 0; i < active.size(); ++i) {
                    g[i] = current.gradient[active[i]];
                    for (size_t j = 0; j < active.size(); ++j) h(i, j) = information(active[i], active[j]);
                }
                Eigen::LLT<Eigen::MatrixXd> decomposition(h);
                Vector step = Vector::Zero();
                if (decomposition.info() == Eigen::Success) {
                    const Eigen::VectorXd solution = decomposition.solve(-g);
                    for (size_t i = 0; i < active.size(); ++i) step[active[i]] = solution[i];
                }
                return step;
            };
            Vector step = direction(current.hessian);
            if (!step.allFinite() || !(step.dot(current.gradient) < 0.0)) step = direction(current.fisher);
            if (!step.allFinite() || !(step.dot(current.gradient) < 0.0)) {
                for (int i : active) step[i] = -current.gradient[i] / std::max(current.fisher(i, i), 1e-8);
            }
            const double maximum_step = step.cwiseAbs().maxCoeff();
            if (maximum_step > 5.0) step *= 5.0 / maximum_step;
            bool accepted = false;
            Vector next = q;
            Evaluation trial;
            for (int retry = 0; retry < 2 && !accepted; ++retry) {
                double length = 1.0;
                for (int backtrack = 0; backtrack < 45; ++backtrack) {
                    next = q + length * step;
                    for (int i = 0; i < 3; ++i)
                        next[i] = face[i] ? std::max(next[i], i == 2 ? residual_floor : 0.0) : 0.0;
                    const double descent = current.gradient.dot(next - q);
                    trial = problem.evaluate(next);
                    const double roundoff = 1e-13 * std::max(1.0, std::abs(current.nll));
                    if (descent <= 0.0 && trial.valid &&
                        trial.nll <= current.nll + 1e-4 * descent + roundoff) {
                        accepted = true;
                        break;
                    }
                    length *= 0.5;
                }
                if (!accepted) {
                    step.setZero();
                    for (int i : active) step[i] = -current.gradient[i] / std::max(current.fisher(i, i), 1e-8);
                    const double size = step.cwiseAbs().maxCoeff();
                    if (size > 1.0) step /= size;
                }
            }
            if (!accepted) break;
            delta = std::abs(current.nll - trial.nll);
            q = next;
            current = std::move(trial);
            result.iterations = iteration + 1;
        }
    }
    const Vector theta = q.cwiseProduct(problem.conversion);
    result.sigma_a2 = theta[0];
    result.sigma_d2 = theta[1];
    result.sigma_e2 = theta[2];
    result.log_likelihood = -current.nll;
    result.likelihood_delta = delta;
    result.boundary_status = boundaryStatus(q);
    result.gradient_norm = problem.projectedNorm(q, current.gradient, face);
    result.kkt_residual = problem.projectedNorm(q, current.gradient, {true, true, true});
    result.raw_gradient_norm = current.gradient.cwiseQuotient(problem.conversion).cwiseAbs().maxCoeff();
    const bool likelihood_ok = delta <= problem.likelihood_tolerance * std::max(1.0, std::abs(current.nll));
    const bool residual_interior = q[2] > residual_floor * 1.01;
    result.face_converged = max_iterations > 0 && current.valid && residual_interior && likelihood_ok &&
                            result.gradient_norm <= problem.gradient_tolerance;
    result.kkt_valid = current.valid && residual_interior &&
                       result.kkt_residual <= problem.gradient_tolerance;
    return result;
}

} // namespace

BgenAdNullResult fitBgenAdReml(const BgenAdKernelPair& kernels,
                             const Eigen::VectorXd& y,
                             const Eigen::MatrixXf& covariates,
                             const EigenVceAD::Options& options) {
    const Eigen::Index residual_dimension = y.size() - covariates.cols();
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(covariates.cast<double>());
    const Eigen::MatrixXd orthogonal = qr.householderQ() * Eigen::MatrixXd::Identity(y.size(), y.size());
    const Eigen::MatrixXd contrast = orthogonal.rightCols(residual_dimension);
    Problem problem;
    problem.y = contrast.transpose() * y;
    problem.kernel[0] = contrast.transpose() * kernels.additive * contrast;
    problem.kernel[1] = contrast.transpose() * kernels.dominance * contrast;
    problem.kernel[2] = Eigen::MatrixXd::Identity(residual_dimension, residual_dimension);
    const double scale = problem.y.squaredNorm() / static_cast<double>(residual_dimension);
    if (!(scale > 1e-20) || !std::isfinite(scale)) {
        throw std::invalid_argument("BGEN AD REML: degenerate phenotype residual variance");
    }
    // Exact change of optimization coordinates, NOT a change of KA/KD or
    // their scientific normalization. theta_i = conversion_i * q_i.
    for (int i = 0; i < 3; ++i) {
        problem.kernel[i] = (0.5 * (problem.kernel[i] + problem.kernel[i].transpose())).eval();
        const double mean_diagonal = problem.kernel[i].trace() / static_cast<double>(residual_dimension);
        if (mean_diagonal < -1e-12) throw std::invalid_argument("BGEN AD REML: negative kernel trace");
        problem.conversion[i] = scale / std::max(mean_diagonal, 1e-12);
        problem.kernel[i] *= problem.conversion[i];
    }
    problem.gradient_tolerance = std::min(options.tol, 1e-7);
    problem.likelihood_tolerance = std::min(options.tol, 1e-9);

    // Original optimizer is only a warm start. Its convergence Boolean cannot
    // certify a BGEN null; every refined candidate is rechecked independently.
    EigenVceAD::Options legacy_options = options;
    legacy_options.verbose = false;
    EigenVceAD legacy(legacy_options);
    legacy.prepare(y, covariates, kernels.additive, kernels.dominance);
    if (!legacy.runNullModel()) throw std::runtime_error("BGEN AD REML: legacy initialization failed");
    const Vector legacy_theta(legacy.getSigma2a(), legacy.getSigma2d(), legacy.getSigma2e());
    const Vector legacy_q = legacy_theta.cwiseQuotient(problem.conversion);
    const std::array<std::pair<std::string, Vector>, 5> starts = {{
        {"balanced", Vector(0.333, 0.333, 0.334)},
        {"additive_dominated", Vector(0.85, 0.05, 0.10)},
        {"dominance_dominated", Vector(0.05, 0.85, 0.10)},
        {"near_boundary", Vector(1e-4, 1e-4, 0.9998)},
        {"legacy_refined", legacy_q}
    }};
    BgenAdNullResult result;
    result.gradient_tolerance = problem.gradient_tolerance;
    result.likelihood_tolerance = problem.likelihood_tolerance;
    const std::array<std::pair<std::string, Face>, 3> faces = {{
        {"FULL", {true, true, true}},
        {"A_ZERO", {false, true, true}},
        {"D_ZERO", {true, false, true}}
    }};
    for (const auto& model : faces) {
        for (const auto& start : starts) {
            Vector initial = start.second;
            for (int i = 0; i < 3; ++i) if (!model.second[i]) initial[i] = 0.0;
            if (model.first != "FULL") {
                const double sum = initial.sum();
                if (sum > 0.0) initial /= sum;
            }
            result.candidates.push_back(optimize(problem, model.second, initial,
                                                model.first, start.first, options.max_em_iter));
        }
    }
    result.candidates.push_back(optimize(problem, {false, false, true}, Vector(0.0, 0.0, 1.0),
                                        "AD_ZERO", "analytic_residual", options.max_em_iter));

    const BgenAdRemlCandidate* best_valid = nullptr;
    const BgenAdRemlCandidate* best_any = nullptr;
    for (const auto& candidate : result.candidates) {
        if (!std::isfinite(candidate.log_likelihood)) continue;
        if (!best_any || candidate.log_likelihood > best_any->log_likelihood) best_any = &candidate;
        if (candidate.face_converged && candidate.kkt_valid &&
            (!best_valid || candidate.log_likelihood > best_valid->log_likelihood)) best_valid = &candidate;
    }
    const auto* selected = best_valid ? best_valid : best_any;
    if (!selected) throw std::runtime_error("BGEN AD REML: no finite feasible null candidate");
    result.sigma_a2 = selected->sigma_a2;
    result.sigma_d2 = selected->sigma_d2;
    result.sigma_e2 = selected->sigma_e2;
    result.log_likelihood = selected->log_likelihood;
    result.iterations = selected->iterations;
    result.likelihood_delta = selected->likelihood_delta;
    result.gradient_norm = selected->gradient_norm;
    result.kkt_residual = selected->kkt_residual;
    result.raw_gradient_norm = selected->raw_gradient_norm;
    result.boundary_status = selected->boundary_status;
    result.selected_start = selected->start;
    result.converged = best_valid != nullptr;
    // A materially higher uncertified candidate prevents a false positive.
    if (best_valid && best_any->log_likelihood > best_valid->log_likelihood +
        problem.likelihood_tolerance * std::max(1.0, std::abs(best_valid->log_likelihood))) result.converged = false;

    const Vector fitted_q = Vector(result.sigma_a2, result.sigma_d2, result.sigma_e2).cwiseQuotient(problem.conversion);
    const Evaluation fitted = problem.evaluate(fitted_q);
    Eigen::SelfAdjointEigenSolver<Matrix> information(fitted.fisher);
    if (information.info() == Eigen::Success) {
        const double cutoff = 1e-10 * std::max(1.0, information.eigenvalues().maxCoeff());
        result.information_rank = (information.eigenvalues().array() > cutoff).count();
        if (result.information_rank == 3) {
            const Matrix inverse_information = information.eigenvectors() *
                information.eigenvalues().cwiseInverse().asDiagonal() * information.eigenvectors().transpose();
            result.sigma_d2_local_se = problem.conversion[1] * std::sqrt(std::max(0.0, inverse_information(1, 1)));
        }
    }
    return result;
}

} // namespace cosmic::detail
