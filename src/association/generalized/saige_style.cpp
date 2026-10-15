#include "association/generalized/saige_style.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
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

double normalP(double score, double variance) {
    if (!(variance > 0.0) || !std::isfinite(score) || !std::isfinite(variance)) return 1.0;
    return std::erfc(std::abs(score) / std::sqrt(2.0 * variance));
}

} // namespace

SaigeStyleNullModel::SaigeStyleNullModel(SaigeStyleOptions options) : options_(options) {
    if (options_.max_iterations < 1 || !std::isfinite(options_.tolerance) ||
        options_.tolerance <= 0.0 || !std::isfinite(options_.maximum_tau) ||
        options_.maximum_tau < 0.0)
        throw std::invalid_argument("Invalid SAIGE-style null-model options");
}

double SaigeStyleNullModel::estimateTau(const Eigen::VectorXd& pearson) const {
    double numerator = 0.0, denominator = 0.0;
    for (int col = 0; col < kinship_.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator entry(kinship_, col); entry; ++entry) {
            if (entry.row() <= entry.col()) continue;
            numerator += 2.0 * entry.value() * pearson[entry.row()] * pearson[entry.col()];
            denominator += 2.0 * entry.value() * entry.value();
        }
    }
    if (!(denominator > 1e-12) || !std::isfinite(numerator)) return 0.0;
    return std::min(options_.maximum_tau, std::max(0.0, numerator / denominator));
}

void SaigeStyleNullModel::factorizeWorkingCovariance() {
    working_covariance_ = tau_ * kinship_;
    working_covariance_.makeCompressed();
    for (int i = 0; i < working_covariance_.rows(); ++i)
        working_covariance_.coeffRef(i, i) += 1.0 / weight_[i];
    working_covariance_.makeCompressed();
    factor_.compute(working_covariance_);
    if (factor_.info() != Eigen::Success)
        throw std::runtime_error("SAIGE-style sparse working covariance factorization failed");
    a_inverse_x_ = factor_.solve(X_);
    Eigen::MatrixXd information = X_.transpose() * a_inverse_x_;
    Eigen::LDLT<Eigen::MatrixXd> fixed_solver(information);
    if (fixed_solver.info() != Eigen::Success || !fixed_solver.isPositive())
        throw std::runtime_error("SAIGE-style fixed-effect information matrix is singular");
    xt_a_inverse_x_inverse_ = fixed_solver.solve(
        Eigen::MatrixXd::Identity(information.rows(), information.cols()));
}

void SaigeStyleNullModel::fit(const Eigen::VectorXd& phenotype,
                              const Eigen::MatrixXd& covariates,
                              const Eigen::SparseMatrix<double>& kinship) {
    if (phenotype.size() < 10 || covariates.rows() != phenotype.size() || covariates.cols() < 1 ||
        kinship.rows() != phenotype.size() || kinship.cols() != phenotype.size() ||
        !phenotype.allFinite() || !covariates.allFinite())
        throw std::invalid_argument("SAIGE-style null model requires aligned finite phenotype, covariates, and GRM");
    int cases = 0;
    for (int i = 0; i < phenotype.size(); ++i) {
        if (phenotype[i] != 0.0 && phenotype[i] != 1.0)
            throw std::invalid_argument("SAIGE-style phenotype must be coded 0/1");
        cases += phenotype[i] == 1.0;
    }
    if (cases == 0 || cases == phenotype.size())
        throw std::invalid_argument("SAIGE-style null model requires cases and controls");
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> rank(covariates);
    if (rank.rank() != covariates.cols())
        throw std::invalid_argument("SAIGE-style covariate matrix is rank deficient");

    y_ = phenotype;
    X_ = covariates;
    kinship_ = kinship;
    kinship_.makeCompressed();
    beta_ = Eigen::VectorXd::Zero(X_.cols());
    beta_[0] = std::log((cases + 0.5) / (phenotype.size() - cases + 0.5));
    Eigen::VectorXd eta = X_ * beta_;
    tau_ = 0.0;
    metadata_ = {};
    metadata_.sparse_nnz = static_cast<int>(kinship_.nonZeros());

    for (int iteration = 1; iteration <= options_.max_iterations; ++iteration) {
        mu_.resize(y_.size());
        for (int i = 0; i < eta.size(); ++i)
            mu_[i] = std::min(1.0 - 1e-6, std::max(1e-6, logistic(eta[i])));
        weight_ = (mu_.array() * (1.0 - mu_.array())).matrix();
        const Eigen::VectorXd pearson = (y_ - mu_).array() / weight_.array().sqrt();
        const double tau_candidate = estimateTau(pearson);
        const double previous_tau = tau_;
        const double next_tau = iteration == 1 ? tau_candidate : 0.5 * (tau_ + tau_candidate);
        const Eigen::VectorXd working_response = eta +
            ((y_ - mu_).array() / weight_.array()).matrix();
        tau_ = next_tau;
        factorizeWorkingCovariance();
        const Eigen::VectorXd a_inverse_z = factor_.solve(working_response);
        const Eigen::VectorXd next_beta = xt_a_inverse_x_inverse_ * X_.transpose() * a_inverse_z;
        const Eigen::VectorXd working_residual = working_response - X_ * next_beta;
        const Eigen::VectorXd random_effect = tau_ * (kinship_ * factor_.solve(working_residual));
        const Eigen::VectorXd next_eta = X_ * next_beta + random_effect;
        const double change = std::max({(next_eta - eta).cwiseAbs().maxCoeff(),
                                        (next_beta - beta_).cwiseAbs().maxCoeff(),
                                        std::abs(next_tau - previous_tau)});
        eta = next_eta;
        beta_ = next_beta;
        metadata_.iterations = iteration;
        if (change < options_.tolerance) {
            metadata_.converged = true;
            break;
        }
    }
    if (!metadata_.converged)
        throw std::runtime_error("SAIGE-style binary mixed null model did not converge");

    for (int i = 0; i < eta.size(); ++i)
        mu_[i] = std::min(1.0 - 1e-6, std::max(1e-6, logistic(eta[i])));
    weight_ = (mu_.array() * (1.0 - mu_.array())).matrix();
    residual_ = y_ - mu_;
    factorizeWorkingCovariance();
    const Eigen::MatrixXd xtwx = X_.transpose() * weight_.asDiagonal() * X_;
Eigen::LDLT<Eigen::MatrixXd> weighted_fixed(xtwx);
        if (weighted_fixed.info() != Eigen::Success || !weighted_fixed.isPositive())
        throw std::runtime_error("SAIGE-style weighted covariate projection is singular");
    xt_w_x_inverse_ = weighted_fixed.solve(Eigen::MatrixXd::Identity(X_.cols(), X_.cols()));
    if (!mu_.allFinite() || !beta_.allFinite())
        throw std::runtime_error("SAIGE-style null model produced non-finite fitted parameters");
}

double SaigeStyleNullModel::naiveVariance(const Eigen::VectorXd& marker,
                                          Eigen::VectorXd* adjusted) const {
    const Eigen::VectorXd coefficient = xt_w_x_inverse_ * X_.transpose() *
                                        (weight_.array() * marker.array()).matrix();
    const Eigen::VectorXd projected = marker - X_ * coefficient;
    if (adjusted) *adjusted = projected;
    return (weight_.array() * projected.array().square()).sum();
}

double SaigeStyleNullModel::exactVariance(const Eigen::VectorXd& marker) const {
    const Eigen::VectorXd a_inverse_marker = factor_.solve(marker);
    if (factor_.info() != Eigen::Success) return std::numeric_limits<double>::quiet_NaN();
    const Eigen::VectorXd fixed_score = X_.transpose() * a_inverse_marker;
    return marker.dot(a_inverse_marker) -
           fixed_score.dot(xt_a_inverse_x_inverse_ * fixed_score);
}

// MAC bin helpers
namespace {
int computeMac(const Eigen::VectorXd& marker, int n_samples) {
    double sum = marker.sum();
    int n_observed = 0;
    for (int i = 0; i < marker.size(); ++i) {
        if (std::isfinite(marker[i]) && marker[i] >= 0 && marker[i] <= 2) {
            ++n_observed;
        }
    }
    if (n_observed < 3) return 0;
    double freq = sum / (2.0 * n_observed);
    double maf = std::min(freq, 1.0 - freq);
    return static_cast<int>(2.0 * n_observed * maf);
}

int macBinIndex(int mac) {
    if (mac <= 10) return 0;
    if (mac <= 20) return 1;
    if (mac <= 50) return 2;
    return 3;
}
} // namespace

double SaigeStyleNullModel::calibrateVarianceRatio(const std::vector<Eigen::VectorXd>& markers) {
    // MAC bins: BIN1: MAC<=10, BIN2: 10<MAC<=20, BIN3: 20<MAC<=50, BIN4: MAC>50
    struct BinAccumulator {
        std::vector<double> ratios;
        int n_markers = 0;
        int min_mac = 0;
        int max_mac = 0;
    };
    std::array<BinAccumulator, 4> bins;
    bins[0].min_mac = 0; bins[0].max_mac = 10;
    bins[1].min_mac = 11; bins[1].max_mac = 20;
    bins[2].min_mac = 21; bins[2].max_mac = 50;
    bins[3].min_mac = 51; bins[3].max_mac = 0; // 0 = unlimited

    // Global fallback
    std::vector<double> global_ratios;
    global_ratios.reserve(markers.size());

    for (const auto& marker : markers) {
        if (marker.size() != y_.size() || !marker.allFinite()) continue;
        const double naive = naiveVariance(marker, nullptr);
        const double exact = exactVariance(marker);
        if (naive > 1e-10 && exact > 1e-10 && std::isfinite(exact)) {
            double ratio = exact / naive;
            int mac = computeMac(marker, static_cast<int>(y_.size()));
            int bin_idx = 0;
            if (mac <= 10) bin_idx = 0;
            else if (mac <= 20) bin_idx = 1;
            else if (mac <= 50) bin_idx = 2;
            else bin_idx = 3;

            bins[bin_idx].ratios.push_back(exact / naive);
            bins[bin_idx].n_markers++;
            global_ratios.push_back(exact / naive);
        }
    }

    if (global_ratios.empty()) {
        throw std::runtime_error("SAIGE-style variance-ratio calibration found no testable markers");
    }

    // Global median as fallback
    const size_t middle = global_ratios.size() / 2;
    std::nth_element(global_ratios.begin(), global_ratios.begin() + middle, global_ratios.end());
    double global_ratio = global_ratios[middle];
    if (global_ratios.size() % 2 == 0) {
        const auto lower = std::max_element(global_ratios.begin(), global_ratios.begin() + middle);
        global_ratio = 0.5 * (global_ratio + *lower);
    }
    if (!(global_ratio > 0.0) || !std::isfinite(global_ratio))
        throw std::runtime_error("SAIGE-style variance-ratio calibration was non-finite");

    // Compute per-bin ratios with fallback
    mac_bins_.clear();
    mac_bins_.reserve(4);
    for (int i = 0; i < 4; ++i) {
        MacVarianceRatioBin bin;
        bin.min_mac = bins[i].min_mac;
        bin.max_mac = bins[i].max_mac;
        bin.n_calibration_markers = bins[i].n_markers;
        bin.fallback = false;
        bin.fallback_source = "";

        if (bins[i].n_markers > 0) {
            std::vector<double>& ratios = bins[i].ratios;
            const size_t middle = ratios.size() / 2;
            std::nth_element(ratios.begin(), ratios.begin() + middle, ratios.end());
            bin.variance_ratio = ratios[middle];
            if (ratios.size() % 2 == 0) {
                const auto lower = std::max_element(ratios.begin(), ratios.begin() + middle);
                bin.variance_ratio = 0.5 * (bin.variance_ratio + *lower);
            }
            bin.fallback = false;
            bin.fallback_source = "";
        } else {
            // Empty bin: fallback to global
            bin.variance_ratio = global_ratio;
            bin.fallback = true;
            bin.fallback_source = "global";
        }
        if (!(bin.variance_ratio > 0.0) || !std::isfinite(bin.variance_ratio)) {
            throw std::runtime_error("SAIGE-style variance-ratio calibration was non-finite");
        }
        mac_bins_.push_back(bin);
    }

    // Set global ratio as the median of bin ratios (or global if no bins populated)
    variance_ratio_ = global_ratio;
    calibration_markers_ = static_cast<int>(global_ratios.size());
return variance_ratio_;
}

double SaigeStyleNullModel::varianceRatioForMac(int mac) const {
    if (mac_bins_.empty()) return variance_ratio_;
    int idx = 0;
    if (mac <= 10) idx = 0;
    else if (mac <= 20) idx = 1;
    else if (mac <= 50) idx = 2;
    else idx = 3;
    return mac_bins_[idx].variance_ratio;
}

SaigeStyleResult SaigeStyleNullModel::testMarker(const Eigen::VectorXd& marker, SPATest* spa) const {
    SaigeStyleResult result;
    if (marker.size() != y_.size() || !marker.allFinite()) return result;
    Eigen::VectorXd adjusted;
    const double naive = naiveVariance(marker, &adjusted);
    
    // Compute MAC and select appropriate variance ratio
    int mac = 0;
    {
        double sum = marker.sum();
        int n_observed = 0;
        for (int i = 0; i < marker.size(); ++i) {
            if (std::isfinite(marker[i]) && marker[i] >= 0 && marker[i] <= 2) {
                ++n_observed;
            }
        }
        if (n_observed >= 3) {
            double freq = marker.sum() / (2.0 * n_observed);
            double maf = std::min(freq, 1.0 - freq);
            mac = static_cast<int>(2.0 * n_observed * maf);
        }
    }
    double vr = varianceRatioForMac(mac);
    const double variance = vr * naive;
    if (!(variance > 1e-10) || !std::isfinite(variance)) return result;
    result.score = marker.dot(residual_);
    result.score_variance = variance;
    result.beta = result.score / variance;
    result.se = 1.0 / std::sqrt(variance);
    result.normal_p_value = normalP(result.score, variance);
    result.p_value = result.normal_p_value;
    if (spa && result.normal_p_value <= 0.05) {
        result.p_value = spa->compute_pvalue(adjusted, result.score / std::sqrt(vr), naive);
        result.spa_used = std::isfinite(result.p_value);
    }
    result.tested = std::isfinite(result.beta) && std::isfinite(result.se) && result.se > 0.0 &&
                    std::isfinite(result.p_value) && result.p_value >= 0.0 && result.p_value <= 1.0;
    return result;
}

SaigeStyleResult SaigeStyleNullModel::testMarkerExact(const Eigen::VectorXd& marker) const {
    SaigeStyleResult result;
    if (marker.size() != y_.size() || !marker.allFinite()) return result;
    const double variance = exactVariance(marker);
    if (!(variance > 1e-10) || !std::isfinite(variance)) return result;
    result.score = marker.dot(residual_);
    result.score_variance = variance;
    result.beta = result.score / variance;
    result.se = 1.0 / std::sqrt(variance);
    result.normal_p_value = normalP(result.score, variance);
    result.p_value = result.normal_p_value;
    result.tested = std::isfinite(result.beta) && std::isfinite(result.se) && result.se > 0.0 &&
                    std::isfinite(result.p_value);
    return result;
}

SaigeStyleTest::SaigeStyleTest(SaigeStyleOptions options) : null_model_(options) {}

SaigeStyleTest::SaigeStyleTest(Eigen::SparseMatrix<double> kinship, SaigeStyleOptions options)
    : kinship_(std::move(kinship)), null_model_(options) {}

void SaigeStyleTest::prepare(const ScanContext& context) {
    if (context.observation_to_sample.size() != static_cast<size_t>(context.phenotypes.rows()))
        throw std::invalid_argument("SAIGE-style observation-to-sample map is misaligned");
    if (kinship_.rows() == 0) {
        kinship_.resize(context.phenotypes.rows(), context.phenotypes.rows());
        kinship_.setIdentity();
    }
    observation_to_sample_ = context.observation_to_sample;
    null_model_.fit(context.phenotypes.col(0), context.covariates, kinship_);
}

std::vector<AssociationResult> SaigeStyleTest::scanBlock(const GenotypeBlock& block) {
    if (block.dosages.cols() != static_cast<int>(block.variant_ids.size()))
        throw std::invalid_argument("Invalid SAIGE-style genotype block");
    results_.clear();
    std::vector<AssociationResult> output;
    output.reserve(block.dosages.cols());
    results_.reserve(block.dosages.cols());
    for (int marker = 0; marker < block.dosages.cols(); ++marker) {
        Eigen::VectorXd dosage(observation_to_sample_.size());
        for (int i = 0; i < dosage.size(); ++i) {
            if (observation_to_sample_[i] < 0 || observation_to_sample_[i] >= block.dosages.rows())
                throw std::invalid_argument("SAIGE-style genotype sample index is out of range");
            dosage[i] = block.dosages(observation_to_sample_[i], marker);
        }
        SaigeStyleResult tested = null_model_.testMarker(dosage);
        AssociationResult generic;
        generic.test_id = block.variant_ids[marker];
        generic.n_obs = dosage.size();
        if (tested.tested) {
            generic.status = AssociationStatus::SUCCESS;
            generic.effects = Eigen::VectorXd::Constant(1, tested.beta);
            generic.effect_covariance = Eigen::MatrixXd::Constant(1, 1, tested.se * tested.se);
            generic.statistic = tested.score * tested.score / tested.score_variance;
            generic.degrees_of_freedom = 1.0;
            generic.p_value = tested.p_value;
        } else {
            generic.status = AssociationStatus::FILTERED;
            generic.message = "Marker is constant after null-model projection";
        }
        results_.push_back(tested);
        output.push_back(std::move(generic));
    }
    return output;
}

} // namespace cosmic::association
