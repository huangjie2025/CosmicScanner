#include "bgen_ad_scan_engine.h"
#include "bgen_dominance_encoding.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace cosmic {

BgenAdScanEngine::BgenAdScanEngine() : BgenAdScanEngine(Options{}) {}

BgenAdScanEngine::BgenAdScanEngine(const Options& options) : options_(options) {
    if (!std::isfinite(options.min_maf) || options.min_maf < 0.0 || options.min_maf > 0.5 ||
        !std::isfinite(options.max_missing) || options.max_missing < 0.0 || options.max_missing > 1.0 ||
        !std::isfinite(options.min_info_a) || options.min_info_a < 0.0 ||
        !std::isfinite(options.min_info_d) || options.min_info_d < 0.0 || options.block_size == 0) {
        throw std::invalid_argument("BGEN AD scan: invalid options");
    }
}

BgenAdPreparedMarker BgenAdScanEngine::encode(const BgenGenotypeProbs& probabilities,
                                            const std::vector<int>& sample_indices) {
    if (probabilities.p0.size() != probabilities.p1.size() ||
        probabilities.p1.size() != probabilities.p2.size() ||
        probabilities.p2.size() != probabilities.missing.size()) {
        throw std::invalid_argument("BGEN AD marker: probability vector dimensions");
    }
    std::vector<int> selected = sample_indices;
    if (selected.empty()) {
        selected.resize(probabilities.size());
        std::iota(selected.begin(), selected.end(), 0);
    }
    BgenAdPreparedMarker result;
    result.additive = Eigen::VectorXd::Zero(selected.size());
    result.dominance = Eigen::VectorXd::Zero(selected.size());
    double sum_x = 0.0, sum_h = 0.0;
    for (size_t i = 0; i < selected.size(); ++i) {
        const int source = selected[i];
        if (source < 0 || source >= probabilities.size()) throw std::out_of_range("BGEN AD marker sample index");
        const float missing = probabilities.missing[source];
        if (missing == 1.0f) continue;
        if (missing != 0.0f) throw std::invalid_argument("BGEN AD marker invalid missing mask");
        const double p0 = probabilities.p0[source], p1 = probabilities.p1[source], p2 = probabilities.p2[source];
        if (!std::isfinite(p0) || !std::isfinite(p1) || !std::isfinite(p2) ||
            p0 < 0.0 || p1 < 0.0 || p2 < 0.0 || p0 > 1.0 || p1 > 1.0 || p2 > 1.0 ||
            std::abs(p0 + p1 + p2 - 1.0) > 1e-6) {
            throw std::invalid_argument("BGEN AD marker invalid posterior triplet");
        }
        result.additive[i] = p1 + 2.0 * p2;
        result.dominance[i] = p1;
        sum_x += result.additive[i];
        sum_h += p1;
        ++result.observed;
    }
    if (selected.empty() || result.observed == 0) return result;
    result.missing_rate = static_cast<double>(selected.size() - result.observed) / selected.size();
    const double mean_x = sum_x / result.observed;
    const double mean_h = sum_h / result.observed;
    result.frequency = mean_x / 2.0;
    double additive_sq = 0.0, mu_sum = 0.0, mu_sq = 0.0;
    for (size_t i = 0; i < selected.size(); ++i) {
        const int source = selected[i];
        if (probabilities.missing[source] != 0.0f) continue;
        result.additive[i] -= mean_x;
        result.dominance[i] -= mean_h;
        additive_sq += result.additive[i] * result.additive[i];
        const double mu = bgen::posteriorMeanVitezica(probabilities.p0[source], probabilities.p1[source],
                                                     probabilities.p2[source], result.frequency);
        mu_sum += mu;
        mu_sq += mu * mu;
    }
    const double denominator_a = 2.0 * result.frequency * (1.0 - result.frequency);
    const double denominator_d = bgen::vitezicaDenominator(result.frequency);
    if (denominator_a > 0.0) result.info_a = additive_sq / (result.observed * denominator_a);
    if (denominator_d > 0.0) {
        const double mean_mu = mu_sum / result.observed;
        result.info_d = std::max(0.0, mu_sq / result.observed - mean_mu * mean_mu) / denominator_d;
    }
    return result;
}

void BgenAdScanEngine::fitNull(const BgenAdKernelResult& kernels, const Eigen::VectorXd& y,
                              const Eigen::MatrixXf& covariates, const EigenVceAD::Options& options,
                              const std::string& excluded_chromosome) {
    prepared_ = false;
    null_result_.reset();
    const auto pair = excluded_chromosome.empty() ? kernels.global() : kernels.loco(excluded_chromosome);
    const auto fit = BgenAdNullBackend::fit(pair, y, covariates, options);
    ++null_fit_count_;
    prepareCertified(kernels, y, covariates, fit, excluded_chromosome);
}

void BgenAdScanEngine::prepareCertified(const BgenAdKernelResult& kernels, const Eigen::VectorXd& y,
                                       const Eigen::MatrixXf& covariates, const BgenAdNullResult& fit,
                                       const std::string& excluded_chromosome) {
    prepared_ = false;
    null_result_.reset();
    if (!fit.converged || !std::isfinite(fit.kkt_residual) || !std::isfinite(fit.gradient_norm) ||
        !std::isfinite(fit.likelihood_delta) || !std::isfinite(fit.log_likelihood) ||
        !std::isfinite(fit.gradient_tolerance) || !std::isfinite(fit.likelihood_tolerance) ||
        fit.gradient_tolerance <= 0.0 || fit.likelihood_tolerance <= 0.0 ||
        fit.kkt_residual > fit.gradient_tolerance ||
        fit.gradient_norm > fit.gradient_tolerance ||
        fit.likelihood_delta > fit.likelihood_tolerance * std::max(1.0, std::abs(fit.log_likelihood))) {
        throw std::runtime_error("BGEN AD scan requires a certified robust null");
    }
    prepareFixed(kernels, y, covariates, fit.sigma_a2, fit.sigma_d2, fit.sigma_e2, excluded_chromosome);
    null_result_ = fit;
}

void BgenAdScanEngine::prepareFixed(const BgenAdKernelResult& kernels, const Eigen::VectorXd& y,
                                   const Eigen::MatrixXf& covariates, double sigma_a2,
                                   double sigma_d2, double sigma_e2,
                                   const std::string& excluded_chromosome) {
    prepared_ = false;
    null_result_.reset();
    if (!std::isfinite(sigma_a2) || !std::isfinite(sigma_d2) || !std::isfinite(sigma_e2) ||
        sigma_a2 < 0.0 || sigma_d2 < 0.0 || sigma_e2 <= 0.0 ||
        !y.allFinite() || !covariates.allFinite() || covariates.rows() != y.size() ||
        covariates.cols() < 1 || covariates.cols() >= y.size() ||
        kernels.sample_indices.size() != static_cast<size_t>(y.size()) ||
        kernels.sample_ids.size() != kernels.sample_indices.size()) {
        throw std::invalid_argument("BGEN AD fixed null contract");
    }
    const Eigen::MatrixXd design = covariates.cast<double>();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);
    if (qr.rank() != design.cols() ||
        (Eigen::VectorXd::Ones(y.size()) - design * qr.solve(Eigen::VectorXd::Ones(y.size()))).norm() > 1e-8) {
        throw std::invalid_argument("BGEN AD scan requires full-rank covariates with an intercept");
    }
    std::unordered_set<int> unique;
    for (int sample : kernels.sample_indices) {
        if (sample < 0 || !unique.insert(sample).second) throw std::invalid_argument("BGEN AD scan sample indices");
    }
    const auto pair = excluded_chromosome.empty() ? kernels.global() : kernels.loco(excluded_chromosome);
    EigenVceAD::Options fixed_options;
    fixed_options.verbose = false;
    EigenVceAD fixed(fixed_options);
    fixed.prepare(y, covariates, pair.additive, pair.dominance);
    fixed.fixedComponents(sigma_a2, sigma_d2, sigma_e2);
    rotation_ = fixed.getProjectedW();
    variance_eigenvalues_ = fixed.getEigenvalues();
    if (!rotation_.allFinite() || !variance_eigenvalues_.allFinite() ||
        rotation_.rows() != y.size() || rotation_.cols() != variance_eigenvalues_.size() ||
        (variance_eigenvalues_.array() <= 0.0).any()) {
        throw std::runtime_error("BGEN AD scan invalid covariance rotation");
    }
    inverse_variance_ = variance_eigenvalues_.cwiseInverse();
    weighted_y_ = fixed.getUty().cwiseProduct(inverse_variance_);
    sample_indices_ = kernels.sample_indices;
    sample_ids_ = kernels.sample_ids;
    prepared_ = true;
}

void BgenAdScanEngine::updateFixedPhenotype(const Eigen::VectorXd& y) {
    if (!prepared_ || null_result_) throw std::logic_error("Only an explicit fixed-covariance validation null may be rescored");
    if (y.size() != rotation_.rows() || !y.allFinite()) throw std::invalid_argument("BGEN AD fixed phenotype dimensions");
    weighted_y_ = (rotation_.transpose() * y).cwiseProduct(inverse_variance_);
}

void BgenAdScanEngine::evaluateBlock(const Eigen::MatrixXd& additive, const Eigen::MatrixXd& dominance,
                                    const std::vector<int>& indices, std::vector<AdResult>& results) {
    const Eigen::Index count = static_cast<Eigen::Index>(indices.size());
    if (count == 0) return;
    const Eigen::MatrixXd ga = rotation_.transpose() * additive.leftCols(count);
    const Eigen::MatrixXd gd = rotation_.transpose() * dominance.leftCols(count);
    const Eigen::VectorXd score_a = ga.transpose() * weighted_y_;
    const Eigen::VectorXd score_d = gd.transpose() * weighted_y_;
    for (Eigen::Index col = 0; col < count; ++col) {
        const double aa = (ga.col(col).array().square() * inverse_variance_.array()).sum();
        const double dd = (gd.col(col).array().square() * inverse_variance_.array()).sum();
        const double ad = (ga.col(col).array() * gd.col(col).array() * inverse_variance_.array()).sum();
        const double determinant = aa * dd - ad * ad;
        auto& result = results[indices[col]];
        if (!(aa > 1e-12) || !(dd > 1e-12) || !(determinant > 1e-12 * aa * dd) ||
            !std::isfinite(determinant) || !std::isfinite(score_a[col]) || !std::isfinite(score_d[col])) {
            ++statistics_.rank_deficient;
            continue;
        }
        Eigen::Matrix2d information;
        information << aa, ad, ad, dd;
        Eigen::LLT<Eigen::Matrix2d> factor(information);
        if (factor.info() != Eigen::Success) { ++statistics_.rank_deficient; continue; }
        const Eigen::Vector2d score(score_a[col], score_d[col]);
        const Eigen::Vector2d beta = factor.solve(score);
        const Eigen::Matrix2d covariance = factor.solve(Eigen::Matrix2d::Identity());
        const Eigen::Vector2d whitened = factor.matrixL().solve(score);
        if (!beta.allFinite() || !covariance.allFinite() || covariance(0, 0) <= 0 || covariance(1, 1) <= 0) {
            ++statistics_.rank_deficient;
            continue;
        }
        result.beta_add = beta[0];
        result.beta_dom = beta[1];
        result.se_add = std::sqrt(covariance(0, 0));
        result.se_dom = std::sqrt(covariance(1, 1));
        result.p_add = std::erfc(std::abs(beta[0] / result.se_add) / std::sqrt(2.0));
        result.p_dom = std::erfc(std::abs(beta[1] / result.se_dom) / std::sqrt(2.0));
        result.joint_chisq = whitened.squaredNorm();
        result.p_joint = std::exp(-0.5 * result.joint_chisq);
        result.log10_p_dom = result.p_dom > 0 ? -std::log10(result.p_dom) : 300.0;
        result.log10_p_joint = result.joint_chisq / (2.0 * std::log(10.0));
        result.rank_deficient = false;
        ++statistics_.tested;
    }
}

void BgenAdScanEngine::scan(BgenReader& reader, std::vector<AdResult>& results) {
    if (!prepared_) throw std::logic_error("BGEN AD scan requires a prepared null");
    if (reader.getCompressionType() == 0) throw std::invalid_argument("BGEN AD scan requires compressed probability input");
    for (size_t i = 0; i < sample_indices_.size(); ++i) {
        const int source = sample_indices_[i];
        if (static_cast<uint32_t>(source) >= reader.getNumSamples() ||
            reader.getSampleIds().at(source) != sample_ids_[i]) {
            throw std::invalid_argument("BGEN AD scan genotype/null sample order mismatch");
        }
    }
    statistics_ = Statistics{};
    results.assign(reader.getNumVariants(), AdResult{});
    Eigen::MatrixXd ga(sample_indices_.size(), options_.block_size);
    Eigen::MatrixXd gd(sample_indices_.size(), options_.block_size);
    std::vector<int> block;
    const auto flush = [&]() { evaluateBlock(ga, gd, block, results); block.clear(); };
    reader.reset();
    for (uint32_t j = 0; j < reader.getNumVariants(); ++j) {
        std::string chrom, id, returned_a1, returned_a2;
        uint32_t position = 0;
        BgenGenotypeProbs probabilities;
        if (!reader.readVariantProbs(chrom, id, position, returned_a1, returned_a2, probabilities)) {
            throw std::runtime_error("BGEN AD scan truncated variant stream");
        }
        auto& result = results[j];
        result.chrom = chrom; result.snp_id = id; result.pos = position;
        result.a1 = returned_a2; result.a2 = returned_a1;
        result.n_obs = static_cast<int>(sample_indices_.size());
        if ((!options_.chromosome.empty() && chrom != options_.chromosome) ||
            (options_.from_bp >= 0 && position < options_.from_bp) ||
            (options_.to_bp >= 0 && position > options_.to_bp) ||
            (!options_.extract.empty() && !options_.extract.count(id)) || options_.exclude.count(id)) {
            ++statistics_.target_filtered;
            continue;
        }
        const auto marker = encode(probabilities, sample_indices_);
        result.freq = marker.frequency;
        const double maf = std::min(marker.frequency, 1.0 - marker.frequency);
        if (marker.observed < 3 || marker.missing_rate > options_.max_missing || maf < options_.min_maf ||
            maf < 1e-8 || marker.info_a < options_.min_info_a || marker.info_d < options_.min_info_d) {
            ++statistics_.qc_filtered;
            continue;
        }
        const Eigen::Index column = static_cast<Eigen::Index>(block.size());
        ga.col(column) = marker.additive;
        gd.col(column) = marker.dominance;
        block.push_back(j);
        if (block.size() == options_.block_size) flush();
    }
    flush();
}

} // namespace cosmic
