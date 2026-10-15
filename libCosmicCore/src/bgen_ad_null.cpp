#include "bgen_ad_null.h"
#include "bgen_dominance_encoding.h"
#include "bgen_ad_reml.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace cosmic {
namespace {

BgenAdRawPartition emptyPartition(Eigen::Index n) {
    BgenAdRawPartition result;
    result.additive = Eigen::MatrixXd::Zero(n, n);
    result.dominance = Eigen::MatrixXd::Zero(n, n);
    return result;
}

BgenAdKernelPair normalize(const BgenAdRawPartition& raw) {
    if (raw.markers <= 0 || !(raw.denominator_a > 1e-12) ||
        !(raw.denominator_d > 1e-12) || !std::isfinite(raw.denominator_a) ||
        !std::isfinite(raw.denominator_d)) {
        throw std::runtime_error("BGEN AD null: no valid shared markers/normalization");
    }
    BgenAdKernelPair result;
    result.additive = raw.additive / raw.denominator_a;
    result.dominance = raw.dominance / raw.denominator_d;
    result.denominator_a = raw.denominator_a;
    result.denominator_d = raw.denominator_d;
    result.markers = raw.markers;
    return result;
}

} // namespace

BgenAdKernelPair BgenAdKernelResult::global() const {
    return normalize(global_raw);
}

BgenAdKernelPair BgenAdKernelResult::loco(const std::string& chromosome) const {
    const auto found = chromosome_raw.find(chromosome);
    if (found == chromosome_raw.end()) {
        throw std::invalid_argument("BGEN AD null: unknown LOCO chromosome " + chromosome);
    }
    BgenAdRawPartition outside;
    outside.additive = global_raw.additive - found->second.additive;
    outside.dominance = global_raw.dominance - found->second.dominance;
    outside.denominator_a = global_raw.denominator_a - found->second.denominator_a;
    outside.denominator_d = global_raw.denominator_d - found->second.denominator_d;
    outside.markers = global_raw.markers - found->second.markers;
    return normalize(outside);
}

BgenAdKernelResult BgenAdNullBackend::build(BgenReader& reader,
                                          const BgenAdKernelOptions& options,
                                          const std::vector<int>& sample_indices) {
    if (reader.getCompressionType() == 0) {
        throw std::invalid_argument("BGEN AD null requires compressed Layout-2 probabilities");
    }
    if (!std::isfinite(options.min_maf) || options.min_maf < 0.0 || options.min_maf > 0.5 ||
        !std::isfinite(options.max_missing) || options.max_missing < 0.0 || options.max_missing > 1.0 ||
        !std::isfinite(options.min_info_a) || options.min_info_a < 0.0 ||
        !std::isfinite(options.min_info_d) || options.min_info_d < 0.0 || options.block_size == 0) {
        throw std::invalid_argument("BGEN AD null: invalid kernel options");
    }

    BgenAdKernelResult result;
    result.sample_indices = sample_indices;
    if (result.sample_indices.empty()) {
        result.sample_indices.resize(reader.getNumSamples());
        std::iota(result.sample_indices.begin(), result.sample_indices.end(), 0);
    }
    if (result.sample_indices.empty()) throw std::invalid_argument("BGEN AD null: no samples");
    std::unordered_set<int> unique_indices;
    for (int index : result.sample_indices) {
        if (index < 0 || static_cast<uint32_t>(index) >= reader.getNumSamples() ||
            !unique_indices.insert(index).second) {
            throw std::invalid_argument("BGEN AD null: duplicate/out-of-range sample index");
        }
        result.sample_ids.push_back(reader.getSampleIds().at(static_cast<size_t>(index)));
    }
    const Eigen::Index n = static_cast<Eigen::Index>(result.sample_indices.size());
    result.global_raw = emptyPartition(n);
    Eigen::MatrixXd za(n, options.block_size);
    Eigen::MatrixXd zd(n, options.block_size);
    std::vector<std::string> block_chromosomes;
    const auto flush = [&]() {
        const Eigen::Index count = static_cast<Eigen::Index>(block_chromosomes.size());
        if (count == 0) return;
        result.global_raw.additive.noalias() += za.leftCols(count) * za.leftCols(count).transpose();
        result.global_raw.dominance.noalias() += zd.leftCols(count) * zd.leftCols(count).transpose();
        std::map<std::string, std::vector<Eigen::Index>> columns;
        for (Eigen::Index j = 0; j < count; ++j) columns[block_chromosomes[j]].push_back(j);
        for (const auto& entry : columns) {
            Eigen::MatrixXd a(n, entry.second.size());
            Eigen::MatrixXd d(n, entry.second.size());
            for (size_t j = 0; j < entry.second.size(); ++j) {
                a.col(j) = za.col(entry.second[j]);
                d.col(j) = zd.col(entry.second[j]);
            }
            auto& partition = result.chromosome_raw.at(entry.first);
            partition.additive.noalias() += a * a.transpose();
            partition.dominance.noalias() += d * d.transpose();
        }
        block_chromosomes.clear();
    };

    reader.reset();
    const uint32_t variants = reader.getNumVariants();
    for (uint32_t j = 0; j < variants; ++j) {
        std::string chrom, id, returned_a1, returned_a2;
        uint32_t position = 0;
        BgenGenotypeProbs probabilities;
        if (!reader.readVariantProbs(chrom, id, position, returned_a1, returned_a2, probabilities)) {
            throw std::runtime_error("BGEN AD null: truncated variant stream");
        }
        if ((!options.include_chromosomes.empty() && !options.include_chromosomes.count(chrom)) ||
            options.exclude_chromosomes.count(chrom) ||
            (!options.extract.empty() && !options.extract.count(id)) || options.exclude.count(id)) continue;
        ++result.candidates;
        if (!result.chromosome_raw.count(chrom)) result.chromosome_raw.emplace(chrom, emptyPartition(n));
        const auto qc = bgen::dominanceQcMetrics(probabilities.p0, probabilities.p1,
                                                probabilities.p2, probabilities.missing,
                                                result.sample_indices);
        if (qc.observed == 0 || qc.missing_rate > options.max_missing) {
            ++result.rejected_missing;
            continue;
        }
        const double p = qc.frequency;
        const double maf = std::min(p, 1.0 - p);
        if (maf < options.min_maf || maf < 1e-8) {
            ++result.rejected_maf;
            continue;
        }
        const double denominator_a = 2.0 * p * (1.0 - p);
        const double denominator_d = bgen::vitezicaDenominator(p);
        double additive_sq_sum = 0.0;
        const Eigen::Index column = static_cast<Eigen::Index>(block_chromosomes.size());
        for (Eigen::Index i = 0; i < n; ++i) {
            const int source = result.sample_indices[i];
            if (probabilities.missing[source] != 0.0f) {
                za(i, column) = zd(i, column) = 0.0;
            } else {
                const double x = static_cast<double>(probabilities.p1[source]) +
                    2.0 * static_cast<double>(probabilities.p2[source]);
                za(i, column) = x - 2.0 * p;
                zd(i, column) = bgen::posteriorMeanVitezica(
                    probabilities.p0[source], probabilities.p1[source], probabilities.p2[source], p);
                additive_sq_sum += za(i, column) * za(i, column);
            }
        }
        const double info_a = additive_sq_sum / (static_cast<double>(qc.observed) * denominator_a);
        if (info_a < options.min_info_a || qc.info_d < options.min_info_d) {
            ++result.rejected_info;
            continue;
        }
        result.accepted_variant_ids.push_back(id);
        ++result.global_raw.markers;
        result.global_raw.denominator_a += denominator_a;
        result.global_raw.denominator_d += denominator_d;
        auto& partition = result.chromosome_raw.at(chrom);
        ++partition.markers;
        partition.denominator_a += denominator_a;
        partition.denominator_d += denominator_d;
        block_chromosomes.push_back(chrom);
        if (block_chromosomes.size() == options.block_size) flush();
    }
    flush();
    // Validate normalization before returning any usable pair.
    if (result.global_raw.markers == 0 || result.global_raw.denominator_a <= 1e-12 ||
        result.global_raw.denominator_d <= 1e-12) {
        throw std::runtime_error("BGEN AD null: no valid shared markers/normalization");
    }
    return result;
}

BgenAdNullResult BgenAdNullBackend::fit(const BgenAdKernelPair& kernels,
                                      const Eigen::VectorXd& y,
                                      const Eigen::MatrixXf& covariates,
                                      const EigenVceAD::Options& options) {
    if (!options.use_covariate_projection || covariates.cols() == 0 ||
        covariates.rows() != y.size() || covariates.cols() >= y.size() ||
        kernels.additive.rows() != y.size() || kernels.additive.cols() != y.size() ||
        kernels.dominance.rows() != y.size() || kernels.dominance.cols() != y.size() ||
        options.max_em_iter < 0 || !std::isfinite(options.tol) || options.tol <= 0.0 ||
        !y.allFinite() || !covariates.allFinite() ||
        !kernels.additive.allFinite() || !kernels.dominance.allFinite()) {
        throw std::invalid_argument("BGEN AD null: invalid phenotype/covariate/kernel contract");
    }
    if (!kernels.additive.isApprox(kernels.additive.transpose(), 1e-10) ||
        !kernels.dominance.isApprox(kernels.dominance.transpose(), 1e-10)) {
        throw std::invalid_argument("BGEN AD null: kernels must be symmetric");
    }
    const Eigen::MatrixXd design = covariates.cast<double>();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);
    if (qr.rank() != design.cols() ||
        (Eigen::VectorXd::Ones(y.size()) - design * qr.solve(Eigen::VectorXd::Ones(y.size()))).norm() > 1e-8) {
        throw std::invalid_argument("BGEN AD null: covariates must be full rank and include an intercept");
    }
    BgenAdNullResult result = detail::fitBgenAdReml(kernels, y, covariates, options);
    if (!std::isfinite(result.sigma_a2) || !std::isfinite(result.sigma_d2) ||
        !std::isfinite(result.sigma_e2) || !std::isfinite(result.log_likelihood) ||
        result.sigma_a2 < 0.0 || result.sigma_d2 < 0.0 || !(result.sigma_e2 > 0.0)) {
        throw std::runtime_error("BGEN AD null: invalid fitted variance components");
    }
    return result;
}

} // namespace cosmic
