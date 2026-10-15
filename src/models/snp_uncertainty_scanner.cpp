#include "snp_uncertainty_scanner.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace cosmic {
namespace {

using Clock = std::chrono::steady_clock;

double secondsSince(const Clock::time_point& start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

void validate(const SnpUncertaintyScanInput& input) {
    const int m = input.reader.getNumSnps();
    const int k = static_cast<int>(input.core_genotyped_indices.size());
    if (m <= 0 || k <= 0 || input.effects.size() != m || input.denominator <= 1e-12) {
        throw std::runtime_error("SnpUncertaintyScanner invalid dimensions or denominator");
    }
    if (input.gw_cc.rows() != k || input.gw_cc.cols() != k ||
        input.pec_cc.rows() != k || input.pec_cc.cols() != k ||
        (input.predictor_covariance != nullptr &&
         (input.predictor_covariance->rows() != k || input.predictor_covariance->cols() != k))) {
        throw std::runtime_error("SnpUncertaintyScanner core matrix dimension mismatch");
    }
    if (input.block_size <= 0 || input.n_threads <= 0) {
        throw std::runtime_error("SnpUncertaintyScanner block size and thread count must be positive");
    }
    for (int marker : input.snp_indices) {
        if (marker < 0 || marker >= m) {
            throw std::runtime_error("SnpUncertaintyScanner marker index is out of range");
        }
    }
    const int n_samples = input.reader.getNumSamples();
    for (int sample : input.genotyped_keep) {
        if (sample < 0 || sample >= n_samples) {
            throw std::runtime_error("SnpUncertaintyScanner sample index is out of range");
        }
    }
    for (int core_index : input.core_genotyped_indices) {
        if (core_index < 0 || core_index >= static_cast<int>(input.genotyped_keep.size())) {
            throw std::runtime_error("SnpUncertaintyScanner core index is out of range");
        }
    }
}

SnpUncertaintyScanResult initializeResult(int m, const std::string& mode, int block_size) {
    SnpUncertaintyScanResult out;
    out.variance.assign(m, 0.0);
    out.se.assign(m, 0.0);
    out.z.assign(m, 0.0);
    out.p.assign(m, 1.0);
    out.log10p.assign(m, 0.0);
    out.status.assign(m, "NA");
    out.audit.mode = mode;
    out.audit.block_size = block_size;
    return out;
}

void finishMarker(int j, double variance, double variance_scale,
                  SnpUncertaintyScanResult& out, const Eigen::VectorXd& effects) {
    const double tolerance = 1e-10 * std::max(1.0, std::fabs(variance_scale));
    if (variance < 0.0) {
        if (variance >= -tolerance) {
            variance = 0.0;
            ++out.clamped;
            out.status[j] = "NEGATIVE_VARIANCE_CLAMPED";
        } else {
            ++out.failures;
            out.status[j] = "NUMERICAL_FAILURE";
            out.variance[j] = variance;
            out.se[j] = out.z[j] = out.p[j] = out.log10p[j] =
                std::numeric_limits<double>::quiet_NaN();
            return;
        }
    }
    out.variance[j] = variance;
    const double se = std::sqrt(variance);
    const double effect = effects(j);
    const double z = effect / se;
    double p = std::erfc(std::fabs(z) / std::sqrt(2.0));
    if (!(p > 0.0)) p = 0.0;
    out.se[j] = se;
    out.z[j] = z;
    out.p[j] = p;
    out.log10p[j] = p > 0.0 ? -std::log10(p) : 320.0;
    if (se <= 1e-15 * std::max(1.0, std::fabs(effect))) {
        ++out.zero_se;
        out.status[j] = "ZERO_SE";
    } else if (out.status[j] == "NA") {
        out.status[j] = "OK";
    }
}

long long coreBytesPerSnp(const std::vector<int>& core_keep) {
    std::unordered_set<int> bytes;
    bytes.reserve(core_keep.size());
    for (int sample : core_keep) bytes.insert(sample >> 2);
    return static_cast<long long>(bytes.size());
}

void finalizeAudit(SnpUncertaintyScanAudit& audit) {
    const double measured = audit.time_setup + audit.time_onehot + audit.time_genotype +
        audit.time_core_gather + audit.time_gcc_solve + audit.time_gw_multiply +
        audit.time_pec_multiply + audit.time_quadratic + audit.time_statistics;
    audit.time_other = std::max(0.0, audit.time_total - measured);
}

void multiplyCoreBlock(const Eigen::MatrixXd& left, const Eigen::MatrixXd& right,
                       int columns, Eigen::MatrixXd& output) {
    const int rows = static_cast<int>(left.rows());
    for (int column = 0; column < columns; ++column) {
        for (int row = 0; row < rows; ++row) {
            double sum = 0.0;
            for (int inner = 0; inner < rows; ++inner) {
                sum += left(row, inner) * right(inner, column);
            }
            output(row, column) = sum;
        }
    }
}

} // namespace

SnpUncertaintyScanResult ScalarReferenceSnpUncertaintyScanner::scan(
    const SnpUncertaintyScanInput& input) const {
    validate(input);
    const auto total_start = Clock::now();
    const int m = input.reader.getNumSnps();
    const int n = static_cast<int>(input.genotyped_keep.size());
    const int k = static_cast<int>(input.core_genotyped_indices.size());
    const auto& inv_stds = input.reader.getInvStds();
    const double denominator2 = input.denominator * input.denominator;

    auto start = Clock::now();
    SnpUncertaintyScanResult out = initializeResult(m, "scalar", 1);
    Eigen::VectorXd onehot = Eigen::VectorXd::Zero(m);
    Eigen::VectorXd genotype(n), core(k), weight(k), gw_product(k), pec_product(k);
    out.audit.heap_allocations = 6;
    out.audit.peak_transient_bytes = static_cast<long long>(sizeof(double)) * (m + n + 4LL * k);
    out.audit.time_setup = secondsSince(start);

    for (int j : input.snp_indices) {
        const double inv_std = j < static_cast<int>(inv_stds.size()) ? inv_stds[j] : 0.0;
        if (!(inv_std > 0.0)) {
            out.status[j] = "NA_MONOMORPHIC";
            continue;
        }
        start = Clock::now();
        onehot.setZero();
        onehot(j) = 1.0;
        out.audit.time_onehot += secondsSince(start);

        start = Clock::now();
        input.reader.multiply_Z_v(onehot, genotype, input.genotyped_keep, input.n_threads);
        out.audit.time_genotype += secondsSince(start);
        ++out.audit.multiply_z_v_calls;
        ++out.audit.genotype_passes;
        out.audit.bed_bytes_read_estimate +=
            static_cast<long long>(input.reader.getBytesPerSnp()) * m;

        start = Clock::now();
        for (int i = 0; i < k; ++i) core(i) = genotype(input.core_genotyped_indices[i]) / inv_std;
        out.audit.time_core_gather += secondsSince(start);

        start = Clock::now();
        weight = input.core_inverse.applyCoreInverse(core);
        out.audit.time_gcc_solve += secondsSince(start);
        ++out.audit.gcc_solve_calls;
        ++out.audit.gcc_rhs_columns;

        start = Clock::now();
        gw_product.noalias() = input.gw_cc * weight;
        out.audit.time_gw_multiply += secondsSince(start);
        ++out.audit.gw_gemm_calls;

        double variance;
        start = Clock::now();
        if (input.predictor_covariance != nullptr) {
            pec_product.noalias() = *input.predictor_covariance * weight;
        } else {
            pec_product.noalias() = input.pec_cc * weight;
        }
        out.audit.time_pec_multiply += secondsSince(start);
        ++out.audit.pec_gemm_calls;

        start = Clock::now();
        const double gquad = weight.dot(gw_product);
        if (input.predictor_covariance != nullptr) {
            variance = weight.dot(pec_product) / denominator2;
        } else {
            variance = (input.sigma2_a * gquad - input.sigma2_e * weight.dot(pec_product)) /
                       denominator2;
        }
        const double scale = input.sigma2_a * gquad / denominator2;
        out.audit.time_quadratic += secondsSince(start);

        start = Clock::now();
        finishMarker(j, variance, scale, out, input.effects);
        out.audit.time_statistics += secondsSince(start);
        ++out.audit.variants_tested;
        ++out.audit.variants_decoded;
    }
    out.audit.blocks = out.audit.variants_decoded;
    out.audit.time_total = secondsSince(total_start);
    finalizeAudit(out.audit);
    return out;
}

SnpUncertaintyScanResult BatchedOperatorSnpUncertaintyScanner::scan(
    const SnpUncertaintyScanInput& input) const {
    validate(input);
    const auto total_start = Clock::now();
    const int m = input.reader.getNumSnps();
    const int k = static_cast<int>(input.core_genotyped_indices.size());
    const auto& inv_stds = input.reader.getInvStds();
    const auto& means = input.reader.getMeans();
    constexpr long long kTransientMemoryCap = 256LL * 1024LL * 1024LL;
    const long long bytes_per_column = std::max(1, k) * 6LL * sizeof(double);
    const int memory_limited_capacity = static_cast<int>(std::max(
        1LL, std::min(static_cast<long long>(std::numeric_limits<int>::max()),
                      kTransientMemoryCap / bytes_per_column)));
    const int capacity = std::min({input.block_size, memory_limited_capacity,
                                  std::max(1, static_cast<int>(input.snp_indices.size()))});
    const double denominator2 = input.denominator * input.denominator;

    auto start = Clock::now();
    SnpUncertaintyScanResult out = initializeResult(m, "batched", capacity);
    std::vector<int> core_keep(k), block_indices(capacity);
    for (int i = 0; i < k; ++i) core_keep[i] = input.genotyped_keep[input.core_genotyped_indices[i]];
    Eigen::MatrixXd centered(k, capacity), weights(k, capacity);
    Eigen::MatrixXd gw_products(k, capacity), pec_products(k, capacity);
    out.audit.heap_allocations = 6;
    out.audit.peak_transient_bytes = bytes_per_column * capacity;
    out.audit.time_setup = secondsSince(start);

    const uint8_t* bed = input.reader.getData() + 3;
    const size_t bytes_per_snp = input.reader.getBytesPerSnp();
    const long long bytes_per_marker = coreBytesPerSnp(core_keep);
    for (int begin = 0; begin < static_cast<int>(input.snp_indices.size()); begin += capacity) {
        const int width = std::min(capacity,
            static_cast<int>(input.snp_indices.size()) - begin);
        int active = 0;
        for (int offset = 0; offset < width; ++offset) {
            const int j = input.snp_indices[begin + offset];
            if (!(j < static_cast<int>(inv_stds.size()) && inv_stds[j] > 0.0)) {
                out.status[j] = "NA_MONOMORPHIC";
                continue;
            }
            block_indices[active++] = j;
        }
        if (active == 0) continue;

        start = Clock::now();
        for (int column = 0; column < active; ++column) {
            const int j = block_indices[column];
            const uint8_t* snp = bed + static_cast<size_t>(j) * bytes_per_snp;
            const double mean = means[j];
            for (int i = 0; i < k; ++i) {
                const int sample = core_keep[i];
                const int code = (snp[sample >> 2] >> ((sample & 3) * 2)) & 3;
                if (code == 0) centered(i, column) = 2.0 - mean;
                else if (code == 2) centered(i, column) = 1.0 - mean;
                else if (code == 3) centered(i, column) = -mean;
                else {
                    centered(i, column) = 0.0;
                    ++out.audit.missing_to_zero;
                }
            }
        }
        out.audit.time_genotype += secondsSince(start);
        out.audit.variants_decoded += active;
        out.audit.bed_bytes_read_estimate += bytes_per_marker * active;

        start = Clock::now();
        const Eigen::MatrixXd rhs = centered.leftCols(active);
        weights.leftCols(active) = input.core_inverse.applyCoreInverse(rhs);
        out.audit.time_gcc_solve += secondsSince(start);
        ++out.audit.gcc_solve_calls;
        out.audit.gcc_rhs_columns += active;

        start = Clock::now();
        multiplyCoreBlock(input.gw_cc, weights, active, gw_products);
        out.audit.time_gw_multiply += secondsSince(start);
        ++out.audit.gw_gemm_calls;

        start = Clock::now();
        if (input.predictor_covariance != nullptr) {
            multiplyCoreBlock(*input.predictor_covariance, weights, active, pec_products);
        } else {
            multiplyCoreBlock(input.pec_cc, weights, active, pec_products);
        }
        out.audit.time_pec_multiply += secondsSince(start);
        ++out.audit.pec_gemm_calls;

        start = Clock::now();
        const double statistics_before = out.audit.time_statistics;
        for (int column = 0; column < active; ++column) {
            const int j = block_indices[column];
            const double gquad = weights.col(column).dot(gw_products.col(column));
            double variance;
            if (input.predictor_covariance != nullptr) {
                variance = weights.col(column).dot(pec_products.col(column)) / denominator2;
            } else {
                const double pquad = weights.col(column).dot(pec_products.col(column));
                variance = (input.sigma2_a * gquad - input.sigma2_e * pquad) / denominator2;
            }
            const double scale = input.sigma2_a * gquad / denominator2;
            const auto statistics_start = Clock::now();
            finishMarker(j, variance, scale, out, input.effects);
            out.audit.time_statistics += secondsSince(statistics_start);
            ++out.audit.variants_tested;
        }
        out.audit.time_quadratic +=
            secondsSince(start) - (out.audit.time_statistics - statistics_before);
        ++out.audit.blocks;
    }
    out.audit.genotype_passes = out.audit.variants_decoded > 0 ? 1 : 0;
    out.audit.time_total = secondsSince(total_start);
    finalizeAudit(out.audit);
    return out;
}

SnpUncertaintyScanResult BatchedExactSnpUncertaintyScanner::scan(
    const ExactSnpUncertaintyScanInput& input) const {
    const auto total_start = Clock::now();
    const int m = input.reader.getNumSnps();
    const int d = static_cast<int>(input.sample_keep.size());
    if (d <= 0 || m <= 0 || input.effects.size() != m ||
        !(input.denominator > 1e-12) || !std::isfinite(input.denominator) ||
        !(input.sigma2_a > 0.0) || input.block_size <= 0 ||
        (input.dense_effect_inverse == nullptr) == (input.core_inverse == nullptr)) {
        throw std::runtime_error("Batched exact SNP scan: invalid dimensions/inverse target");
    }
    for (const auto* matrix : {&input.evaluation_relationship, &input.predictor_covariance}) {
        if (matrix->rows() != d || matrix->cols() != d || !matrix->allFinite()) {
            throw std::runtime_error("Batched exact SNP scan: covariance/prior dimension mismatch");
        }
    }
    if (input.dense_effect_inverse && (input.dense_effect_inverse->rows() != d ||
        input.dense_effect_inverse->cols() != d || !input.dense_effect_inverse->allFinite())) {
        throw std::runtime_error("Batched exact SNP scan: dense inverse dimension mismatch");
    }
    std::unordered_set<int> unique_samples;
    for (int sample : input.sample_keep) {
        if (sample < 0 || sample >= input.reader.getNumSamples() || !unique_samples.insert(sample).second) {
            throw std::runtime_error("Batched exact SNP scan: invalid/repeated BED sample index");
        }
    }
    for (int j : input.snp_indices) {
        if (j < 0 || j >= m) throw std::runtime_error("Batched exact SNP scan: SNP index out of range");
    }
    const auto& means = input.reader.getMeans();
    const auto& inv_stds = input.reader.getInvStds();
    if (means.size() < static_cast<size_t>(m) || inv_stds.size() < static_cast<size_t>(m)) {
        throw std::runtime_error("Batched exact SNP scan requires prepared reader statistics");
    }
    constexpr long long cap_bytes = 256LL * 1024 * 1024;
    const long long bytes_per_column = 6LL * sizeof(double) * d;
    if (bytes_per_column > cap_bytes) {
        throw std::runtime_error("Batched exact SNP scan: one RHS exceeds transient memory budget");
    }
    const int capacity = std::min({input.block_size,
        static_cast<int>(cap_bytes / bytes_per_column),
        std::max(1, static_cast<int>(input.snp_indices.size()))});
    auto start = Clock::now();
    auto out = initializeResult(m, "batched_exact_selected", capacity);
    Eigen::MatrixXd centered(d, capacity), weights(d, capacity);
    Eigen::MatrixXd prior_product(d, capacity), covariance_product(d, capacity);
    std::vector<int> indices(capacity);
    out.audit.heap_allocations = 5; // named work buffers only, not BLAS/backend scratch
    out.audit.peak_transient_bytes = bytes_per_column * capacity;
    out.audit.time_setup = secondsSince(start);
    const uint8_t* bed = input.reader.getData() + 3;
    const size_t stride = input.reader.getBytesPerSnp();
    const long long bytes_per_marker = coreBytesPerSnp(input.sample_keep);
    const double S2 = input.denominator * input.denominator;
    for (int begin = 0; begin < static_cast<int>(input.snp_indices.size()); begin += capacity) {
        const int width = std::min(capacity, static_cast<int>(input.snp_indices.size()) - begin);
        int active = 0;
        for (int col = 0; col < width; ++col) {
            const int j = input.snp_indices[begin + col];
            if (!(inv_stds[j] > 0.0)) {
                out.status[j] = "NA_MONOMORPHIC";
            } else indices[active++] = j;
        }
        if (!active) continue;
        start = Clock::now();
        for (int col = 0; col < active; ++col) {
            const int j = indices[col];
            const uint8_t* snp = bed + static_cast<size_t>(j) * stride;
            for (int row = 0; row < d; ++row) {
                const int sample = input.sample_keep[row];
                const int code = (snp[sample >> 2] >> ((sample & 3) * 2)) & 3;
                if (code == 0) centered(row, col) = 2.0 - means[j];
                else if (code == 2) centered(row, col) = 1.0 - means[j];
                else if (code == 3) centered(row, col) = -means[j];
                else { centered(row, col) = 0.0; ++out.audit.missing_to_zero; }
            }
        }
        out.audit.time_genotype += secondsSince(start);
        out.audit.variants_decoded += active;
        out.audit.bed_bytes_read_estimate += bytes_per_marker * active;
        start = Clock::now();
        if (input.dense_effect_inverse) {
            weights.leftCols(active).noalias() = *input.dense_effect_inverse * centered.leftCols(active);
        } else {
            // Gcc only. Full APY apply(Matrix) may rescan the panel per RHS.
            const Eigen::MatrixXd rhs = centered.leftCols(active);
            weights.leftCols(active) = input.core_inverse->applyCoreInverse(rhs);
            out.audit.heap_allocations += 2; // explicit RHS and returned matrix
        }
        out.audit.time_gcc_solve += secondsSince(start);
        ++out.audit.gcc_solve_calls;
        out.audit.gcc_rhs_columns += active;
        if (!weights.leftCols(active).allFinite()) {
            throw std::runtime_error("Batched exact SNP scan: nonfinite inverse action");
        }
        start = Clock::now();
        prior_product.leftCols(active).noalias() = input.evaluation_relationship * weights.leftCols(active);
        out.audit.time_gw_multiply += secondsSince(start);
        ++out.audit.gw_gemm_calls;
        start = Clock::now();
        covariance_product.leftCols(active).noalias() = input.predictor_covariance * weights.leftCols(active);
        out.audit.time_pec_multiply += secondsSince(start);
        ++out.audit.pec_gemm_calls;
        for (int col = 0; col < active; ++col) {
            start = Clock::now();
            const int j = indices[col];
            double variance = weights.col(col).dot(covariance_product.col(col)) / S2;
            const double scale = input.sigma2_a * weights.col(col).dot(prior_product.col(col)) / S2;
            out.audit.time_quadratic += secondsSince(start);
            start = Clock::now();
            const double tolerance = 1e-10 * std::max(1.0, std::fabs(scale));
            if (variance < 0.0) {
                if (variance >= -tolerance) {
                    variance = 0.0; ++out.clamped; out.status[j] = "NEGATIVE_VARIANCE_CLAMPED";
                } else {
                    ++out.failures; out.status[j] = "NUMERICAL_FAILURE";
                    out.variance[j] = variance;
                    out.se[j] = out.z[j] = out.p[j] = out.log10p[j] = std::numeric_limits<double>::quiet_NaN();
                    out.audit.time_statistics += secondsSince(start);
                    ++out.audit.variants_tested;
                    continue;
                }
            }
            out.variance[j] = variance;
            const double se = std::sqrt(variance);
            out.se[j] = se;
            if (se <= 1e-15 * std::max(1.0, std::fabs(input.effects(j)))) {
                ++out.zero_se; out.status[j] = "ZERO_SE";
                out.z[j] = out.p[j] = out.log10p[j] = std::numeric_limits<double>::quiet_NaN();
            } else {
                const double z = input.effects(j) / se;
                double p = std::erfc(std::fabs(z) / std::sqrt(2.0));
                if (!(p > 0.0)) p = 0.0;
                out.z[j] = z; out.p[j] = p; out.log10p[j] = p > 0.0 ? -std::log10(p) : 320.0;
                if (out.status[j] == "NA") out.status[j] = "OK";
            }
            out.audit.time_statistics += secondsSince(start);
            ++out.audit.variants_tested;
        }
        ++out.audit.blocks;
    }
    out.audit.genotype_passes = out.audit.variants_decoded > 0 ? 1 : 0;
    out.audit.time_total = secondsSince(total_start);
    finalizeAudit(out.audit);
    return out;
}

} // namespace cosmic
