#include "ld_calculator.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace cosmic {
namespace genetics {
namespace {

void normalize(Eigen::MatrixXd& dosages) {
    if (dosages.rows() < 2 || dosages.cols() > std::numeric_limits<int>::max()) {
        throw std::invalid_argument("LD: expected >= 2 samples and an int-sized variant count");
    }
    for (Eigen::Index j = 0; j < dosages.cols(); ++j) {
        double sum = 0.0, minimum = 2.0, maximum = 0.0;
        Eigen::Index count = 0;
        for (Eigen::Index i = 0; i < dosages.rows(); ++i) {
            const double value = dosages(i, j);
            if (std::isnan(value)) continue;
            if (!std::isfinite(value) || value > 2.0) {
                throw std::invalid_argument("LD: invalid dosage");
            }
            if (value < 0.0) continue;
            sum += value;
            ++count;
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
        }
        if (count == 0 || minimum == maximum) {
            throw std::invalid_argument("LD: all-missing or monomorphic variant");
        }
        const double mean = sum / static_cast<double>(count);
        for (Eigen::Index i = 0; i < dosages.rows(); ++i) {
            double& value = dosages(i, j);
            value = std::isnan(value) || value < 0.0 ? 0.0 : value - mean;
        }
        dosages.col(j) /= dosages.col(j).stableNorm();
    }
}

void validateWindow(Eigen::Index m, const std::vector<std::string>& chrom,
                    const std::vector<long long>& pos,
                    int window_bp, int step, double r2) {
    if (chrom.size() != static_cast<std::size_t>(m) || pos.size() != chrom.size() ||
        window_bp < 0 || step < 1 || !std::isfinite(r2) || r2 < 0.0 || r2 > 1.0) {
        throw std::invalid_argument("LD: invalid window parameters or metadata dimensions");
    }
    for (std::size_t j = 0; j < chrom.size(); ++j) {
        if (chrom[j].empty() || pos[j] < 0 ||
            (j > 0 && (chrom[j] < chrom[j - 1] ||
                       (chrom[j] == chrom[j - 1] && pos[j] < pos[j - 1])))) {
            throw std::invalid_argument("LD: metadata must be sorted by (chromosome, position)");
        }
    }
}

double correlation(const Eigen::MatrixXd& normalized, int i, int j) {
    double r = normalized.col(i).dot(normalized.col(j));
    // Preserve perfect-LD boundary behavior despite normalization roundoff.
    if (std::abs(r) >= 1.0 - 16.0 * std::numeric_limits<double>::epsilon()) {
        r = std::copysign(1.0, r);
    }
    return r;
}

} // namespace

Eigen::MatrixXd LDCalculator::computePairwiseLD(Eigen::MatrixXd dosages) {
    normalize(dosages);
    const int m = static_cast<int>(dosages.cols());
    Eigen::MatrixXd result = Eigen::MatrixXd::Identity(m, m);
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            result(i, j) = result(j, i) = correlation(dosages, i, j);
        }
    }
    return result;
}

std::vector<int> LDPruner::prune(Eigen::MatrixXd dosages,
                                const std::vector<std::string>& chrom,
                                const std::vector<long long>& pos,
                                int window_bp, int step, double r2) {
    validateWindow(dosages.cols(), chrom, pos, window_bp, step, r2);
    normalize(dosages);
    const int m = static_cast<int>(dosages.cols());
    std::vector<bool> removed(m, false);
    for (int begin = 0; begin < m;) {
        int end = begin + 1;
        while (end < m && chrom[end] == chrom[begin]) ++end;
        for (int start = begin; start < end; start += std::min(step, end - start)) {
            int stop = start + 1;
            while (stop < end && pos[stop] - pos[start] <= window_bp) ++stop;
            for (int i = start; i < stop; ++i) {
                if (removed[i]) continue;
                for (int j = i + 1; j < stop; ++j) {
                    if (removed[j]) continue;
                    const double r = correlation(dosages, i, j);
                    if (r * r >= r2) removed[j] = true;
                }
            }
        }
        begin = end;
    }
    std::vector<int> kept;
    for (int j = 0; j < m; ++j) {
        if (!removed[j]) kept.push_back(j);
    }
    return kept;
}

std::vector<int> LDClumper::clump(Eigen::MatrixXd dosages,
                                const std::vector<std::string>& chrom,
                                const std::vector<long long>& pos,
                                int window_bp, int step, double r2,
                                const Eigen::VectorXd& p, double threshold) {
    validateWindow(dosages.cols(), chrom, pos, window_bp, step, r2);
    if (p.size() != dosages.cols() || !p.allFinite() ||
        (p.array() < 0.0).any() || (p.array() > 1.0).any() ||
        !std::isfinite(threshold) || threshold < 0.0 || threshold > 1.0) {
        throw std::invalid_argument("LD: invalid p-values or lead threshold");
    }
    normalize(dosages);
    const int m = static_cast<int>(dosages.cols());
    std::vector<int> order(m), lead(m, -1);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return p[a] < p[b] || (p[a] == p[b] && a < b);
    });
    for (int index : order) {
        if (p[index] > threshold) break;
        if (lead[index] != -1) continue;
        lead[index] = index;
        int begin = index, end = index + 1;
        while (begin > 0 && chrom[begin - 1] == chrom[index] &&
               pos[index] - pos[begin - 1] <= window_bp) --begin;
        while (end < m && chrom[end] == chrom[index] &&
               pos[end] - pos[index] <= window_bp) ++end;
        for (int batch = begin; batch < end;) {
            const int stop = batch + std::min(step, end - batch);
            for (int j = batch; j < stop; ++j) {
                if (lead[j] != -1) continue;
                const double r = correlation(dosages, index, j);
                if (r * r >= r2) lead[j] = index;
            }
            batch = stop;
        }
    }
    return lead;
}

} // namespace genetics
} // namespace cosmic
