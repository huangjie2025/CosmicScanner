#include "cross_genomic_operator.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace cosmic {

namespace {

// Inverse read back of the reader's per-SNP standard deviation: recovers the
// centered (unstandardized) SNP column from the standardized Z-column.  A
// monomorphic SNP (inv_std == 0) contributes zero in both directions, so any
// positive scale is irrelevant there (matches GccOperator / Lanczos E1).
double stdScale(double inv_std) {
    return inv_std > 0.0 ? 1.0 / inv_std : 1.0;
}

// Z * (diag(1/inv_std) * v) for a subset of rows (used to turn the reader's
// standardized product into the centered product).  Output size == keep.size().
Eigen::VectorXd multiplyCentered(PlinkReader& reader,
                                 const Eigen::VectorXd& weighted_m,
                                 const std::vector<int>& keep,
                                 int n_threads) {
    Eigen::VectorXd y;
    reader.multiply_Z_v(weighted_m, y, keep, n_threads);
    return y;
}

} // namespace

CrossGenomicOperator::CrossGenomicOperator(
    PlinkReader& reader, std::vector<int> core_keep,
    std::vector<int> noncore_keep, double S_2pq, double tune_a, double tune_b,
    double blend_alpha, GcnBlockView view, const Eigen::MatrixXd* a22_full,
    std::vector<int> core_idx, std::vector<int> noncore_idx, int n_threads,
    const std::vector<char>* snp_valid, const std::vector<double>* snp_freq)
    : reader_(&reader), core_keep_(std::move(core_keep)),
      noncore_keep_(std::move(noncore_keep)), S_2pq_(S_2pq), tune_a_(tune_a),
      tune_b_(tune_b), blend_alpha_(blend_alpha), view_(view),
      n_threads_(n_threads > 0 ? n_threads : 1) {
    if (core_keep_.empty()) {
        throw std::runtime_error("CrossGenomicOperator: empty core block");
    }
    if (a22_full != nullptr && blend_alpha_ > 0.0) {
        const int nc = static_cast<int>(core_keep_.size());
        const int nn = static_cast<int>(noncore_keep_.size());
        if (nn == 0) {
            throw std::runtime_error("CrossGenomicOperator: empty noncore block");
        }
        if (static_cast<int>(core_idx.size()) != nc ||
            static_cast<int>(noncore_idx.size()) != nn) {
            throw std::runtime_error(
                "CrossGenomicOperator: genotyped-order partition index mismatch for A22");
        }
        if (a22_full->rows() != a22_full->cols() || a22_full->cols() <= 0) {
            throw std::runtime_error("CrossGenomicOperator: A22 target must be square");
        }
        const int a22_n = static_cast<int>(a22_full->cols());
        for (int a = 0; a < nc; ++a) {
            if (core_idx[a] < 0 || core_idx[a] >= a22_n) {
                throw std::runtime_error("CrossGenomicOperator: core index out of A22 range");
            }
        }
        for (int b = 0; b < nn; ++b) {
            if (noncore_idx[b] < 0 || noncore_idx[b] >= a22_n) {
                throw std::runtime_error("CrossGenomicOperator: noncore index out of A22 range");
            }
        }
        a22_cn_.resize(nc, nn);
        for (int a = 0; a < nc; ++a) {
            const int ca = core_idx[a];
            for (int b = 0; b < nn; ++b) {
                a22_cn_(a, b) = (*a22_full)(ca, noncore_idx[b]);
            }
        }
    }
    // Pre-compute the allele stats once so concurrent first applies cannot race
    // the lazily-computed reader stats.
    reader_->computeStats(n_threads_);
    if (snp_valid != nullptr) {
        if (snp_valid->size() != static_cast<size_t>(reader_->getNumSnps())) {
            throw std::runtime_error("CrossGenomicOperator: SNP validity mask size mismatch");
        }
        snp_valid_ = *snp_valid;
    }
    if (snp_freq != nullptr) {
        if (snp_freq->size() != static_cast<size_t>(reader_->getNumSnps())) {
            throw std::runtime_error("CrossGenomicOperator: SNP frequency vector size mismatch");
        }
        snp_freq_ = *snp_freq;
    }
    if (snp_valid_.empty() != snp_freq_.empty()) {
        throw std::runtime_error(
            "CrossGenomicOperator: SNP validity and frequency state must be supplied together");
    }
}

Eigen::VectorXd CrossGenomicOperator::applyRawCross(
    const Eigen::VectorXd& x, const std::vector<int>& input_keep,
    const std::vector<int>& output_keep) const {
    Eigen::VectorXd y = Eigen::VectorXd::Zero(output_keep.size());
    const uint8_t* bed = reader_->getData() + 3;
    const size_t bytes_per_snp = reader_->getBytesPerSnp();
    const int m = reader_->getNumSnps();

#pragma omp parallel num_threads(n_threads_)
    {
        Eigen::VectorXd local = Eigen::VectorXd::Zero(output_keep.size());
#pragma omp for
        for (int j = 0; j < m; ++j) {
            if (!snp_valid_[j]) continue;
            const uint8_t* snp = bed + static_cast<size_t>(j) * bytes_per_snp;
            const double mean = static_cast<double>(static_cast<float>(2.0 * snp_freq_[j]));
            double dot = 0.0;
            for (int i = 0; i < static_cast<int>(input_keep.size()); ++i) {
                const int sample = input_keep[i];
                const int code = (snp[sample >> 2] >> ((sample & 3) * 2)) & 3;
                if (code == 0) dot += (2.0 - mean) * x(i);
                else if (code == 2) dot += (1.0 - mean) * x(i);
                else if (code == 3) dot += -mean * x(i);
            }
            for (int i = 0; i < static_cast<int>(output_keep.size()); ++i) {
                const int sample = output_keep[i];
                const int code = (snp[sample >> 2] >> ((sample & 3) * 2)) & 3;
                if (code == 0) local(i) += (2.0 - mean) * dot;
                else if (code == 2) local(i) += (1.0 - mean) * dot;
                else if (code == 3) local(i) -= mean * dot;
            }
        }
#pragma omp critical
        y += local;
    }
    if (S_2pq_ > 0.0) y /= S_2pq_;
    return y;
}

Eigen::VectorXd CrossGenomicOperator::applyCoreFromNoncore(
    const Eigen::VectorXd& x_n) const {
    if (x_n.size() != static_cast<Eigen::Index>(noncore_keep_.size())) {
        throw std::runtime_error(
            "CrossGenomicOperator::applyCoreFromNoncore: size mismatch");
    }
    Eigen::VectorXd y;
    if (!snp_valid_.empty()) {
        y = applyRawCross(x_n, noncore_keep_, core_keep_);
    } else {
        Eigen::VectorXd zt_std;
        reader_->multiply_Zt_v(x_n, zt_std, noncore_keep_, n_threads_);
        const auto& inv_stds = reader_->getInvStds();
        const int m = reader_->getNumSnps();
        Eigen::VectorXd weighted(m);
        for (int j = 0; j < m; ++j) {
            const double std_j = stdScale(inv_stds[j]);
            weighted(j) = std_j * std_j * zt_std(j);
        }
        y = multiplyCentered(*reader_, weighted, core_keep_, n_threads_);
        if (S_2pq_ > 0.0) y /= S_2pq_;
    }
    if (view_ != GcnBlockView::Raw) {
        if (tune_b_ != 1.0) y *= tune_b_;
        if (tune_a_ != 0.0) y.array() += tune_a_ * x_n.sum();
    }
    if (view_ == GcnBlockView::Final) {
        y *= (1.0 - blend_alpha_);
        if (blend_alpha_ > 0.0 && a22_cn_.size() != 0) {
            y.noalias() += blend_alpha_ * (a22_cn_ * x_n);
        }
    }
    return y;
}

Eigen::VectorXd CrossGenomicOperator::applyNoncoreFromCore(
    const Eigen::VectorXd& x_c) const {
    if (x_c.size() != static_cast<Eigen::Index>(core_keep_.size())) {
        throw std::runtime_error(
            "CrossGenomicOperator::applyNoncoreFromCore: size mismatch");
    }
    Eigen::VectorXd y;
    if (!snp_valid_.empty()) {
        y = applyRawCross(x_c, core_keep_, noncore_keep_);
    } else {
        Eigen::VectorXd zt_std;
        reader_->multiply_Zt_v(x_c, zt_std, core_keep_, n_threads_);
        const auto& inv_stds = reader_->getInvStds();
        const int m = reader_->getNumSnps();
        Eigen::VectorXd weighted(m);
        for (int j = 0; j < m; ++j) {
            const double std_j = stdScale(inv_stds[j]);
            weighted(j) = std_j * std_j * zt_std(j);
        }
        y = multiplyCentered(*reader_, weighted, noncore_keep_, n_threads_);
        if (S_2pq_ > 0.0) y /= S_2pq_;
    }
    if (view_ != GcnBlockView::Raw) {
        if (tune_b_ != 1.0) y *= tune_b_;
        if (tune_a_ != 0.0) y.array() += tune_a_ * x_c.sum();
    }
    if (view_ == GcnBlockView::Final) {
        y *= (1.0 - blend_alpha_);
        if (blend_alpha_ > 0.0 && a22_cn_.size() != 0) {
            y.noalias() += blend_alpha_ * (a22_cn_.transpose() * x_c);
        }
    }
    return y;
}

} // namespace cosmic
