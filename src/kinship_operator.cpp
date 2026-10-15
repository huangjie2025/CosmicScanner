#include "kinship_operator.h"
#include "plink_reader.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstring>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace cosmic {

namespace {

// PLINK 2-bit code -> genotype dosage.  Must match PlinkReader::lookup_table.
constexpr float kDosage[4] = {2.0f, 0.0f, 1.0f, 0.0f};
constexpr double kObserved[4] = {1.0, 0.0, 1.0, 1.0};
constexpr uint8_t kMissing = 1;

inline int codeAt(const uint8_t* col, int sample) {
    return (col[sample >> 2] >> ((sample & 3) * 2)) & 0x03;
}

} // namespace

PlinkKinshipOperator::PlinkKinshipOperator(PlinkReader* reader,
                                           std::vector<int> keep_indices,
                                           std::vector<int> variant_indices)
    : reader_(reader),
      keep_(std::move(keep_indices)),
      variants_(std::move(variant_indices)) {
    if (reader_ == nullptr) throw std::runtime_error("PlinkKinshipOperator: null reader");

    const int total_samples = static_cast<int>(reader_->getNumSamples());
    const int total_snps = static_cast<int>(reader_->getNumSnps());

    if (!keep_.empty()) {
        for (int s : keep_) {
            if (s < 0 || s >= total_samples) {
                throw std::runtime_error("PlinkKinshipOperator: sample index out of range");
            }
        }
        identity_samples_ = (static_cast<int>(keep_.size()) == total_samples);
        for (int i = 0; identity_samples_ && i < total_samples; ++i) {
            identity_samples_ = keep_[i] == i;
        }
    } else {
        identity_samples_ = true;
    }

    if (!variants_.empty()) {
        for (int j : variants_) {
            if (j < 0 || j >= total_snps) {
                throw std::runtime_error("PlinkKinshipOperator: variant index out of range");
            }
        }
    }

    num_samples_ = keep_.empty() ? total_samples : static_cast<int>(keep_.size());
    num_variants_ = variants_.empty() ? total_snps : static_cast<int>(variants_.size());
    if (num_variants_ == 0) throw std::runtime_error("PlinkKinshipOperator: no variants selected");
    scaling_ = 1.0 / static_cast<double>(num_variants_);

    // Genotype statistics (means / inverse standard deviations) are required.
    reader_->computeStats(1);
}

int PlinkKinshipOperator::rows() const { return num_samples_; }
int PlinkKinshipOperator::cols() const { return num_samples_; }
int PlinkKinshipOperator::numVariants() const { return num_variants_; }

void PlinkKinshipOperator::multiplyZt(const Eigen::VectorXd& v, Eigen::VectorXd& u,
                                      int n_threads) const {
    if (n_threads < 1) throw std::invalid_argument("Kinship vector threads must be positive");
    if (static_cast<int>(v.size()) != num_samples_) {
        throw std::runtime_error("PlinkKinshipOperator::multiplyZt: dimension mismatch");
    }
    const int m = num_variants_;
    u.setZero(m);

    const uint8_t* base = reader_->getData() + 3;  // skip BED magic bytes
    const size_t bytes_per_snp = reader_->getBytesPerSnp();
    const std::vector<double>& means = reader_->getMeans();
    const std::vector<double>& inv_stds = reader_->getInvStds();
    const bool identity = identity_samples_;
    const std::vector<int>& keep = keep_;

    const int variant_begin = variants_.empty() ? 0 : 0;
    (void)variant_begin;

    #pragma omp parallel for num_threads(n_threads) schedule(static) if (n_threads > 1)
    for (int t = 0; t < m; ++t) {
        const int j = variants_.empty() ? t : variants_[static_cast<std::size_t>(t)];
        const double inv_std = inv_stds[static_cast<std::size_t>(j)];
        if (inv_std == 0.0) continue;  // monomorphic or fully missing
        const double mean = means[static_cast<std::size_t>(j)];
        const uint8_t* col = base + static_cast<std::size_t>(j) * bytes_per_snp;

        double sum_xv = 0.0;
        double sum_v = 0.0;
        if (identity) {
            // Byte-unrolled, branch-free.  Missing genotypes contribute zero to
            // both accumulators, so they are excluded from sum_xv - mean*sum_v.
            const int full_bytes = num_samples_ >> 2;
            const int tail = num_samples_ & 3;
            for (int b = 0; b < full_bytes; ++b) {
                const uint8_t byte = col[b];
                const int r = b << 2;
                const double v0 = v[r];
                const double v1 = v[r + 1];
                const double v2 = v[r + 2];
                const double v3 = v[r + 3];
                sum_xv += v0 * kDosage[byte & 3] + v1 * kDosage[(byte >> 2) & 3] +
                          v2 * kDosage[(byte >> 4) & 3] + v3 * kDosage[(byte >> 6) & 3];
                sum_v += v0 * kObserved[byte & 3] + v1 * kObserved[(byte >> 2) & 3] +
                         v2 * kObserved[(byte >> 4) & 3] + v3 * kObserved[(byte >> 6) & 3];
            }
            for (int r = full_bytes << 2; r < num_samples_; ++r) {
                const int code = codeAt(col, r);
                if (code == kMissing) continue;
                sum_xv += static_cast<double>(kDosage[code]) * v[r];
                sum_v += v[r];
            }
            (void)tail;
        } else {
            for (int r = 0; r < num_samples_; ++r) {
                const int code = codeAt(col, keep[static_cast<std::size_t>(r)]);
                if (code == kMissing) continue;
                sum_xv += static_cast<double>(kDosage[code]) * v[r];
                sum_v += v[r];
            }
        }
        u[t] = (sum_xv - mean * sum_v) * inv_std;
    }
}

void PlinkKinshipOperator::multiplyZ(const Eigen::VectorXd& u, Eigen::VectorXd& y,
                                     int n_threads) const {
    if (n_threads < 1) throw std::invalid_argument("Kinship vector threads must be positive");
    if (static_cast<int>(u.size()) != num_variants_) {
        throw std::runtime_error("PlinkKinshipOperator::multiplyZ: dimension mismatch");
    }
    const int m = num_variants_;
    y.setZero(num_samples_);

    const uint8_t* base = reader_->getData() + 3;
    const size_t bytes_per_snp = reader_->getBytesPerSnp();
    const std::vector<double>& means = reader_->getMeans();
    const std::vector<double>& inv_stds = reader_->getInvStds();
    const bool identity = identity_samples_;
    const std::vector<int>& keep = keep_;

    // Each worker owns an output vector. Reduce in worker-index order after
    // the parallel region; never update the shared output from variant workers.
    const auto accumulate = [&](int t, double* dest) {
        const int j = variants_.empty() ? t : variants_[static_cast<std::size_t>(t)];
        const double inv_std = inv_stds[static_cast<std::size_t>(j)];
        if (inv_std == 0.0) return;
        const double mean = means[static_cast<std::size_t>(j)];
        const double weight = u[t] * inv_std;
        const uint8_t* col = base + static_cast<std::size_t>(j) * bytes_per_snp;

        if (identity) {
            // Branch-free byte-unrolled accumulation; missing genotypes map to 0.
            const double vmap[4] = {(2.0 - mean) * weight, 0.0,
                                    (1.0 - mean) * weight, (0.0 - mean) * weight};
            const int full_bytes = num_samples_ >> 2;
            for (int b = 0; b < full_bytes; ++b) {
                const uint8_t byte = col[b];
                const int r = b << 2;
                dest[r]     += vmap[byte & 3];
                dest[r + 1] += vmap[(byte >> 2) & 3];
                dest[r + 2] += vmap[(byte >> 4) & 3];
                dest[r + 3] += vmap[(byte >> 6) & 3];
            }
            for (int r = full_bytes << 2; r < num_samples_; ++r) {
                dest[r] += vmap[codeAt(col, r)];
            }
        } else {
            for (int r = 0; r < num_samples_; ++r) {
                const int code = codeAt(col, keep[static_cast<std::size_t>(r)]);
                if (code == kMissing) continue;
                dest[r] += (static_cast<double>(kDosage[code]) - mean) * weight;
            }
        }
    };
    int workers = 1;
#ifdef _OPENMP
    workers = std::min(n_threads, m);
#endif
    if (workers == 1) {
        for (int t = 0; t < m; ++t) accumulate(t, y.data());
        return;
    }
    Eigen::MatrixXd local = Eigen::MatrixXd::Zero(num_samples_, workers);
    #pragma omp parallel num_threads(workers)
    {
        int worker = 0;
#ifdef _OPENMP
        worker = omp_get_thread_num();
#endif
        double* dest = local.col(worker).data();
        #pragma omp for schedule(static)
        for (int t = 0; t < m; ++t) accumulate(t, dest);
    }
    for (int worker = 0; worker < workers; ++worker) y += local.col(worker);
}

void PlinkKinshipOperator::apply(const Eigen::VectorXd& v, Eigen::VectorXd& out,
                                 int n_threads) const {
    Eigen::VectorXd u;
    multiplyZt(v, u, n_threads);
    Eigen::VectorXd w;
    multiplyZ(u, w, n_threads);
    out.noalias() = scaling_ * w;
}

void PlinkKinshipOperator::apply(const Eigen::MatrixXd& X, Eigen::MatrixXd& out,
                                 int n_threads) const {
    apply_block(X, out, n_threads);
}

void PlinkKinshipOperator::apply_block(const Eigen::MatrixXd& X, Eigen::MatrixXd& out,
                                      int n_threads, BlockProfile* profile) const {
    if (X.rows() != num_samples_) {
        throw std::runtime_error("PlinkKinshipOperator::apply: dimension mismatch");
    }
    if (n_threads < 1) throw std::invalid_argument("Kinship block threads must be positive");
    if (profile) *profile = BlockProfile{};
    // Keep out untouched until completion, including when X aliases out.
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(num_samples_, X.cols());
    if (X.cols() == 0) {
        out.swap(result);
        return;
    }
    constexpr int block_size = 512;
    const auto* bed = reader_->getData() + 3;
    const size_t stride = reader_->getBytesPerSnp();
    const auto& means = reader_->getMeans();
    const auto& inv_stds = reader_->getInvStds();
    std::vector<uint8_t> packed(stride * std::min(block_size, num_variants_));
    Eigen::MatrixXd z, projected;
    using Clock = std::chrono::steady_clock;
    const auto elapsed = [](Clock::time_point start) {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };
    for (int first = 0; first < num_variants_; first += block_size) {
        const int count = std::min(block_size, num_variants_ - first);
        z.resize(num_samples_, count);
        // Copy selected packed columns once. Decoding then reuses the bounded
        // staging buffer; mmap access and decoding can be timed independently.
        auto start = profile ? Clock::now() : Clock::time_point{};
        for (int j = 0; j < count; ++j) {
            const int snp = variants_.empty() ? first+j : variants_[first+j];
            std::memcpy(packed.data() + size_t(j)*stride, bed + size_t(snp)*stride, stride);
        }
        if (profile) {
            profile->packed_read_seconds += elapsed(start);
            profile->packed_bytes += size_t(count)*stride;
            ++profile->blocks;
        }
        start = profile ? Clock::now() : Clock::time_point{};
        #pragma omp parallel for num_threads(n_threads) schedule(static) if (n_threads > 1)
        for (int j = 0; j < count; ++j) {
            const int snp = variants_.empty() ? first+j : variants_[first+j];
            const double mean = means[snp], inv_std = inv_stds[snp];
            double* dest = z.col(j).data();
            const auto* col = packed.data() + size_t(j)*stride;
            if (inv_std == 0.0) {
                std::fill(dest, dest + num_samples_, 0.0);
                continue;
            }
            const double values[4] = {(2.0-mean)*inv_std, 0.0,
                                      (1.0-mean)*inv_std, -mean*inv_std};
            for (int r = 0; r < num_samples_; ++r) {
                // Match the existing vector operator's selection semantics.
                const int sample = identity_samples_ ? r : keep_[r];
                dest[r] = values[codeAt(col, sample)];
            }
        }
        if (profile) {
            profile->decode_seconds += elapsed(start);
            for (int j = 0; j < count; ++j) {
                const int snp = variants_.empty() ? first+j : variants_[first+j];
                if (inv_stds[snp] != 0.0) profile->decoded_genotypes += num_samples_;
            }
        }
        start = profile ? Clock::now() : Clock::time_point{};
        projected.noalias() = z.transpose() * X;
        if (profile) profile->transpose_product_seconds += elapsed(start);
        start = profile ? Clock::now() : Clock::time_point{};
        result.noalias() += z * projected;
        if (profile) profile->forward_product_seconds += elapsed(start);
    }
    result *= scaling_;
    out.swap(result);
}

void PlinkKinshipOperator::applyTranspose(const Eigen::VectorXd& v, Eigen::VectorXd& out,
                                          int n_threads) const {
    apply(v, out, n_threads);
}

double PlinkKinshipOperator::trace() const {
    const int m = num_variants_;
    const uint8_t* base = reader_->getData() + 3;
    const size_t bytes_per_snp = reader_->getBytesPerSnp();
    const std::vector<double>& means = reader_->getMeans();
    const std::vector<double>& inv_stds = reader_->getInvStds();
    const bool identity = identity_samples_;
    const std::vector<int>& keep = keep_;

    double total = 0.0;
    #pragma omp parallel for num_threads(1) schedule(static)
    for (int t = 0; t < m; ++t) {
        const int j = variants_.empty() ? t : variants_[static_cast<std::size_t>(t)];
        const double inv_std = inv_stds[static_cast<std::size_t>(j)];
        if (inv_std == 0.0) continue;
        const double mean = means[static_cast<std::size_t>(j)];
        const uint8_t* col = base + static_cast<std::size_t>(j) * bytes_per_snp;

        double sum_sq = 0.0;
        if (identity) {
            for (int r = 0; r < num_samples_; ++r) {
                const int code = codeAt(col, r);
                if (code == kMissing) continue;
                const double z = (static_cast<double>(kDosage[code]) - mean) * inv_std;
                sum_sq += z * z;
            }
        } else {
            for (int r = 0; r < num_samples_; ++r) {
                const int code = codeAt(col, keep[static_cast<std::size_t>(r)]);
                if (code == kMissing) continue;
                const double z = (static_cast<double>(kDosage[code]) - mean) * inv_std;
                sum_sq += z * z;
            }
        }
        #pragma omp atomic
        total += sum_sq;
    }
    return scaling_ * total;
}

void PlinkKinshipOperator::getDiagonal(Eigen::VectorXd& diag) const {
    diag.setZero(num_samples_);
    const int m = num_variants_;
    const uint8_t* base = reader_->getData() + 3;
    const size_t bytes_per_snp = reader_->getBytesPerSnp();
    const std::vector<double>& means = reader_->getMeans();
    const std::vector<double>& inv_stds = reader_->getInvStds();

    for (int t = 0; t < m; ++t) {
        const int j = variants_.empty() ? t : variants_[static_cast<std::size_t>(t)];
        const double inv_std = inv_stds[static_cast<std::size_t>(j)];
        if (inv_std == 0.0) continue;
        const double mean = means[static_cast<std::size_t>(j)];
        const uint8_t* col = base + static_cast<std::size_t>(j) * bytes_per_snp;
        for (int r = 0; r < num_samples_; ++r) {
            const int sample = keep_.empty() ? r : keep_[static_cast<std::size_t>(r)];
            const int code = codeAt(col, sample);
            if (code == kMissing) continue;
            const double z = (static_cast<double>(kDosage[code]) - mean) * inv_std;
            diag[r] += scaling_ * z * z;
        }
    }
}

std::unique_ptr<KinshipOperator> createKinshipOperator(PlinkReader* reader,
                                                       const std::vector<int>& keep_indices,
                                                       const std::vector<int>& variant_indices) {
    return std::make_unique<PlinkKinshipOperator>(reader, keep_indices, variant_indices);
}

} // namespace cosmic
