#include "models/block_streamed_grm.h"
#include "qc_utils.h"
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace cosmic {

namespace {

double centeredDosage(const uint8_t* snp_ptr, int sample_index, double mean_dosage) {
    const uint8_t byte = snp_ptr[static_cast<size_t>(sample_index >> 2)];
    const int code = static_cast<int>((byte >> ((sample_index & 3) * 2)) & 3);
    if (code == 0) return 2.0 - mean_dosage;
    if (code == 2) return 1.0 - mean_dosage;
    if (code == 3) return -mean_dosage;
    return 0.0;
}

} // namespace

BlockStreamedGrmBuilder::BlockStreamedGrmBuilder(const PlinkReader& reader,
                                                  const std::vector<int>& keep_indices,
                                                  const BlockStreamedGrmBuilder::Options& opts)
    : reader_(reader), keep_indices_(keep_indices), opts_(opts) {
    n_samples_ = (int)reader_.getNumSamples();
    n_snps_ = (int)reader_.getNumSnps();
    n_geno_ = (int)keep_indices_.size();
    if (n_samples_ == 0 || n_snps_ == 0) {
        throw std::runtime_error("BlockStreamedGrmBuilder: empty dataset");
    }
    sample_is_geno_.assign(n_samples_, 0);
    for (int gp = 0; gp < n_geno_; ++gp) {
        const int s = keep_indices_[gp];
        if (s < 0 || s >= n_samples_) {
            throw std::runtime_error("BlockStreamedGrmBuilder: keep index out of range");
        }
        sample_is_geno_[s] = 1;
    }
}

void BlockStreamedGrmBuilder::computeStats() {
    auto t_start = std::chrono::steady_clock::now();
    const uint8_t* raw_data = reader_.getData();
    const size_t bytes_per_snp = reader_.getBytesPerSnp();
    const auto& bim = reader_.getBimInfo();
    const int m = n_snps_;
    const int n = n_samples_;

    audit_.n_samples = n_samples_;
    audit_.n_snps = n_snps_;
    audit_.n_genotyped = n_geno_;

    snp_freq_.assign(m, 0.0);
    snp_valid_.assign(m, 0);
    snp_2pq_.assign(m, 0.0);
    // Per-SNP partial scalar sums over the genotyped rows (deterministic: the
    // block pass accumulates the same per-SNP partials in ascending SNP order).
    std::vector<double> snp_gcs(m, 0.0), snp_gsq(m, 0.0);

    const int block = std::max(1, opts_.snp_block_size);
    const long long num_blocks = (long long)((m + block - 1) / block);
    long long n_valid = 0;

#pragma omp parallel num_threads(opts_.n_threads) reduction(+ : n_valid)
    {
#pragma omp for schedule(static)
        for (long long bi = 0; bi < num_blocks; ++bi) {
            const int start = (int)(bi * (long long)block);
            const int end = std::min(start + block, m);
            for (int j = start; j < end; ++j) {
                if (!opts_.extract_set.empty() && !opts_.extract_set.count(bim[j].id)) continue;
                if (!opts_.exclude_set.empty() && opts_.exclude_set.count(bim[j].id)) continue;

                const uint8_t* snp_ptr = raw_data + 3 + (size_t)j * bytes_per_snp;
                long long global_sum_x = 0, global_n_obs = 0;
                long long n_hom1 = 0, n_het = 0, n_hom2 = 0;  // 00->2, 10->1, 11->0
                long long i = 0;
                for (; i + 3 < n; i += 4) {
                    const uint8_t byte = snp_ptr[(size_t)(i >> 2)];
                    const uint8_t c0 = byte & 3, c1 = (byte >> 2) & 3, c2 = (byte >> 4) & 3,
                                  c3 = byte >> 6;
                    int codes[4] = {(int)c0, (int)c1, (int)c2, (int)c3};
                    for (int k = 0; k < 4; ++k) {
                        const int code = codes[k];
                        if (code == 1) continue;
                        ++global_n_obs;
                        if (code == 0) { global_sum_x += 2; ++n_hom1; }
                        else if (code == 2) { global_sum_x += 1; ++n_het; }
                        else { ++n_hom2; }
                    }
                }
                for (; i < n; ++i) {
                    const uint8_t byte = snp_ptr[(size_t)(i >> 2)];
                    const int code = (int)((byte >> ((i & 3) * 2)) & 3);
                    if (code == 1) continue;
                    ++global_n_obs;
                    if (code == 0) { global_sum_x += 2; ++n_hom1; }
                    else if (code == 2) { global_sum_x += 1; ++n_het; }
                    else { ++n_hom2; }
                }

                if (global_n_obs == 0) continue;

                const double global_freq = (double)global_sum_x / (2.0 * (double)global_n_obs);
                const double global_maf = std::min(global_freq, 1.0 - global_freq);
                if (global_maf < opts_.min_maf || std::abs(global_maf) < 1e-8) continue;
                if (opts_.min_hwe > 0.0f) {
                    const double p_hwe = calculateHWE(n_hom1, n_het, n_hom2);
                    if (p_hwe < opts_.min_hwe) continue;
                }

                snp_freq_[j] = global_freq;
                snp_valid_[j] = 1;
                // VanRaden divisor, stored float exactly as the dense builder:
                //   snp_2pq[j] = (float)(weight * 2.0 * freq * (1.0 - freq))
                snp_2pq_[j] = (double)(float)(2.0 * global_freq * (1.0 - global_freq));

                // z = dosage - 2*freq (global freq, all samples), missing -> 0;
                // the centering replicates the dense Z: mean = (float)(2*freq),
                // code0 -> 2-mean, code2 -> 1-mean, code3 -> 0-mean.
                const double m2 = (double)(float)(2.0 * global_freq);
                double gcs = 0.0, gsq = 0.0;
                for (int gp = 0; gp < n_geno_; ++gp) {
                    const int idx = keep_indices_[gp];
                    const uint8_t byte = snp_ptr[(size_t)(idx >> 2)];
                    const int code = (int)((byte >> ((idx & 3) * 2)) & 3);
                    double zv = 0.0;
                    if (code == 0) zv = 2.0 - m2;
                    else if (code == 2) zv = 1.0 - m2;
                    else if (code == 3) zv = 0.0 - m2;
                    gcs += zv;
                    gsq += zv * zv;
                }
                snp_gcs[j] = gcs;
                snp_gsq[j] = gsq;
                ++n_valid;
            }
        }
    }

    // E3-D3 audit: computeStats = 1 full-sample sweep + 1 genotyped-row sweep.
    ++audit_.full_passes;
    ++audit_.genotyped_passes;
    audit_.snp_bytes_read += (long long)bytes_per_snp * (long long)m;

    S_2pq_ = 0.0;
    raw_diag_sum_ = 0.0;
    raw_off_sum_ = 0.0;
    for (int j = 0; j < m; ++j) {
        if (!snp_valid_[j]) continue;
        S_2pq_ += snp_2pq_[j];
        raw_diag_sum_ += snp_gsq[j];
        raw_off_sum_ += snp_gcs[j] * snp_gcs[j] - snp_gsq[j];
    }
    // Same convention as the production dense G (G = ZZ'/S): the scalar trace /
    // off-diagonal sums are reported in the ZZ'/S units of the G blocks.
    if (S_2pq_ > 0.0) {
        raw_diag_sum_ /= S_2pq_;
        raw_off_sum_ /= S_2pq_;
    }
    std::cout << "Block-streamed GRM stats: " << n_valid << " valid SNPs, "
              << "S = sum 2p(1-p) = " << S_2pq_ << "\n";
    audit_.compute_stats_time_s =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
}

void BlockStreamedGrmBuilder::setTune(double a, double b) {
    tune_a_ = a;
    tune_b_ = b;
}

void BlockStreamedGrmBuilder::computeBlocks(const std::vector<int>& core_idx,
                                           const std::vector<int>& noncore_idx) {
    auto t_start = std::chrono::steady_clock::now();
    ++audit_.compute_blocks_calls;
    n_core_ = (int)core_idx.size();
    n_noncore_ = (int)noncore_idx.size();
    if (n_core_ + n_noncore_ != n_geno_) {
        throw std::runtime_error("BlockStreamedGrmBuilder: core + noncore != n_genotyped");
    }
    core_idx_ = core_idx;
    noncore_idx_ = noncore_idx;
    std::vector<int> core_pos(n_geno_, -1), noncore_pos(n_geno_, -1);
    for (int a = 0; a < n_core_; ++a) {
        if (core_idx[a] < 0 || core_idx[a] >= n_geno_ || core_pos[core_idx[a]] != -1) {
            throw std::runtime_error("BlockStreamedGrmBuilder: invalid / duplicate core index");
        }
        core_pos[core_idx[a]] = a;
    }
    for (int b = 0; b < n_noncore_; ++b) {
        if (noncore_idx[b] < 0 || noncore_idx[b] >= n_geno_ ||
            noncore_pos[noncore_idx[b]] != -1 || core_pos[noncore_idx[b]] != -1) {
            throw std::runtime_error("BlockStreamedGrmBuilder: invalid / duplicate noncore index");
        }
        noncore_pos[noncore_idx[b]] = b;
    }

    gcc_raw_.setZero(n_core_, n_core_);
    if (opts_.materialize_gcn) {
        gcn_raw_.setZero(n_core_, n_noncore_);
    } else {
        gcn_raw_.resize(0, 0);
        gcn_tuned_.resize(0, 0);
        gcn_final_.resize(0, 0);
    }
    gnn_diag_raw_.setZero(n_noncore_);

    const uint8_t* raw_data = reader_.getData();
    const size_t bytes_per_snp = reader_.getBytesPerSnp();
    const int chunk = std::max(1, opts_.gcn_chunk_size);
    const int nchunks = (n_noncore_ + chunk - 1) / chunk;
    std::vector<double> z(n_geno_);

    for (int j = 0; j < n_snps_; ++j) {
        if (!snp_valid_[j]) continue;
        const double m2 = (double)(float)(2.0 * snp_freq_[j]);
        const uint8_t* snp_ptr = raw_data + 3 + (size_t)j * bytes_per_snp;
        for (int gp = 0; gp < n_geno_; ++gp) {
            z[gp] = centeredDosage(snp_ptr, keep_indices_[gp], m2);
        }

        // Gcc rank-1 (n_core^2 per SNP).
        for (int a = 0; a < n_core_; ++a) {
            const double za = z[core_idx[a]];
            for (int b = 0; b < n_core_; ++b) gcc_raw_(a, b) += za * z[core_idx[b]];
        }

        // Gcn + diag(Gnn) rank-1 over disjoint noncore chunks (each element
        // written by exactly one thread in SNP order -> thread-count stable).
        if (n_noncore_ > 0) {
#pragma omp parallel num_threads(opts_.n_threads)
            {
#pragma omp for schedule(static)
                for (int c = 0; c < nchunks; ++c) {
                    const int cs = c * chunk;
                    const int ce = std::min(cs + chunk, n_noncore_);
                    for (int b = cs; b < ce; ++b) {
                        const double zb = z[noncore_idx[b]];
                        if (opts_.materialize_gcn) {
                            double* col = &gcn_raw_(0, b);
                            for (int a = 0; a < n_core_; ++a) col[a] += z[core_idx[a]] * zb;
                        }
                        gnn_diag_raw_(b) += zb * zb;
                    }
                }
            }
        }
    }

    if (S_2pq_ > 0.0) {
        gcc_raw_ /= S_2pq_;
        if (opts_.materialize_gcn) gcn_raw_ /= S_2pq_;
        gnn_diag_raw_ /= S_2pq_;
    }
    // E3-D3 audit: computeBlocks = 1 genotyped-row sweep (no full-sample pass);
    // genotype bytes touched == the full BED arena (each SNP column touched).
    ++audit_.genotyped_passes;
    audit_.snp_bytes_read += (long long)bytes_per_snp * (long long)n_snps_;
    audit_.n_core = n_core_;
    audit_.n_noncore = n_noncore_;
    audit_.compute_blocks_time_s =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    blocks_computed_ = true;
    finalized_ = false;
}

void BlockStreamedGrmBuilder::finalize(double blend_alpha, const Eigen::MatrixXd* a22_full) {
    const int ng = n_geno_;
    const int nc = n_core_;
    const int nn = n_noncore_;

    if (!blocks_computed_) {
        throw std::runtime_error("BlockStreamedGrmBuilder: computeBlocks must precede finalize");
    }
    blend_alpha_ = blend_alpha;
    blend_target_ = a22_full;

    // Tuned blocks: G_tuned = b*G_raw + a*J (elementwise +a).
    gcc_tuned_ = tune_b_ * gcc_raw_;
    gcc_tuned_.array() += tune_a_;
    if (opts_.materialize_gcn) {
        gcn_tuned_ = tune_b_ * gcn_raw_;
        gcn_tuned_.array() += tune_a_;
    }
    gnn_diag_tuned_ = tune_b_ * gnn_diag_raw_;
    gnn_diag_tuned_.array() += tune_a_;

    // Tuned scalar sums over the genotyped block (i != j off-diagonals of the
    // full ng x ng genotyped block, which includes Gnn).
    const double tune_diag_sum = tune_b_ * raw_diag_sum_ + tune_a_ * (double)ng;
    const double tune_off_sum = tune_b_ * raw_off_sum_ + tune_a_ * (double)ng * (double)(ng - 1);

    if (blend_alpha > 0.0 && a22_full != nullptr) {
        const int n22 = (int)a22_full->cols();
        if (a22_full->rows() != n22 || n22 != ng) {
            throw std::runtime_error("BlockStreamedGrmBuilder: A22 must be n_g x n_g");
        }
        Eigen::MatrixXd a22_cc(nc, nc);
        Eigen::VectorXd a22_nn_diag(nn);
        for (int a = 0; a < nc; ++a) {
            const int ca = core_idx_[a];
            for (int b = 0; b < nc; ++b) a22_cc(a, b) = (*a22_full)(ca, core_idx_[b]);
        }
        double sum_A_off = 0.0;
        for (int b = 0; b < nn; ++b) {
            const int nb = noncore_idx_[b];
            a22_nn_diag(b) = (*a22_full)(nb, nb);
        }
        for (int i = 0; i < n22; ++i) {
            for (int j = 0; j < n22; ++j) {
                if (i != j) sum_A_off += (*a22_full)(i, j);
            }
        }
        const double cm = 1.0 - blend_alpha;
        gcc_final_ = cm * gcc_tuned_ + blend_alpha * a22_cc;
        if (opts_.materialize_gcn) {
            gcn_final_ = cm * gcn_tuned_;
            for (int a = 0; a < nc; ++a) {
                const int ca = core_idx_[a];
                for (int b = 0; b < nn; ++b) {
                    gcn_final_(a, b) += blend_alpha * (*a22_full)(ca, noncore_idx_[b]);
                }
            }
        }
        gnn_diag_final_ = cm * gnn_diag_tuned_ + blend_alpha * a22_nn_diag;
        final_diag_sum_ = cm * tune_diag_sum + blend_alpha * a22_full->trace();
        final_off_sum_ = cm * tune_off_sum + blend_alpha * sum_A_off;
    } else {
        // Identity blend: G_final = (1-alpha)*G_tuned + alpha*I.
        const double cm = 1.0 - blend_alpha;
        gcc_final_ = cm * gcc_tuned_;
        gcc_final_.diagonal().array() += blend_alpha;
        if (opts_.materialize_gcn) gcn_final_ = cm * gcn_tuned_;
        gnn_diag_final_ = cm * gnn_diag_tuned_;
        gnn_diag_final_.array() += blend_alpha;
        final_diag_sum_ = cm * tune_diag_sum + blend_alpha * (double)ng;
        final_off_sum_ = cm * tune_off_sum;
    }
    finalized_ = true;
}

void BlockStreamedGrmBuilder::generateGcnBlock(int start_col, int end_col,
                                                GcnBlockView view,
                                                Eigen::MatrixXd& block) const {
    if (!blocks_computed_ || start_col < 0 || end_col < start_col || end_col > n_noncore_) {
        throw std::runtime_error("BlockStreamedGrmBuilder: invalid genotype-backed Gcn block request");
    }
    if (view != GcnBlockView::Raw && !finalized_) {
        throw std::runtime_error("BlockStreamedGrmBuilder: finalize must precede tuned/final Gcn iteration");
    }
    const int cols = end_col - start_col;
    block.setZero(n_core_, cols);
    if (cols == 0) return;
    ++gcn_regen_calls_;

    const uint8_t* raw_data = reader_.getData();
    const size_t bytes_per_snp = reader_.getBytesPerSnp();
    std::vector<double> z_core(n_core_);
    for (int j = 0; j < n_snps_; ++j) {
        if (!snp_valid_[j]) continue;
        const double m2 = static_cast<double>(static_cast<float>(2.0 * snp_freq_[j]));
        const uint8_t* snp_ptr = raw_data + 3 + static_cast<size_t>(j) * bytes_per_snp;
        for (int a = 0; a < n_core_; ++a) {
            z_core[a] = centeredDosage(
                snp_ptr, keep_indices_[core_idx_[a]], m2);
        }
        for (int b = start_col; b < end_col; ++b) {
            const double zb = centeredDosage(
                snp_ptr, keep_indices_[noncore_idx_[b]], m2);
            double* col = &block(0, b - start_col);
            for (int a = 0; a < n_core_; ++a) col[a] += z_core[a] * zb;
        }
    }
    if (S_2pq_ > 0.0) block /= S_2pq_;
    if (view == GcnBlockView::Raw) return;

    block *= tune_b_;
    block.array() += tune_a_;
    if (view == GcnBlockView::Tuned) return;

    block *= (1.0 - blend_alpha_);
    if (blend_alpha_ > 0.0 && blend_target_ != nullptr) {
        for (int a = 0; a < n_core_; ++a) {
            const int ca = core_idx_[a];
            for (int b = start_col; b < end_col; ++b) {
                block(a, b - start_col) +=
                    blend_alpha_ * (*blend_target_)(ca, noncore_idx_[b]);
            }
        }
    }
}

GcnBlockIteratorImpl::GcnBlockIteratorImpl(const BlockStreamedGrmBuilder* builder,
                                           int block_size, GcnBlockView view)
    : builder_(builder),
      block_size_(block_size > 0 ? block_size : (builder ? builder->opts_.gcn_chunk_size : 4096)),
      view_(view) {
    if (builder_) n_chunks_ = (builder_->nNoncore() + block_size_ - 1) / block_size_;
}

bool GcnBlockIteratorImpl::next(GcnBlock& out_block) {
    if (current_chunk_ >= n_chunks_ || !builder_) return false;
    const int start = current_chunk_ * block_size_;
    const int end = std::min(start + block_size_, builder_->nNoncore());
    out_block.start_col = start;
    out_block.end_col = end;
    out_block.noncore_indices.assign(builder_->noncore_idx_.begin() + start,
                                     builder_->noncore_idx_.begin() + end);
    if (builder_->opts_.materialize_gcn) {
        const Eigen::MatrixXd* source = &builder_->gcn_raw_;
        if (view_ == GcnBlockView::Tuned) source = &builder_->gcn_tuned_;
        if (view_ == GcnBlockView::Final) source = &builder_->gcn_final_;
        out_block.block = source->middleCols(start, end - start);
    } else {
        builder_->generateGcnBlock(start, end, view_, out_block.block);
    }
    ++current_chunk_;
    return true;
}

GcnBlockIteratorPtr BlockStreamedGrmBuilder::makeGcnIterator(int block_size,
                                                              GcnBlockView view) const {
    if (!blocks_computed_) {
        throw std::runtime_error("BlockStreamedGrmBuilder: computeBlocks must precede Gcn iteration");
    }
    if (view != GcnBlockView::Raw && !finalized_) {
        throw std::runtime_error("BlockStreamedGrmBuilder: finalize must precede tuned/final Gcn iteration");
    }
    return std::make_unique<GcnBlockIteratorImpl>(this, block_size, view);
}

std::string BlockStreamedGrmBuilder::streamAuditSummary() const {
    std::ostringstream os;
    os << "COMPUTE_BLOCKS_CALLS=" << audit_.compute_blocks_calls
       << " COMPUTE_BLOCKS_TIME=" << audit_.compute_blocks_time_s << "s"
       << " COMPUTE_STATS_TIME=" << audit_.compute_stats_time_s << "s"
       << " GENOTYPE_PASSES=" << (audit_.full_passes + audit_.genotyped_passes)
       << " (full=" << audit_.full_passes << ", genotyped=" << audit_.genotyped_passes << ")"
       << " SNP_BYTES_TOUCHED=" << audit_.snp_bytes_read
       << " n_samples=" << audit_.n_samples
       << " n_snps=" << audit_.n_snps
       << " n_genotyped=" << audit_.n_genotyped
       << " n_core=" << audit_.n_core
       << " n_noncore=" << audit_.n_noncore;
    return os.str();
}

} // namespace cosmic
