#pragma once

#include <Eigen/Dense>

#include <vector>

#include "plink_reader.h"
#include "block_streamed_grm.h"

namespace cosmic {

// Matrix-free cross-block genomic operator for the APY streamed path
// (Phase E3-D2D).
//
// Reproduces the K x Nn cross block Gcn = Zc Zn' / S of the block-streamed
// Tuned/Final G directly from the PLINK reader primitives
// (O(M * (K + Nn)) per apply) instead of regenerating the full K x Nn block on
// every APY MME / PEV operator apply:
//   Raw   : y = G_raw * v = Zc diag(std) Zn' v / S  (Gcn)          [K output]
//           y = G_raw' * v = Zn diag(std) Zc' v / S (Gcn')          [Nn output]
//   Tuned : y = b * raw + a * sum(v) * ones   (preGSf90 tuneG, J block)
//   Final : y = (1 - alpha) * tuned; when blending with A22 the dense target
//           contributes  alpha * A22_cn * v (identity blend contributes zero
//           to the cross block, matching BlockStreamedGrmBuilder::finalize).
//
// keep arrays are SAMPLE indices (see BlockStreamedGrmBuilder::keep_indices_);
// core_idx / noncore_idx are genotyped-order indices used to slice the dense
// A22 target (kept as a pre-sliced K x Nn matrix; the A22 sub-block matvec is a
// documented bounded fallback until a scalable A22 operator is added).
class CrossGenomicOperator {
public:
    CrossGenomicOperator(PlinkReader& reader,
                         std::vector<int> core_keep,
                         std::vector<int> noncore_keep,
                         double S_2pq,
                         double tune_a,
                         double tune_b,
                         double blend_alpha,
                         GcnBlockView view,
                         const Eigen::MatrixXd* a22_full,
                         std::vector<int> core_idx,
                         std::vector<int> noncore_idx,
                         int n_threads = 1,
                         const std::vector<char>* snp_valid = nullptr,
                         const std::vector<double>* snp_freq = nullptr);

    CrossGenomicOperator(const CrossGenomicOperator&) = delete;
    CrossGenomicOperator& operator=(const CrossGenomicOperator&) = delete;

    // y_core = Gcn * x_n (K entries).
    Eigen::VectorXd applyCoreFromNoncore(const Eigen::VectorXd& x_n) const;
    // y_noncore = Gcn' * x_c (Nn entries).
    Eigen::VectorXd applyNoncoreFromCore(const Eigen::VectorXd& x_c) const;

    int coreRows() const { return static_cast<int>(core_keep_.size()); }
    int noncoreRows() const { return static_cast<int>(noncore_keep_.size()); }

private:
    PlinkReader* reader_ = nullptr;
    std::vector<int> core_keep_;
    std::vector<int> noncore_keep_;
    double S_2pq_ = 0.0;
    double tune_a_ = 0.0;
    double tune_b_ = 1.0;
    double blend_alpha_ = 0.0;
    GcnBlockView view_ = GcnBlockView::Raw;
    Eigen::MatrixXd a22_cn_;   // K x Nn pre-sliced A22 target (empty when unused)
    int n_threads_ = 1;
    std::vector<char> snp_valid_;
    std::vector<double> snp_freq_;

    Eigen::VectorXd applyRawCross(const Eigen::VectorXd& x,
                                  const std::vector<int>& input_keep,
                                  const std::vector<int>& output_keep) const;
};

} // namespace cosmic
