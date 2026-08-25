#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <string>
#include <vector>
#include "matrix_adapter.h"

namespace cosmic {

// ---------------------------------------------------------------------------
// APY (Algorithm for Proven and Young) block inverse of a genomic relationship
// matrix  G  (n_g x n_g) over a core/noncore partition.
//
// Reference: Leite et al. (2024) GSE 56:59 (PMC11340074), Eqs. 5-6:
//
//   G_APY^-1 = [ Gcc^-1 + Q Mnn^-1 Q'   -Q Mnn^-1 ]
//              [ -Mnn^-1 Q'             Mnn^-1    ]
//
//   Q   = Gcc^-1 Gcn          (core x noncore)
//   Mnn = diag(Gnn - Gnc Gcc^-1 Gcn)  (noncore conditional variances, kept
//         diagonal; the published APY approximation)
//
// G_APY^-1 is the EXACT inverse of the modified matrix
//   G_apy = [ Gcc  Gcn ; Gnc  Mnn + Q' Gcn ] ,
// i.e. G with the noncore off-diagonal block replaced by its conditional
// prediction given the core (Q = Gcc^-1 Gcn); it equals G itself when the
// noncore residuals  Gnn - Gnc Gcc^-1 Gcn  are diagonal (the APY conditional-
// independence assumption).  With the FULL core (n_core = n_genotyped)
// G_apy = G and G_APY^-1 = G^-1 (exact degeneracy; C2 gate G2).
//
// Storage is O(n_core^2 + n_core*n_noncore): no dense n_g^2 allocation of the
// inverse is performed.  The dense Gcc^-1 and the core-block
//   [Gcc^-1 + Q Mnn^-1 Q']                                 (core x core)
// are materialized on build; the MME/applications use the matrix-free
// applyInverse (O(n_core^2 + n_core*n_noncore) per matvec).
// ---------------------------------------------------------------------------
class ApyGinv {
public:
    struct Info {
        int n_core = 0;
        int n_noncore = 0;
        int n_mnn_nonpositive = 0;  // conditional-variance violations (build fails)
        double mnn_min = 0.0;       // smallest noncore conditional variance
        double mnn_max = 0.0;
        double gcc_min_diag_d = 0.0;  // LDLT pivot diagnostics of Gcc
        double gcc_max_diag_d = 0.0;
        size_t dense_nnz = 0;       // nonzero count of the full G_APY^-1 pattern
        bool core_is_all = false;   // n_core == n_genotyped (exact degeneracy)
    };

    // Core indices: unique ascending values in [0, n_g); n_core may equal n_g
    // (exact degeneracy; noncore empty).  Throws on invalid partitions or when
    // Gcc is not SPD or any noncore conditional variance is <= 0.
    void build(const Eigen::MatrixXd& G, const std::vector<int>& core_idx,
               const std::string& tag);

    int n() const { return n_core_ + n_noncore_; }
    int nCore() const { return n_core_; }
    int nNoncore() const { return n_noncore_; }
    const std::vector<int>& coreIndices() const { return core_idx_; }
    const std::vector<int>& noncoreIndices() const { return noncore_idx_; }
    const Info& info() const { return info_; }

    // y = G_APY^-1 x  (x in genotyped order).
    Eigen::VectorXd applyInverse(const Eigen::VectorXd& x) const;

    // y = Gcc^-1 x_c  where x_c is a core-ordered (n_core) vector; used by the
    // per-SNP core weights  w_c = Gcc^-1 M_cj  (Leite Eq. 7 / Bermann 2022a
    // Eq. 21-22 form).
    Eigen::VectorXd applyCoreInverse(const Eigen::VectorXd& x_core) const;

    // Dense core x core block of G_APY^-1:  Gcc^-1 + Q Mnn^-1 Q'.
    const Eigen::MatrixXd& coreBlock() const { return core_block_; }

    // Phase D (Eq. 8 Schur form): Gcc^-1, Q = Gcc^-1 Gcn, and the noncore
    // conditional variances m_jj = g_jj - g_jc' Gcc^-1 g_cj (diag of Mnn).
    const Eigen::MatrixXd& gccInv() const { return gcc_inv_; }
    const Eigen::MatrixXd& q() const { return q_; }
    const Eigen::VectorXd& mnn() const { return mnn_; }

    // Full pattern of G_APY^-1 (both orientations, symmetric values).
    void visitTriplets(const std::function<void(int, int, double)>& cb) const;

    // Diagonals of G_APY^-1 into diag (alpha-scaled; genotyped order).
    void addDiagonalTo(Eigen::VectorXd& diag, double alpha) const;

    const std::string& tag() const { return tag_; }

private:
    std::string tag_;
    std::vector<int> core_idx_;       // ascending
    std::vector<int> noncore_idx_;    // ascending
    int n_core_ = 0;
    int n_noncore_ = 0;
    Eigen::MatrixXd gcc_inv_;     // n_core^2  dense inverse of Gcc
    Eigen::MatrixXd q_;           // n_core x n_noncore, Q = Gcc^-1 Gcn
    Eigen::VectorXd mnn_;         // n_noncore, conditional variances m_jj
    Eigen::VectorXd mnn_inv_;     // n_noncore, 1 / m_jj
    Eigen::MatrixXd core_block_;  // n_core^2 = Gcc_inv + Q Mnn^-1 Q'
    Info info_;
};

// AbstractMatrix adapter of G_APY^-1 for the MME (SplitMatrixAdapter, PCG,
// LHS builder).  scale = tau (h_tau) of the Hinv genotyped block:
// maps v -> scale * G_APY^-1 v.
class ApyGinvAdapter : public AbstractMatrix {
public:
    ApyGinvAdapter(const ApyGinv* g, double scale) : g_(g), scale_(scale) {}

    int rows() const override { return g_->n(); }
    int cols() const override { return g_->n(); }
    size_t nonZeros() const override { return g_->info().dense_nnz; }

    Eigen::VectorXd operator*(const Eigen::VectorXd& x) const override {
        return scale_ * g_->applyInverse(x);
    }

    void add_diagonal_to(Eigen::VectorXd& diag, double alpha) const override {
        g_->addDiagonalTo(diag, alpha * scale_);
    }

    Eigen::MatrixXd getBlock(int start_row, int size) const override {
        (void)start_row; (void)size;
        throw std::runtime_error(
            "ApyGinvAdapter::getBlock is not supported (use applyInverse / visitTriplets)");
    }

    void visit_triplets(std::function<void(int, int, double)> callback) const override {
        g_->visitTriplets([&](int r, int c, double v) { callback(r, c, scale_ * v); });
    }

    Eigen::MatrixXd toDense() const override {
        throw std::runtime_error(
            "ApyGinvAdapter::toDense is not supported (no dense n_g^2 materialization); "
            "use applyInverse / visitTriplets");
    }

private:
    const ApyGinv* g_;
    double scale_;
};

} // namespace cosmic