#pragma once
#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <vector>
#include "plink_reader.h"

namespace cosmic {

// Abstract interface for a matrix-free kinship operator.
//
// Represents a genetic relationship matrix
//     K = (1/M) * Z * Z'
// over a selected set of samples and variants, without ever materializing the
// N x N matrix.  Z holds mean-imputed, variance-standardized genotypes
//     Z_ij = (X_ij - mean_j) * inv_std_j
// with missing genotypes contributing zero.
class KinshipOperator {
public:
    virtual ~KinshipOperator() = default;

    // Number of samples (N) represented by the operator
    virtual int rows() const = 0;
    virtual int cols() const = 0;

    // Number of variants (M) contributing to the GRM
    virtual int numVariants() const = 0;

    // Core operation: out = K * v
    virtual void apply(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1) const = 0;

    // Apply multiple vectors at once: out = K * X (X is rows() x k)
    virtual void apply(const Eigen::MatrixXd& X, Eigen::MatrixXd& out, int n_threads = 1) const = 0;

    // K is symmetric, so K' == K
    virtual void applyTranspose(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1) const = 0;

    // trace(K)
    virtual double trace() const = 0;

    // Diagonal of K
    virtual void getDiagonal(Eigen::VectorXd& diag) const = 0;

    // True when the operator is backed by a materialized dense matrix
    virtual bool isDense() const = 0;
};

// PLINK BED backed matrix-free operator.
//
// Both the sample set and the variant set are configurable, which is required
// for LOCO (exclude one chromosome) and for honoring sample/variant QC.  When
// a selection list is empty the full reader dimension is used.
class PlinkKinshipOperator : public KinshipOperator {
public:
    struct BlockProfile {
        std::uint64_t blocks = 0;
        std::uint64_t packed_bytes = 0;
        std::uint64_t decoded_genotypes = 0;
        double packed_read_seconds = 0.0;
        double decode_seconds = 0.0;
        double transpose_product_seconds = 0.0;
        double forward_product_seconds = 0.0;
    };
    PlinkKinshipOperator(PlinkReader* reader,
                         std::vector<int> keep_indices = {},
                         std::vector<int> variant_indices = {});

    int rows() const override;
    int cols() const override;
    int numVariants() const override;
    void apply(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1) const override;
    void apply(const Eigen::MatrixXd& X, Eigen::MatrixXd& out, int n_threads = 1) const override;
    // One decode per selected SNP, reused across every RHS and both products.
    // Profile is optional and local to the call (no output or global state).
    void apply_block(const Eigen::MatrixXd& X, Eigen::MatrixXd& out,
                     int n_threads = 1, BlockProfile* profile = nullptr) const;
    void applyTranspose(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1) const override;
    double trace() const override;
    void getDiagonal(Eigen::VectorXd& diag) const override;
    bool isDense() const override { return false; }

private:
    PlinkReader* reader_;
    std::vector<int> keep_;        // sample indices; empty means "all samples"
    std::vector<int> variants_;    // SNP indices; empty means "all SNPs"
    bool identity_samples_ = true;
    int num_samples_ = 0;
    int num_variants_ = 0;
    double scaling_ = 1.0;         // 1 / M

    // y = Z' * v over the selected variants (size num_variants_)
    void multiplyZt(const Eigen::VectorXd& v, Eigen::VectorXd& u, int n_threads) const;
    // y = Z * v over the selected samples (size rows())
    void multiplyZ(const Eigen::VectorXd& u, Eigen::VectorXd& y, int n_threads) const;
};

// Factory for the matrix-free PLINK operator.  The other readers do not yet
// provide the streaming Z-multiplies required by the matrix-free path.
std::unique_ptr<KinshipOperator> createKinshipOperator(PlinkReader* reader,
                                                       const std::vector<int>& keep_indices = {},
                                                       const std::vector<int>& variant_indices = {});

} // namespace cosmic
