#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <functional>
#include <string>
#include <vector>

namespace cosmic {

struct PEVSolveInfo {
    bool converged = false;
    int iterations = 0;
    double residual_norm = 0.0;
    double relative_residual = 0.0;
    std::string termination;
};

// Applies the unscaled inverse of a lambda-form coefficient matrix:
// apply(x) = C_lambda^-1 x. Callers retain responsibility for sigma_e^2 scaling.
class PEVOperator {
public:
    virtual ~PEVOperator() = default;

    virtual int rows() const = 0;
    virtual std::string backendName() const = 0;
    virtual Eigen::VectorXd apply(const Eigen::VectorXd& x) const = 0;
    virtual Eigen::MatrixXd apply(const Eigen::MatrixXd& x) const;

    double quadraticForm(const Eigen::VectorXd& x) const;
};

// Materialized-MME reference backend. The name identifies the reference PEV
// route; its coefficient matrix and factorization remain sparse.
class DensePEVOperator final : public PEVOperator {
public:
    explicit DensePEVOperator(const Eigen::SparseMatrix<double>& coefficient);
    explicit DensePEVOperator(const Eigen::MatrixXd& coefficient);

    int rows() const override { return static_cast<int>(coefficient_.rows()); }
    std::string backendName() const override { return "dense-reference/sparse-ldlt"; }
    Eigen::VectorXd apply(const Eigen::VectorXd& x) const override;
    Eigen::MatrixXd apply(const Eigen::MatrixXd& x) const override;

private:
    Eigen::SparseMatrix<double> coefficient_;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factor_;

    void factorize();
};

class MatrixFreePEVOperator final : public PEVOperator {
public:
    using ApplyFunction =
        std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>;
    using PreconditionerFunction =
        std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>;

    struct Options {
        double tolerance = 1e-10;
        int max_iterations = 2000;
    };

    MatrixFreePEVOperator(int dimension, ApplyFunction apply,
                          PreconditionerFunction preconditioner = {});
    MatrixFreePEVOperator(int dimension, ApplyFunction apply,
                          PreconditionerFunction preconditioner,
                          const Options& options);

    int rows() const override { return dimension_; }
    std::string backendName() const override { return "matrix-free-pcg"; }
    Eigen::VectorXd apply(const Eigen::VectorXd& x) const override;
    Eigen::MatrixXd apply(const Eigen::MatrixXd& x) const override;

    const PEVSolveInfo& lastSolveInfo() const { return last_solve_info_; }
    const std::vector<PEVSolveInfo>& lastBatchSolveInfo() const {
        return last_batch_solve_info_;
    }

private:
    int dimension_;
    ApplyFunction apply_;
    PreconditionerFunction preconditioner_;
    Options options_;
    mutable PEVSolveInfo last_solve_info_;
    mutable std::vector<PEVSolveInfo> last_batch_solve_info_;

    Eigen::VectorXd solve(const Eigen::VectorXd& rhs, PEVSolveInfo& info) const;
};

} // namespace cosmic
