#pragma once

#include <Eigen/Dense>

#include <memory>
#include <string>

#include "apy_ginv.h"
#include "block_streamed_grm.h"
#include "gcc_solver.h"
#include "pev_operator.h"

namespace cosmic {

class CrossGenomicOperator;

// Inverse application for the APY Eq. 8 core precision
//   Kc = Dc + lambda Gcc^-1 + Q diag(s) Q'.
class APYCorePEVOperator : public PEVOperator {
public:
    virtual void applyCoefficient(const Eigen::VectorXd& x,
                                  Eigen::VectorXd& y) const = 0;
};

class DenseAPYCorePEVOperator final : public APYCorePEVOperator {
public:
    DenseAPYCorePEVOperator(const ApyGinv& apy,
                            const Eigen::VectorXd& dc,
                            const Eigen::VectorXd& dn,
                            double lambda);

    int rows() const override { return static_cast<int>(coefficient_.rows()); }
    std::string backendName() const override { return "dense-apy-core-pev"; }
    Eigen::VectorXd apply(const Eigen::VectorXd& x) const override;
    Eigen::MatrixXd apply(const Eigen::MatrixXd& x) const override;
    void applyCoefficient(const Eigen::VectorXd& x,
                          Eigen::VectorXd& y) const override;

    const Eigen::MatrixXd& coefficient() const { return coefficient_; }
    const Eigen::VectorXd& schurWeights() const { return schur_weights_; }

private:
    Eigen::MatrixXd coefficient_;
    Eigen::VectorXd schur_weights_;
    std::unique_ptr<DensePEVOperator> inverse_;
};

class MatrixFreeAPYCorePEVOperator final : public APYCorePEVOperator {
public:
    MatrixFreeAPYCorePEVOperator(
        GccSolverPtr gcc_solver,
        GcnBlockIteratorFactory gcn_factory,
        const Eigen::VectorXd& mnn,
        const Eigen::VectorXd& dc,
        const Eigen::VectorXd& dn,
        double lambda,
        const MatrixFreePEVOperator::Options& options = {},
        const CrossGenomicOperator* cross_operator = nullptr);

    MatrixFreeAPYCorePEVOperator(const MatrixFreeAPYCorePEVOperator&) = delete;
    MatrixFreeAPYCorePEVOperator& operator=(const MatrixFreeAPYCorePEVOperator&) = delete;
    MatrixFreeAPYCorePEVOperator(MatrixFreeAPYCorePEVOperator&&) = delete;
    MatrixFreeAPYCorePEVOperator& operator=(MatrixFreeAPYCorePEVOperator&&) = delete;

    int rows() const override { return static_cast<int>(dc_.size()); }
    std::string backendName() const override { return "matrix-free-apy-core-pev"; }
    Eigen::VectorXd apply(const Eigen::VectorXd& x) const override;
    Eigen::MatrixXd apply(const Eigen::MatrixXd& x) const override;
    void applyCoefficient(const Eigen::VectorXd& x,
                          Eigen::VectorXd& y) const override;

    const Eigen::VectorXd& schurWeights() const { return schur_weights_; }
    const Eigen::VectorXd& coefficientDiagonal() const { return diagonal_; }
    const PEVSolveInfo& lastSolveInfo() const;
    const std::vector<PEVSolveInfo>& lastBatchSolveInfo() const;

private:
    GccSolverPtr gcc_solver_;
    GcnBlockIteratorFactory gcn_factory_;
    Eigen::VectorXd dc_;
    Eigen::VectorXd schur_weights_;
    Eigen::VectorXd diagonal_;
    double lambda_ = 0.0;
    int n_noncore_ = 0;
    const CrossGenomicOperator* cross_operator_ = nullptr;
    std::unique_ptr<MatrixFreePEVOperator> inverse_;

    void validateStream() const;
};

} // namespace cosmic
