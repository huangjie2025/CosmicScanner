#pragma once

#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "apy_ginv.h"
#include "block_streamed_grm.h"
#include "gcc_solver.h"

namespace cosmic {

class CrossGenomicOperator;

class ApyInverseOperator {
public:
    virtual ~ApyInverseOperator() = default;

    virtual Eigen::VectorXd apply(const Eigen::VectorXd& x) const = 0;
    virtual Eigen::MatrixXd apply(const Eigen::MatrixXd& x) const = 0;
    virtual Eigen::VectorXd applyCoreInverse(const Eigen::VectorXd& x) const = 0;
    virtual Eigen::MatrixXd applyCoreInverse(const Eigen::MatrixXd& x) const = 0;
    virtual const Eigen::VectorXd& diagonal() const = 0;
    virtual int rows() const = 0;
    virtual std::string backendName() const = 0;
};

class DenseApyInverseOperator final : public ApyInverseOperator {
public:
    explicit DenseApyInverseOperator(const ApyGinv& inverse) : inverse_(inverse) {}

    Eigen::VectorXd apply(const Eigen::VectorXd& x) const override;
    Eigen::MatrixXd apply(const Eigen::MatrixXd& x) const override;
    Eigen::VectorXd applyCoreInverse(const Eigen::VectorXd& x) const override;
    Eigen::MatrixXd applyCoreInverse(const Eigen::MatrixXd& x) const override;
    const Eigen::VectorXd& diagonal() const override;
    int rows() const override { return inverse_.n(); }
    std::string backendName() const override { return "dense-apy"; }

private:
    const ApyGinv& inverse_;
    mutable Eigen::VectorXd diagonal_;
};

GcnBlockIteratorFactory makeDenseGcnBlockIteratorFactory(
    const Eigen::MatrixXd& g,
    const std::vector<int>& core_indices,
    const std::vector<int>& noncore_indices,
    int block_size = 4096);

class MatrixFreeApyInverseOperator final : public ApyInverseOperator {
public:
    MatrixFreeApyInverseOperator(
        GccSolverPtr gcc_solver,
        GcnBlockIteratorFactory gcn_factory,
        const Eigen::VectorXd& gnn_diag,
        int ng,
        const std::vector<int>& core_indices,
        const std::string& tag,
        const CrossGenomicOperator* cross_operator = nullptr);

    Eigen::VectorXd apply(const Eigen::VectorXd& x) const override;
    Eigen::MatrixXd apply(const Eigen::MatrixXd& x) const override;
    Eigen::VectorXd applyCoreInverse(const Eigen::VectorXd& x) const override;
    Eigen::MatrixXd applyCoreInverse(const Eigen::MatrixXd& x) const override;
    const Eigen::VectorXd& diagonal() const override { return diagonal_; }
    int rows() const override { return n_core_ + n_noncore_; }
    std::string backendName() const override;

    const Eigen::VectorXd& mnn() const { return mnn_; }
    const std::vector<int>& coreIndices() const { return core_indices_; }
    const std::vector<int>& noncoreIndices() const { return noncore_indices_; }
    GccSolverPtr gccSolverPtr() const { return gcc_solver_; }
    size_t persistentBytes() const;

private:
    GccSolverPtr gcc_solver_;
    GcnBlockIteratorFactory gcn_factory_;
    std::vector<int> core_indices_;
    std::vector<int> noncore_indices_;
    Eigen::VectorXd mnn_;
    Eigen::VectorXd mnn_inv_;
    Eigen::VectorXd diagonal_;
    int n_core_ = 0;
    int n_noncore_ = 0;
    std::string tag_;
    const CrossGenomicOperator* cross_operator_ = nullptr;

    void validateAndBuildMnn(const Eigen::VectorXd& gnn_diag);
};

} // namespace cosmic
