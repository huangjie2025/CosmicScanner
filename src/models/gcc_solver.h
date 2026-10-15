#ifndef COSMIC_GCC_SOLVER_H
#define COSMIC_GCC_SOLVER_H

#include <Eigen/Dense>
#include <memory>
#include <string>

namespace cosmic {

class GccSolver {
public:
    virtual ~GccSolver() = default;

    virtual Eigen::VectorXd solve(const Eigen::VectorXd& rhs) const = 0;
    virtual Eigen::MatrixXd solve(const Eigen::MatrixXd& rhs) const = 0;

    virtual int rows() const = 0;
    virtual int cols() const = 0;
    virtual bool isValid() const = 0;
    virtual std::string backendName() const = 0;

    virtual Eigen::VectorXd diagonal() const = 0;
    virtual Eigen::VectorXd vectorD() const = 0;

    virtual double minDiagD() const = 0;
    virtual double maxDiagD() const = 0;
};

class GccSolverDense : public GccSolver {
public:
    GccSolverDense() = default;
    explicit GccSolverDense(const Eigen::MatrixXd& Gcc);

    void compute(const Eigen::MatrixXd& Gcc);

    Eigen::VectorXd solve(const Eigen::VectorXd& rhs) const override;
    Eigen::MatrixXd solve(const Eigen::MatrixXd& rhs) const override;

    int rows() const override { return ldlt_.rows(); }
    int cols() const override { return ldlt_.cols(); }
    bool isValid() const override { return computed_; }
    std::string backendName() const override { return "DenseLDLT"; }

    Eigen::VectorXd diagonal() const override;
    Eigen::VectorXd vectorD() const override { return ldlt_.vectorD(); }
    double minDiagD() const override { return min_diag_d_; }
    double maxDiagD() const override { return max_diag_d_; }

    const Eigen::LDLT<Eigen::MatrixXd>& ldlt() const { return ldlt_; }

private:
    Eigen::LDLT<Eigen::MatrixXd> ldlt_;
    bool computed_ = false;
    double min_diag_d_ = 0.0;
    double max_diag_d_ = 0.0;
};

// Forward declarations for PCG backend
class GccSolverPCG;
using GccSolverPCGPtr = std::shared_ptr<GccSolverPCG>;

using GccSolverPtr = std::shared_ptr<GccSolver>;

// Factory functions
GccSolverPtr makeDenseSolver();
GccSolverPtr makeDenseSolver(const Eigen::MatrixXd& Gcc);
GccSolverPCGPtr makePCGSolver();

} // namespace cosmic

#endif // COSMIC_GCC_SOLVER_H