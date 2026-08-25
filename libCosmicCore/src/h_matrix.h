#pragma once
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace cosmic {

using Matrix = Eigen::MatrixXd;

struct NamedMatrix {
  Matrix data;
  std::vector<std::string> rows;
  std::vector<std::string> cols;
};

NamedMatrix Hmatrix(const NamedMatrix& A_in,
                    const NamedMatrix& G_in,
                    const Matrix* markers,
                    double c,
                    const std::string& method,
                    double tau,
                    double omega,
                    double missingValue,
                    double maf,
                    int ploidy,
                    int roundVar,
                    bool ASVflag);

NamedMatrix HInverseMatrix(const NamedMatrix& A22,
                           const Eigen::SparseMatrix<double>& Ainv_sparse,
                           const std::vector<std::string>& a_ids,
                           const NamedMatrix& Ginv_in,
                           double tau,
                           double omega);

} // namespace cosmic
