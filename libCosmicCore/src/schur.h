#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include "design.h"
#include "matrix_adapter.h"

namespace cosmic {

Eigen::MatrixXd computeSchurFullSinv(const Eigen::MatrixXd& Xfull,
                                     const std::vector<GenRecord>& recs,
                                     const AbstractMatrix* Qinv,
                                     double sigma2_e,
                                     double sigma2_u);

}