#pragma once
#include "bgen_ad_null.h"

namespace cosmic::detail {
BgenAdNullResult fitBgenAdReml(const BgenAdKernelPair& kernels,
                             const Eigen::VectorXd& y,
                             const Eigen::MatrixXf& covariates,
                             const EigenVceAD::Options& options);
}
