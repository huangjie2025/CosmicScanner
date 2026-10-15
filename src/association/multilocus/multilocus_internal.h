#pragma once

#include <Eigen/Dense>
#include <vector>

namespace cosmic::association::detail {

struct MarkerFit {
    double beta = 0.0;
    double se = 0.0;
    double p = 1.0;
    bool tested = false;
};

MarkerFit fitMarker(const Eigen::VectorXd& y, const Eigen::MatrixXd& base,
                    const Eigen::MatrixXd& genotypes,
                    const std::vector<int>& conditioning, int marker);

std::vector<int> ldPruneRanked(const Eigen::MatrixXd& genotypes,
                               const std::vector<int>& ranked,
                               int maximum, double r2);

std::vector<int> selectByBic(const Eigen::VectorXd& y,
                             const Eigen::MatrixXd& base,
                             const Eigen::MatrixXd& genotypes,
                             const std::vector<int>& candidates,
                             int maximum);

} // namespace cosmic::association::detail
