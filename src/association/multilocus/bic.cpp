#include "association/multilocus/multilocus_internal.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cosmic::association::detail {

std::vector<int> selectByBic(const Eigen::VectorXd& y,
                             const Eigen::MatrixXd& base,
                             const Eigen::MatrixXd& genotypes,
                             const std::vector<int>& candidates,
                             int maximum) {
    if (y.size() != base.rows() || y.size() != genotypes.rows())
        throw std::invalid_argument("BIC dimensions do not match");
    const int residual_df_limit = std::max(0, static_cast<int>(y.size()) -
                                               static_cast<int>(base.cols()) - 2);
    const int limit = std::min({maximum, static_cast<int>(candidates.size()), residual_df_limit});
    double best = std::numeric_limits<double>::infinity();
    int best_count = 0;
    for (int count = 0; count <= limit; ++count) {
        Eigen::MatrixXd design(base.rows(), base.cols() + count);
        design.leftCols(base.cols()) = base;
        for (int j = 0; j < count; ++j) design.col(base.cols() + j) = genotypes.col(candidates[j]);
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);
        if (qr.rank() != design.cols()) continue;
        const Eigen::VectorXd residual = y - design * qr.solve(y);
        const double rss = residual.squaredNorm();
        if (!(rss > 0.0) || !std::isfinite(rss)) continue;
        const double bic = static_cast<double>(y.size()) * std::log(rss / y.size()) +
                           static_cast<double>(design.cols()) * std::log(static_cast<double>(y.size()));
        if (bic < best) {
            best = bic;
            best_count = count;
        }
    }
    return {candidates.begin(), candidates.begin() + best_count};
}

} // namespace cosmic::association::detail
