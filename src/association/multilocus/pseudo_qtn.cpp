#include "association/multilocus/multilocus_internal.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cosmic::association::detail {
namespace {

double betaFraction(double a, double b, double x) {
    constexpr int max_iterations = 200;
    constexpr double epsilon = 3e-14;
    constexpr double floor = 1e-300;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::abs(d) < floor) d = floor;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= max_iterations; ++m) {
        const int m2 = 2 * m;
        double aa = static_cast<double>(m) * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < floor) d = floor;
        c = 1.0 + aa / c;
        if (std::abs(c) < floor) c = floor;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < floor) d = floor;
        c = 1.0 + aa / c;
        if (std::abs(c) < floor) c = floor;
        d = 1.0 / d;
        const double delta = d * c;
        h *= delta;
        if (std::abs(delta - 1.0) <= epsilon) break;
    }
    return h;
}

double fSurvival(double f, double denominator_df) {
    if (!(f > 0.0) || !std::isfinite(f)) return 1.0;
    const double x = denominator_df / (denominator_df + f);
    const double a = denominator_df / 2.0, b = 0.5;
    const double log_beta = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                            a * std::log(x) + b * std::log1p(-x);
    const double value = x < (a + 1.0) / (a + b + 2.0)
        ? std::exp(log_beta) * betaFraction(a, b, x) / a
        : 1.0 - std::exp(log_beta) * betaFraction(b, a, 1.0 - x) / b;
    return std::clamp(value, 0.0, 1.0);
}

} // namespace

MarkerFit fitMarker(const Eigen::VectorXd& y, const Eigen::MatrixXd& base,
                    const Eigen::MatrixXd& genotypes,
                    const std::vector<int>& conditioning, int marker) {
    if (marker < 0 || marker >= genotypes.cols() || y.size() != genotypes.rows() ||
        y.size() != base.rows()) throw std::invalid_argument("Marker fit dimensions do not match");
    int count = 0;
    for (int index : conditioning) if (index != marker) ++count;
    Eigen::MatrixXd design(base.rows(), base.cols() + count);
    design.leftCols(base.cols()) = base;
    int column = base.cols();
    for (int index : conditioning) {
        if (index == marker) continue;
        design.col(column++) = genotypes.col(index);
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);
    if (qr.rank() != design.cols()) return {};
    const int df = static_cast<int>(y.size()) - qr.rank() - 1;
    if (df <= 0) return {};
    const Eigen::VectorXd y_residual = y - design * qr.solve(y);
    const Eigen::VectorXd genotype = genotypes.col(marker);
    const Eigen::VectorXd g_residual = genotype - design * qr.solve(genotype);
    const double denominator = g_residual.squaredNorm();
    if (!(denominator > 1e-12) || !std::isfinite(denominator)) return {};
    MarkerFit result;
    result.beta = g_residual.dot(y_residual) / denominator;
    const double rss = std::max(y_residual.squaredNorm() - result.beta * result.beta * denominator,
                                std::numeric_limits<double>::min());
    result.se = std::sqrt(rss / (static_cast<double>(df) * denominator));
    const double f = result.beta * result.beta / (result.se * result.se);
    result.p = fSurvival(f, df);
    result.tested = std::isfinite(result.beta) && std::isfinite(result.se) && result.se > 0.0 &&
                    std::isfinite(result.p) && result.p >= 0.0 && result.p <= 1.0;
    return result;
}

} // namespace cosmic::association::detail
