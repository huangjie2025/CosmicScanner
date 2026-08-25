#include "design.h"
#include <set>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <cctype>
#include <stdexcept>

namespace cosmic {


FixedDesignG buildFixedDesignGeneric(const std::vector<GenRecord>& recs,
                                     bool legacy_names,
                                     const std::vector<std::string>& cat_names,
                                     const std::vector<std::string>& num_names) {
    int n = (int)recs.size();
    std::vector<std::vector<std::string>> cat_levels;
    if (!recs.empty()) {
        int cnum = (int)recs.front().cats.size();
        cat_levels.resize(cnum);
        for (int j = 0; j < cnum; ++j) {
            std::set<std::string> lv;
            for (const auto& r : recs) {
                if (!r.cats[j].empty()) lv.insert(r.cats[j]);
            }
            // Try to sort numerically if all look like integers
            bool all_digits = false; // Disable numerical sorting to match HIBLUP's string sort
            /*
            bool all_digits = true;
            for (const auto& s : lv) {
                if (s.empty()) continue;
                for (char c : s) if (!isdigit(c)) { all_digits = false; break; }
                if (!all_digits) break;
            }
            */
            if (all_digits && !lv.empty()) {
                std::vector<int> ints;
                for (const auto& s : lv) ints.push_back(std::stoi(s));
                std::sort(ints.begin(), ints.end());
                for (int v : ints) cat_levels[j].push_back(std::to_string(v));
            } else {
                for (auto& s : lv) cat_levels[j].push_back(s);
            }
        }
    }
    int p = 1;
    std::vector<std::string> names; names.push_back("mu");
    std::map<std::string, std::vector<int>> factor_cols;
    std::vector<std::string> factor_names;
    
    factor_cols["mu"].push_back(0);
    factor_names.push_back("mu");
    
    for (size_t j = 0; j < cat_levels.size(); ++j) {
        auto& lv = cat_levels[j];
        if (lv.empty()) continue;
        std::string base = lv.front();
        std::string prefix;
        if (legacy_names && cat_levels.size() == 4) {
            static const std::vector<std::string> legacy = {"parity","herd","ys","sex"};
            prefix = legacy[j];
        } else if (!cat_names.empty() && j < cat_names.size()) {
            prefix = cat_names[j];
        } else {
            prefix = std::string("factor") + std::to_string(j+1);
        }
        factor_names.push_back(prefix);

        for (size_t k = 0; k < lv.size(); ++k) {
            if (lv[k] == base) continue;
            names.push_back(prefix + std::string("_") + lv[k]);
            factor_cols[prefix].push_back(p);
            ++p;
        }
    }
    if (!recs.empty()) {
        int qn = (int)recs.front().nums.size();
        for (int j = 0; j < qn; ++j) {
            std::string name;
            if (legacy_names && qn == 1) name = "bw";
            else if (!num_names.empty() && j < num_names.size()) name = num_names[j];
            else name = std::string("cov") + std::to_string(j+1);
            names.push_back(name);
            factor_cols[name].push_back(p);
            factor_names.push_back(name);
            ++p;
        }
    }
    
    FixedDesignG fd; fd.p = p; fd.names = names; fd.rows.resize(n); fd.factor_cols = factor_cols; fd.factor_names = factor_names;
    for (int i = 0; i < n; ++i) {
        std::vector<std::pair<int,double>> nz; nz.emplace_back(0, 1.0);
        int col = 1;
        for (size_t j = 0; j < cat_levels.size(); ++j) {
            auto& lv = cat_levels[j]; if (lv.empty()) continue; std::string base = lv.front(); std::string val = recs[i].cats[j];
            for (size_t k = 0; k < lv.size(); ++k) {
                if (lv[k] == base) continue;
                if (val == lv[k]) nz.emplace_back(col, 1.0);
                ++col;
            }
        }
        if (!recs[i].nums.empty()) {
            for (size_t j = 0; j < recs[i].nums.size(); ++j) { nz.emplace_back(col, recs[i].nums[j]); ++col; }
        }
        fd.rows[i] = std::move(nz);
    }
    return fd;
}

// Helper for factorial
static double factorial(int n) {
    if (n <= 1) return 1.0;
    double res = 1.0;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}

Eigen::MatrixXd buildLegendreMatrix(const Eigen::VectorXd& time, double tmin, double tmax, int order) {
    int n = time.size();
    Eigen::MatrixXd pmat(n, order + 1);
    
    // Scale time to [-1, 1]
    Eigen::VectorXd tvec = 2.0 * (time.array() - tmin) / (tmax - tmin) - 1.0;
    
    for (int k = 0; k <= order; ++k) {
        int c = k / 2;
        int j = k;
        Eigen::VectorXd p = Eigen::VectorXd::Zero(n);
        
        for (int r = 0; r <= c; ++r) {
            double coeff = std::sqrt((2.0 * j + 1.0) / 2.0) * std::pow(0.5, j) * 
                           (std::pow(-1.0, r) * factorial(2 * j - 2 * r) / 
                           (factorial(r) * factorial(j - r) * factorial(j - 2 * r)));
            
            p.array() += coeff * tvec.array().pow(j - 2 * r);
        }
        pmat.col(k) = p;
    }
    return pmat;
}

Eigen::MatrixXd buildBSplineMatrix(const Eigen::VectorXd& time, double tmin, double tmax, int order) {
    int n = time.size();
    int n_basis = order + 1;
    if (n_basis < 2) n_basis = 2;
    int degree = std::min(3, n_basis - 1);
    int n_knots = n_basis + degree + 1;
    std::vector<double> knots(n_knots, 0.0);

    for (int i = 0; i <= degree; ++i) knots[i] = 0.0;
    for (int i = n_basis; i < n_knots; ++i) knots[i] = 1.0;
    int interior = n_basis - degree - 1;
    for (int j = 1; j <= interior; ++j) knots[degree + j] = (double)j / (interior + 1.0);

    double denom_time = tmax - tmin;
    if (!(denom_time > 0.0) || !std::isfinite(denom_time)) denom_time = 1.0;

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n, n_basis);
    for (int r = 0; r < n; ++r) {
        double x = (time(r) - tmin) / denom_time;
        if (!std::isfinite(x)) x = 0.0;
        x = std::max(0.0, std::min(1.0, x));

        std::vector<double> prev(n_basis + degree, 0.0), curr(n_basis + degree, 0.0);
        for (int i = 0; i < n_basis + degree; ++i) {
            bool in_interval = (knots[i] <= x && x < knots[i + 1]);
            bool at_right_edge = (x == 1.0 && knots[i] <= x && x <= knots[i + 1] && knots[i + 1] == 1.0);
            prev[i] = (in_interval || at_right_edge) ? 1.0 : 0.0;
        }
        for (int d = 1; d <= degree; ++d) {
            std::fill(curr.begin(), curr.end(), 0.0);
            for (int i = 0; i < n_basis + degree - d; ++i) {
                double left = 0.0, right = 0.0;
                double denom1 = knots[i + d] - knots[i];
                if (denom1 > 0.0) left = (x - knots[i]) / denom1 * prev[i];
                double denom2 = knots[i + d + 1] - knots[i + 1];
                if (denom2 > 0.0) right = (knots[i + d + 1] - x) / denom2 * prev[i + 1];
                curr[i] = left + right;
            }
            prev.swap(curr);
        }
        double row_sum = 0.0;
        for (int i = 0; i < n_basis; ++i) row_sum += prev[i];
        if (row_sum > 0.0) {
            for (int i = 0; i < n_basis; ++i) B(r, i) = prev[i] / row_sum;
        } else {
            B(r, std::min(n_basis - 1, std::max(0, (int)std::floor(x * n_basis)))) = 1.0;
        }
    }
    return B;
}

Eigen::MatrixXd buildRrmBasisMatrix(const Eigen::VectorXd& time, double tmin, double tmax,
                                   int order, const std::string& basis) {
    if (basis == "legendre" || basis.empty()) return buildLegendreMatrix(time, tmin, tmax, order);
    if (basis == "bspline" || basis == "b-spline" || basis == "spline") return buildBSplineMatrix(time, tmin, tmax, order);
    throw std::runtime_error("Unknown --rrm-basis: " + basis + " (expected legendre or bspline)");
}

}
