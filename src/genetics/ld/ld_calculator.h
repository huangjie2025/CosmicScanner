#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace cosmic {
namespace genetics {

// All dosage matrices have samples in rows (>= 2), variants in columns.
// NaN and finite negative dosages are missing and are column-mean imputed.
// Other dosages must be finite in [0,2]; constant/all-missing columns throw.
class LDCalculator {
public:
    // Dense signed Pearson r, not r-squared. Zero variants yields a 0x0 matrix.
    static Eigen::MatrixXd computePairwiseLD(Eigen::MatrixXd dosages);
};

// Metadata must match columns and be lexicographically sorted by (chrom,pos),
// with nonempty chromosome labels and nonnegative positions; ties are allowed.
// window_bp >= 0, step >= 1, finite r2 in [0,1]. Comparisons use r*r >= r2.
class LDPruner {
public:
    // Returns retained original column indices in input order. Each chromosome
    // starts windows at variant offsets 0,step,2*step,..., including positions
    // through start+window_bp. Within a window, earlier retained variants win.
    static std::vector<int> prune(Eigen::MatrixXd dosages,
                                 const std::vector<std::string>& chrom,
                                 const std::vector<long long>& pos,
                                 int window_bp, int step, double r2);
};

class LDClumper {
public:
    // Leads are unassigned variants with p <= threshold, ordered by (p,index).
    // Each lead claims unassigned same-chromosome variants within +/-window_bp
    // meeting r2, regardless of their p. Leads map to themselves; others not
    // claimed map to -1. p and threshold must be finite in [0,1]. step is the
    // neighbor scan batch size; it never skips variants or changes membership.
    static std::vector<int> clump(Eigen::MatrixXd dosages,
                                 const std::vector<std::string>& chrom,
                                 const std::vector<long long>& pos,
                                 int window_bp, int step, double r2,
                                 const Eigen::VectorXd& p, double threshold);
};

} // namespace genetics
} // namespace cosmic
