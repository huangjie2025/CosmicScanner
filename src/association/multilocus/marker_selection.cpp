#include "association/multilocus/multilocus_internal.h"
#include "genetics/ld/ld_calculator.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cosmic::association::detail {

std::vector<int> ldPruneRanked(const Eigen::MatrixXd& genotypes,
                               const std::vector<int>& ranked,
                               int maximum, double r2) {
    if (maximum < 0 || !std::isfinite(r2) || r2 < 0.0 || r2 > 1.0)
        throw std::invalid_argument("Invalid multi-locus marker-selection options");
    const int candidate_count = std::min<int>(ranked.size(), std::min(2000, std::max(maximum, maximum * 4)));
    Eigen::MatrixXd candidates(genotypes.rows(), candidate_count);
    for (int j = 0; j < candidate_count; ++j) {
        const int marker = ranked[j];
        if (marker < 0 || marker >= genotypes.cols())
            throw std::invalid_argument("Invalid ranked marker index");
        candidates.col(j) = genotypes.col(marker);
    }
    const Eigen::MatrixXd ld = genetics::LDCalculator::computePairwiseLD(std::move(candidates));
    std::vector<int> selected;
    std::vector<int> selected_rank_positions;
    for (int j = 0; j < candidate_count; ++j) {
        bool linked = false;
        for (int retained : selected_rank_positions) {
            if (ld(j, retained) * ld(j, retained) >= r2) {
                linked = true;
                break;
            }
        }
        if (!linked) {
            selected.push_back(ranked[j]);
            selected_rank_positions.push_back(j);
        }
        if (static_cast<int>(selected.size()) >= maximum) break;
    }
    return selected;
}

} // namespace cosmic::association::detail
