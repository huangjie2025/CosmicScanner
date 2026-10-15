#include "repeatability_helpers.h"

#include <algorithm>

namespace cosmic {

std::vector<std::string> sortChroms(const std::set<std::string>& chroms_set) {
    std::vector<std::string> chroms(chroms_set.begin(), chroms_set.end());
    std::sort(chroms.begin(), chroms.end(), [](const std::string& a, const std::string& b) {
        try { return std::stoi(a) < std::stoi(b); } catch (...) { return a < b; }
    });
    return chroms;
}

bool uses_per_snp_scaling(int alg) {
    return alg == 0 || alg == 3 || alg == 5;
}

} // namespace cosmic
