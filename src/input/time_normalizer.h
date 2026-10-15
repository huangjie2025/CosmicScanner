#pragma once

#include <string>
#include <vector>

namespace cosmic {

std::vector<double> normalizeBaseTimes(const std::vector<std::string>& raw_times,
                                      const std::string& parser);
int resolveTimeMultiplier(const std::string& value, int max_duplicate_count);

} // namespace cosmic
