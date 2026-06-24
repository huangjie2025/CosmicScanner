#include "input/time_normalizer.h"

#include "input/input_common.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace cosmic {
namespace {

int daysFromCivil(int y, unsigned m, unsigned d) {
    y -= m <= 2;
    const int era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = static_cast<unsigned>(y - era * 400);
    const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + static_cast<int>(doe) - 719468;
}

} // namespace

std::vector<double> normalizeBaseTimes(const std::vector<std::string>& raw_times,
                                      const std::string& parser) {
    std::vector<double> bases(raw_times.size(), NAN);
    bool all_numeric = true;
    for (const auto& t : raw_times) {
        double v = 0.0;
        if (!parseDouble(t, v)) { all_numeric = false; break; }
    }
    bool use_date = parser == "date" || (parser == "auto" && !all_numeric);
    if (parser == "integer" || (parser == "auto" && all_numeric)) {
        for (size_t i = 0; i < raw_times.size(); ++i) {
            double v = 0.0;
            if (!parseDouble(raw_times[i], v)) throw std::runtime_error("TIME_NOT_PARSEABLE: " + raw_times[i]);
            bases[i] = v;
        }
        return bases;
    }
    if (use_date) {
        std::vector<int> days(raw_times.size());
        int min_day = 0;
        bool min_set = false;
        for (size_t i = 0; i < raw_times.size(); ++i) {
            int y = 0, m = 0, d = 0;
            if (!parseDate(raw_times[i], y, m, d)) throw std::runtime_error("TIME_NOT_PARSEABLE: " + raw_times[i]);
            days[i] = daysFromCivil(y, static_cast<unsigned>(m), static_cast<unsigned>(d));
            if (!min_set || days[i] < min_day) { min_day = days[i]; min_set = true; }
        }
        for (size_t i = 0; i < raw_times.size(); ++i) bases[i] = days[i] - min_day;
        return bases;
    }
    std::map<std::string, int> batch_to_index;
    for (size_t i = 0; i < raw_times.size(); ++i) {
        auto it = batch_to_index.find(raw_times[i]);
        if (it == batch_to_index.end()) {
            int idx = static_cast<int>(batch_to_index.size());
            it = batch_to_index.emplace(raw_times[i], idx).first;
        }
        bases[i] = it->second;
    }
    return bases;
}

int resolveTimeMultiplier(const std::string& value, int max_duplicate_count) {
    if (value.empty() || value == "auto") {
        int multiplier = 1;
        while (multiplier <= std::max(1, max_duplicate_count)) multiplier *= 10;
        return multiplier;
    }
    try {
        int multiplier = std::stoi(value);
        if (multiplier <= max_duplicate_count) {
            throw std::runtime_error("INVALID_TIME_MULTIPLIER: --time-index-multiplier must be greater than max duplicate records per ID-time");
        }
        return multiplier;
    } catch (const std::invalid_argument&) {
        throw std::runtime_error("INVALID_TIME_MULTIPLIER: " + value);
    }
}

} // namespace cosmic
