#include "input/wide_to_long.h"

#include <cctype>
#include <stdexcept>

namespace cosmic {
namespace {

bool allDigits(const std::string& s) {
    if (s.empty()) return false;
    for (char c : s) if (!std::isdigit(static_cast<unsigned char>(c))) return false;
    return true;
}

bool parseTimeToken(const std::string& token, std::string& out_time) {
    if (token.empty()) return false;
    if (allDigits(token)) { out_time = token; return true; }
    if ((token[0] == 'd' || token[0] == 'D' || token[0] == 't' || token[0] == 'T') && allDigits(token.substr(1))) {
        out_time = token.substr(1);
        return true;
    }
    const std::string day = "day";
    if (token.size() > day.size()) {
        std::string prefix = token.substr(0, day.size());
        for (char& c : prefix) c = std::tolower(static_cast<unsigned char>(c));
        if (prefix == day && allDigits(token.substr(day.size()))) {
            out_time = token.substr(day.size());
            return true;
        }
    }
    if (token.size() == 8 && allDigits(token)) { out_time = token; return true; }
    return false;
}

} // namespace

WideSpec parseWideSpec(const std::string& column) {
    size_t pos = column.find_last_of('_');
    if (pos == std::string::npos || pos == 0 || pos + 1 >= column.size()) {
        throw std::runtime_error("WIDE_COLUMN_PARSE_FAILED: Cannot parse wide trait-time column: " + column + " (expected trait_1, trait_d1, day1_trait, or trait_day1)");
    }
    std::string left = column.substr(0, pos);
    std::string right = column.substr(pos + 1);
    std::string trait;
    std::string raw_time;
    if (parseTimeToken(right, raw_time)) {
        trait = left;
    } else if (parseTimeToken(left, raw_time)) {
        trait = right;
    } else {
        throw std::runtime_error("WIDE_COLUMN_PARSE_FAILED: Cannot parse wide time token from column: " + column);
    }
    if (trait.empty() || raw_time.empty()) throw std::runtime_error("WIDE_COLUMN_PARSE_FAILED: Cannot parse wide column: " + column);
    return {column, trait, raw_time};
}

} // namespace cosmic
