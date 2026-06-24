#include "input/input_common.h"

#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace cosmic {

std::string trim(const std::string& s) {
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

bool isMissing(const std::string& v) {
    std::string s = trim(v);
    std::string l;
    l.reserve(s.size());
    for (char c : s) l.push_back(std::tolower(static_cast<unsigned char>(c)));
    return s.empty() || l == "na" || l == "nan" || l == "null" || l == "." || l == "-9";
}

bool parseDouble(const std::string& s, double& out) {
    try {
        size_t pos = 0;
        std::string t = trim(s);
        out = std::stod(t, &pos);
        return pos == t.size() && std::isfinite(out);
    } catch (...) {
        return false;
    }
}

bool parseDate(const std::string& s, int& y, int& m, int& d) {
    std::string t = trim(s);
    if (t.size() < 8) return false;
    for (char& c : t) if (c == '/' || c == '.') c = '-';
    char sep1 = 0, sep2 = 0;
    std::stringstream ss(t);
    if ((ss >> y >> sep1 >> m >> sep2 >> d) && sep1 == '-' && sep2 == '-') return true;
    return false;
}

std::string numberToString(double x) {
    std::ostringstream os;
    os << std::setprecision(15) << x;
    return os.str();
}

std::string key2(const std::string& fid, const std::string& iid) {
    return fid + ":" + iid;
}

std::string key3(const std::string& fid, const std::string& iid, const std::string& time) {
    return fid + ":" + iid + ":" + time;
}

} // namespace cosmic
