#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace cosmic {

struct Table {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    std::unordered_map<std::string, int> col_index;
};

struct PreparedRecord {
    std::string fid;
    std::string iid;
    std::string raw_id;
    std::string raw_time;
    std::string time;
    std::string time_base;
    int within_time_index = 0;
    int raw_row_index = 0;
    std::unordered_map<std::string, std::string> values;
};

struct WideSpec {
    std::string column;
    std::string trait;
    std::string raw_time;
};

std::string trim(const std::string& s);
bool isMissing(const std::string& v);
bool parseDouble(const std::string& s, double& out);
bool parseDate(const std::string& s, int& y, int& m, int& d);
std::string numberToString(double x);
std::string key2(const std::string& fid, const std::string& iid);
std::string key3(const std::string& fid, const std::string& iid, const std::string& time);

} // namespace cosmic
