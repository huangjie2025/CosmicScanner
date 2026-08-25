#pragma once

#include "input/input_common.h"

namespace cosmic {

char detectDelimiter(const std::string& path, const std::string& format, const std::string& first_line);
std::vector<std::string> parseDelimitedLine(const std::string& line, char delim);
Table readTable(const std::string& path, const std::string& format);
void requireColumn(const Table& table, const std::string& col, const std::string& arg);
std::string getCell(const Table& table, const std::vector<std::string>& row, const std::string& col);

} // namespace cosmic
