#include "input/input_source_reader.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace cosmic {

char detectDelimiter(const std::string& path, const std::string& format, const std::string& first_line) {
    if (format == "csv") return ',';
    if (format == "tsv") return '\t';
    if (format == "txt") return '\t';
    size_t tabs = std::count(first_line.begin(), first_line.end(), '\t');
    size_t commas = std::count(first_line.begin(), first_line.end(), ',');
    if (commas > tabs) return ',';
    if (tabs > 0) return '\t';
    if (path.size() >= 4 && path.substr(path.size() - 4) == ".csv") return ',';
    return '\t';
}

std::vector<std::string> parseDelimitedLine(const std::string& line, char delim) {
    std::vector<std::string> fields;
    std::string cur;
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (c == '"') {
            if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (c == delim && !in_quotes) {
            fields.push_back(trim(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    fields.push_back(trim(cur));
    return fields;
}

Table readTable(const std::string& path, const std::string& format) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("INPUT_FILE_NOT_READABLE: Cannot open --input-file: " + path);
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("INPUT_FILE_EMPTY: " + path);
    if (!line.empty() && line.back() == '\r') line.pop_back();
    char delim = detectDelimiter(path, format, line);
    Table table;
    table.header = parseDelimitedLine(line, delim);
    for (size_t i = 0; i < table.header.size(); ++i) {
        std::string name = trim(table.header[i]);
        if (name.empty()) throw std::runtime_error("MISSING_REQUIRED_COLUMN: Input header contains an empty column name");
        if (table.col_index.count(name)) throw std::runtime_error("DUPLICATE_INPUT_COLUMN: " + name);
        table.header[i] = name;
        table.col_index[name] = static_cast<int>(i);
    }
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (trim(line).empty()) continue;
        auto row = parseDelimitedLine(line, delim);
        row.resize(table.header.size());
        table.rows.push_back(row);
    }
    return table;
}

void requireColumn(const Table& table, const std::string& col, const std::string& arg) {
    if (col.empty()) throw std::runtime_error("MISSING_REQUIRED_COLUMN: " + arg + " is required for --input-file MVP");
    if (!table.col_index.count(col)) throw std::runtime_error("MISSING_REQUIRED_COLUMN: " + arg + " column not found: " + col);
}

std::string getCell(const Table& table, const std::vector<std::string>& row, const std::string& col) {
    auto it = table.col_index.find(col);
    if (it == table.col_index.end()) return "";
    return it->second < (int)row.size() ? row[it->second] : "";
}

} // namespace cosmic
