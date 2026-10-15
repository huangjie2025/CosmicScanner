#include "covar_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <iomanip>

namespace cosmic {

CovarReader::CovarReader() : num_covars(0), missing_value_code(-9.0) {}

std::string CovarReader::makeTimeKey(const std::string& fid, const std::string& iid, double time) const {
    std::ostringstream os;
    os << std::setprecision(15) << time;
    return fid + ":" + iid + ":" + os.str();
}

void CovarReader::load(const std::string& file_path, double missing_val) {
    missing_value_code = missing_val;
    covar_map.clear();
    num_covars = 0;
    
    std::ifstream ifs(file_path);
    if (!ifs) throw std::runtime_error("Cannot open covariate file: " + file_path);
    
    std::string line;
    int line_num = 0;
    
    while (std::getline(ifs, line)) {
        line_num++;
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string fid, iid;
        ss >> fid >> iid;
        
        if (fid.empty() || iid.empty()) continue;
        
        std::vector<float> vals;
        float val;
        while (ss >> val) {
            // Check missing
            if (std::abs(val - missing_value_code) < 1e-8) {
                vals.push_back(NAN);
            } else {
                vals.push_back(val);
            }
        }
        
        if (num_covars == 0) {
            num_covars = vals.size();
        } else if (vals.size() != (size_t)num_covars) {
            // Warn or skip? GCTA usually strictly checks column count.
            // Let's warn and skip for now.
            std::cerr << "Warning: Line " << line_num << " in " << file_path 
                      << " has " << vals.size() << " covariates, expected " << num_covars << ". Skipping.\n";
            continue;
        }
        
        covar_map[makeKey(fid, iid)] = vals;
    }
    
    std::cout << "Loaded " << covar_map.size() << " covariate records (p=" << num_covars << ") from " << file_path << std::endl;
}

void CovarReader::loadCategorical(const std::string& file_path, const std::string& missing_val) {
    covar_map.clear();
    num_covars = 0;
    
    std::ifstream ifs(file_path);
    if (!ifs) throw std::runtime_error("Cannot open covariate file: " + file_path);
    
    std::string line;
    int line_num = 0;
    
    // 1. First pass: Collect all unique levels for each column
    std::vector<std::map<std::string, int>> levels_per_col;
    std::vector<std::vector<std::string>> raw_data_rows;
    std::vector<std::pair<std::string, std::string>> row_ids;

    while (std::getline(ifs, line)) {
        line_num++;
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string fid, iid;
        ss >> fid >> iid;
        
        if (fid.empty() || iid.empty()) continue;
        
        std::vector<std::string> row_vals;
        std::string val;
        while (ss >> val) {
            row_vals.push_back(val);
        }
        
        if (levels_per_col.empty()) {
            levels_per_col.resize(row_vals.size());
        } else if (row_vals.size() != levels_per_col.size()) {
             std::cerr << "Warning: Line " << line_num << " in " << file_path 
                      << " has " << row_vals.size() << " covariates, expected " << levels_per_col.size() << ". Skipping.\n";
            continue;
        }

        raw_data_rows.push_back(row_vals);
        row_ids.push_back({fid, iid});

        for(size_t i=0; i<row_vals.size(); ++i) {
            if (row_vals[i] != missing_val) {
                levels_per_col[i][row_vals[i]]++; // Count occurrences (could be used for frequency check)
            }
        }
    }

    // 2. Define Encoding (Reference level = first encountered or alphabetically first?)
    // Let's use alphabetically first as reference (dropped), others get a dummy column.
    // If only 2 levels, 1 column. If k levels, k-1 columns.
    
    int total_dummy_cols = 0;
    
    // We need offsets per column
    std::vector<int> col_offsets(levels_per_col.size(), 0);
    int accum = 0;
    for(size_t i=0; i<levels_per_col.size(); ++i) {
        col_offsets[i] = accum;
        int k = levels_per_col[i].size();
        if (k > 1) accum += (k - 1);
    }
    total_dummy_cols = accum;
    
    num_covars = total_dummy_cols;
    
    // 3. Encode Data
    for(size_t r=0; r<raw_data_rows.size(); ++r) {
        std::vector<float> encoded(num_covars, 0.0f);
        const auto& row = raw_data_rows[r];
        
        for(size_t i=0; i<row.size(); ++i) {
            int k = levels_per_col[i].size();
            if (k <= 1) continue;
            
            std::string val = row[i];
            if (val == missing_val) {
                 for(int d=0; d < k-1; ++d) {
                    encoded[col_offsets[i] + d] = NAN;
                }
            } else {
                // Check if it is a reference level (first in map)
                auto it = levels_per_col[i].begin();
                if (val != it->first) {
                    // Find index among non-refs
                    int non_ref_idx = 0;
                    bool found = false;
                    it++; // Skip ref
                    while(it != levels_per_col[i].end()) {
                        if (it->first == val) {
                            found = true;
                            break;
                        }
                        non_ref_idx++;
                        it++;
                    }
                    if (found) {
                        encoded[col_offsets[i] + non_ref_idx] = 1.0f;
                    }
                }
                // If it IS reference, all dummies remain 0.
            }
        }
        covar_map[makeKey(row_ids[r].first, row_ids[r].second)] = encoded;
    }

    std::cout << "Loaded " << covar_map.size() << " categorical records. Converted " 
              << levels_per_col.size() << " vars to " << num_covars << " dummy columns.\n";
}

void CovarReader::loadCategoricalWithTime(const std::string& file_path, int time_col_idx, const std::string& missing_val) {
    covar_map.clear();
    num_covars = 0;

    std::ifstream ifs(file_path);
    if (!ifs) throw std::runtime_error("Cannot open covariate file: " + file_path);

    std::string line;
    int line_num = 0;
    std::vector<std::map<std::string, int>> levels_per_col;
    std::vector<std::vector<std::string>> raw_data_rows;
    struct RowId { std::string fid, iid; double time; };
    std::vector<RowId> row_ids;

    while (std::getline(ifs, line)) {
        line_num++;
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string fid, iid;
        ss >> fid >> iid;
        if (fid.empty() || iid.empty()) continue;
        if (line_num == 1 && (fid == "FID" || fid == "fid")) continue;

        std::vector<std::string> tokens;
        std::string token;
        while (ss >> token) tokens.push_back(token);
        if (time_col_idx < 1 || time_col_idx > static_cast<int>(tokens.size())) {
            std::cerr << "Warning: Line " << line_num << " in " << file_path
                      << " lacks TIME column. Skipping.\n";
            continue;
        }

        double time_val = NAN;
        try { time_val = std::stod(tokens[time_col_idx - 1]); }
        catch (...) { continue; }

        std::vector<std::string> row_vals;
        for (int i = 0; i < static_cast<int>(tokens.size()); ++i) {
            if (i == time_col_idx - 1) continue;
            row_vals.push_back(tokens[i]);
        }

        if (levels_per_col.empty()) {
            levels_per_col.resize(row_vals.size());
        } else if (row_vals.size() != levels_per_col.size()) {
            std::cerr << "Warning: Line " << line_num << " in " << file_path
                      << " has " << row_vals.size() << " covariates, expected "
                      << levels_per_col.size() << ". Skipping.\n";
            continue;
        }

        raw_data_rows.push_back(row_vals);
        row_ids.push_back({fid, iid, time_val});

        for (size_t i = 0; i < row_vals.size(); ++i) {
            if (row_vals[i] != missing_val) levels_per_col[i][row_vals[i]]++;
        }
    }

    std::vector<int> col_offsets(levels_per_col.size(), 0);
    int accum = 0;
    for (size_t i = 0; i < levels_per_col.size(); ++i) {
        col_offsets[i] = accum;
        int k = levels_per_col[i].size();
        if (k > 1) accum += (k - 1);
    }
    num_covars = accum;

    for (size_t r = 0; r < raw_data_rows.size(); ++r) {
        std::vector<float> encoded(num_covars, 0.0f);
        const auto& row = raw_data_rows[r];
        for (size_t i = 0; i < row.size(); ++i) {
            int k = levels_per_col[i].size();
            if (k <= 1) continue;

            const std::string& val = row[i];
            if (val == missing_val) {
                for (int d = 0; d < k - 1; ++d) encoded[col_offsets[i] + d] = NAN;
            } else {
                auto it = levels_per_col[i].begin();
                if (val != it->first) {
                    int non_ref_idx = 0;
                    bool found = false;
                    ++it;
                    while (it != levels_per_col[i].end()) {
                        if (it->first == val) { found = true; break; }
                        non_ref_idx++;
                        ++it;
                    }
                    if (found) encoded[col_offsets[i] + non_ref_idx] = 1.0f;
                }
            }
        }
        covar_map[makeTimeKey(row_ids[r].fid, row_ids[r].iid, row_ids[r].time)] = encoded;
    }

    std::cout << "Loaded " << covar_map.size()
              << " record-level categorical records. Converted "
              << levels_per_col.size() << " vars to " << num_covars
              << " dummy columns.\n";
}

void CovarReader::loadWithTime(const std::string& file_path, int time_col_idx, double missing_val) {
    missing_value_code = missing_val;
    covar_map.clear();
    num_covars = 0;

    std::ifstream ifs(file_path);
    if (!ifs) throw std::runtime_error("Cannot open covariate file: " + file_path);

    std::string line;
    int line_num = 0;
    while (std::getline(ifs, line)) {
        line_num++;
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string fid, iid;
        ss >> fid >> iid;
        if (fid.empty() || iid.empty()) continue;
        if (line_num == 1 && (fid == "FID" || fid == "fid")) continue;

        std::vector<std::string> tokens;
        std::string token;
        while (ss >> token) tokens.push_back(token);
        if (time_col_idx < 1 || time_col_idx > static_cast<int>(tokens.size())) continue;

        double time_val = NAN;
        try { time_val = std::stod(tokens[time_col_idx - 1]); }
        catch (...) { continue; }

        std::vector<float> vals;
        for (int i = 0; i < static_cast<int>(tokens.size()); ++i) {
            if (i == time_col_idx - 1) continue;
            try {
                float val = std::stof(tokens[i]);
                vals.push_back(std::abs(val - missing_value_code) < 1e-8 ? NAN : val);
            } catch (...) {
                vals.push_back(NAN);
            }
        }

        if (num_covars == 0) num_covars = vals.size();
        else if (vals.size() != static_cast<size_t>(num_covars)) {
            std::cerr << "Warning: Line " << line_num << " in " << file_path
                      << " has " << vals.size() << " covariates, expected "
                      << num_covars << ". Skipping.\n";
            continue;
        }
        covar_map[makeTimeKey(fid, iid, time_val)] = vals;
    }

    std::cout << "Loaded " << covar_map.size()
              << " record-level covariate records (p=" << num_covars
              << ") from " << file_path << std::endl;
}

void CovarReader::loadMixed(const std::string& file_path, const std::string& missing_val) {
    covar_map.clear();
    num_covars = 0;
    
    std::ifstream ifs(file_path);
    if (!ifs) throw std::runtime_error("Cannot open covariate file: " + file_path);
    
    std::string line;
    std::vector<std::vector<std::string>> raw_data_rows;
    std::vector<std::pair<std::string, std::string>> row_ids;

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string fid, iid;
        ss >> fid >> iid;
        if (fid.empty() || iid.empty()) continue;
        
        std::vector<std::string> row_vals;
        std::string val;
        while (ss >> val) {
            row_vals.push_back(val);
        }
        
        if (!raw_data_rows.empty() && row_vals.size() != raw_data_rows[0].size()) {
            std::cerr << "Warning: row has different number of columns. Skipping.\n";
            continue;
        }

        raw_data_rows.push_back(row_vals);
        row_ids.push_back({fid, iid});
    }

    if (raw_data_rows.empty()) return;

    int n_cols = raw_data_rows[0].size();
    std::vector<bool> is_numeric(n_cols, true);
    std::vector<std::map<std::string, int>> levels_per_col(n_cols);

    for (int c = 0; c < n_cols; ++c) {
        for (size_t r = 0; r < raw_data_rows.size(); ++r) {
            std::string val = raw_data_rows[r][c];
            if (val == missing_val || val == "NA" || val == "NaN") continue;
            
            if (is_numeric[c]) {
                try {
                    size_t pos;
                    std::stod(val, &pos);
                    if (pos != val.length()) {
                        is_numeric[c] = false;
                    }
                } catch (...) {
                    is_numeric[c] = false;
                }
            }
            if (!is_numeric[c]) {
                levels_per_col[c][val]++;
            }
        }
        
        if (!is_numeric[c]) {
            levels_per_col[c].clear();
            for (size_t r = 0; r < raw_data_rows.size(); ++r) {
                std::string val = raw_data_rows[r][c];
                if (val != missing_val && val != "NA" && val != "NaN") {
                    levels_per_col[c][val]++;
                }
            }
        }
    }

    std::vector<int> col_offsets(n_cols, 0);
    int accum = 0;
    for(int c = 0; c < n_cols; ++c) {
        col_offsets[c] = accum;
        if (is_numeric[c]) {
            accum += 1;
        } else {
            int k = levels_per_col[c].size();
            if (k > 1) accum += (k - 1);
            else accum += 1; 
        }
    }
    num_covars = accum;

    for(size_t r=0; r<raw_data_rows.size(); ++r) {
        std::vector<float> encoded(num_covars, 0.0f);
        const auto& row = raw_data_rows[r];
        
        for(int c = 0; c < n_cols; ++c) {
            std::string val = row[c];
            bool is_miss = (val == missing_val || val == "NA" || val == "NaN");

            if (is_numeric[c]) {
                if (is_miss) {
                    encoded[col_offsets[c]] = NAN;
                } else {
                    encoded[col_offsets[c]] = std::stof(val);
                }
            } else {
                int k = levels_per_col[c].size();
                if (k <= 1) {
                    encoded[col_offsets[c]] = is_miss ? NAN : 0.0f;
                    continue;
                }
                
                if (is_miss) {
                    for(int d=0; d < k-1; ++d) encoded[col_offsets[c] + d] = NAN;
                } else {
                    auto it = levels_per_col[c].begin();
                    if (val != it->first) {
                        int non_ref_idx = 0;
                        bool found = false;
                        it++;
                        while(it != levels_per_col[c].end()) {
                            if (it->first == val) { found = true; break; }
                            non_ref_idx++;
                            it++;
                        }
                        if (found) {
                            encoded[col_offsets[c] + non_ref_idx] = 1.0f;
                        }
                    }
                }
            }
        }
        covar_map[makeKey(row_ids[r].first, row_ids[r].second)] = encoded;
    }

    std::cout << "Loaded " << covar_map.size() << " mixed covariates. Final dimension p=" << num_covars << ".\n";
}

size_t CovarReader::align(const std::vector<std::string>& target_fids, 
                          const std::vector<std::string>& target_iids, 
                          Eigen::MatrixXf& out_X) const {
    if (target_fids.size() != target_iids.size()) {
        throw std::runtime_error("Target FID/IID size mismatch in align");
    }
    
    size_t n = target_fids.size();
    
    if (num_covars == 0) {
        // No covariates loaded, return empty matrix? 
        // Or should we throw?
        // Let's resize to n x 0.
        out_X.resize(n, 0);
        return 0;
    }

    out_X.resize(n, num_covars);
    out_X.setConstant(NAN);
    
    size_t matched_count = 0;
    
    for (size_t i = 0; i < n; ++i) {
        std::string key = makeKey(target_fids[i], target_iids[i]);
        auto it = covar_map.find(key);
        if (it != covar_map.end()) {
            bool has_nan = false;
            for(int j=0; j<num_covars; ++j) {
                out_X(i, j) = it->second[j];
                if (std::isnan(out_X(i, j))) has_nan = true;
            }
            if (!has_nan) matched_count++;
        }
    }
    
    return matched_count;
}

size_t CovarReader::alignRecords(const std::vector<std::string>& target_fids,
                                 const std::vector<std::string>& target_iids,
                                 const std::vector<double>& target_times,
                                 Eigen::MatrixXf& out_X) const {
    if (target_fids.size() != target_iids.size() || target_fids.size() != target_times.size()) {
        throw std::runtime_error("Target FID/IID/TIME size mismatch in alignRecords");
    }

    size_t n = target_fids.size();
    if (num_covars == 0) {
        out_X.resize(n, 0);
        return 0;
    }

    out_X.resize(n, num_covars);
    out_X.setConstant(NAN);

    size_t matched_count = 0;
    for (size_t i = 0; i < n; ++i) {
        if (std::isnan(target_times[i])) continue;
        std::string key = makeTimeKey(target_fids[i], target_iids[i], target_times[i]);
        auto it = covar_map.find(key);
        if (it != covar_map.end()) {
            bool has_nan = false;
            for (int j = 0; j < num_covars; ++j) {
                out_X(i, j) = it->second[j];
                if (std::isnan(out_X(i, j))) has_nan = true;
            }
            if (!has_nan) matched_count++;
        }
    }
    return matched_count;
}

} // namespace cosmic
