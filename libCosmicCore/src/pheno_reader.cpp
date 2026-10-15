#include "pheno_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cosmic {

PhenoReader::PhenoReader() : missing_value_code(-9.0) {}

void PhenoReader::load(const std::string& file_path, int col_idx, double missing_val) {
    missing_value_code = missing_val;
    pheno_map.clear();
    
    std::ifstream ifs(file_path);
    if (!ifs) throw std::runtime_error("Cannot open phenotype file: " + file_path);
    
    std::string line;
    int line_num = 0;
    while (std::getline(ifs, line)) {
        line_num++;
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string fid, iid;
        ss >> fid >> iid;
        
        if (fid.empty() || iid.empty()) continue;
        
        std::string val_str;
        double val = NAN;
        bool found = false;
        
        for (int i = 0; i < col_idx; ++i) {
            if (!(ss >> val_str)) {
                found = false; 
                break;
            }
            found = true;
        }
        
        if (found) {
            if (val_str == "NA" || val_str == "na" || val_str == "NaN" || val_str == "nan" || val_str == "NAN") {
                val = NAN;
            } else {
                try {
                    val = std::stod(val_str);
                } catch (...) {
                    val = NAN;
                }
            }
            if (std::isnan(val) || std::abs(val - missing_value_code) < 1e-6) {
                // It is missing
                val = NAN;
            }
            pheno_map[makeKey(fid, iid)] = val;
        }
    }
    
    std::cout << "Loaded " << pheno_map.size() << " phenotype records from " << file_path << std::endl;
}

void PhenoReader::loadMatrix(const std::string& file_path, const std::vector<int>& cols, double missing_val) {
    missing_value_code = missing_val;
    pheno_matrix_map.clear();
    num_traits = 0;
    
    std::ifstream ifs(file_path);
    if (!ifs) throw std::runtime_error("Cannot open phenotype file: " + file_path);
    
    std::string line;
    int line_num = 0;
    
    while (std::getline(ifs, line)) {
        line_num++;
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string fid, iid;
        ss >> fid >> iid;
        
        if (fid.empty() || iid.empty()) continue;
        
        std::vector<double> vals;
        double val;
        
        // Read all values
        std::vector<double> row_vals;
        while (ss >> val) {
            if (std::abs(val - missing_value_code) < 1e-8) {
                row_vals.push_back(NAN);
            } else {
                row_vals.push_back(val);
            }
        }
        
        // Select columns
        if (cols.empty()) {
            // Use all
            vals = row_vals;
        } else {
            for(int c : cols) {
                int idx = c - 1; // 1-based to 0-based
                if (idx >= 0 && idx < (int)row_vals.size()) {
                    vals.push_back(row_vals[idx]);
                } else {
                    vals.push_back(NAN); // Column missing
                }
            }
        }
        
        if (num_traits == 0) {
            num_traits = vals.size();
        } else if ((int)vals.size() != num_traits) {
            std::cerr << "Warning: Line " << line_num << " has different column count. Skipping.\n";
            continue;
        }
        
        Eigen::VectorXd vec(num_traits);
        for(int i=0; i<num_traits; ++i) vec[i] = vals[i];
        
        pheno_matrix_map[makeKey(fid, iid)] = vec;
    }
    
    std::cout << "Loaded " << pheno_matrix_map.size() << " multivariate phenotype records (" << num_traits << " traits) from " << file_path << std::endl;
}

size_t PhenoReader::align(const std::vector<std::string>& target_fids, 
                          const std::vector<std::string>& target_iids, 
                          Eigen::VectorXd& out_y) const {
    if (target_fids.size() != target_iids.size()) {
        throw std::runtime_error("Target FID/IID size mismatch in align");
    }
    
    size_t n = target_fids.size();
    out_y.resize(n);
    out_y.setConstant(NAN);
    
    size_t matched_count = 0;
    
    if (n > 0 && pheno_map.size() > 0) {
        std::cout << "[DEBUG] align: Target[0] = '" << makeKey(target_fids[0], target_iids[0]) 
                  << "', PhenoMap first key = '" << pheno_map.begin()->first << "'" << std::endl;
    }

    for (size_t i = 0; i < n; ++i) {
        std::string key = makeKey(target_fids[i], target_iids[i]);
        auto it = pheno_map.find(key);
        if (it != pheno_map.end()) {
            out_y[i] = it->second;
            if (!std::isnan(out_y[i])) {
                matched_count++;
            }
        }
    }

    // A trait with no variation carries no information: every marker would receive a
    // degenerate estimate and an ordinary-looking P-value. Refuse instead of emitting
    // significance-looking output for a constant trait. Non-finite values are excluded
    // here: they are reported by the input-validation stage (Validation : FAIL).
    double sum = 0.0;
    double sum_sq = 0.0;
    size_t finite_count = 0;
    for (size_t i = 0; i < n; ++i) {
        if (!std::isfinite(out_y[i])) continue;
        sum += out_y[i];
        sum_sq += out_y[i] * out_y[i];
        ++finite_count;
    }
    if (finite_count >= 2) {
        const double mean = sum / static_cast<double>(finite_count);
        const double mean_sq = sum_sq / static_cast<double>(finite_count);
        if (!(mean_sq - mean * mean > 1e-12)) {
            throw std::runtime_error(
                "Phenotype has zero variance across " + std::to_string(finite_count) +
                " finite matched samples (every observed value is identical to " + std::to_string(mean) +
                "). Association tests are not estimable for a constant trait; check the "
                "phenotype column (--mpheno) and the sample matching.");
        }
    }

    return matched_count;
}

size_t PhenoReader::alignMatrix(const std::vector<std::string>& target_fids, 
                                const std::vector<std::string>& target_iids, 
                                Eigen::MatrixXd& out_Y) const {
    if (target_fids.size() != target_iids.size()) {
        throw std::runtime_error("Target FID/IID size mismatch in align");
    }
    
    size_t n = target_fids.size();
    out_Y.resize(n, num_traits);
    out_Y.setConstant(NAN);
    
    size_t matched_count = 0;
    
    for (size_t i = 0; i < n; ++i) {
        std::string key = makeKey(target_fids[i], target_iids[i]);
        auto it = pheno_matrix_map.find(key);
        if (it != pheno_matrix_map.end()) {
            out_Y.row(i) = it->second;
            // Check if full row is valid? Or partial?
            // Usually for multivariate, we require full phenotype or handle missingness.
            // Current MvLmmEngine assumes complete data (or handled before).
            // Let's count if at least one trait is not NAN.
            if (!it->second.hasNaN()) {
                matched_count++;
            }
        }
    }
    
    return matched_count;
}

std::vector<PhenoReader::PhenoRecord> PhenoReader::loadAllRecords(const std::string& file_path, int col_idx, double missing_val) {
    std::vector<PhenoRecord> records;
    missing_value_code = missing_val;
    
    std::ifstream ifs(file_path);
    if (!ifs) throw std::runtime_error("Cannot open phenotype file: " + file_path);
    
    std::string line;
    int line_num = 0;
    while (std::getline(ifs, line)) {
        line_num++;
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string fid, iid;
        
        if (!(ss >> fid >> iid)) continue;
        
        // Simple header check
        if (line_num == 1 && (fid == "FID" || fid == "fid")) continue;
        
        double val = NAN;
        bool found = false;
        
        for (int i = 0; i < col_idx; ++i) {
            if (!(ss >> val)) {
                found = false; 
                break;
            }
            found = true;
        }
        
        if (found) {
            if (std::abs(val - missing_value_code) < 1e-8) {
                val = NAN;
            }
            // For now we don't filter NAN here, but usually we skip them?
            // load() keeps NAN in map but align() might use it.
            // Let's keep NAN to maintain record count if needed, but usually we drop missing for VCE.
            // Let's keep it consistent with load().
            records.push_back({fid, iid, val, line_num});
        }
    }
    
    std::cout << "Loaded " << records.size() << " records from " << file_path << " (Long Format)" << std::endl;
    return records;
}

std::vector<PhenoReader::PhenoRecord> PhenoReader::loadLongitudinalRecords(const std::string& file_path, int t_col_idx, int y_col_idx, double missing_val) {
    std::vector<PhenoRecord> records;
    missing_value_code = missing_val;
    
    std::ifstream ifs(file_path);
    if (!ifs) throw std::runtime_error("Cannot open phenotype file: " + file_path);
    
    std::string line;
    int line_num = 0;
    while (std::getline(ifs, line)) {
        line_num++;
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string fid, iid;
        
        if (!(ss >> fid >> iid)) continue;
        
        // Simple header check
        if (line_num == 1 && (fid == "FID" || fid == "fid")) continue;
        
        // We need to extract columns at t_col_idx and y_col_idx (1-based)
        int max_col = std::max(t_col_idx, y_col_idx);
        std::vector<std::string> tokens;
        std::string token;
        while (ss >> token) {
            tokens.push_back(token);
        }
        
        if (tokens.size() < max_col) continue;
        
        try {
            double time_val = std::stod(tokens[t_col_idx - 1]);
            double pheno_val = std::stod(tokens[y_col_idx - 1]);
            
            if (std::abs(pheno_val - missing_value_code) < 1e-8) {
                pheno_val = NAN;
            }
            
            records.push_back({fid, iid, pheno_val, time_val, line_num});
        } catch (...) {
            // Ignore parsing errors (e.g. non-numeric values)
            continue;
        }
    }
    
    std::cout << "Loaded " << records.size() << " longitudinal records from " << file_path << std::endl;
    return records;
}

} // namespace cosmic
