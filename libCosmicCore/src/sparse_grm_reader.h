#pragma once
#include <string>
#include <vector>
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>

namespace cosmic {

class SparseGrmReader {
public:
    struct Meta {
        int n_samples;
        float threshold;
        std::vector<std::string> fids;
        std::vector<std::string> iids;
    };

    static Meta readMeta(const std::string& prefix);
    
    static void read(const std::string& prefix, Eigen::SparseMatrix<float>& G);
    static void read(const std::string& prefix, Eigen::SparseMatrix<double>& G);
};

} // namespace cosmic
