#pragma once
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "plink_reader.h"

namespace cosmic {

// Helper class to manage subsetting of samples based on phenotype/covariate availability
class DataFilter {
public:
    struct KeepInfo {
        std::vector<int> keep_indices; // Indices in the original PLINK file
        size_t n_keep;
    };

    // Filter samples based on missing phenotype and covariates
    // y: Phenotype vector (n_total) - NAN indicates missing
    // X: Covariate matrix (n_total x p) - NAN indicates missing
    // Returns KeepInfo containing valid indices
    static KeepInfo filterSamples(const Eigen::VectorXd& y, const Eigen::MatrixXf& X) {
        KeepInfo info;
        size_t n = y.size();
        if (X.rows() != n) throw std::runtime_error("DataFilter: X and y dimension mismatch");

        info.keep_indices.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            bool keep = true;
            if (std::isnan(y[i])) keep = false;
            else {
                for (int j = 0; j < X.cols(); ++j) {
                    if (std::isnan(X(i, j))) {
                        keep = false;
                        break;
                    }
                }
            }
            if (keep) {
                info.keep_indices.push_back(i);
            }
        }
        info.n_keep = info.keep_indices.size();
        return info;
    }

    // Subset VectorXd
    static Eigen::VectorXd subset(const Eigen::VectorXd& v, const KeepInfo& info) {
        Eigen::VectorXd res(info.n_keep);
        for (size_t i = 0; i < info.n_keep; ++i) {
            res[i] = v[info.keep_indices[i]];
        }
        return res;
    }

    // Subset MatrixXf (Rows)
    static Eigen::MatrixXf subsetRows(const Eigen::MatrixXf& M, const KeepInfo& info) {
        Eigen::MatrixXf res(info.n_keep, M.cols());
        for (size_t i = 0; i < info.n_keep; ++i) {
            res.row(i) = M.row(info.keep_indices[i]);
        }
        return res;
    }

    // Subset MatrixXd (Rows and Cols) - for Double Precision GRM
    static Eigen::MatrixXd subsetSquare(const Eigen::MatrixXd& M, const KeepInfo& info) {
        if (M.rows() != M.cols()) throw std::runtime_error("subsetSquare: Matrix not square");
        Eigen::MatrixXd res(info.n_keep, info.n_keep);
        for (size_t i = 0; i < info.n_keep; ++i) {
            for (size_t j = 0; j < info.n_keep; ++j) {
                res(i, j) = M(info.keep_indices[i], info.keep_indices[j]);
            }
        }
        return res;
    }

    // Subset MatrixXf (Rows and Cols) - for GRM
    static Eigen::MatrixXf subsetSquare(const Eigen::MatrixXf& M, const KeepInfo& info) {
        if (M.rows() != M.cols()) throw std::runtime_error("subsetSquare: Matrix not square");
        Eigen::MatrixXf res(info.n_keep, info.n_keep);
        // This is slow for large matrices, but necessary if we computed full GRM.
        // Better to compute GRM only for subset if possible, but PlinkReader works on full file.
        // We can optimize GRM builder to take a subset list later.
        for (size_t i = 0; i < info.n_keep; ++i) {
            for (size_t j = 0; j < info.n_keep; ++j) {
                res(i, j) = M(info.keep_indices[i], info.keep_indices[j]);
            }
        }
        return res;
    }
};

} // namespace cosmic
