#include "ridge_regenie.h"
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <omp.h>

namespace cosmic {

RegenieStep1::RegenieStep1(const Options& opts) : options_(opts) {
    loco_predictions_.resize(22, Eigen::VectorXd());
}

bool RegenieStep1::add_block(const Eigen::Ref<const Eigen::MatrixXd>& X_block_in, const Eigen::VectorXd& Y, int block_chrom) {
    int n = X_block_in.rows();
    int p = X_block_in.cols();
    
    if (n == 0 || p == 0) return false;
    
    // Center and scale X_block
    Eigen::MatrixXd X_block = X_block_in;
    for (int j = 0; j < p; ++j) {
        double mean = X_block.col(j).mean();
        X_block.col(j).array() -= mean;
        double sd = std::sqrt(X_block.col(j).squaredNorm() / (n - 1));
        if (sd > 1e-8) {
            X_block.col(j).array() /= sd;
        }
    }
    
    Eigen::MatrixXd oof_preds;
    fit_block_ridge(X_block, Y, oof_preds);
    
    level0_oof_preds_.push_back(oof_preds);
    block_chroms_.push_back(block_chrom);
    
    return true;
}

void RegenieStep1::fit_block_ridge(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, Eigen::MatrixXd& oof_preds) {
    int n = X.rows();
    int p = X.cols();
    int k_folds = options_.cv_folds;
    int num_lambdas = options_.ridge_lambdas.size();
    double lambda_scale = options_.n_model_snps > 0 ? (double)options_.n_model_snps : (double)p;
    
    oof_preds = Eigen::MatrixXd::Zero(n, num_lambdas);
    Eigen::MatrixXd XtX_full = X.transpose() * X;
    Eigen::VectorXd XtY_full = X.transpose() * Y;
    
    int fold_size = n / k_folds;
    
    // K-fold CV for Out-Of-Fold predictions
    for (int k = 0; k < k_folds; ++k) {
        int start = k * fold_size;
        int end = (k == k_folds - 1) ? n : start + fold_size;
        int test_n = end - start;
        auto X_test = X.middleRows(start, test_n);
        auto Y_test = Y.segment(start, test_n);
        Eigen::MatrixXd XtX = XtX_full - X_test.transpose() * X_test;
        Eigen::VectorXd XtY = XtY_full - X_test.transpose() * Y_test;
        
        for (int i = 0; i < num_lambdas; ++i) {
            Eigen::MatrixXd XtX_reg = XtX;
            double h2 = options_.ridge_lambdas[i];
            double lambda = lambda_scale * (1.0 - h2) / h2;
            XtX_reg.diagonal().array() += lambda;
            Eigen::VectorXd beta = XtX_reg.ldlt().solve(XtY);
            oof_preds.block(start, i, test_n, 1) = X_test * beta;
        }
    }
}

void RegenieStep1::fit_level1_stacking(const Eigen::VectorXd& Y) {
    int n = Y.size();
    int m_blocks = level0_oof_preds_.size();
    int num_lambdas = options_.ridge_lambdas.size();
    int num_features = m_blocks * num_lambdas;
    
    if (m_blocks == 0) return;
    
    std::cout << "  REGENIE Level-1 stacking: " << m_blocks << " blocks, " << num_features << " features." << std::endl;
    
    // Construct X_level1_oof
    Eigen::MatrixXd X_level1_oof(n, num_features);
    for (int b = 0; b < m_blocks; ++b) {
        X_level1_oof.middleCols(b * num_lambdas, num_lambdas) = level0_oof_preds_[b];
    }
    
    // Scale X_level1_oof
    Eigen::VectorXd scale_factors(num_features);
    for (int j = 0; j < num_features; ++j) {
        double mean = X_level1_oof.col(j).mean();
        X_level1_oof.col(j).array() -= mean;
        double sd = std::sqrt(X_level1_oof.col(j).squaredNorm() / (n - 1));
        if (sd < 1e-8) sd = 1.0;
        X_level1_oof.col(j).array() /= sd;
        scale_factors(j) = sd;
    }
    
    std::vector<double> l1_h2 = options_.ridge_lambdas;
    int k_folds = options_.cv_folds;
    int fold_size = n / k_folds;
    std::vector<double> cv_errors(l1_h2.size(), 0.0);
    std::vector<std::vector<Eigen::VectorXd>> fold_betas(k_folds, std::vector<Eigen::VectorXd>(l1_h2.size()));
    
    for (int k = 0; k < k_folds; ++k) {
        int start = k * fold_size;
        int end = (k == k_folds - 1) ? n : start + fold_size;
        int test_n = end - start;
        int train_n = n - test_n;
        
        Eigen::MatrixXd X_train(train_n, num_features);
        Eigen::VectorXd Y_train(train_n);
        Eigen::MatrixXd X_test = X_level1_oof.middleRows(start, test_n);
        Eigen::VectorXd Y_test = Y.segment(start, test_n);
        
        if (start > 0) {
            X_train.topRows(start) = X_level1_oof.topRows(start);
            Y_train.head(start) = Y.head(start);
        }
        if (end < n) {
            X_train.bottomRows(n - end) = X_level1_oof.bottomRows(n - end);
            Y_train.tail(n - end) = Y.tail(n - end);
        }
        
        Eigen::MatrixXd XtX = X_train.transpose() * X_train;
        Eigen::VectorXd XtY = X_train.transpose() * Y_train;
        
        for (size_t i = 0; i < l1_h2.size(); ++i) {
            Eigen::MatrixXd XtX_reg = XtX;
            double h2 = l1_h2[i];
            double tau = (double)num_features * (1.0 - h2) / h2;
            XtX_reg.diagonal().array() += tau;
            Eigen::VectorXd alpha = XtX_reg.ldlt().solve(XtY);
            fold_betas[k][i] = alpha;
            Eigen::VectorXd preds = X_test * alpha;
            cv_errors[i] += (preds - Y_test).squaredNorm();
        }
    }
    
    auto min_it = std::min_element(cv_errors.begin(), cv_errors.end());
    size_t best_idx = (size_t)std::distance(cv_errors.begin(), min_it);
    
    // Generate LOCO predictions
    if (options_.use_loco) {
        Eigen::MatrixXd chrom_contrib = Eigen::MatrixXd::Zero(n, 23);
        for (int k = 0; k < k_folds; ++k) {
            int start = k * fold_size;
            int end = (k == k_folds - 1) ? n : start + fold_size;
            int test_n = end - start;
            const Eigen::VectorXd& alpha = fold_betas[k][best_idx];
            for (int b = 0; b < m_blocks; ++b) {
                int chr = block_chroms_[b];
                if (chr <= 0 || chr > 22) continue;
                int feat_start = b * num_lambdas;
                chrom_contrib.block(start, chr, test_n, 1) +=
                    X_level1_oof.block(start, feat_start, test_n, num_lambdas) *
                    alpha.segment(feat_start, num_lambdas);
            }
        }
        Eigen::VectorXd total_pred = chrom_contrib.rowwise().sum();
        for (int chr = 1; chr <= 22; ++chr) {
            loco_predictions_[chr - 1] = total_pred - chrom_contrib.col(chr);
        }
    } else {
        Eigen::VectorXd global_pred = Eigen::VectorXd::Zero(n);
        for (int k = 0; k < k_folds; ++k) {
            int start = k * fold_size;
            int end = (k == k_folds - 1) ? n : start + fold_size;
            int test_n = end - start;
            global_pred.segment(start, test_n) = X_level1_oof.middleRows(start, test_n) * fold_betas[k][best_idx];
        }
        loco_predictions_[0] = global_pred;
    }
}

Eigen::VectorXd RegenieStep1::get_loco_prediction(int exclude_chr) const {
    if (exclude_chr >= 0 && exclude_chr <= 22) {
        if (exclude_chr == 0) return loco_predictions_[0]; // Non-LOCO
        if (loco_predictions_[exclude_chr - 1].size() > 0) {
            return loco_predictions_[exclude_chr - 1];
        }
    }
    return Eigen::VectorXd::Zero(1); // Should not happen if correctly populated
}

} // namespace cosmic
