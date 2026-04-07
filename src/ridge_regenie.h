#ifndef COSMIC_SCANNER_RIDGE_REGENIE_H
#define COSMIC_SCANNER_RIDGE_REGENIE_H

#include <Eigen/Dense>
#include <vector>
#include <string>

namespace cosmic {

/**
 * @brief Regenie Step 1 Skeleton: Block-wise Ridge Regression
 * Reference: REGENIE (Mbatchou et al., 2021)
 * 
 * Core Logic:
 * 1. Divide genome into blocks of SNPs (e.g. 1000 SNPs/block)
 * 2. Fit Ridge Regression with multiple lambda penalties on each block.
 * 3. Use cross-validation (K-fold CV) to select optimal lambda or combine predictions.
 * 4. Level 1: Stack predictions from blocks and run a second level Ridge Regression.
 * 5. Calculate LOCO (Leave-One-Chromosome-Out) predictions to use as covariates in Step 2.
 */
class RegenieStep1 {
public:
    struct Options {
        int block_size = 1000;
        int cv_folds = 5;
        std::vector<double> ridge_lambdas = {0.01, 0.1, 1.0, 10.0, 100.0};
        bool use_loco = true;
        int n_threads = 1;
    };

    RegenieStep1(const Options& opts);
    ~RegenieStep1() = default;

    // Add a block of SNPs and fit Level 0 Ridge Regression
    // Returns true on success
    bool add_block(const Eigen::MatrixXd& X_block, const Eigen::VectorXd& Y, int block_chrom);

    // Level 1: Combine block predictions using out-of-fold predictions
    void fit_level1_stacking(const Eigen::VectorXd& Y);

    // Return the final LOCO prediction for a given chromosome
    Eigen::VectorXd get_loco_prediction(int exclude_chr) const;

private:
    Options options_;
    
    // Internal data structures to store Level 0 predictions
    std::vector<Eigen::MatrixXd> level0_oof_preds_; // Out-of-fold predictions for each block
    std::vector<Eigen::MatrixXd> level0_insample_preds_; // In-sample predictions for each block
    std::vector<int> block_chroms_;
    
    std::vector<Eigen::VectorXd> loco_predictions_; 
    
    // Fit ridge on a single block (computes OOF and In-sample preds for all lambdas)
    void fit_block_ridge(const Eigen::MatrixXd& X_block, const Eigen::VectorXd& Y, Eigen::MatrixXd& oof_preds, Eigen::MatrixXd& insample_preds);
};

} // namespace cosmic

#endif // COSMIC_SCANNER_RIDGE_REGENIE_H