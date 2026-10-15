#include "ridge_regenie.h"
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

using namespace cosmic;
using namespace Eigen;

void test_regenie_step1() {
    std::cout << "--- Testing REGENIE Step 1 (Ridge + LOCO) ---" << std::endl;
    
    RegenieStep1::Options opts;
    opts.block_size = 50;
    opts.cv_folds = 2; // Use 2 folds for small dataset
    opts.ridge_lambdas = {0.1, 1.0, 10.0};
    opts.use_loco = true;
    
    RegenieStep1 regenie(opts);
    
    int n = 100;
    int m = 150; // 3 blocks of 50 SNPs
    
    // Create synthetic genotype matrix
    MatrixXd X = MatrixXd::Random(n, m);
    
    // Create true effect for SNPs 0-49 (Chrom 1), 50-99 (Chrom 2), 100-149 (Chrom 3)
    VectorXd beta = VectorXd::Zero(m);
    for (int i = 0; i < 10; ++i) beta(i) = 1.0;         // Chr 1 signal
    for (int i = 50; i < 60; ++i) beta(i) = 0.5;        // Chr 2 signal
    for (int i = 100; i < 110; ++i) beta(i) = -0.8;     // Chr 3 signal
    
    VectorXd Y = X * beta + VectorXd::Random(n) * 0.1; // Small noise
    
    // Add blocks
    bool ok1 = regenie.add_block(X.middleCols(0, 50), Y, 1);
    bool ok2 = regenie.add_block(X.middleCols(50, 50), Y, 2);
    bool ok3 = regenie.add_block(X.middleCols(100, 50), Y, 3);
    
    assert(ok1 && ok2 && ok3);
    
    // Fit Level 1
    regenie.fit_level1_stacking(Y);
    
    // Check LOCO predictions
    VectorXd loco_chr1 = regenie.get_loco_prediction(1);
    VectorXd loco_chr2 = regenie.get_loco_prediction(2);
    VectorXd loco_chr3 = regenie.get_loco_prediction(3);
    
    std::cout << "LOCO Chr1 norm: " << loco_chr1.norm() << std::endl;
    std::cout << "LOCO Chr2 norm: " << loco_chr2.norm() << std::endl;
    std::cout << "LOCO Chr3 norm: " << loco_chr3.norm() << std::endl;
    
    assert(loco_chr1.size() == n);
    assert(loco_chr2.size() == n);
    assert(loco_chr3.size() == n);
    
    // Non-LOCO prediction
    VectorXd global_pred = regenie.get_loco_prediction(0);
    std::cout << "Global Pred norm: " << global_pred.norm() << std::endl;
    
    std::cout << "PASS: REGENIE Step 1 completed successfully without crashing." << std::endl;
}

int main() {
    std::cout << "=======================================" << std::endl;
    std::cout << "   CosmicScanner REGENIE Unit Test     " << std::endl;
    std::cout << "=======================================" << std::endl;
    
    test_regenie_step1();
    
    std::cout << "ALL REGENIE TESTS PASSED!" << std::endl;
    return 0;
}
