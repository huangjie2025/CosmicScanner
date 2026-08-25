#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "eigen_vce.h"

using namespace cosmic;
using namespace Eigen;

int main() {
    std::cout << "Testing Sparse REML support in EigenVce..." << std::endl;
    
    int n = 100;
    VectorXd y = VectorXd::Random(n);
    MatrixXf X = MatrixXf::Ones(n, 1);
    
    // Create a sparse diagonal-dominant matrix
    SparseMatrix<double> G_sparse(n, n);
    std::vector<Triplet<double>> triplets;
    for(int i = 0; i < n; ++i) {
        triplets.push_back(Triplet<double>(i, i, 1.0 + 0.1 * (rand() % 100) / 100.0));
        if (i < n - 1) {
            triplets.push_back(Triplet<double>(i, i + 1, 0.1));
            triplets.push_back(Triplet<double>(i + 1, i, 0.1));
        }
    }
    G_sparse.setFromTriplets(triplets.begin(), triplets.end());
    
    EigenVce::Options opts;
    opts.use_covariate_projection = true;
    EigenVce vce(opts);
    
    try {
        vce.prepareSparse(y, X, G_sparse);
        vce.runNullModel();
        
        std::cout << "Sparse VCE successful!" << std::endl;
        std::cout << "Vg = " << vce.getSigma2g() << std::endl;
        std::cout << "Ve = " << vce.getSigma2e() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error during Sparse VCE: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
