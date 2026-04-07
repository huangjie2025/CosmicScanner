#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "../../Common/src/eigen_vce.h"

int main() {
    try {
        int n = 100;
        int d = 2; // 2 traits
        int p = 1; // 1 covariate

        // Random Data
        Eigen::MatrixXd Y = Eigen::MatrixXd::Random(n, d);
        Eigen::MatrixXd X = Eigen::MatrixXd::Random(n, p);
        Eigen::MatrixXd G = Eigen::MatrixXd::Random(n, n);
        G = G * G.transpose(); // SPD

        // Setup Engine
        cosmic::MvEigenVce::Options opts;
        opts.verbose = false;
        cosmic::MvEigenVce vce(opts);
        
        std::cout << "Preparing..." << std::endl;
        vce.prepare(Y, X, G);

        std::cout << "Running Null Model..." << std::endl;
        vce.runNullModel();
        
        std::cout << "SUCCESS: Null Model converged." << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
