#include <iostream>
#include <Eigen/Dense>

int main() {
    std::cout << "Allocating 5000x5000..." << std::endl;
    Eigen::MatrixXd G = Eigen::MatrixXd::Random(5000, 5000);
    G = G + G.transpose().eval(); // make symmetric
    std::cout << "Solving..." << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(G);
    std::cout << "Done! First eigenvalue: " << es.eigenvalues()(0) << std::endl;
    return 0;
}
