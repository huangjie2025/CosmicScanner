#include "schur.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <Eigen/SparseLU>
#include "pcg_solver.h"

namespace cosmic {
using namespace Eigen; using namespace std;

// Old cg_solve_sparse removed as we use PCGSolver



MatrixXd computeSchurFullSinv(const MatrixXd& Xfull,
                              const vector<GenRecord>& recs,
                              const AbstractMatrix* Qinv,
                              double sigma2_e,
                              double sigma2_u) {
    int na = Qinv->rows();
    double lambda = sigma2_e / sigma2_u;
    MatrixXd XfTXf = Xfull.transpose() * Xfull;
    vector<int> obs_cnt(na, 0);
    for (size_t i = 0; i < recs.size(); ++i) {
        if (recs[i].aid - 1 >= na || recs[i].aid - 1 < 0) {
            std::cerr << "Error: Individual AID " << recs[i].aid << " is out of bounds (1.." << na << ")." << std::endl;
            continue; 
        }
        obs_cnt[recs[i].aid - 1]++;
    }
    vector<Triplet<double>> tA22; tA22.reserve(Qinv->nonZeros() + na);
    
    Qinv->visit_triplets([&](int r, int c, double v) {
        tA22.emplace_back(r, c, lambda * v);
    });
    
    for (int a = 0; a < na; ++a) if (obs_cnt[a] > 0) tA22.emplace_back(a, a, (double)obs_cnt[a]);
    sort(tA22.begin(), tA22.end(), [](const Triplet<double>& a, const Triplet<double>& b){ return (a.row() < b.row()) || (a.row() == b.row() && a.col() < b.col()); });
    SparseMatrix<double> A22(na, na); A22.setFromTriplets(tA22.begin(), tA22.end()); A22.makeCompressed();
    
    // Convert to RowMajor once to avoid per-thread copy in PCGSolver
    SparseMatrix<double, RowMajor> A22_row = A22;

    // Use PCG instead of SparseLU for speed and robustness against dense blocks in SSGBLUP
    int ccols = Xfull.cols();
    MatrixXd Sfull = MatrixXd::Zero(ccols, ccols);
    
    // Parallelize the loop over columns using OpenMP
    // Each thread solves A22 * y = b_j using its own PCGSolver instance.
    // NOTE: OpenMP disabled for robustness (stack/memory issues). Serial execution is fast enough for small p_fixed.
    {
        // Instantiate solver once per thread, referencing the shared RowMajor matrix
        PCGSolver<SparseMatrix<double, RowMajor>> pcg(A22_row, VectorXd::Zero(na), "diag");
        pcg.setQuiet(true);
        
        for (int j = 0; j < ccols; ++j) {
            VectorXd b_j = VectorXd::Zero(na); 
            for (int i = 0; i < Xfull.rows(); ++i) {
                if (recs[i].aid - 1 >= 0 && recs[i].aid - 1 < na)
                    b_j(recs[i].aid - 1) += Xfull(i, j);
            }
            
            // Use PCG to solve A22 * y = b
            // Reuse solver
            VectorXd yj = pcg.solve_with_rhs(b_j, 1e-12, 3000, 3000); 
            
            for (int k = 0; k < ccols; ++k) {
                double acc = 0.0;
                for (int i = 0; i < Xfull.rows(); ++i) {
                    if (recs[i].aid - 1 >= 0 && recs[i].aid - 1 < na)
                        acc += Xfull(i, k) * yj(recs[i].aid - 1);
                }
                // Critical section needed? Sfull is shared.
                // But we write to Sfull(k, j). Different j means different columns.
                // However, k loop iterates all rows for column j.
                // Sfull(k, j) is touched only by thread handling j. So no race condition on Sfull elements.
                // Wait, Sfull(k, j) = ...
                // Yes, j is the outer loop index. Each thread writes to a unique column j. Safe.
                Sfull(k, j) = XfTXf(k, j) - acc;
            }
        }
    }
    // Always use Generalized Inverse via Eigendecomposition to handle rank deficiency
    SelfAdjointEigenSolver<MatrixXd> es(Sfull);
    VectorXd eval = es.eigenvalues(); MatrixXd U = es.eigenvectors(); MatrixXd Dinv = MatrixXd::Zero(ccols, ccols);
    // Use a relative tolerance for eigenvalues
    double max_ev = eval.maxCoeff();
    double tol_ev = max_ev * 1e-9;
    if (tol_ev < 1e-14) tol_ev = 1e-14;

    for (int i = 0; i < ccols; ++i) { 
        double e = eval(i); 
        if (e > tol_ev) Dinv(i, i) = 1.0 / e; 
        else Dinv(i, i) = 0.0;
    }
    MatrixXd Sinv_full = U * Dinv * U.transpose();
    
    return Sinv_full;
}

}
