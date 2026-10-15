#include "models/pedigree_relationship.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <queue>
#include "a_matrix.h"

namespace cosmic {

namespace {

std::string normalizeParent(const std::string& v) {
    return (v == "NA" || v == ".") ? "0" : v;
}

} // anonymous namespace

void PedigreeRelationship::load(const std::string& ped_file) {
    std::ifstream ifs(ped_file);
    if (!ifs.is_open()) throw std::runtime_error("Cannot open pedigree file: " + ped_file);
    std::vector<PedigreeRecord> records;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string id, sire, dam;
        if (!(iss >> id >> sire >> dam)) continue;
        records.push_back({id, normalizeParent(sire), normalizeParent(dam)});
    }
    if (records.empty()) throw std::runtime_error("Pedigree file is empty: " + ped_file);
    load(records);
    std::cout << "Loaded pedigree with " << records_.size() << " individuals.\n";
}

void PedigreeRelationship::load(const std::vector<PedigreeRecord>& records) {
    records_ = records;
    id_to_idx_.clear();
    id_to_idx_.reserve(records_.size());
    for (int i = 0; i < (int)records_.size(); ++i) {
        const auto& r = records_[i];
        if (r.id.empty()) {
            throw std::runtime_error("Pedigree record with empty ID at row " + std::to_string(i + 1));
        }
        if (id_to_idx_.count(r.id)) {
            throw std::runtime_error("Duplicate pedigree ID: " + r.id);
        }
        id_to_idx_[r.id] = i;
    }
    for (int i = 0; i < (int)records_.size(); ++i) {
        const auto& r = records_[i];
        if (r.sire != "0" && r.sire == r.id) {
            throw std::runtime_error("Pedigree record lists itself as sire: " + r.id);
        }
        if (r.dam != "0" && r.dam == r.id) {
            throw std::runtime_error("Pedigree record lists itself as dam: " + r.id);
        }
        for (const std::string& p : {r.sire, r.dam}) {
            if (p == "0") continue;
            if (id_to_idx_.find(p) == id_to_idx_.end()) {
                throw std::runtime_error("Unknown parent ID in pedigree: " + p);
            }
        }
    }

    // Topological order (parents before offspring), stable with respect to the
    // original row order (min-heap by original index). This guarantees that
    // the inbreeding computation in invertA_henderson_fixed is correct even
    // when the input pedigree rows are shuffled.
    const int n = (int)records_.size();
    std::vector<int> indeg(n, 0);
    std::vector<std::vector<int>> children(n);
    for (int i = 0; i < n; ++i) {
        const auto& r = records_[i];
        std::vector<int> pars;
        for (const std::string& p : {r.sire, r.dam}) {
            if (p == "0") continue;
            auto it = id_to_idx_.find(p);
            if (it != id_to_idx_.end()) pars.push_back(it->second);
        }
        std::sort(pars.begin(), pars.end());
        pars.erase(std::unique(pars.begin(), pars.end()), pars.end());
        for (int p : pars) {
            indeg[i]++;
            children[p].push_back(i);
        }
    }
    std::priority_queue<int, std::vector<int>, std::greater<int>> ready;
    for (int i = 0; i < n; ++i) {
        if (indeg[i] == 0) ready.push(i);
    }
    std::vector<int> order;
    order.reserve(n);
    while (!ready.empty()) {
        int i = ready.top();
        ready.pop();
        order.push_back(i);
        for (int c : children[i]) {
            if (--indeg[c] == 0) ready.push(c);
        }
    }
    if ((int)order.size() != n) {
        throw std::runtime_error(
            "Pedigree contains a cycle (an individual is its own ancestor); cannot build A.");
    }
    order_ = order;
}

const Eigen::SparseMatrix<double>& PedigreeRelationship::buildAInverse() {
    if (a_inv_built_) return a_inv_;
    const int n = (int)records_.size();
    if (n == 0) throw std::runtime_error("Pedigree is empty; cannot build A-inverse.");

    std::vector<int> pos_of(n);
    for (int pos = 0; pos < n; ++pos) pos_of[order_[pos]] = pos;

    // dam/sire index vectors in topological order (-1 = unknown parent).
    std::vector<int> dam(n, -1), sire(n, -1);
    for (int pos = 0; pos < n; ++pos) {
        int orig = order_[pos];
        const auto& r = records_[orig];
        if (r.sire != "0") sire[pos] = pos_of[id_to_idx_.at(r.sire)];
        if (r.dam != "0") dam[pos] = pos_of[id_to_idx_.at(r.dam)];
    }

    Eigen::SparseMatrix<double> Ainv_new = genmatrix::invertA_henderson_fixed(dam, sire);
    if (Ainv_new.rows() != n) {
        throw std::runtime_error("A-inverse construction returned a matrix of unexpected size.");
    }

    // Permute back to original input order.
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(Ainv_new.nonZeros());
    for (int k = 0; k < Ainv_new.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Ainv_new, k); it; ++it) {
            trips.emplace_back(order_[it.row()], order_[it.col()], it.value());
        }
    }
    Eigen::SparseMatrix<double> Ainv_orig(n, n);
    Ainv_orig.setFromTriplets(trips.begin(), trips.end());
    a_inv_ = Ainv_orig;
    a_inv_built_ = true;
    return a_inv_;
}

Eigen::MatrixXd PedigreeRelationship::buildA22(const std::vector<int>& analysis_indices) {
    const Eigen::SparseMatrix<double>& Ainv = buildAInverse();
    const int n = (int)Ainv.rows();
    const int n2 = (int)analysis_indices.size();
    if (n2 == 0) {
        throw std::runtime_error("Pedigree A22 requires at least one analysis individual.");
    }
    std::vector<char> is_analysis(n, 0);
    for (int t : analysis_indices) {
        if (t < 0 || t >= n) {
            throw std::runtime_error("Analysis individual index is outside the pedigree matrix.");
        }
        is_analysis[t] = 1;
    }
    std::vector<int> untyped;
    untyped.reserve(n - n2);
    for (int i = 0; i < n; ++i) {
        if (!is_analysis[i]) untyped.push_back(i);
    }
    const int n1 = (int)untyped.size();

    std::vector<int> untyped_pos(n, -1), typed_pos(n, -1);
    for (int i = 0; i < n1; ++i) untyped_pos[untyped[i]] = i;
    for (int i = 0; i < n2; ++i) typed_pos[analysis_indices[i]] = i;

    // Extract A^-1 blocks: [A11 A12; A21 A22] (1 = untyped, 2 = analysis).
    Eigen::MatrixXd A11 = Eigen::MatrixXd::Zero(n1, n1);
    Eigen::MatrixXd A12 = Eigen::MatrixXd::Zero(n1, n2);
    Eigen::MatrixXd Ainv22 = Eigen::MatrixXd::Zero(n2, n2);
    for (int k = 0; k < Ainv.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Ainv, k); it; ++it) {
            const int r = it.row();
            const int c = it.col();
            const bool r_t = is_analysis[r];
            const bool c_t = is_analysis[c];
            if (r_t && c_t) {
                Ainv22(typed_pos[r], typed_pos[c]) += it.value();
            } else if (!r_t && !c_t) {
                A11(untyped_pos[r], untyped_pos[c]) += it.value();
            } else if (!r_t) {
                A12(untyped_pos[r], typed_pos[c]) += it.value();
            }
        }
    }

    // A22 = (Ainv22 - Ainv21 * Ainv11^-1 * Ainv12)^-1 via the Schur complement
    // of the untyped block in A^-1 (standard block-inverse identity). A^-1 is
    // symmetric, so Ainv21 = Ainv12'.
    if (n1 > 0) {
        Eigen::SparseMatrix<double> A11_sp(n1, n1);
        std::vector<Eigen::Triplet<double>> t;
        t.reserve(A11.rows() * A11.cols());
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n1; ++j) {
                if (A11(i, j) != 0.0) t.emplace_back(i, j, A11(i, j));
            }
        }
        A11_sp.setFromTriplets(t.begin(), t.end());
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A11_sp);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Failed to factorize the untyped block of pedigree A-inverse.");
        }
        Eigen::MatrixXd X = solver.solve(A12);  // Ainv11^-1 * Ainv12 (n1 x n2)
        Eigen::MatrixXd S = Ainv22 - A12.transpose() * X;
        S = 0.5 * (S + S.transpose());
        return invertSpdWithRidge(S, "pedigree A22");
    }
    return invertSpdWithRidge(Ainv22, "pedigree A22");
}

Eigen::MatrixXd invertSpdWithRidge(const Eigen::MatrixXd& matrix,
                                   const std::string& label,
                                   double ridge) {
    if (matrix.rows() != matrix.cols()) {
        throw std::runtime_error(label + " must be square.");
    }
    Eigen::MatrixXd sym = 0.5 * (matrix + matrix.transpose());
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(sym.rows(), sym.cols());
    // Attempt the EXACT inverse first (ridge = 0); the geometric ridge is only
    // a fallback for numerically non-positive-definite inputs (this keeps the
    // inverses bit-consistent with e.g. a dense exact reference).
    double r = 0.0;
    for (int attempt = 0; attempt < 6; ++attempt) {
        Eigen::LDLT<Eigen::MatrixXd> ldlt(sym + r * I);
        if (ldlt.info() == Eigen::Success && ldlt.isPositive()) {
            return ldlt.solve(I);
        }
        r = (attempt == 0) ? ridge : r * 10.0;
    }
    Eigen::FullPivLU<Eigen::MatrixXd> lu(sym + r * I);
    if (!lu.isInvertible()) {
        throw std::runtime_error("Cannot invert " + label + ".");
    }
    return lu.inverse();
}

} // namespace cosmic
