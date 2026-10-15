#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace cosmic {

struct PedigreeRecord {
    std::string id;
    std::string sire;  // "0" if unknown
    std::string dam;   // "0" if unknown
};

/// Shared pedigree A / A22 builder for the pedigree-based Simple Animal Model
/// and for ssGWAS.
///
/// - Parses `ID SIRE DAM` (no header), normalizes "NA"/"." to "0".
/// - Topologically orders the pedigree (parents before offspring) so that
///   `genmatrix::invertA_henderson_fixed` is numerically correct regardless of
///   the row order in the input file (shuffle-invariant by construction).
/// - Returns A^-1 as a sparse matrix in ORIGINAL input row order.
/// - Extracts the A22 block of an analysis (typed) subset via the Schur
///   complement of the analysis block within A^-1, avoiding a full dense A.
class PedigreeRelationship {
public:
    void load(const std::string& ped_file);
    void load(const std::vector<PedigreeRecord>& records);

    int nPedigree() const { return (int)records_.size(); }
    const std::vector<PedigreeRecord>& records() const { return records_; }
    const std::unordered_map<std::string, int>& idToIndex() const { return id_to_idx_; }

    /// Sparse A^-1 in original input order (rows/cols = original indices).
    const Eigen::SparseMatrix<double>& buildAInverse();

    /// Dense A22 in the order given by `analysis_indices` (original-order
    /// indices of the analysis individuals). Uses the Schur complement of the
    /// analysis block within A^-1; no full dense A is materialized.
    Eigen::MatrixXd buildA22(const std::vector<int>& analysis_indices);

private:
    std::vector<PedigreeRecord> records_;
    std::unordered_map<std::string, int> id_to_idx_;
    std::vector<int> order_;  // topological order: position -> original index
    mutable bool a_inv_built_ = false;
    mutable Eigen::SparseMatrix<double> a_inv_;
};

/// Invert an SPD-ish dense matrix with a ridge that grows geometrically until
/// an LDLT with positive diagonal succeeds; falls back to FullPivLU.
Eigen::MatrixXd invertSpdWithRidge(const Eigen::MatrixXd& matrix,
                                   const std::string& label,
                                   double ridge = 1e-8);

} // namespace cosmic
