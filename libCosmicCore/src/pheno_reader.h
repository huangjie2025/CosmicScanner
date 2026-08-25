#pragma once
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>

namespace cosmic {

// Handles reading phenotype files
// Supported formats:
// 1. GCTA/PLINK .pheno (FID IID Pheno1 [Pheno2 ...])
// 2. .fam file (FID IID ... Pheno)
class PhenoReader {
public:
    struct Sample {
        std::string fid;
        std::string iid;
    };

    PhenoReader();

    // Read phenotype from file
    // file_path: Path to phenotype file
    // col_idx: 1-based index of phenotype column (default 1). 
    //          If reading from .fam, col_idx is ignored (always last column).
    // missing_val: Value representing missing data (default -9)
    void load(const std::string& file_path, int col_idx = 1, double missing_val = -9.0);

    // Load multiple columns (Multivariate)
    // cols: List of 1-based column indices. If empty, load ALL columns after FID IID.
    void loadMatrix(const std::string& file_path, const std::vector<int>& cols = {}, double missing_val = -9.0);

    // Align phenotype with genotype samples
    // target_fids/iids: The order of samples in the genotype file (e.g. from PlinkReader)
    // Returns aligned phenotype vector (same size as target samples)
    // Missing samples are set to NAN
    // Returns number of matched samples
    size_t align(const std::vector<std::string>& target_fids, 
                 const std::vector<std::string>& target_iids, 
                 Eigen::VectorXd& out_y) const;

    // Align matrix phenotype
    size_t alignMatrix(const std::vector<std::string>& target_fids, 
                       const std::vector<std::string>& target_iids, 
                       Eigen::MatrixXd& out_Y) const;

    // Structure for raw phenotype record (Long Format)
    struct PhenoRecord {
        std::string fid;
        std::string iid;
        double value;
        double time; // Time variable for longitudinal data (Random Regression Model)
        int line_num;
    };

    // Load all records (allows duplicates for repeated measures)
    // Returns a vector of all valid records found in the file.
    std::vector<PhenoRecord> loadAllRecords(const std::string& file_path, int col_idx = 1, double missing_val = -9.0);

    // Load longitudinal records (FID IID Time Pheno)
    // t_col_idx: 1-based index of time column
    // y_col_idx: 1-based index of phenotype column
    std::vector<PhenoRecord> loadLongitudinalRecords(const std::string& file_path, int t_col_idx, int y_col_idx, double missing_val = -9.0);

private:
    std::map<std::string, double> pheno_map; // Key: "FID:IID", Value: Phenotype
    std::map<std::string, Eigen::VectorXd> pheno_matrix_map; // Key: "FID:IID", Value: Phenotype Vector
    double missing_value_code;
    int num_traits = 0;
    
    std::string makeKey(const std::string& fid, const std::string& iid) const {
        return fid + ":" + iid;
    }
};

} // namespace cosmic
