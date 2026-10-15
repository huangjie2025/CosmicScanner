#pragma once
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>

namespace cosmic {

// Handles reading covariate files
// Supported format:
// GCTA/PLINK .covar/.qcovar (FID IID Cov1 [Cov2 ...])
// Note: Intercept is NOT added by default, caller should add it.
class CovarReader {
public:
    CovarReader();

    // Read covariate from file
    // file_path: Path to covariate file
    // missing_val: Value representing missing data (default -9)
    // is_categorical: If true, treats all columns as strings and performs one-hot encoding.
    //                 If false, treats as numeric.
    void load(const std::string& file_path, double missing_val = -9.0);
    
    // Load categorical covariates (e.g. Sex, Location) and dummy encode them.
    // Each unique level (except first/reference) gets a column.
    void loadCategorical(const std::string& file_path, const std::string& missing_val = "-9");

    // Load categorical covariates with a record-level TIME column.
    // Columns after FID/IID are 1-based; time_col_idx identifies the TIME column
    // and is excluded from encoded covariates. Records are keyed by FID:IID:TIME.
    void loadCategoricalWithTime(const std::string& file_path, int time_col_idx = 1,
                                 const std::string& missing_val = "-9");

    // Load mixed covariates. Auto-detects column type (numeric vs categorical).
    // If a column has any non-numeric value (other than missing_val), it is treated as categorical
    // and dummy-encoded. Otherwise, it is treated as continuous.
    void loadMixed(const std::string& file_path, const std::string& missing_val = "-9");

    // Numeric covariates with record-level TIME key; TIME is excluded from covariates.
    void loadWithTime(const std::string& file_path, int time_col_idx = 1, double missing_val = -9.0);

    // Align covariates with genotype samples
    // target_fids/iids: The order of samples in the genotype file
    // out_X: Output matrix (n x p), where p is number of covariates loaded.
    //        If some samples are missing in covar file, their rows will contain NAN.
    // Returns number of samples that were found in the covariate file.
    size_t align(const std::vector<std::string>& target_fids, 
                 const std::vector<std::string>& target_iids, 
                 Eigen::MatrixXf& out_X) const;

    // Align record-level covariates by FID/IID/TIME.
    size_t alignRecords(const std::vector<std::string>& target_fids,
                        const std::vector<std::string>& target_iids,
                        const std::vector<double>& target_times,
                        Eigen::MatrixXf& out_X) const;

    int getNumCovars() const { return num_covars; }

private:
    // Key: "FID:IID", Value: Vector of covariates
    std::map<std::string, std::vector<float>> covar_map; 
    int num_covars = 0;
    double missing_value_code;
    
    std::string makeKey(const std::string& fid, const std::string& iid) const {
        return fid + ":" + iid;
    }

    std::string makeTimeKey(const std::string& fid, const std::string& iid, double time) const;
};

} // namespace cosmic
