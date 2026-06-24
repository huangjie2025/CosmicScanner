#pragma once

#include <string>
#include <vector>

namespace cosmic {

struct InputPrepareOptions {
    std::string input_file;
    std::string input_format = "auto";
    std::string input_model = "repeatability";
    std::string output_dir;
    std::string id_col;
    std::string iid_col;
    std::string time_col;
    std::vector<std::string> trait_cols;
    std::string input_trait;
    std::vector<std::string> record_qcovars;
    std::vector<std::string> record_ccovars;
    std::vector<std::string> subject_qcovars;
    std::vector<std::string> subject_ccovars;
    std::string time_parser = "auto";
    std::string time_unique_policy = "within-time-index";
    std::string time_index_multiplier = "auto";
    bool preview_schema = false;
    bool infer_columns = false;
    bool wide_to_long = false;
};

struct PreparedInputResult {
    std::string output_dir;
    std::string pheno_file;
    std::string rrm_pheno_file;
    std::string subject_pheno_file;
    std::string multivariate_pheno_file;
    std::string qcovar_file;
    std::string covar_file;
    std::string subject_qcovar_file;
    std::string subject_covar_file;
    std::string selected_trait;
    int raw_rows = 0;
    int prepared_records = 0;
    int matched_records = 0;
    int unmatched_records = 0;
};

class InputPreparer {
public:
    static PreparedInputResult prepare(const InputPrepareOptions& options,
                                       const std::vector<std::string>& genotype_fids,
                                       const std::vector<std::string>& genotype_iids);

    static PreparedInputResult loadPrepared(const std::string& prepared_dir,
                                            const std::string& input_trait,
                                            const std::string& input_model = "repeatability");
};

} // namespace cosmic
