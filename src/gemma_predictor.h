#pragma once

#include <string>
#include <vector>

namespace cosmic {

class PlinkReader;
class BgenReader;
class PgenReader;

struct GemmaPredictOptions {
    std::string effect_file;
    std::string output_file;
    std::string snp_col;
    std::string effect_col;
};

class GemmaPredictor {
public:
    explicit GemmaPredictor(GemmaPredictOptions options);

    void run(const std::vector<std::string>& fids,
             const std::vector<std::string>& iids,
             PlinkReader& plink_reader,
             BgenReader& bgen_reader,
             PgenReader& pgen_reader,
             bool is_bgen,
             bool is_pgen);

private:
    GemmaPredictOptions options_;
};

} // namespace cosmic
