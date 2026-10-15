#pragma once

#include <string>
#include <vector>

namespace cosmic {
namespace genetics {

struct QCOptions {
    double geno = 1;
    double mind = 1;
    double maf = .01;
    double hwe = 0;
    std::string keep_file;
    std::string remove_file;
    std::string extract_file;
    std::string exclude_file;
};

struct QCResult {
    std::string prefix;
    std::vector<int> samples;
    std::vector<int> variants;
};

class QCFilter {
public:
    // BED hardcalls only. Indices refer to the input, in original order.
    // geno/mind are maximum missing fractions; maf/hwe are inclusive minima.
    // Variants with no surviving called genotypes are always rejected.
    // Writes output_prefix.{bed,bim,fam,qc_report.txt}; never overwrites.
    static QCResult preprocess(std::string input_prefix, std::string output_prefix,
                               QCOptions options);
};

} // namespace genetics
} // namespace cosmic
