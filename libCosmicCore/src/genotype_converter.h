#pragma once
#include <string>
#include <cmath>

namespace cosmic {

class GenotypeConverter {
public:
    // Convert HapMap string genotype (e.g., "AA", "AG", "NN") to dosage (0, 1, 2)
    // ref: Reference allele character
    // alt: Alternate allele character
    // Returns: 
    //   0.0 if Homozygous REF (e.g., "AA" and ref='A')
    //   1.0 if Heterozygous (e.g., "AG" and ref='A', alt='G')
    //   2.0 if Homozygous ALT (e.g., "GG" and alt='G')
    //   NAN if missing ("NN", "N", "--") or unrecognized
    static float hapmapToDosage(const std::string& gt, char ref, char alt);

    // Convert standard VCF GT string (e.g., "0/0", "0/1", "1/1", "./.") to dosage
    // Returns number of ALT alleles (0, 1, 2) or NAN for missing
    static float vcfGtToDosage(const std::string& gt);
};

} // namespace cosmic
