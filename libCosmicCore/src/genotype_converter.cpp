#include "genotype_converter.h"

namespace cosmic {

float GenotypeConverter::hapmapToDosage(const std::string& gt, char ref, char alt) {
    if (gt.empty() || gt == "NN" || gt == "N" || gt == "--") {
        return NAN;
    }

    if (gt.length() != 2) {
        return NAN;
    }

    char a1 = gt[0];
    char a2 = gt[1];

    if (a1 == ref && a2 == ref) {
        return 0.0f; // Homozygous REF
    }
    if (a1 == alt && a2 == alt) {
        return 2.0f; // Homozygous ALT
    }
    if ((a1 == ref && a2 == alt) || (a1 == alt && a2 == ref)) {
        return 1.0f; // Heterozygous
    }

    // Any other ambiguity or multiallelic state
    return NAN;
}

float GenotypeConverter::vcfGtToDosage(const std::string& gt) {
    if (gt.empty() || gt == "." || gt == "./." || gt == ".|.") {
        return NAN;
    }

    // Simple parser for standard "0/0", "0/1", "1/1", "0|1", etc.
    float dosage = 0.0f;
    bool valid = false;

    for (char c : gt) {
        if (c >= '0' && c <= '9') {
            // In a biallelic context, we count the number of non-zero alleles (ALT)
            // Typically '1' is ALT. If it's a multiallelic '2', it still counts as an ALT allele 
            // relative to REF. Here we just sum the values if they are 1.
            int allele_idx = c - '0';
            if (allele_idx > 0) {
                dosage += 1.0f;
            }
            valid = true;
        } else if (c == '.') {
            return NAN; // If any allele is missing, treat whole genotype as missing
        }
    }

    return valid ? dosage : NAN;
}

} // namespace cosmic
