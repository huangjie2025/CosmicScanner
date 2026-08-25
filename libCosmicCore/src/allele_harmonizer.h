#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cctype>

namespace cosmic {

class AlleleHarmonizer {
public:
    // Flips a single allele character (A->T, C->G, G->C, T->A)
    static char complement(char allele);

    // Flips a string of alleles
    static std::string complementString(const std::string& alleles);

    // Checks if the REF/ALT pair is ambiguous (A/T or C/G)
    static bool isAmbiguous(const std::string& ref, const std::string& alt);

    // Harmonize test alleles to reference alleles.
    // Returns an integer vector mapping test allele index to ref allele index.
    // Returns empty vector if harmonization is not possible.
    // Example: ref=["A", "C"], test=["T", "G"] -> returns [0, 1] (after complement)
    //          ref=["A", "C"], test=["C", "A"] -> returns [1, 0] (direct match, swapped)
    static std::vector<int> findHarmonizationMap(const std::vector<std::string>& ref_alleles, 
                                                 const std::vector<std::string>& test_alleles,
                                                 bool& out_used_complement);

    // Remap a genotype dosage value (e.g. 0, 1, 2 for number of ALT alleles)
    // based on the harmonization map.
    // Assuming biallelic: map size is 2. map[0] is new REF idx, map[1] is new ALT idx.
    static float remapDosage(float dosage, const std::vector<int>& index_map);
};

} // namespace cosmic
