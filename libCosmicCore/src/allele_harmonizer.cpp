#include "allele_harmonizer.h"

namespace cosmic {

char AlleleHarmonizer::complement(char allele) {
    switch (std::toupper(allele)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return allele;
    }
}

std::string AlleleHarmonizer::complementString(const std::string& alleles) {
    std::string result = alleles;
    for (char& c : result) {
        c = complement(c);
    }
    return result;
}

bool AlleleHarmonizer::isAmbiguous(const std::string& ref, const std::string& alt) {
    if (ref.length() != 1 || alt.length() != 1) return false;
    char r = std::toupper(ref[0]);
    char a = std::toupper(alt[0]);
    return (r == 'A' && a == 'T') || (r == 'T' && a == 'A') ||
           (r == 'C' && a == 'G') || (r == 'G' && a == 'C');
}

std::vector<int> AlleleHarmonizer::findHarmonizationMap(const std::vector<std::string>& ref_alleles, 
                                             const std::vector<std::string>& test_alleles,
                                             bool& out_used_complement) {
    out_used_complement = false;
    if (ref_alleles.empty() || test_alleles.empty()) return {};

    auto build_map = [&](const std::vector<std::string>& src, const std::vector<std::string>& tgt) -> std::vector<int> {
        std::unordered_map<std::string, int> tgt_idx;
        for (size_t i = 0; i < tgt.size(); ++i) {
            tgt_idx[tgt[i]] = i;
        }
        std::vector<int> map(src.size(), -1);
        for (size_t i = 0; i < src.size(); ++i) {
            auto it = tgt_idx.find(src[i]);
            if (it == tgt_idx.end()) return {}; // Not all src alleles found in tgt
            map[i] = it->second;
        }
        return map;
    };

    // Try direct mapping first
    std::vector<int> map = build_map(test_alleles, ref_alleles);
    if (!map.empty()) return map;

    // Try complement mapping
    std::vector<std::string> test_alleles_comp;
    for (const auto& a : test_alleles) {
        test_alleles_comp.push_back(complementString(a));
    }
    map = build_map(test_alleles_comp, ref_alleles);
    if (!map.empty()) {
        out_used_complement = true;
        return map;
    }

    return {};
}

float AlleleHarmonizer::remapDosage(float dosage, const std::vector<int>& index_map) {
    // If missing or not biallelic index map, return as is (or missing)
    if (dosage != dosage || index_map.size() != 2) return dosage;

    // Assuming dosage is number of ALT alleles (test_alleles[1])
    // index_map[0] -> where the old REF points in the new alleles
    // index_map[1] -> where the old ALT points in the new alleles

    // If index_map is [0, 1] -> no change
    if (index_map[0] == 0 && index_map[1] == 1) {
        return dosage;
    }
    // If index_map is [1, 0] -> alleles are swapped (ALT becomes REF)
    else if (index_map[0] == 1 && index_map[1] == 0) {
        // Dosage of ALT becomes (2 - dosage)
        return 2.0f - dosage;
    }
    
    return dosage; // Default fallback
}

} // namespace cosmic
