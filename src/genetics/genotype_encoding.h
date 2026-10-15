#pragma once
#include <cmath>
#include <stdexcept>
#include <string>
#include <Eigen/Dense>
#include "qc_utils.h"

namespace cosmic { namespace genetics {
// Nonlinear coding is defined on hardcalls, not on expected dosages.
inline double encodeGenotype(double g, const std::string& model) {
    if (model != "additive" && model != "dominant" && model != "recessive" && model != "heterozygote")
        throw std::invalid_argument("Unknown SNP model");
    if (std::isnan(g)) return g;
    if (!std::isfinite(g) || g < 0 || g > 2) throw std::invalid_argument("Invalid genotype");
    if (model == "additive") return g;
    const double hardcall = std::round(g);
    if (std::abs(g - hardcall) > 1e-6)
        throw std::invalid_argument("Non-additive encoding requires hardcalls");
    if (model == "dominant") return hardcall > 0 ? 1.0 : 0.0;
    if (model == "recessive") return hardcall == 2 ? 1.0 : 0.0;
    return hardcall == 1 ? 1.0 : 0.0;
}

struct PreparedMarker {
    Eigen::VectorXd values;
    double allele_frequency = 0.0;
    int observed = 0;
    bool hardcalls = true;
    bool passes_maf = false;
    bool passes_hwe = true;
};

inline PreparedMarker prepareAssociationMarker(const Eigen::VectorXd& raw,
                                               const std::string& model,
                                               double min_maf,
                                               double min_hwe,
                                               bool exact_hwe_available = true,
                                               bool center = true) {
    PreparedMarker out;
    out.values.resize(raw.size());
    double dosage_sum = 0.0;
    long long hom0 = 0, het = 0, hom2 = 0;
    for (int i = 0; i < raw.size(); ++i) {
        const double value = raw[i];
        if (std::isnan(value)) continue;
        if (!std::isfinite(value) || value < 0.0 || value > 2.0)
            throw std::invalid_argument("Invalid genotype");
        dosage_sum += value;
        ++out.observed;
        const double rounded = std::round(value);
        if (std::abs(value - rounded) > 1e-6) {
            out.hardcalls = false;
        } else if (rounded == 0.0) {
            ++hom0;
        } else if (rounded == 1.0) {
            ++het;
        } else {
            ++hom2;
        }
    }
    if (out.observed == 0) return out;
    out.allele_frequency = dosage_sum / (2.0 * out.observed);
    out.passes_maf = out.observed >= 3 &&
        std::min(out.allele_frequency, 1.0 - out.allele_frequency) >= min_maf;
    if (min_hwe > 0.0 && exact_hwe_available && out.hardcalls)
        out.passes_hwe = calculateHWE(hom0, het, hom2) >= min_hwe;
    if (model != "additive" && !out.hardcalls)
        throw std::invalid_argument("Non-additive encoding requires hardcalls");

    double encoded_sum = 0.0;
    for (int i = 0; i < raw.size(); ++i) {
        out.values[i] = encodeGenotype(raw[i], model);
        if (std::isfinite(out.values[i])) encoded_sum += out.values[i];
    }
    const double encoded_mean = encoded_sum / out.observed;
    for (int i = 0; i < out.values.size(); ++i) {
        if (!std::isfinite(out.values[i])) out.values[i] = encoded_mean;
        if (center) out.values[i] -= encoded_mean;
    }
    return out;
}
} }
