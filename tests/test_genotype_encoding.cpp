#include "genetics/genotype_encoding.h"

#include <cmath>
#include <stdexcept>

using cosmic::genetics::prepareAssociationMarker;

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    Eigen::VectorXd raw(7);
    raw << 0.0, 1.0, 2.0, 0.0, 1.0, 2.0,
        std::numeric_limits<double>::quiet_NaN();
    for (const auto& model : {"additive", "dominant", "recessive", "heterozygote"}) {
        const auto marker = prepareAssociationMarker(raw, model, 0.0, 0.0);
        check(marker.observed == 6 && marker.hardcalls && marker.passes_maf,
              "Hard-call marker preparation failed");
        check(std::abs(marker.allele_frequency - 0.5) < 1e-12 &&
              std::abs(marker.values.mean()) < 1e-12,
              "Frequency, imputation, or centering is incorrect");
    }

    Eigen::VectorXd dosage(4);
    dosage << 0.1, 0.9, 1.4, 1.8;
    const auto additive = prepareAssociationMarker(dosage, "additive", 0.0, 0.5, false);
    check(!additive.hardcalls && additive.passes_hwe,
          "Continuous dosage HWE must be explicitly skipped");
    bool rejected = false;
    try {
        (void)prepareAssociationMarker(dosage, "dominant", 0.0, 0.0);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, "Non-additive coding accepted continuous dosage");
    Eigen::VectorXd near_integer(3);
    near_integer << 0.0, 1.0000005, 2.0;
    const auto tolerant = prepareAssociationMarker(near_integer, "heterozygote", 0.0, 0.0);
    check(tolerant.hardcalls, "Near-integer hardcall tolerance is inconsistent");

    Eigen::VectorXd disequilibrium = Eigen::VectorXd::Zero(40);
    disequilibrium.tail(20).setConstant(2.0);
    const auto hwe = prepareAssociationMarker(disequilibrium, "additive", 0.0, 0.05);
    check(!hwe.passes_hwe, "Exact HWE filter did not reject a disequilibrium marker");
}
