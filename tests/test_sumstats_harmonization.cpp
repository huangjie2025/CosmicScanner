#include "association/sumstats/sumstats.h"

#include <cmath>
#include <stdexcept>
#include <string>

using cosmic::association::sumstats::AlleleMatchStatus;
using cosmic::association::sumstats::SumstatsRecord;
using cosmic::association::sumstats::harmonizeAlleles;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void close(double observed, double expected, const std::string& message) {
    require(std::abs(observed - expected) < 1e-12, message);
}

SumstatsRecord record(const std::string& a1, const std::string& a2) {
    SumstatsRecord value;
    value.snp_id = "rs1";
    value.a1 = a1;
    value.a2 = a2;
    value.beta = 0.4;
    value.se = 0.2;
    value.freq = 0.3;
    value.z_score = 999.0;
    return value;
}

void checkAccepted(const std::string& a1, const std::string& a2,
                   AlleleMatchStatus status, double beta, double frequency) {
    const auto result = harmonizeAlleles(record(a1, a2), "A", "C");
    require(result.harmonized && result.status == status, "Unexpected accepted orientation state");
    require(result.record.a1 == "A" && result.record.a2 == "C", "Alleles were not canonicalized");
    close(result.record.beta, beta, "Incorrect harmonized beta");
    close(result.record.freq, frequency, "Incorrect harmonized frequency");
    close(result.record.z_score, beta / 0.2, "Harmonized Z was not recomputed");
    close(result.record.weight, 25.0, "Harmonized weight was not recomputed");
}

} // namespace

int main() {
    checkAccepted("A", "C", AlleleMatchStatus::MATCH, 0.4, 0.3);
    checkAccepted("C", "A", AlleleMatchStatus::SWAPPED, -0.4, 0.7);
    checkAccepted("T", "G", AlleleMatchStatus::STRAND_FLIP, 0.4, 0.3);
    checkAccepted("G", "T", AlleleMatchStatus::STRAND_FLIP_SWAPPED, -0.4, 0.7);
    checkAccepted("a", "c", AlleleMatchStatus::MATCH, 0.4, 0.3);

    for (const auto& alleles : {std::pair{"A", "T"}, std::pair{"T", "A"},
                                std::pair{"C", "G"}, std::pair{"G", "C"}}) {
        const auto ambiguous = harmonizeAlleles(record(alleles.first, alleles.second),
                                                alleles.first, alleles.second);
        require(!ambiguous.harmonized && ambiguous.status == AlleleMatchStatus::AMBIGUOUS,
                "Ambiguous SNP was silently accepted");
    }

    const auto mismatch = harmonizeAlleles(record("A", "G"), "A", "C");
    require(!mismatch.harmonized && mismatch.status == AlleleMatchStatus::MISMATCH,
            "Mismatched SNP was not rejected");
    return 0;
}
