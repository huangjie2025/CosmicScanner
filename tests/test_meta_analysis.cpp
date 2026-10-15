#include "association/meta/meta_engine.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cosmic::association::meta;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void close(double observed, double expected, double tolerance, const std::string& message) {
    require(std::isfinite(observed) && std::abs(observed - expected) <= tolerance,
            message + ": observed=" + std::to_string(observed) + " expected=" + std::to_string(expected));
}

StudySummary summary(const std::string& id, double beta, double se,
                     const std::string& a1 = "A", const std::string& a2 = "C",
                     double frequency = 0.3) {
    StudySummary value;
    value.snp_id = id;
    value.chrom = "1";
    value.position = id == "late" ? 2000 : 1000;
    value.a1 = a1;
    value.a2 = a2;
    value.beta = beta;
    value.se = se;
    value.p_value = std::erfc(std::abs(beta / se) / std::sqrt(2.0));
    value.n = 1000;
    value.freq = frequency;
    return value;
}

MetaResult runHeterogeneous(MetaAnalysisMethod method) {
    MetaAnalysisOptions options;
    options.method = method;
    MetaAnalysisEngine engine(options);
    const std::vector<double> betas = {0.2, 0.8, -0.1, 1.1};
    const std::vector<double> ses = {0.1, 0.2, 0.15, 0.25};
    for (int i = 0; i < 4; ++i)
        engine.addStudy("study" + std::to_string(i + 1), {summary("heterogeneous", betas[i], ses[i])});
    const auto results = engine.run();
    require(results.size() == 1, "Heterogeneous meta-analysis did not return one result");
    return results.front();
}

void checkNumericalOracles() {
    const auto fixed = runHeterogeneous(MetaAnalysisMethod::FIXED_EFFECT);
    close(fixed.beta, 0.28663870581186346, 1e-12, "FE beta");
    close(fixed.se, 0.0734333069472097, 1e-12, "FE SE");
    close(fixed.p_value, 9.485536611702767e-05, 1e-15, "FE P");
    close(fixed.heterogeneity.q_statistic, 24.568004793289397, 1e-11, "Q");
    close(fixed.heterogeneity.q_df, 3.0, 0.0, "Q df");
    close(fixed.heterogeneity.q_p_value, 1.9008881099771135e-05, 1e-14, "Q P");
    close(fixed.heterogeneity.i2, 87.78899619549313, 1e-11, "I2");
    close(fixed.heterogeneity.tau2, 0.18574303405572756, 1e-11, "DL tau2 reported by FE");

    const auto dl = runHeterogeneous(MetaAnalysisMethod::RANDOM_EFFECT_DL);
    close(dl.heterogeneity.tau2, 0.18574303405572756, 1e-11, "DL tau2");
    close(dl.beta, 0.46363756979265613, 1e-11, "DL beta");
    close(dl.se, 0.23332317709249825, 1e-11, "DL SE");
    close(dl.p_value, 0.04691080123260502, 1e-12, "DL P");

    const auto reml = runHeterogeneous(MetaAnalysisMethod::RANDOM_EFFECT_REML);
    close(reml.heterogeneity.tau2, 0.255847204023964, 1e-9, "REML tau2");
    close(reml.beta, 0.4722333464454319, 1e-9, "REML beta");
    close(reml.se, 0.2684571774130413, 1e-9, "REML SE");
    close(reml.p_value, 0.07856664494357501, 1e-9, "REML P");
}

void checkHarmonizationAndOrdering() {
    MetaAnalysisEngine engine(MetaAnalysisOptions{});
    engine.addStudy("z-first", {
        summary("orientation", 0.2, 0.1),
        summary("ambiguous", 0.3, 0.1, "A", "T"),
        summary("mismatch", 0.4, 0.1)});
    engine.addStudy("a-second", {
        summary("orientation", -0.2, 0.1, "C", "A", 0.7),
        summary("ambiguous", 0.3, 0.1, "A", "T"),
        summary("mismatch", 0.4, 0.1, "A", "G"),
        summary("late", 0.5, 0.2)});
    engine.addStudy("m-third", {
        summary("orientation", 0.2, 0.1, "T", "G"),
        summary("ambiguous", 0.3, 0.1, "T", "A"),
        summary("late", 0.5, 0.2)});
    engine.addStudy("b-fourth", {
        summary("orientation", -0.2, 0.1, "G", "T", 0.7),
        summary("late", 0.5, 0.2)});

    const auto results = engine.run();
    require(results.size() == 2, "Ambiguous/mismatched variants were not rejected");
    require(results[0].snp_id == "orientation" && results[1].snp_id == "late",
            "Meta result order did not preserve deterministic first occurrence");
    require(results[0].a1 == "A" && results[0].a2 == "C",
            "Reference alleles did not come from the first inserted study");
    close(results[0].beta, 0.2, 1e-14, "Harmonized identical effect");
    close(results[0].freq, 0.3, 1e-14, "Harmonized frequency");
    close(results[0].heterogeneity.q_statistic, 0.0, 1e-14, "Identical-effect Q");
    close(results[1].beta, 0.5, 1e-14, "SNP absent from first study");

    const auto repeated = engine.run();
    require(repeated.size() == results.size(), "Repeated deterministic run changed result count");
    for (size_t i = 0; i < results.size(); ++i) {
        require(repeated[i].snp_id == results[i].snp_id, "Repeated run changed result order");
        close(repeated[i].beta, results[i].beta, 0.0, "Repeated run changed beta");
    }
}

} // namespace

int main() {
    checkNumericalOracles();
    checkHarmonizationAndOrdering();
    return 0;
}
