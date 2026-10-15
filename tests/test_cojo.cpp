#include "association/cojo/cojo_engine.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

using namespace cosmic::association::cojo;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void close(double observed, double expected, double tolerance, const std::string& message) {
    require(std::isfinite(observed) && std::abs(observed - expected) <= tolerance,
            message + ": observed=" + std::to_string(observed) + " expected=" + std::to_string(expected));
}

SummaryRecord summary(const std::string& id, long long position, double z, double se) {
    SummaryRecord value;
    value.snp_id = id;
    value.chrom = "1";
    value.position = position;
    value.a1 = "A";
    value.a2 = "C";
    value.beta = z * se;
    value.se = se;
    value.p_value = std::erfc(std::abs(z) / std::sqrt(2.0));
    value.freq = 0.3;
    value.n = 1000;
    value.z_score = z;
    value.weight = 1.0 / (se * se);
    return value;
}

LDReference reference(const std::vector<std::string>& ids, const Eigen::MatrixXd& matrix) {
    LDReference value;
    value.snp_ids = ids;
    value.ld_matrix = matrix;
    for (size_t i = 0; i < ids.size(); ++i) {
        value.chroms.push_back("1");
        value.positions.push_back(1000 + static_cast<long long>(i) * 1000);
        value.a1s.push_back("A");
        value.a2s.push_back("C");
        value.snp_to_idx[ids[i]] = static_cast<int>(i);
    }
    return value;
}

std::unordered_map<std::string, CoJoResult> byId(const std::vector<CoJoResult>& rows) {
    std::unordered_map<std::string, CoJoResult> result;
    for (const auto& row : rows) result.emplace(row.snp_id, row);
    return result;
}

void checkMappingConditionalAndJoint(const LDReference& ld) {
    const std::vector<SummaryRecord> summaries = {
        summary("D", 4000, 1.5, 0.5), summary("A", 1000, 2.0, 0.2),
        summary("C", 3000, -1.0, 0.4)};

    CoJoOptions conditional_options;
    conditional_options.mode = CoJoMode::CONDITIONAL;
    conditional_options.maf_filter = 0.0;
    CoJoEngine conditional(conditional_options);
    conditional.setSummaryStatsForTesting(summaries);
    conditional.setLDReferenceForTesting(ld);
    conditional.setConditionSnpsForTesting({"A", "C"});
    const auto conditional_rows = byId(conditional.run());

    Eigen::Matrix3d selected_ld;
    selected_ld << 1.0, 0.4, -0.1,
                   0.4, 1.0, 0.25,
                  -0.1, 0.25, 1.0;
    Eigen::Vector3d z(2.0, -1.0, 1.5);
    const Eigen::Matrix3d inverse = selected_ld.inverse();
    const Eigen::Vector3d solved = inverse * z;
    const auto& candidate = conditional_rows.at("D");
    close(candidate.conditional_beta, 0.5 * solved[2], 1e-11, "Conditional candidate beta oracle");
    close(candidate.conditional_se, 0.5 * std::sqrt(inverse(2, 2)), 1e-11,
          "Conditional candidate SE oracle");

    CoJoOptions joint_options = conditional_options;
    joint_options.mode = CoJoMode::JOINT;
    CoJoEngine joint(joint_options);
    joint.setSummaryStatsForTesting(summaries);
    joint.setLDReferenceForTesting(ld);
    joint.setConditionSnpsForTesting({"A", "C", "D"});
    const auto joint_rows = byId(joint.run());
    const std::vector<std::string> ids = {"A", "C", "D"};
    const std::vector<double> ses = {0.2, 0.4, 0.5};
    for (int i = 0; i < 3; ++i) {
        const auto& row = joint_rows.at(ids[i]);
        require(row.selected, "Joint SNP was not selected");
        close(row.conditional_beta, ses[i] * solved[i], 1e-11, "Joint beta oracle");
        close(row.conditional_se, ses[i] * std::sqrt(inverse(i, i)), 1e-11, "Joint SE oracle");
    }
}

void checkCausalProxyAndIndependentSignals() {
    Eigen::Matrix3d matrix;
    matrix << 1.0, 0.95, 0.0,
              0.95, 1.0, 0.0,
              0.0, 0.0, 1.0;
    const auto ld = reference({"causal", "proxy", "independent"}, matrix);
    const Eigen::Vector3d truth(5.0, 0.0, 4.0);
    const Eigen::Vector3d marginal = matrix * truth;
    const std::vector<SummaryRecord> summaries = {
        summary("causal", 1000, marginal[0], 0.1),
        summary("proxy", 2000, marginal[1], 0.1),
        summary("independent", 3000, marginal[2], 0.1)};

    CoJoOptions conditional_options;
    conditional_options.mode = CoJoMode::CONDITIONAL;
    conditional_options.maf_filter = 0.0;
    CoJoEngine conditional(conditional_options);
    conditional.setSummaryStatsForTesting(summaries);
    conditional.setLDReferenceForTesting(ld);
    conditional.setConditionSnpsForTesting({"causal"});
    const auto conditioned = byId(conditional.run());
    require(conditioned.at("proxy").conditional_p > 0.99, "High-LD proxy was not attenuated");
    require(conditioned.at("causal").selected, "Causal condition SNP was not retained");

    CoJoOptions stepwise_options;
    stepwise_options.mode = CoJoMode::STEPWISE;
    stepwise_options.threshold = 1e-4;
    stepwise_options.collinearity = 0.99;
    stepwise_options.maf_filter = 0.0;
    CoJoEngine stepwise(stepwise_options);
    stepwise.setSummaryStatsForTesting(summaries);
    stepwise.setLDReferenceForTesting(ld);
    const auto selected = byId(stepwise.run());
    require(selected.at("causal").selected && selected.at("independent").selected,
            "Two independent causal SNPs were not retained");
    require(!selected.at("proxy").selected, "High-LD null proxy was selected after conditioning");
}

void checkExternalReference(const LDReference& expected) {
    const auto path = std::filesystem::temp_directory_path() / "cosmicscanner_cojo_dense_ld.tsv";
    {
        std::ofstream out(path);
        out << "SNP CHR BP A1 A2 A B C D\n";
        for (int i = 0; i < 4; ++i) {
            out << expected.snp_ids[i] << " 1 " << expected.positions[i] << " A C";
            for (int j = 0; j < 4; ++j) out << ' ' << expected.ld_matrix(i, j);
            out << '\n';
        }
    }
    CoJoEngine loaded(CoJoOptions{});
    loaded.loadLDReference(path.string());
    const auto& observed = loaded.getLDReference();
    require(observed.snp_ids == expected.snp_ids, "External LD SNP order changed");
    require((observed.ld_matrix - expected.ld_matrix).cwiseAbs().maxCoeff() < 1e-12,
            "External LD matrix values changed");
    checkMappingConditionalAndJoint(observed);
    std::filesystem::remove(path);
}

} // namespace

int main() {
    CoJoOptions defaults;
    require(defaults.mode == CoJoMode::CONDITIONAL && defaults.threshold == 5e-8,
            "COJO defaults changed");

    Eigen::Matrix4d matrix;
    matrix << 1.0, 0.2, 0.4, -0.1,
              0.2, 1.0, 0.3, 0.05,
              0.4, 0.3, 1.0, 0.25,
             -0.1, 0.05, 0.25, 1.0;
    const auto ld = reference({"A", "B", "C", "D"}, matrix);
    checkMappingConditionalAndJoint(ld);
    checkExternalReference(ld);
    checkCausalProxyAndIndependentSignals();
    return 0;
}
