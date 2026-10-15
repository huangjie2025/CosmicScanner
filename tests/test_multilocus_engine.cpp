#include "association/multilocus/multilocus_engine.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cosmic;
using namespace cosmic::association;

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

int minimumP(const std::vector<MultiLocusResult>& results) {
    int best = 0;
    for (int j = 1; j < static_cast<int>(results.size()); ++j)
        if (results[j].p_value < results[best].p_value) best = j;
    return best;
}
} // namespace

int main() {
    constexpr int n = 80, m = 6;
    Eigen::MatrixXd genotypes(n, m);
    Eigen::MatrixXd phenotype(n, 1);
    Eigen::MatrixXd covariates = Eigen::MatrixXd::Ones(n, 1);
    for (int i = 0; i < n; ++i) {
        genotypes(i, 0) = i % 3;
        genotypes(i, 1) = i % 3;
        genotypes(i, 2) = (i / 3) % 3;
        genotypes(i, 3) = (i / 7) % 3;
        genotypes(i, 4) = (i * 5 + 1) % 3;
        genotypes(i, 5) = (i * 7 + 2) % 3;
        phenotype(i, 0) = 1.5 + 1.8 * genotypes(i, 0) + 0.7 * genotypes(i, 3) +
                          0.15 * std::sin(static_cast<double>(i));
    }
    const std::vector<Eigen::Index> identity = [] {
        std::vector<Eigen::Index> values(n);
        for (int i = 0; i < n; ++i) values[i] = i;
        return values;
    }();
    const ScanContext context{phenotype, covariates, identity};
    GenotypeBlock block{genotypes, {"lead", "lead_proxy", "v2", "secondary", "v4", "v5"}};

    for (auto options : {farmCpuOptions(5, 3), blinkOptions(3)}) {
        MultiLocusEngine engine(options);
        engine.prepare(context);
        const auto generic = engine.scanBlock(block);
        check(generic.size() == m && engine.results().size() == m, "Missing multi-locus results");
        check(minimumP(engine.results()) == 0 || minimumP(engine.results()) == 1,
              "Multi-locus trend missed the causal lead/proxy");
        check(!engine.selectedIndices().empty(), "BIC did not select a pseudo-QTN");
        for (const auto& result : engine.results())
            check(std::isfinite(result.p_value) && result.p_value >= 0 && result.p_value <= 1,
                  "Invalid multi-locus P-value");
    }

    MultiLocusEngine conditional(conditionalOptions({"lead"}));
    conditional.prepare(context);
    const auto generic = conditional.scanBlock(block);
    check(generic[1].status == AssociationStatus::FILTERED,
          "Perfect lead proxy must be collinear after conditioning");
    check(conditional.results()[3].p_value < 1e-8,
          "Conditional scan did not retain the independent secondary signal");
    check(conditional.selectedIndices().size() == 1 && conditional.selectedIndices()[0] == 0,
          "Conditional lead resolution failed");

    MultiLocusEngine joint(conditionalJointOptions({"lead", "lead_proxy", "secondary"}));
    joint.prepare(context);
    joint.scanBlock(block);
    check(joint.selectedIndices().size() == 2 && joint.selectedIndices()[0] == 0 &&
          joint.selectedIndices()[1] == 3,
          "Conditional joint LD/BIC selection did not recover independent leads");
    check(joint.results()[0].selected_pseudo_qtn && joint.results()[3].selected_pseudo_qtn &&
          std::isfinite(joint.results()[0].p_value) && std::isfinite(joint.results()[3].p_value),
          "Conditional joint selected-marker estimates are missing");

    // Independent deterministic polygenic simulation, separate from production selection code.
    constexpr int ns = 240, ms = 40;
    Eigen::MatrixXd simulated(ns, ms);
    Eigen::MatrixXd simulated_y(ns, 1);
    for (int i = 0; i < ns; ++i) {
        for (int j = 0; j < ms; ++j)
            simulated(i, j) = (i * (2 * j + 3) + i / (j + 2) + j * j) % 3;
        simulated(i, 8) = simulated(i, 7);
        simulated_y(i, 0) = 1.4 * simulated(i, 7) - 1.1 * simulated(i, 29) +
                            0.25 * std::sin(0.37 * i) + 0.15 * std::cos(0.11 * i);
    }
    std::vector<std::string> ids(ms);
    for (int j = 0; j < ms; ++j) ids[j] = "sim" + std::to_string(j);
    std::vector<Eigen::Index> simulated_identity(ns);
    for (int i = 0; i < ns; ++i) simulated_identity[i] = i;
    const ScanContext simulated_context{
        simulated_y, Eigen::MatrixXd::Ones(ns, 1), simulated_identity};
    const GenotypeBlock simulated_block{simulated, ids};
    for (auto options : {farmCpuOptions(8, 8), blinkOptions(8)}) {
        MultiLocusEngine engine(options);
        engine.prepare(simulated_context);
        engine.scanBlock(simulated_block);
        std::vector<int> order(ms);
        for (int j = 0; j < ms; ++j) order[j] = j;
        std::sort(order.begin(), order.end(), [&](int left, int right) {
            return engine.results()[left].p_value < engine.results()[right].p_value;
        });
        bool lead_recovered = false, secondary_recovered = false;
        for (int k = 0; k < 10; ++k) {
            lead_recovered = lead_recovered || order[k] == 7 || order[k] == 8;
            secondary_recovered = secondary_recovered || order[k] == 29;
        }
        check(lead_recovered && secondary_recovered,
              "FarmCPU/BLINK-style scan failed deterministic two-signal recovery");
    }
}
