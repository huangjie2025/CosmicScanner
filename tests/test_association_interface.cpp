#include "../src/association/association_test.h"
#include "../src/association/multilocus/multilocus_engine.h"

#include <cmath>
#include <stdexcept>
#include <type_traits>

using namespace cosmic;

static_assert(std::is_abstract<AssociationTest>::value, "Interface must remain abstract");
static_assert(std::has_virtual_destructor<AssociationTest>::value, "Polymorphic destruction");
static_assert(std::is_same<decltype(ScanContext::phenotypes), const Eigen::MatrixXd&>::value,
              "Phenotypes must be borrowed matrices");
static_assert(std::is_same<decltype(ScanContext::covariates), const Eigen::MatrixXd&>::value,
              "Covariates must be borrowed");
static_assert(std::is_same<decltype(ScanContext::observation_to_sample),
                           const std::vector<Eigen::Index>&>::value,
              "Observation mapping must be borrowed");

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

// Compile contract only, not a numerical adapter or an advertised implementation.
class ContractTest final : public AssociationTest {
public:
    AssociationModel getModelType() const override { return AssociationModel::FARMCPU; }
    AssociationCapabilities capabilities() const override { return {}; }
    void prepare(const ScanContext&) override {}
};

template <typename Function>
void rejectsUnsupported(Function function) {
    try {
        function();
    } catch (const std::runtime_error& error) {
        check(std::string(error.what()).find("unsupported") != std::string::npos,
              "Rejection must explicitly identify unsupported operation");
        return;
    }
    throw std::runtime_error("Unsupported operation returned instead of throwing");
}
} // namespace

int main() {
    Eigen::MatrixXd phenotypes = Eigen::MatrixXd::Zero(3, 2);
    Eigen::MatrixXd covariates = Eigen::MatrixXd::Ones(3, 1);
    const std::vector<Eigen::Index> mapping{1, 0, 1};
    const ScanContext context{phenotypes, covariates, mapping};
    check(&context.phenotypes == &phenotypes && &context.covariates == &covariates &&
          &context.observation_to_sample == &mapping, "Context must not copy inputs");
    check(context.phenotypes.cols() == 2 && context.observation_to_sample[2] == 1,
          "Matrix phenotypes and repeated observations must be representable");

    AssociationResult result;
    check(result.status == AssociationStatus::NOT_TESTED && std::isnan(result.p_value) &&
          std::isnan(result.statistic) && std::isnan(result.degrees_of_freedom),
          "Default result must not imply successful inference");
    result.effects.resize(2);
    result.effect_covariance.resize(2, 2);
    result.degrees_of_freedom = 2;
    check(result.effects.size() == 2 && result.effect_covariance.cols() == 2 &&
          result.degrees_of_freedom == 2, "Generic multi-df result contract");

    ContractTest test;
    AssociationTest& interface = test;
    interface.prepare(context);
    const auto capabilities = interface.capabilities();
    check(!capabilities.block_scan && !capabilities.iterative_scan && !capabilities.set_test,
          "Default capabilities must not promise any execution workflow");
    AssociationCapabilities iterative{false, true, false};
    AssociationCapabilities set{false, false, true};
    check(iterative.iterative_scan && !iterative.block_scan && set.set_test && !set.block_scan,
          "Iterative and set capabilities must be distinct from block scanning");
    rejectsUnsupported([&] { interface.scanBlock(GenotypeBlock{}); });

    for (const auto model : {
             AssociationModel::LMM, AssociationModel::MVLMM, AssociationModel::RRM,
             AssociationModel::REPEATABILITY, AssociationModel::MATERNAL,
             AssociationModel::SINGLE_STEP, AssociationModel::ANIMAL,
             AssociationModel::THRESHOLD, AssociationModel::COUNT, AssociationModel::ORDINAL,
             AssociationModel::SURVIVAL, AssociationModel::GXE, AssociationModel::COMMON_ENV,
             AssociationModel::SKAT, AssociationModel::BURDEN, AssociationModel::CONDITIONAL,
             static_cast<AssociationModel>(-1), static_cast<AssociationModel>(999)}) {
        rejectsUnsupported([&] { AssociationTestFactory::create(model); });
    }
    const auto farmcpu = AssociationTestFactory::create(AssociationModel::FARMCPU);
    const auto blink = AssociationTestFactory::create(AssociationModel::BLINK);
    const auto firth = AssociationTestFactory::create(AssociationModel::FIRTH_LOGISTIC);
    const auto saige = AssociationTestFactory::create(AssociationModel::SAIGE_STYLE);
    check(farmcpu->getModelType() == AssociationModel::FARMCPU &&
          blink->getModelType() == AssociationModel::BLINK &&
          firth->getModelType() == AssociationModel::FIRTH_LOGISTIC &&
          saige->getModelType() == AssociationModel::SAIGE_STYLE,
          "Factory must expose native association adapters");
}
