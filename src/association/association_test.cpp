#include "association_test.h"
#include "association/multilocus/multilocus_engine.h"
#include "association/generalized/firth_logistic.h"
#include "association/generalized/saige_style.h"
#include "association/generalized/logistic.h"

#include <stdexcept>

namespace cosmic {

std::vector<AssociationResult> AssociationTest::scanBlock(const GenotypeBlock&) {
    throw std::runtime_error("AssociationTest: block scanning is unsupported");
}

AssociationTestPtr AssociationTestFactory::create(AssociationModel model) {
    if (model == AssociationModel::FARMCPU)
        return std::make_unique<association::MultiLocusEngine>(association::farmCpuOptions(5, 10));
    if (model == AssociationModel::BLINK)
        return std::make_unique<association::MultiLocusEngine>(association::blinkOptions());
    if (model == AssociationModel::FIRTH_LOGISTIC)
        return std::make_unique<association::FirthLogisticTest>();
    if (model == AssociationModel::SAIGE_STYLE)
        return std::make_unique<association::SaigeStyleTest>();
    if (model == AssociationModel::LOGISTIC)
        return std::make_unique<association::LogisticTest>();
    throw std::runtime_error(
        "AssociationTestFactory: unsupported model " + std::to_string(static_cast<int>(model)) +
        "; no numerical adapter is registered for this model");
}

} // namespace cosmic
