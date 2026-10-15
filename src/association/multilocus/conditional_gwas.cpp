#include "association/multilocus/multilocus_engine.h"

#include <utility>

namespace cosmic::association {

MultiLocusOptions conditionalOptions(std::vector<std::string> condition_ids) {
    MultiLocusOptions options;
    options.method = MultiLocusMethod::CONDITIONAL;
    options.max_iterations = 1;
    options.max_pseudo_qtn = static_cast<int>(condition_ids.size());
    options.condition_ids = std::move(condition_ids);
    return options;
}

MultiLocusOptions conditionalJointOptions(std::vector<std::string> condition_ids) {
    MultiLocusOptions options = conditionalOptions(std::move(condition_ids));
    options.method = MultiLocusMethod::CONDITIONAL_JOINT;
    return options;
}

} // namespace cosmic::association
