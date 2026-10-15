#include "association/multilocus/multilocus_engine.h"

namespace cosmic::association {

MultiLocusOptions blinkOptions(int max_pseudo_qtn) {
    MultiLocusOptions options;
    options.method = MultiLocusMethod::BLINK;
    options.max_iterations = 5;
    options.max_pseudo_qtn = max_pseudo_qtn;
    options.ld_r2 = 0.2;
    return options;
}

} // namespace cosmic::association
