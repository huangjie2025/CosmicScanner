#pragma once

namespace cosmic {
struct GwasConfig;

// Call before loading genotypes or preparing phenotype input. Throws on invalid
// requests; true means standalone completion, false continues the existing engine.
// Only explicitly requested QC/PCA/LD routes materialize an out_prefix.qc BED.
bool prepareGeneticsWorkflow(GwasConfig& cfg);
} // namespace cosmic
