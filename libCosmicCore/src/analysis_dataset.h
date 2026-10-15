#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <cmath>
#include <sstream>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cosmic {

inline constexpr std::size_t kAnalysisNpos = std::numeric_limits<std::size_t>::max();

struct SampleKey {
    std::string fid;
    std::string iid;

    bool operator==(const SampleKey& other) const {
        return fid == other.fid && iid == other.iid;
    }
    bool operator<(const SampleKey& other) const {
        if (fid != other.fid) return fid < other.fid;
        return iid < other.iid;
    }
    std::string display() const {
        if (!fid.empty() && !iid.empty()) return fid + ":" + iid;
        if (!iid.empty()) return iid;
        return fid;
    }
};

struct SampleKeyHash {
    std::size_t operator()(const SampleKey& key) const {
        const std::size_t h1 = std::hash<std::string>{}(key.fid);
        const std::size_t h2 = std::hash<std::string>{}(key.iid);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

inline SampleKey make_iid_sample_key(const std::string& id) {
    return SampleKey{"", id};
}

inline std::vector<SampleKey> sample_keys_from_fid_iid(const std::vector<std::string>& fids,
                                                       const std::vector<std::string>& iids) {
    if (fids.size() != iids.size()) {
        throw std::runtime_error("FID/IID size mismatch while building sample keys.");
    }
    std::vector<SampleKey> keys;
    keys.reserve(iids.size());
    for (std::size_t i = 0; i < iids.size(); ++i) {
        keys.push_back(SampleKey{fids[i], iids[i]});
    }
    return keys;
}

struct AnalysisSample {
    SampleKey key;
    std::size_t analysis_index = kAnalysisNpos;
    std::size_t genotype_index = kAnalysisNpos;
    std::size_t phenotype_row_index = kAnalysisNpos;
    std::size_t covariate_row_index = kAnalysisNpos;
    std::size_t pedigree_index = kAnalysisNpos;
};

struct SampleAlignmentDiagnostics {
    std::size_t genotype_samples = 0;
    std::size_t phenotype_rows = 0;
    std::size_t covariate_rows = 0;
    std::size_t pedigree_rows = 0;
    std::size_t analysis_samples = 0;
    std::size_t matched_phenotypes = 0;
    std::size_t matched_covariates = 0;
    std::size_t matched_pedigree = 0;
    std::size_t duplicate_genotype_ids = 0;
    std::size_t duplicate_phenotype_ids = 0;
    std::size_t duplicate_covariate_ids = 0;
    std::size_t duplicate_pedigree_ids = 0;
    std::size_t genotype_without_phenotype = 0;
    std::size_t genotype_without_covariate = 0;
    std::size_t genotype_without_pedigree = 0;
    std::size_t phenotype_without_genotype = 0;
    std::size_t covariate_without_genotype = 0;
    std::size_t pedigree_without_genotype = 0;
};

class SampleRegistry {
public:
    static SampleRegistry fromGenotypeOrder(const std::vector<SampleKey>& genotype_samples,
                                            const std::vector<SampleKey>& phenotype_samples = {},
                                            const std::vector<SampleKey>& covariate_samples = {},
                                            const std::vector<SampleKey>& pedigree_samples = {}) {
        SampleRegistry registry;
        registry.diagnostics_.genotype_samples = genotype_samples.size();
        registry.diagnostics_.phenotype_rows = phenotype_samples.size();
        registry.diagnostics_.covariate_rows = covariate_samples.size();
        registry.diagnostics_.pedigree_rows = pedigree_samples.size();

        auto phenotype_map = first_index_map(phenotype_samples, registry.diagnostics_.duplicate_phenotype_ids);
        auto covariate_map = first_index_map(covariate_samples, registry.diagnostics_.duplicate_covariate_ids);
        auto pedigree_map = first_index_map(pedigree_samples, registry.diagnostics_.duplicate_pedigree_ids);

        std::unordered_map<SampleKey, std::size_t, SampleKeyHash> genotype_seen;
        registry.samples_.reserve(genotype_samples.size());
        for (std::size_t i = 0; i < genotype_samples.size(); ++i) {
            const auto& key = genotype_samples[i];
            if (!genotype_seen.emplace(key, i).second) {
                ++registry.diagnostics_.duplicate_genotype_ids;
            }

            AnalysisSample sample;
            sample.key = key;
            sample.analysis_index = registry.samples_.size();
            sample.genotype_index = i;
            auto pit = phenotype_map.find(key);
            if (pit != phenotype_map.end()) {
                sample.phenotype_row_index = pit->second;
                ++registry.diagnostics_.matched_phenotypes;
            } else if (!phenotype_samples.empty()) {
                ++registry.diagnostics_.genotype_without_phenotype;
            }
            auto cit = covariate_map.find(key);
            if (cit != covariate_map.end()) {
                sample.covariate_row_index = cit->second;
                ++registry.diagnostics_.matched_covariates;
            } else if (!covariate_samples.empty()) {
                ++registry.diagnostics_.genotype_without_covariate;
            }
            auto dit = pedigree_map.find(key);
            if (dit != pedigree_map.end()) {
                sample.pedigree_index = dit->second;
                ++registry.diagnostics_.matched_pedigree;
            } else if (!pedigree_samples.empty()) {
                ++registry.diagnostics_.genotype_without_pedigree;
            }
            registry.analysis_index_by_key_[key] = sample.analysis_index;
            registry.samples_.push_back(std::move(sample));
        }

        registry.diagnostics_.analysis_samples = registry.samples_.size();
        registry.diagnostics_.phenotype_without_genotype =
            count_absent_from(phenotype_samples, genotype_seen);
        registry.diagnostics_.covariate_without_genotype =
            count_absent_from(covariate_samples, genotype_seen);
        registry.diagnostics_.pedigree_without_genotype =
            count_absent_from(pedigree_samples, genotype_seen);
        return registry;
    }

    const std::vector<AnalysisSample>& samples() const { return samples_; }
    const SampleAlignmentDiagnostics& diagnostics() const { return diagnostics_; }
    std::size_t size() const { return samples_.size(); }
    bool empty() const { return samples_.empty(); }

    const AnalysisSample* find(const SampleKey& key) const {
        auto it = analysis_index_by_key_.find(key);
        if (it == analysis_index_by_key_.end()) return nullptr;
        return &samples_[it->second];
    }

private:
    static std::unordered_map<SampleKey, std::size_t, SampleKeyHash>
    first_index_map(const std::vector<SampleKey>& keys, std::size_t& duplicate_count) {
        std::unordered_map<SampleKey, std::size_t, SampleKeyHash> out;
        for (std::size_t i = 0; i < keys.size(); ++i) {
            if (!out.emplace(keys[i], i).second) ++duplicate_count;
        }
        return out;
    }

    static std::size_t count_absent_from(
        const std::vector<SampleKey>& keys,
        const std::unordered_map<SampleKey, std::size_t, SampleKeyHash>& reference) {
        std::size_t n = 0;
        for (const auto& key : keys) {
            if (reference.find(key) == reference.end()) ++n;
        }
        return n;
    }

    std::vector<AnalysisSample> samples_;
    std::unordered_map<SampleKey, std::size_t, SampleKeyHash> analysis_index_by_key_;
    SampleAlignmentDiagnostics diagnostics_;
};

struct PhenotypeTable {
    std::vector<SampleKey> row_samples;
    std::vector<std::string> trait_names;
    std::vector<double> values;       // row-major: row * n_traits + trait
    std::vector<bool> missing_mask;   // true means missing/unusable

    std::size_t rows() const { return row_samples.size(); }
    std::size_t traits() const { return trait_names.size(); }
    double value(std::size_t row, std::size_t trait) const {
        return values[row * traits() + trait];
    }
    bool missing(std::size_t row, std::size_t trait) const {
        return missing_mask[row * traits() + trait];
    }
};

struct ObservationRecord {
    SampleKey sample;
    std::size_t observation_index = kAnalysisNpos;
    double time = 0.0;
    bool has_time = false;
    std::vector<double> trait_values;
    std::vector<bool> missing_mask;
};

struct ObservationTable {
    std::vector<std::string> trait_names;
    std::vector<ObservationRecord> observations;
};

struct CovariateTable {
    std::vector<SampleKey> row_samples;
    std::vector<std::string> numeric_names;
    std::vector<double> numeric_values;      // row-major
    std::vector<bool> numeric_missing_mask;  // true means missing/unusable
    std::vector<std::string> categorical_names;

    std::size_t rows() const { return row_samples.size(); }
    std::size_t numeric_cols() const { return numeric_names.size(); }
    double numeric_value(std::size_t row, std::size_t col) const {
        return numeric_values[row * numeric_cols() + col];
    }
};

struct PedigreeEntry {
    SampleKey id;
    SampleKey sire;
    SampleKey dam;
    bool sire_known = false;
    bool dam_known = false;
};

struct PedigreeTable {
    std::vector<PedigreeEntry> records;
};

enum class GenotypeStorageFormat {
    Unknown,
    PlinkBed,
    PlinkPgen,
    Bgen,
    Vcf,
    ExternalMatrix,
    None
};

inline const char* to_string(GenotypeStorageFormat format) {
    switch (format) {
        case GenotypeStorageFormat::PlinkBed: return "PLINK_BED";
        case GenotypeStorageFormat::PlinkPgen: return "PLINK_PGEN";
        case GenotypeStorageFormat::Bgen: return "BGEN";
        case GenotypeStorageFormat::Vcf: return "VCF";
        case GenotypeStorageFormat::ExternalMatrix: return "EXTERNAL_MATRIX";
        case GenotypeStorageFormat::None: return "NONE";
        case GenotypeStorageFormat::Unknown:
        default: return "UNKNOWN";
    }
}

struct GenotypeSourceCapabilities {
    bool supports_streaming = false;
    bool supports_random_access = false;
    bool supports_prefetch = false;
    bool supports_sample_subset = false;
    bool materializes_full_matrix = false;
};

struct MarkerMetadata {
    std::string chrom;
    std::string id;
    std::string ref;
    std::string alt;
    long long bp_pos = 0;
};

struct GenotypeSourceMetadata {
    GenotypeStorageFormat format = GenotypeStorageFormat::Unknown;
    std::size_t n_samples = 0;
    std::size_t n_markers = 0;
    std::vector<SampleKey> sample_order;
    std::vector<MarkerMetadata> markers;
    GenotypeSourceCapabilities capabilities;
    std::string source_label;
};

class GenotypeSource {
public:
    virtual ~GenotypeSource() = default;
    virtual const GenotypeSourceMetadata& metadata() const = 0;
};

class MetadataOnlyGenotypeSource final : public GenotypeSource {
public:
    explicit MetadataOnlyGenotypeSource(GenotypeSourceMetadata metadata)
        : metadata_(std::move(metadata)) {}

    const GenotypeSourceMetadata& metadata() const override { return metadata_; }

private:
    GenotypeSourceMetadata metadata_;
};

class AnalysisDataset {
public:
    void setSampleRegistry(SampleRegistry registry) { sample_registry_ = std::move(registry); }
    void setPhenotypeTable(PhenotypeTable table) { phenotype_table_ = std::move(table); has_phenotype_ = true; }
    void setObservationTable(ObservationTable table) { observation_table_ = std::move(table); has_observations_ = true; }
    void setCovariateTable(CovariateTable table) { covariate_table_ = std::move(table); has_covariates_ = true; }
    void setPedigreeTable(PedigreeTable table) { pedigree_table_ = std::move(table); has_pedigree_ = true; }
    void setGenotypeSource(std::shared_ptr<const GenotypeSource> source) { genotype_source_ = std::move(source); }

    const SampleRegistry& samples() const { return sample_registry_; }
    const PhenotypeTable* phenotype() const { return has_phenotype_ ? &phenotype_table_ : nullptr; }
    const ObservationTable* observations() const { return has_observations_ ? &observation_table_ : nullptr; }
    const CovariateTable* covariates() const { return has_covariates_ ? &covariate_table_ : nullptr; }
    const PedigreeTable* pedigree() const { return has_pedigree_ ? &pedigree_table_ : nullptr; }
    std::shared_ptr<const GenotypeSource> genotype() const { return genotype_source_; }

    std::string diagnosticSummary(const std::string& mode_label = "shadow mode") const {
        std::ostringstream os;
        os << "AnalysisDataset " << mode_label << "\n";
        os << "Samples              : " << sample_registry_.size() << "\n";
        if (genotype_source_) {
            const auto& meta = genotype_source_->metadata();
            os << "Genotype format      : " << to_string(meta.format) << "\n";
            os << "Genotype samples     : " << meta.n_samples << "\n";
            os << "Markers              : " << meta.n_markers << "\n";
            os << "Genotype materialized: " << (meta.capabilities.materializes_full_matrix ? "YES" : "NO") << "\n";
        } else {
            os << "Genotype format      : NONE\n";
        }
        if (has_phenotype_) {
            os << "Phenotype rows       : " << phenotype_table_.rows() << "\n";
            os << "Traits               : " << phenotype_table_.traits() << "\n";
        }
        if (has_observations_) {
            os << "Observations         : " << observation_table_.observations.size() << "\n";
        }
        if (has_covariates_) {
            os << "Numeric covariates   : " << covariate_table_.numeric_cols() << "\n";
        }
        const auto& diag = sample_registry_.diagnostics();
        os << "Matched phenotypes   : " << diag.matched_phenotypes << "\n";
        os << "Unmatched genotype   : " << diag.genotype_without_phenotype << "\n";
        os << "Unmatched phenotype  : " << diag.phenotype_without_genotype << "\n";
        return os.str();
    }

private:
    SampleRegistry sample_registry_;
    PhenotypeTable phenotype_table_;
    ObservationTable observation_table_;
    CovariateTable covariate_table_;
    PedigreeTable pedigree_table_;
    std::shared_ptr<const GenotypeSource> genotype_source_;
    bool has_phenotype_ = false;
    bool has_observations_ = false;
    bool has_covariates_ = false;
    bool has_pedigree_ = false;
};

inline std::vector<SampleKey> sample_keys_from_idmap_order(const std::map<std::string, int>& idmap) {
    std::vector<SampleKey> keys(idmap.size());
    for (const auto& kv : idmap) {
        if (kv.second >= 1 && kv.second <= static_cast<int>(keys.size())) {
            keys[static_cast<std::size_t>(kv.second - 1)] = make_iid_sample_key(kv.first);
        }
    }
    return keys;
}

inline PhenotypeTable phenotype_table_from_records(const std::vector<SampleKey>& row_samples,
                                                   const std::vector<double>& values,
                                                   const std::vector<std::string>& trait_names,
                                                   const std::vector<bool>& missing_mask = {}) {
    PhenotypeTable table;
    table.row_samples = row_samples;
    table.trait_names = trait_names.empty() ? std::vector<std::string>{"trait"} : trait_names;
    table.values = values;
    table.missing_mask = missing_mask.empty()
        ? std::vector<bool>(values.size(), false)
        : missing_mask;
    return table;
}

inline CovariateTable covariate_table_from_numeric(const std::vector<SampleKey>& row_samples,
                                                   const std::vector<std::string>& names,
                                                   const std::vector<double>& values,
                                                   const std::vector<bool>& missing_mask = {}) {
    CovariateTable table;
    table.row_samples = row_samples;
    table.numeric_names = names;
    table.numeric_values = values;
    table.numeric_missing_mask = missing_mask.empty()
        ? std::vector<bool>(values.size(), false)
        : missing_mask;
    return table;
}

inline std::vector<double> phenotype_values_in_analysis_order(const AnalysisDataset& dataset,
                                                              std::size_t trait_index = 0) {
    const auto* phenotype = dataset.phenotype();
    if (!phenotype) throw std::runtime_error("AnalysisDataset has no phenotype table.");
    if (trait_index >= phenotype->traits()) throw std::runtime_error("Phenotype trait index is out of range.");

    std::vector<double> out(dataset.samples().size(), std::numeric_limits<double>::quiet_NaN());
    for (const auto& sample : dataset.samples().samples()) {
        if (sample.phenotype_row_index == kAnalysisNpos) continue;
        if (sample.phenotype_row_index >= phenotype->rows()) {
            throw std::runtime_error("Phenotype row index is out of range for AnalysisDataset.");
        }
        out[sample.analysis_index] = phenotype->value(sample.phenotype_row_index, trait_index);
    }
    return out;
}

inline std::vector<double> covariate_values_in_analysis_order(const AnalysisDataset& dataset,
                                                              std::size_t* n_columns = nullptr) {
    const auto* covariates = dataset.covariates();
    const std::size_t cols = covariates ? covariates->numeric_cols() : 0;
    if (n_columns) *n_columns = cols;
    if (!covariates || cols == 0) return {};

    std::vector<double> out(dataset.samples().size() * cols,
                            std::numeric_limits<double>::quiet_NaN());
    for (const auto& sample : dataset.samples().samples()) {
        if (sample.covariate_row_index == kAnalysisNpos) continue;
        if (sample.covariate_row_index >= covariates->rows()) {
            throw std::runtime_error("Covariate row index is out of range for AnalysisDataset.");
        }
        for (std::size_t c = 0; c < cols; ++c) {
            out[sample.analysis_index * cols + c] =
                covariates->numeric_value(sample.covariate_row_index, c);
        }
    }
    return out;
}

inline std::shared_ptr<const GenotypeSource>
make_external_matrix_genotype_source(const std::vector<SampleKey>& sample_order,
                                     const std::string& source_label = "") {
    GenotypeSourceMetadata meta;
    meta.format = GenotypeStorageFormat::ExternalMatrix;
    meta.n_samples = sample_order.size();
    meta.n_markers = 0;
    meta.source_label = source_label;
    meta.sample_order = sample_order;
    meta.capabilities.materializes_full_matrix = false;
    return std::make_shared<MetadataOnlyGenotypeSource>(std::move(meta));
}

} // namespace cosmic
