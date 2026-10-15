#pragma once

#include "analysis_dataset.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "plink_reader.h"

#include <memory>
#include <string>
#include <vector>

namespace cosmic {

inline std::vector<SampleKey> sample_keys_from_strings(const std::vector<std::string>& ids) {
    std::vector<SampleKey> keys;
    keys.reserve(ids.size());
    for (const auto& id : ids) keys.push_back(make_iid_sample_key(id));
    return keys;
}

inline std::vector<SampleKey> sample_keys_from_plink_fam(const std::vector<PlinkFamInfo>& fam) {
    std::vector<SampleKey> keys;
    keys.reserve(fam.size());
    for (const auto& row : fam) keys.push_back(SampleKey{row.fid, row.iid});
    return keys;
}

inline std::shared_ptr<const GenotypeSource>
make_plink_bed_genotype_source(const PlinkReader& reader, const std::string& source_label = "") {
    GenotypeSourceMetadata meta;
    meta.format = GenotypeStorageFormat::PlinkBed;
    meta.n_samples = reader.getNumSamples();
    meta.n_markers = reader.getNumSnps();
    meta.source_label = source_label;
    meta.sample_order = sample_keys_from_plink_fam(reader.getFamInfo());
    const auto& bim = reader.getBimInfo();
    meta.markers.reserve(bim.size());
    for (const auto& row : bim) {
        MarkerMetadata marker;
        marker.chrom = row.chrom;
        marker.id = row.id;
        marker.ref = row.ref;
        marker.alt = row.alt;
        marker.bp_pos = row.bp_pos;
        meta.markers.push_back(std::move(marker));
    }
    meta.capabilities.supports_streaming = true;
    meta.capabilities.supports_random_access = true;
    meta.capabilities.supports_sample_subset = true;
    meta.capabilities.materializes_full_matrix = false;
    return std::make_shared<MetadataOnlyGenotypeSource>(std::move(meta));
}

inline std::shared_ptr<const GenotypeSource>
make_plink_pgen_genotype_source(const PgenReader& reader, const std::string& source_label = "") {
    GenotypeSourceMetadata meta;
    meta.format = GenotypeStorageFormat::PlinkPgen;
    meta.n_samples = reader.getNumSamples();
    meta.n_markers = reader.getNumVariants();
    meta.source_label = source_label;
    meta.sample_order = sample_keys_from_strings(reader.getSampleIds());
    const auto& variants = reader.getVariants();
    meta.markers.reserve(variants.size());
    for (const auto& row : variants) {
        MarkerMetadata marker;
        marker.chrom = row.chrom;
        marker.id = row.id;
        marker.ref = row.ref;
        marker.alt = row.alt;
        marker.bp_pos = row.pos;
        meta.markers.push_back(std::move(marker));
    }
    meta.capabilities.supports_streaming = true;
    meta.capabilities.supports_random_access = true;
    meta.capabilities.materializes_full_matrix = false;
    return std::make_shared<MetadataOnlyGenotypeSource>(std::move(meta));
}

inline std::shared_ptr<const GenotypeSource>
make_bgen_genotype_source(const BgenReader& reader, const std::string& source_label = "") {
    GenotypeSourceMetadata meta;
    meta.format = GenotypeStorageFormat::Bgen;
    meta.n_samples = reader.getNumSamples();
    meta.n_markers = reader.getNumVariants();
    meta.source_label = source_label;
    meta.sample_order = sample_keys_from_strings(reader.getSampleIds());
    const auto& variants = reader.getVariants();
    meta.markers.reserve(variants.size());
    for (const auto& row : variants) {
        MarkerMetadata marker;
        marker.chrom = row.chrom;
        marker.id = row.rsid.empty() ? row.id : row.rsid;
        if (!row.alleles.empty()) marker.ref = row.alleles.front();
        if (row.alleles.size() > 1) marker.alt = row.alleles[1];
        marker.bp_pos = row.pos;
        meta.markers.push_back(std::move(marker));
    }
    meta.capabilities.supports_streaming = true;
    meta.capabilities.supports_prefetch = true;
    meta.capabilities.materializes_full_matrix = false;
    return std::make_shared<MetadataOnlyGenotypeSource>(std::move(meta));
}

} // namespace cosmic
