#include "text_geno_reader.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace cosmic {
namespace {

std::vector<std::string> split_fields(const std::string& line) {
    std::istringstream input(line);
    std::vector<std::string> fields;
    std::string field;
    while (input >> field) fields.push_back(field);
    return fields;
}

uint8_t parse_genotype(const std::string& value) {
    if (value == "NA" || value == "na" || value == "NaN" || value == "nan" ||
        value == "." || value == "./." || value == ".|.") return 3;
    if (value == "0/0" || value == "0|0") return 0;
    if (value == "0/1" || value == "1/0" || value == "0|1" || value == "1|0") return 1;
    if (value == "1/1" || value == "1|1") return 2;
    size_t consumed = 0;
    double numeric = 0.0;
    try {
        numeric = std::stod(value, &consumed);
    } catch (...) {
        throw std::runtime_error("Invalid text genotype value: " + value);
    }
    if (consumed != value.size() || !std::isfinite(numeric) ||
        (numeric != 0.0 && numeric != 1.0 && numeric != 2.0)) {
        throw std::runtime_error("Text genotypes must be 0, 1, 2, or NA; found: " + value);
    }
    return static_cast<uint8_t>(numeric);
}

} // namespace

void TextGenoReader::open(const std::string& geno_file,
                          TextGenoOrientation orientation,
                          const std::string& var_info_file,
                          const std::string& sample_info_file) {
    close();
    geno_file_ = geno_file;
    orientation_ = orientation;
    var_info_file_ = var_info_file;
    sample_info_file_ = sample_info_file;

    std::ifstream input(geno_file);
    if (!input) throw std::runtime_error("Cannot open text genotype file: " + geno_file);
    std::vector<std::vector<std::string>> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto fields = split_fields(line);
        if (!fields.empty()) rows.push_back(std::move(fields));
    }
    if (rows.size() < 2 || rows[0].size() < 2) {
        throw std::runtime_error("Text genotype matrix requires a header and at least one data row");
    }

    const size_t expected_width = rows[0].size();
    for (size_t row = 1; row < rows.size(); ++row) {
        if (rows[row].size() != expected_width) {
            throw std::runtime_error("Text genotype row " + std::to_string(row + 1) + " has " +
                                     std::to_string(rows[row].size()) + " columns; expected " +
                                     std::to_string(expected_width));
        }
    }

    if (orientation_ == TextGenoOrientation::VariantsBySamples) {
        num_samples_ = static_cast<uint32_t>(expected_width - 1);
        num_variants_ = static_cast<uint32_t>(rows.size() - 1);
        sample_ids_.assign(rows[0].begin() + 1, rows[0].end());
        genotypes_.assign(num_samples_, std::vector<uint8_t>(num_variants_, 3));
        variants_.reserve(num_variants_);
        for (size_t j = 0; j < num_variants_; ++j) {
            variants_.push_back({rows[j + 1][0], "0", static_cast<uint32_t>(j + 1), "A", "C"});
            for (size_t i = 0; i < num_samples_; ++i) {
                genotypes_[i][j] = parse_genotype(rows[j + 1][i + 1]);
            }
        }
    } else {
        num_variants_ = static_cast<uint32_t>(expected_width - 1);
        num_samples_ = static_cast<uint32_t>(rows.size() - 1);
        genotypes_.assign(num_samples_, std::vector<uint8_t>(num_variants_, 3));
        variants_.reserve(num_variants_);
        for (size_t j = 0; j < num_variants_; ++j) {
            variants_.push_back({rows[0][j + 1], "0", static_cast<uint32_t>(j + 1), "A", "C"});
        }
        sample_ids_.reserve(num_samples_);
        for (size_t i = 0; i < num_samples_; ++i) {
            sample_ids_.push_back(rows[i + 1][0]);
            for (size_t j = 0; j < num_variants_; ++j) {
                genotypes_[i][j] = parse_genotype(rows[i + 1][j + 1]);
            }
        }
    }

    if (std::unordered_set<std::string>(sample_ids_.begin(), sample_ids_.end()).size() != sample_ids_.size()) {
        throw std::runtime_error("Duplicate sample ID in text genotype matrix");
    }
    std::unordered_set<std::string> variant_ids;
    for (const auto& variant : variants_) {
        if (!variant_ids.insert(variant.id).second) {
            throw std::runtime_error("Duplicate variant ID in text genotype matrix: " + variant.id);
        }
    }

    if (!var_info_file.empty()) {
        std::ifstream metadata(var_info_file);
        if (!metadata) throw std::runtime_error("Cannot open variant metadata file: " + var_info_file);
        if (!std::getline(metadata, line)) throw std::runtime_error("Empty variant metadata file: " + var_info_file);
        std::unordered_map<std::string, TextGenoVariant> by_id;
        size_t line_number = 1;
        while (std::getline(metadata, line)) {
            ++line_number;
            if (line.empty() || line[0] == '#') continue;
            const auto fields = split_fields(line);
            if (fields.size() != 5) {
                throw std::runtime_error("Variant metadata row " + std::to_string(line_number) + " must have ID CHROM POS REF ALT");
            }
            TextGenoVariant variant{fields[0], fields[1], 0, fields[3], fields[4]};
            try {
                variant.pos = static_cast<uint32_t>(std::stoul(fields[2]));
            } catch (...) {
                throw std::runtime_error("Invalid variant position at metadata row " + std::to_string(line_number));
            }
            if (!by_id.emplace(variant.id, variant).second) {
                throw std::runtime_error("Duplicate variant ID in metadata: " + variant.id);
            }
        }
        if (by_id.size() != variants_.size()) throw std::runtime_error("Variant metadata count does not match genotype matrix");
        for (auto& variant : variants_) {
            auto found = by_id.find(variant.id);
            if (found == by_id.end()) throw std::runtime_error("Missing variant metadata for: " + variant.id);
            variant = found->second;
        }
    }

    if (!sample_info_file.empty()) {
        std::ifstream metadata(sample_info_file);
        if (!metadata) throw std::runtime_error("Cannot open sample metadata file: " + sample_info_file);
        if (!std::getline(metadata, line)) throw std::runtime_error("Empty sample metadata file: " + sample_info_file);
        std::vector<std::string> mapped_ids;
        while (std::getline(metadata, line)) {
            if (line.empty() || line[0] == '#') continue;
            const auto fields = split_fields(line);
            if (fields.size() == 1) mapped_ids.push_back(fields[0]);
            else if (fields.size() == 2) mapped_ids.push_back(fields[0] + ":" + fields[1]);
            else throw std::runtime_error("Sample metadata rows must contain IID or FID IID");
        }
        if (mapped_ids.size() != sample_ids_.size()) throw std::runtime_error("Sample metadata count does not match genotype matrix");
        sample_ids_ = std::move(mapped_ids);
    }

    std::string last_chrom;
    std::unordered_set<std::string> closed_chroms;
    uint32_t chrom_start = 0;
    for (uint32_t j = 0; j < num_variants_; ++j) {
        const std::string& chrom = variants_[j].chrom;
        if (chrom == last_chrom) continue;
        if (!last_chrom.empty()) {
            chrom_ranges_.emplace_back(chrom_start, j - chrom_start);
            closed_chroms.insert(last_chrom);
        }
        if (closed_chroms.count(chrom)) throw std::runtime_error("Text variants must be grouped by chromosome");
        chrom_labels_.push_back(chrom);
        last_chrom = chrom;
        chrom_start = j;
    }
    chrom_ranges_.emplace_back(chrom_start, num_variants_ - chrom_start);
}

void TextGenoReader::close() {
    num_samples_ = 0;
    num_variants_ = 0;
    sample_ids_.clear();
    variants_.clear();
    chrom_labels_.clear();
    chrom_ranges_.clear();
    genotypes_.clear();
    snp_means_.clear();
    snp_inv_stds_.clear();
    stats_computed_ = false;
    current_variant_idx_ = 0;
}

bool TextGenoReader::readVariant(std::string& chrom, std::string& id, uint32_t& pos,
                                 std::string& ref, std::string& alt,
                                 Eigen::VectorXf& dosage) const {
    if (current_variant_idx_ >= variants_.size()) return false;
    const size_t index = current_variant_idx_++;
    const auto& variant = variants_[index];
    chrom = variant.chrom;
    id = variant.id;
    pos = variant.pos;
    ref = variant.ref;
    alt = variant.alt;
    dosage.resize(num_samples_);
    for (size_t i = 0; i < num_samples_; ++i) {
        dosage[i] = genotypes_[i][index] <= 2
            ? static_cast<float>(genotypes_[i][index])
            : std::numeric_limits<float>::quiet_NaN();
    }
    return true;
}

void TextGenoReader::computeStats(int n_threads) {
    if (stats_computed_) return;
    snp_means_.assign(num_variants_, 0.0);
    snp_inv_stds_.assign(num_variants_, 0.0);
#pragma omp parallel for num_threads(n_threads)
    for (int j = 0; j < static_cast<int>(num_variants_); ++j) {
        double sum = 0.0;
        int count = 0;
        for (size_t i = 0; i < num_samples_; ++i) {
            if (genotypes_[i][j] <= 2) { sum += genotypes_[i][j]; ++count; }
        }
        const double mean = count ? sum / count : 0.0;
        const double p = mean / 2.0;
        const double variance = 2.0 * p * (1.0 - p);
        snp_means_[j] = mean;
        snp_inv_stds_[j] = variance > 1e-12 ? 1.0 / std::sqrt(variance) : 0.0;
    }
    stats_computed_ = true;
}

void TextGenoReader::multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (v.size() != static_cast<int>(num_variants_)) throw std::invalid_argument("Text genotype Z*v dimension mismatch");
    if (!stats_computed_) computeStats(n_threads);
    out.setZero(num_samples_);
#pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < static_cast<int>(num_samples_); ++i) {
        double value = 0.0;
        for (size_t j = 0; j < num_variants_; ++j) {
            const uint8_t gt = genotypes_[i][j];
            if (gt <= 2) value += (gt - snp_means_[j]) * snp_inv_stds_[j] * v[j];
        }
        out[i] = value;
    }
}

void TextGenoReader::multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (v.size() != static_cast<int>(num_samples_)) throw std::invalid_argument("Text genotype Z'*v dimension mismatch");
    if (!stats_computed_) computeStats(n_threads);
    out.setZero(num_variants_);
#pragma omp parallel for num_threads(n_threads)
    for (int j = 0; j < static_cast<int>(num_variants_); ++j) {
        double value = 0.0;
        for (size_t i = 0; i < num_samples_; ++i) {
            const uint8_t gt = genotypes_[i][j];
            if (gt <= 2) value += (gt - snp_means_[j]) * snp_inv_stds_[j] * v[i];
        }
        out[j] = value;
    }
}

void TextGenoReader::seekVariant(uint32_t idx) const {
    current_variant_idx_ = std::min<size_t>(idx, num_variants_);
}

std::pair<uint32_t, uint32_t> TextGenoReader::getChromRange(const std::string& chrom) const {
    for (size_t i = 0; i < chrom_labels_.size(); ++i) {
        if (chrom_labels_[i] == chrom) return chrom_ranges_[i];
    }
    return {0, 0};
}

void TextGenoReader::exportToText(const std::string& geno_file,
                                  const std::string& var_info_file,
                                  const std::string& sample_info_file) const {
    std::ofstream output(geno_file);
    if (!output) throw std::runtime_error("Cannot write genotype file: " + geno_file);
    static const char* codes[] = {"0", "1", "2", "NA"};
    if (orientation_ == TextGenoOrientation::VariantsBySamples) {
        output << "VariantID";
        for (const auto& id : sample_ids_) output << '\t' << id;
        output << '\n';
        for (size_t j = 0; j < num_variants_; ++j) {
            output << variants_[j].id;
            for (size_t i = 0; i < num_samples_; ++i) output << '\t' << codes[genotypes_[i][j]];
            output << '\n';
        }
    } else {
        output << "SampleID";
        for (const auto& variant : variants_) output << '\t' << variant.id;
        output << '\n';
        for (size_t i = 0; i < num_samples_; ++i) {
            output << sample_ids_[i];
            for (size_t j = 0; j < num_variants_; ++j) output << '\t' << codes[genotypes_[i][j]];
            output << '\n';
        }
    }

    if (!var_info_file.empty()) {
        std::ofstream metadata(var_info_file);
        if (!metadata) throw std::runtime_error("Cannot write variant metadata file: " + var_info_file);
        metadata << "ID\tCHROM\tPOS\tREF\tALT\n";
        for (const auto& variant : variants_) {
            metadata << variant.id << '\t' << variant.chrom << '\t' << variant.pos << '\t'
                     << variant.ref << '\t' << variant.alt << '\n';
        }
    }
    if (!sample_info_file.empty()) {
        std::ofstream metadata(sample_info_file);
        if (!metadata) throw std::runtime_error("Cannot write sample metadata file: " + sample_info_file);
        metadata << "FID\tIID\n";
        for (const auto& id : sample_ids_) {
            const size_t separator = id.find(':');
            if (separator == std::string::npos) metadata << id << '\t' << id << '\n';
            else metadata << id.substr(0, separator) << '\t' << id.substr(separator + 1) << '\n';
        }
    }
}

} // namespace cosmic

namespace cosmic {

void TextGenoReader::exportToNpy(const std::string& npy_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream file(npy_file, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open output NPY file: " + npy_file);
    
    uint32_t n_variants = num_variants_;
    uint32_t n_samples = num_samples_;
    
    file.write("\x93NUMPY", 6);
    file.put(1); file.put(0);
    
    size_t header_len_pos = file.tellp();
    file.write("\x00\x00", 2);
    
    if (orientation_ == TextGenoOrientation::VariantsBySamples) {
        for (size_t j = 0; j < num_variants_; ++j) {
            for (size_t i = 0; i < num_samples_; ++i) {
                uint8_t code = genotypes_[i][j];
                float dosage = (code == 3) ? std::numeric_limits<float>::quiet_NaN() : static_cast<float>(code);
                file.write(reinterpret_cast<const char*>(&dosage), 4);
            }
        }
    } else {
        // SamplesByVariants: transpose
        std::vector<float> buffer(n_samples * n_variants);
        for (size_t i = 0; i < n_samples; ++i) {
            for (size_t j = 0; j < n_variants; ++j) {
                uint8_t code = genotypes_[i][j];
                buffer[i * n_variants + j] = (code == 3) ? std::numeric_limits<float>::quiet_NaN() : static_cast<float>(code);
            }
        }
        file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(float));
    }
    
    size_t data_end = file.tellp();
    
    std::ostringstream header;
    header << "{'descr': '<f4', 'fortran_order': False, 'shape': (" << n_variants << ", " << n_samples << "), }";
    std::string header_str = header.str();
    while ((header_str.size() + 10) % 64 != 0) header_str += " ";
    size_t header_len = header_str.size();
    
    file.seekp(header_len_pos);
    uint16_t header_len_le = static_cast<uint16_t>(header_len);
    file.write(reinterpret_cast<const char*>(&header_len_le), 2);
    file.write(header_str.c_str(), header_len);
    file.seekp(data_end);
    
    if (!var_info_file.empty()) {
        std::ofstream metadata(var_info_file);
        metadata << "ID\tCHROM\tPOS\tREF\tALT\n";
        for (const auto& variant : variants_) {
            metadata << variant.id << '\t' << variant.chrom << '\t' << variant.pos << '\t'
                     << variant.ref << '\t' << variant.alt << '\n';
        }
    }
    if (!sample_info_file.empty()) {
        std::ofstream metadata(sample_info_file);
        metadata << "FID\tIID\n";
        for (const auto& id : sample_ids_) {
            const size_t separator = id.find(':');
            if (separator == std::string::npos) metadata << id << '\t' << id << '\n';
            else metadata << id.substr(0, separator) << '\t' << id.substr(separator + 1) << '\n';
        }
    }
}

void TextGenoReader::exportToPedMap(const std::string& ped_file, const std::string& map_file) const {
    std::ofstream map_out(map_file);
    if (!map_out) throw std::runtime_error("Cannot open output MAP file: " + map_file);
    
    for (const auto& variant : variants_) {
        map_out << variant.chrom << " " << variant.id << " 0 " << variant.pos << " " << variant.ref << " " << variant.alt << "\n";
    }
    
    std::ofstream ped_out(ped_file);
    if (!ped_out) throw std::runtime_error("Cannot open output PED file: " + ped_file);
    
    for (size_t i = 0; i < num_samples_; ++i) {
        ped_out << sample_ids_[i] << " " << sample_ids_[i] << " 0 0 1 -9";
        for (size_t j = 0; j < num_variants_; ++j) {
            uint8_t code = genotypes_[i][j];
            if (code == 3) {
                ped_out << " 0 0";
            } else {
                const auto& variant = variants_[j];
                if (code == 0) ped_out << " " << variant.ref << " " << variant.ref;
                else if (code == 1) ped_out << " " << variant.ref << " " << variant.alt;
                else ped_out << " " << variant.alt << " " << variant.alt;
            }
        }
        ped_out << "\n";
    }
}

void TextGenoReader::exportToVcf(const std::string& vcf_file) const {
    std::ofstream vcf_out(vcf_file);
    if (!vcf_out) throw std::runtime_error("Cannot open output VCF file: " + vcf_file);
    
    vcf_out << "##fileformat=VCFv4.2\n";
    vcf_out << "##source=CosmicScanner\n";
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& sid : sample_ids_) {
        vcf_out << "\t" << sid;
    }
    vcf_out << "\n";
    
    for (size_t v = 0; v < variants_.size(); ++v) {
        const auto& variant = variants_[v];
        vcf_out << variant.chrom << "\t" << variant.pos << "\t" << variant.id << "\t"
                << variant.ref << "\t" << variant.alt << "\t.\tPASS\t.\tGT";
        
        for (size_t i = 0; i < num_samples_; ++i) {
            uint8_t code = genotypes_[i][v];
            if (code == 3) {
                vcf_out << "\t./.";
            } else {
                if (code == 0) vcf_out << "\t0/0";
                else if (code == 1) vcf_out << "\t0/1";
                else vcf_out << "\t1/1";
            }
        }
        vcf_out << "\n";
    }
}

} // namespace cosmic
