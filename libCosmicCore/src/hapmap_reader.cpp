#include "hapmap_reader.h"
#include "genotype_converter.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
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

std::pair<std::string, std::string> split_alleles(const std::string& value) {
    const size_t separator = value.find('/');
    if (separator == std::string::npos || separator == 0 || separator + 1 >= value.size()) {
        throw std::runtime_error("Invalid HapMap alleles field: " + value);
    }
    return {value.substr(0, separator), value.substr(separator + 1)};
}

} // namespace

void HapMapReader::open(const std::string& hmp_file) {
    close();
    hmp_file_ = hmp_file;

    std::ifstream input(hmp_file);
    if (!input) throw std::runtime_error("Cannot open HapMap file: " + hmp_file);

    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("Empty HapMap file: " + hmp_file);
    const auto header = split_fields(line);
    if (header.size() < 12 || header[0] != "rs#" || header[1] != "alleles" ||
        header[2] != "chrom" || header[3] != "pos") {
        throw std::runtime_error("HapMap header must contain the standard 11 metadata columns followed by samples");
    }
    sample_ids_.assign(header.begin() + 11, header.end());
    if (std::unordered_set<std::string>(sample_ids_.begin(), sample_ids_.end()).size() != sample_ids_.size()) {
        throw std::runtime_error("Duplicate sample ID in HapMap header");
    }
    num_samples_ = static_cast<uint32_t>(sample_ids_.size());

    std::unordered_set<std::string> variant_ids;
    size_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty()) continue;
        const auto fields = split_fields(line);
        if (fields.size() != 11 + num_samples_) {
            throw std::runtime_error("HapMap row " + std::to_string(line_number) + " has " +
                                     std::to_string(fields.size()) + " columns; expected " +
                                     std::to_string(11 + num_samples_));
        }
        if (!variant_ids.insert(fields[0]).second) {
            throw std::runtime_error("Duplicate HapMap variant ID: " + fields[0]);
        }

        const auto alleles = split_alleles(fields[1]);
        if (alleles.first.size() != 1 || alleles.second.size() != 1) {
            throw std::runtime_error("Only biallelic single-base HapMap variants are supported: " + fields[0]);
        }

        HapMapVariant variant;
        variant.id = fields[0];
        variant.allele1 = alleles.first;
        variant.allele2 = alleles.second;
        variant.chrom = fields[2];
        try {
            variant.pos = static_cast<uint32_t>(std::stoul(fields[3]));
        } catch (...) {
            throw std::runtime_error("Invalid HapMap position at row " + std::to_string(line_number));
        }
        variant.genotypes.resize(num_samples_, 3);
        for (size_t i = 0; i < num_samples_; ++i) {
            const float dosage = GenotypeConverter::hapmapToDosage(
                fields[11 + i], variant.allele1[0], variant.allele2[0]);
            if (std::isfinite(dosage)) variant.genotypes[i] = static_cast<uint8_t>(std::lround(dosage));
        }
        variants_.push_back(std::move(variant));
    }
    num_variants_ = static_cast<uint32_t>(variants_.size());
    if (num_variants_ == 0) throw std::runtime_error("HapMap file has no variants");

    std::string last_chrom;
    std::unordered_set<std::string> closed_chroms;
    uint32_t chrom_start = 0;
    for (uint32_t i = 0; i < num_variants_; ++i) {
        const std::string& chrom = variants_[i].chrom;
        if (chrom == last_chrom) continue;
        if (!last_chrom.empty()) {
            chrom_ranges_.emplace_back(chrom_start, i - chrom_start);
            closed_chroms.insert(last_chrom);
        }
        if (closed_chroms.count(chrom)) {
            throw std::runtime_error("HapMap variants must be grouped by chromosome");
        }
        chrom_labels_.push_back(chrom);
        last_chrom = chrom;
        chrom_start = i;
    }
    chrom_ranges_.emplace_back(chrom_start, num_variants_ - chrom_start);
}

void HapMapReader::close() {
    hmp_file_.clear();
    num_samples_ = 0;
    num_variants_ = 0;
    sample_ids_.clear();
    variants_.clear();
    chrom_labels_.clear();
    chrom_ranges_.clear();
    snp_means_.clear();
    snp_inv_stds_.clear();
    stats_computed_ = false;
    current_variant_idx_ = 0;
}

bool HapMapReader::readVariant(std::string& chrom, std::string& id, uint32_t& pos,
                               std::string& ref, std::string& alt,
                               Eigen::VectorXf& dosage) const {
    if (current_variant_idx_ >= variants_.size()) return false;
    const auto& variant = variants_[current_variant_idx_++];
    chrom = variant.chrom;
    id = variant.id;
    pos = variant.pos;
    ref = variant.allele1;
    alt = variant.allele2;
    dosage.resize(num_samples_);
    for (size_t i = 0; i < num_samples_; ++i) {
        dosage[i] = variant.genotypes[i] <= 2
            ? static_cast<float>(variant.genotypes[i])
            : std::numeric_limits<float>::quiet_NaN();
    }
    return true;
}

void HapMapReader::computeStats(int n_threads) {
    if (stats_computed_) return;
    snp_means_.assign(num_variants_, 0.0);
    snp_inv_stds_.assign(num_variants_, 0.0);
#pragma omp parallel for num_threads(n_threads)
    for (int j = 0; j < static_cast<int>(num_variants_); ++j) {
        double sum = 0.0;
        int count = 0;
        for (uint8_t gt : variants_[j].genotypes) {
            if (gt <= 2) { sum += gt; ++count; }
        }
        const double mean = count ? sum / count : 0.0;
        const double p = mean / 2.0;
        const double variance = 2.0 * p * (1.0 - p);
        snp_means_[j] = mean;
        snp_inv_stds_[j] = variance > 1e-12 ? 1.0 / std::sqrt(variance) : 0.0;
    }
    stats_computed_ = true;
}

void HapMapReader::multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (v.size() != static_cast<int>(num_variants_)) throw std::invalid_argument("HapMap Z*v dimension mismatch");
    if (!stats_computed_) computeStats(n_threads);
    out.setZero(num_samples_);
#pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < static_cast<int>(num_samples_); ++i) {
        double value = 0.0;
        for (size_t j = 0; j < num_variants_; ++j) {
            const uint8_t gt = variants_[j].genotypes[i];
            if (gt <= 2) value += (gt - snp_means_[j]) * snp_inv_stds_[j] * v[j];
        }
        out[i] = value;
    }
}

void HapMapReader::multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (v.size() != static_cast<int>(num_samples_)) throw std::invalid_argument("HapMap Z'*v dimension mismatch");
    if (!stats_computed_) computeStats(n_threads);
    out.setZero(num_variants_);
#pragma omp parallel for num_threads(n_threads)
    for (int j = 0; j < static_cast<int>(num_variants_); ++j) {
        double value = 0.0;
        for (size_t i = 0; i < num_samples_; ++i) {
            const uint8_t gt = variants_[j].genotypes[i];
            if (gt <= 2) value += (gt - snp_means_[j]) * snp_inv_stds_[j] * v[i];
        }
        out[j] = value;
    }
}

void HapMapReader::seekVariant(uint32_t idx) const {
    current_variant_idx_ = std::min<size_t>(idx, num_variants_);
}

std::pair<uint32_t, uint32_t> HapMapReader::getChromRange(const std::string& chrom) const {
    for (size_t i = 0; i < chrom_labels_.size(); ++i) {
        if (chrom_labels_[i] == chrom) return chrom_ranges_[i];
    }
    return {0, 0};
}

} // namespace cosmic

namespace cosmic {

void HapMapReader::exportToText(const std::string& geno_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream geno_out(geno_file);
    if (!geno_out) throw std::runtime_error("Cannot open output geno file: " + geno_file);
    
    geno_out << "VariantID";
    for (const auto& sid : sample_ids_) {
        geno_out << " " << sid;
    }
    geno_out << "\n";
    
    for (const auto& var : variants_) {
        geno_out << var.id;
        for (size_t i = 0; i < num_samples_; ++i) {
            uint8_t code = var.genotypes[i];
            if (code == 3) {
                geno_out << " NA";
            } else {
                geno_out << " " << static_cast<int>(code);
            }
        }
        geno_out << "\n";
    }
    
    std::ofstream var_out(var_info_file);
    var_out << "ID CHROM POS REF ALT\n";
    for (const auto& var : variants_) {
        var_out << var.id << " " << var.chrom << " " << var.pos << " " << var.allele1 << " " << var.allele2 << "\n";
    }
    
    std::ofstream sample_out(sample_info_file);
    sample_out << "FID IID\n";
    for (const auto& sid : sample_ids_) {
        size_t sep = sid.find(':');
        if (sep != std::string::npos) {
            sample_out << sid.substr(0, sep) << " " << sid.substr(sep + 1) << "\n";
        } else {
            sample_out << sid << " " << sid << "\n";
        }
    }
}

void HapMapReader::exportToNpy(const std::string& npy_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream file(npy_file, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open output NPY file: " + npy_file);
    
    uint32_t n_variants = variants_.size();
    uint32_t n_samples = num_samples_;
    
    file.write("\x93NUMPY", 6);
    file.put(1); file.put(0);
    
    size_t header_len_pos = file.tellp();
    file.write("\x00\x00", 2);
    
    for (const auto& var : variants_) {
        for (uint32_t i = 0; i < num_samples_; ++i) {
            uint8_t code = var.genotypes[i];
            float dosage = (code == 3) ? std::numeric_limits<float>::quiet_NaN() : static_cast<float>(code);
            file.write(reinterpret_cast<const char*>(&dosage), 4);
        }
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
    
    std::ofstream var_out(var_info_file);
    var_out << "ID CHROM POS REF ALT\n";
    for (const auto& var : variants_) {
        var_out << var.id << " " << var.chrom << " " << var.pos << " " << var.allele1 << " " << var.allele2 << "\n";
    }
    
    std::ofstream sample_out(sample_info_file);
    sample_out << "FID IID\n";
    for (const auto& sid : sample_ids_) {
        size_t sep = sid.find(':');
        if (sep != std::string::npos) {
            sample_out << sid.substr(0, sep) << " " << sid.substr(sep + 1) << "\n";
        } else {
            sample_out << sid << " " << sid << "\n";
        }
    }
}

void HapMapReader::exportToPedMap(const std::string& ped_file, const std::string& map_file) const {
    std::ofstream map_out(map_file);
    if (!map_out) throw std::runtime_error("Cannot open output MAP file: " + map_file);
    
    for (const auto& var : variants_) {
        map_out << var.chrom << " " << var.id << " 0 " << var.pos << " " << var.allele1 << " " << var.allele2 << "\n";
    }
    
    std::ofstream ped_out(ped_file);
    if (!ped_out) throw std::runtime_error("Cannot open output PED file: " + ped_file);
    
    for (uint32_t i = 0; i < num_samples_; ++i) {
        ped_out << sample_ids_[i] << " " << sample_ids_[i] << " 0 0 1 -9";
    }
    ped_out << "\n";
    
    for (const auto& var : variants_) {
        for (uint32_t i = 0; i < num_samples_; ++i) {
            uint8_t code = var.genotypes[i];
            if (code == 3) {
                ped_out << " 0 0";
            } else {
                if (code == 0) ped_out << " " << var.allele1 << " " << var.allele1;
                else if (code == 1) ped_out << " " << var.allele1 << " " << var.allele2;
                else ped_out << " " << var.allele2 << " " << var.allele2;
            }
        }
        ped_out << "\n";
    }
}

void HapMapReader::exportToVcf(const std::string& vcf_file) const {
    std::ofstream vcf_out(vcf_file);
    if (!vcf_out) throw std::runtime_error("Cannot open output VCF file: " + vcf_file);
    
    vcf_out << "##fileformat=VCFv4.2\n";
    vcf_out << "##source=CosmicScanner\n";
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& sid : sample_ids_) {
        vcf_out << "\t" << sid;
    }
    vcf_out << "\n";
    
    for (const auto& var : variants_) {
        vcf_out << var.chrom << "\t" << var.pos << "\t" << var.id << "\t"
                << var.allele1 << "\t" << var.allele2 << "\t.\tPASS\t.\tGT";
        
        for (uint32_t i = 0; i < num_samples_; ++i) {
            uint8_t code = var.genotypes[i];
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
