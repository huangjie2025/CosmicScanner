#include "npy_reader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cctype>
#include <numeric>

namespace cosmic {

bool NpyReader::parseNpyHeader(std::ifstream& file, std::vector<size_t>& shape, bool& fortran_order, std::string& dtype) {
    char magic[6];
    file.read(magic, 6);
    if (!file || std::strncmp(magic, "\x93NUMPY", 6) != 0) {
        return false;
    }
    
    uint8_t major, minor;
    file.read(reinterpret_cast<char*>(&major), 1);
    file.read(reinterpret_cast<char*>(&minor), 1);
    
    uint16_t header_len;
    file.read(reinterpret_cast<char*>(&header_len), 2);
    
    std::string header(header_len, '\0');
    file.read(&header[0], header_len);
    
    size_t shape_start = header.find('(');
    size_t shape_end = header.find(')');
    if (shape_start == std::string::npos || shape_end == std::string::npos) {
        return false;
    }
    std::string shape_str = header.substr(shape_start + 1, shape_end - shape_start - 1);
    std::stringstream ss(shape_str);
    std::string dim;
    while (std::getline(ss, dim, ',')) {
        if (!dim.empty()) {
            shape.push_back(std::stoull(dim));
        }
    }
    
    fortran_order = header.find("'fortran_order': True") != std::string::npos;
    
    size_t dtype_start = header.find("'descr':");
    if (dtype_start != std::string::npos) {
        size_t dtype_quote1 = header.find('\'', dtype_start + 8);
        size_t dtype_quote2 = header.find('\'', dtype_quote1 + 1);
        if (dtype_quote1 != std::string::npos && dtype_quote2 != std::string::npos) {
            dtype = header.substr(dtype_quote1 + 1, dtype_quote2 - dtype_quote1 - 1);
        }
    }
    
    return true;
}

void NpyReader::loadVarInfo(const std::string& file) {
    if (file.empty()) return;
    std::ifstream ifs(file);
    if (!ifs) {
        std::cerr << "Warning: Cannot open variant info file: " << file << "\n";
        return;
    }
    std::string line;
    std::getline(ifs, line); // skip header if present
    if (line.find("ID") != std::string::npos && line.find("CHROM") != std::string::npos) {
        // has header
    } else {
        ifs.seekg(0);
    }
    
    variant_ids_.clear();
    chrom_labels_.clear();
    ref_alleles_.clear();
    alt_alleles_.clear();
    positions_.clear();
    
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string id, chrom, ref, alt;
        uint32_t pos = 0;
        ss >> id >> chrom >> pos >> ref >> alt;
        if (!id.empty()) {
            variant_ids_.push_back(id);
            chrom_labels_.push_back(chrom);
            ref_alleles_.push_back(ref);
            alt_alleles_.push_back(alt);
            positions_.push_back(pos);
        }
    }
}

void NpyReader::loadSampleInfo(const std::string& file) {
    if (file.empty()) return;
    std::ifstream ifs(file);
    if (!ifs) {
        std::cerr << "Warning: Cannot open sample info file: " << file << "\n";
        return;
    }
    std::string line;
    std::getline(ifs, line); // skip header if present
    if (line.find("FID") != std::string::npos || line.find("IID") != std::string::npos) {
        // has header
    } else {
        ifs.seekg(0);
    }
    
    sample_ids_.clear();
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string fid, iid;
        ss >> fid >> iid;
        if (!iid.empty()) {
            sample_ids_.push_back(fid + ":" + iid);
        } else if (!fid.empty()) {
            sample_ids_.push_back(fid);
        }
    }
}

void NpyReader::open(const std::string& npy_file,
                     NpyOrientation orientation,
                     const std::string& var_info_file,
                     const std::string& sample_info_file) {
    npy_file_ = npy_file;
    orientation_ = orientation;
    var_info_file_ = var_info_file;
    sample_info_file_ = sample_info_file;
    
    std::ifstream file(npy_file, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open NPY file: " + npy_file);
    }
    
    std::vector<size_t> shape;
    bool fortran_order;
    std::string dtype;
    if (!parseNpyHeader(file, shape, fortran_order, dtype)) {
        throw std::runtime_error("Invalid NPY file format: " + npy_file);
    }
    
    if (shape.size() != 2) {
        throw std::runtime_error("NPY genotype matrix must be 2D, got " + std::to_string(shape.size()) + " dimensions");
    }
    
    size_t n_rows = shape[0];
    size_t n_cols = shape[1];
    
    if (orientation == NpyOrientation::VariantsBySamples) {
        num_variants_ = n_rows;
        num_samples_ = n_cols;
    } else {
        num_samples_ = n_rows;
        num_variants_ = n_cols;
    }
    
    loadVarInfo(var_info_file);
    loadSampleInfo(sample_info_file);
    
    if (variant_ids_.empty()) {
        variant_ids_.resize(num_variants_);
        for (size_t i = 0; i < num_variants_; ++i) {
            variant_ids_[i] = "var_" + std::to_string(i);
        }
    }
    if (chrom_labels_.empty()) {
        chrom_labels_.resize(num_variants_, "1");
    }
    if (ref_alleles_.empty()) {
        ref_alleles_.resize(num_variants_, "A");
    }
    if (alt_alleles_.empty()) {
        alt_alleles_.resize(num_variants_, "C");
    }
    if (positions_.empty()) {
        positions_.resize(num_variants_);
        for (size_t i = 0; i < num_variants_; ++i) {
            positions_[i] = i + 1;
        }
    }
    if (sample_ids_.empty()) {
        sample_ids_.resize(num_samples_);
        for (size_t i = 0; i < num_samples_; ++i) {
            sample_ids_[i] = "sample_" + std::to_string(i);
        }
    }
    
    genotype_data_.resize(num_variants_, std::vector<float>(num_samples_));
    
    size_t element_size = 4; // assume float32
    if (dtype.find('f') != std::string::npos) {
        if (dtype.find('8') != std::string::npos || dtype.find("64") != std::string::npos) {
            element_size = 8;
        }
    } else if (dtype.find('i') != std::string::npos || dtype.find('u') != std::string::npos) {
        if (dtype.find('8') != std::string::npos || dtype.find("64") != std::string::npos) {
            element_size = 8;
        } else if (dtype.find('2') != std::string::npos || dtype.find("16") != std::string::npos) {
            element_size = 2;
        } else {
            element_size = 4;
        }
    }
    
    std::vector<char> buffer(element_size);
    for (size_t i = 0; i < num_variants_; ++i) {
        for (size_t j = 0; j < num_samples_; ++j) {
            file.read(buffer.data(), element_size);
            if (!file) {
                throw std::runtime_error("Unexpected end of NPY file");
            }
            float value;
            if (element_size == 4) {
                std::memcpy(&value, buffer.data(), 4);
            } else if (element_size == 8) {
                double dvalue;
                std::memcpy(&dvalue, buffer.data(), 8);
                value = static_cast<float>(dvalue);
            } else if (element_size == 2) {
                int16_t ivalue;
                std::memcpy(&ivalue, buffer.data(), 2);
                value = static_cast<float>(ivalue);
            } else {
                int32_t ivalue;
                std::memcpy(&ivalue, buffer.data(), 4);
                value = static_cast<float>(ivalue);
            }
            
            if (orientation == NpyOrientation::VariantsBySamples) {
                genotype_data_[i][j] = value;
            } else {
                genotype_data_[j][i] = value;
            }
        }
    }
    
    this->buildChromIndex();
}

void NpyReader::buildChromIndex() {
    chrom_ranges_.clear();
    if (chrom_labels_.empty()) return;
    
    std::string current_chrom = chrom_labels_[0];
    uint32_t start_idx = 0;
    for (size_t i = 1; i <= num_variants_; ++i) {
        if (i == num_variants_ || chrom_labels_[i] != current_chrom) {
            chrom_ranges_.push_back({start_idx, i});
            if (i < num_variants_) {
                current_chrom = chrom_labels_[i];
                start_idx = i;
            }
        }
    }
}

void NpyReader::close() {
    genotype_data_.clear();
    variant_ids_.clear();
    sample_ids_.clear();
    chrom_labels_.clear();
    ref_alleles_.clear();
    alt_alleles_.clear();
    positions_.clear();
    snp_means_.clear();
    snp_inv_stds_.clear();
    num_samples_ = 0;
    num_variants_ = 0;
    stats_computed_ = false;
    current_variant_idx_ = 0;
}

bool NpyReader::readVariant(std::string& chrom, std::string& id, uint32_t& pos,
                            std::string& ref, std::string& alt,
                            Eigen::VectorXf& dosage) const {
    if (current_variant_idx_ >= num_variants_) {
        return false;
    }
    
    size_t idx = current_variant_idx_++;
    id = variant_ids_[idx];
    chrom = chrom_labels_[idx];
    pos = positions_[idx];
    ref = ref_alleles_[idx];
    alt = alt_alleles_[idx];
    
    dosage.resize(num_samples_);
    for (size_t i = 0; i < num_samples_; ++i) {
        dosage[i] = genotype_data_[idx][i];
    }
    
    return true;
}

void NpyReader::computeStats(int n_threads) {
    if (stats_computed_) return;
    
    snp_means_.resize(num_variants_);
    snp_inv_stds_.resize(num_variants_);
    
    #pragma omp parallel for num_threads(n_threads)
    for (size_t i = 0; i < num_variants_; ++i) {
        double sum = 0.0;
        double sum_sq = 0.0;
        int count = 0;
        for (size_t j = 0; j < num_samples_; ++j) {
            float val = genotype_data_[i][j];
            if (std::isfinite(val)) {
                sum += val;
                sum_sq += val * val;
                count++;
            }
        }
        if (count > 0) {
            double mean = sum / count;
            double var = (sum_sq / count) - mean * mean;
            snp_means_[i] = mean;
            snp_inv_stds_[i] = (var > 1e-12) ? 1.0 / std::sqrt(var) : 0.0;
        } else {
            snp_means_[i] = 0.0;
            snp_inv_stds_[i] = 0.0;
        }
    }
    stats_computed_ = true;
}

void NpyReader::multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (!stats_computed_) computeStats(n_threads);
    
    out.resize(num_samples_);
    out.setZero();
    
    #pragma omp parallel for num_threads(n_threads)
    for (size_t i = 0; i < num_samples_; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < num_variants_; ++j) {
            float val = genotype_data_[j][i];
            if (std::isfinite(val)) {
                sum += snp_inv_stds_[j] * (val - snp_means_[j]) * v[j];
            }
        }
        out[i] = sum;
    }
}

void NpyReader::multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (!stats_computed_) computeStats(n_threads);
    
    out.resize(num_variants_);
    out.setZero();
    
    #pragma omp parallel for num_threads(n_threads)
    for (size_t j = 0; j < num_variants_; ++j) {
        double sum = 0.0;
        for (size_t i = 0; i < num_samples_; ++i) {
            float val = genotype_data_[j][i];
            if (std::isfinite(val)) {
                sum += snp_inv_stds_[j] * (val - snp_means_[j]) * v[i];
            }
        }
        out[j] = sum;
    }
}

void NpyReader::seekVariant(uint32_t idx) const {
    current_variant_idx_ = idx;
}

std::pair<uint32_t, uint32_t> NpyReader::getChromRange(const std::string& chrom) const {
    for (const auto& range : chrom_ranges_) {
        if (range.first < chrom_labels_.size() && chrom_labels_[range.first] == chrom) {
            return range;
        }
    }
    return {0, 0};
}

void NpyReader::exportToNpy(const std::string& npy_file,
                            const std::string& var_info_file,
                            const std::string& sample_info_file) const {
    std::ofstream file(npy_file, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open output NPY file: " + npy_file);
    }
    
    file.write("\x93NUMPY", 6);
    file.put(1); // major version
    file.put(0); // minor version
    
    size_t n_rows = (orientation_ == NpyOrientation::VariantsBySamples) ? num_variants_ : num_samples_;
    size_t n_cols = (orientation_ == NpyOrientation::VariantsBySamples) ? num_samples_ : num_variants_;
    
    std::ostringstream header;
    header << "{'descr': '<f4', 'fortran_order': False, 'shape': (" << n_rows << ", " << n_cols << "), }";
    std::string header_str = header.str();
    size_t header_len = header_str.size();
    
    while ((header_len + 10) % 64 != 0) {
        header_str += " ";
        header_len = header_str.size();
    }
    
    uint16_t header_len_le = static_cast<uint16_t>(header_len);
    file.write(reinterpret_cast<char*>(&header_len_le), 2);
    file.write(header_str.c_str(), header_len);
    
    for (size_t i = 0; i < num_variants_; ++i) {
        for (size_t j = 0; j < num_samples_; ++j) {
            float val = genotype_data_[i][j];
            file.write(reinterpret_cast<const char*>(&val), 4);
        }
    }
    
    if (!var_info_file.empty()) {
        std::ofstream var_out(var_info_file);
        var_out << "ID CHROM POS REF ALT\n";
        for (size_t i = 0; i < num_variants_; ++i) {
            var_out << variant_ids_[i] << " " << chrom_labels_[i] << " "
                    << positions_[i] << " " << ref_alleles_[i] << " " << alt_alleles_[i] << "\n";
        }
    }
    
    if (!sample_info_file.empty()) {
        std::ofstream sample_out(sample_info_file);
        sample_out << "FID IID\n";
        for (const auto& id : sample_ids_) {
            size_t sep = id.find(':');
            if (sep != std::string::npos) {
                sample_out << id.substr(0, sep) << " " << id.substr(sep + 1) << "\n";
            } else {
                sample_out << id << " " << id << "\n";
            }
        }
    }
}

} // namespace cosmic

namespace cosmic {

void NpyReader::exportToText(const std::string& geno_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream geno_out(geno_file);
    if (!geno_out) throw std::runtime_error("Cannot open output geno file: " + geno_file);
    
    geno_out << "VariantID";
    for (const auto& sid : sample_ids_) {
        geno_out << " " << sid;
    }
    geno_out << "\n";
    
    for (size_t i = 0; i < num_variants_; ++i) {
        geno_out << variant_ids_[i];
        for (size_t j = 0; j < num_samples_; ++j) {
            float val = genotype_data_[i][j];
            if (std::isnan(val)) geno_out << " NA";
            else geno_out << " " << val;
        }
        geno_out << "\n";
    }
    
    std::ofstream var_out(var_info_file);
    var_out << "ID CHROM POS REF ALT\n";
    for (size_t i = 0; i < num_variants_; ++i) {
        var_out << variant_ids_[i] << " " << chrom_labels_[i] << " "
                << positions_[i] << " " << ref_alleles_[i] << " " << alt_alleles_[i] << "\n";
    }
    
    std::ofstream sample_out(sample_info_file);
    sample_out << "FID IID\n";
    for (const auto& id : sample_ids_) {
        size_t sep = id.find(':');
        if (sep != std::string::npos) {
            sample_out << id.substr(0, sep) << " " << id.substr(sep + 1) << "\n";
        } else {
            sample_out << id << " " << id << "\n";
        }
    }
}

void NpyReader::exportToPedMap(const std::string& ped_file, const std::string& map_file) const {
    std::ofstream map_out(map_file);
    if (!map_out) throw std::runtime_error("Cannot open output MAP file: " + map_file);
    
    for (size_t i = 0; i < num_variants_; ++i) {
        map_out << chrom_labels_[i] << " " << variant_ids_[i] << " 0 " << positions_[i] << " " << ref_alleles_[i] << " " << alt_alleles_[i] << "\n";
    }
    
    std::ofstream ped_out(ped_file);
    if (!ped_out) throw std::runtime_error("Cannot open output PED file: " + ped_file);
    
    for (size_t j = 0; j < num_samples_; ++j) {
        ped_out << sample_ids_[j] << " " << sample_ids_[j] << " 0 0 1 -9";
        for (size_t i = 0; i < num_variants_; ++i) {
            float val = genotype_data_[i][j];
            if (std::isnan(val)) {
                ped_out << " 0 0";
            } else {
                if (val < 0.5) ped_out << " 0 0";
                else if (val < 1.5) ped_out << " 1 0";
                else ped_out << " 1 1";
            }
        }
        ped_out << "\n";
    }
}

void NpyReader::exportToVcf(const std::string& vcf_file) const {
    std::ofstream vcf_out(vcf_file);
    if (!vcf_out) throw std::runtime_error("Cannot open output VCF file: " + vcf_file);
    
    vcf_out << "##fileformat=VCFv4.2\n";
    vcf_out << "##source=CosmicScanner\n";
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& sid : sample_ids_) {
        vcf_out << "\t" << sid;
    }
    vcf_out << "\n";
    
    for (size_t i = 0; i < num_variants_; ++i) {
        vcf_out << chrom_labels_[i] << "\t" << positions_[i] << "\t" << variant_ids_[i] << "\t"
                << ref_alleles_[i] << "\t" << alt_alleles_[i] << "\t.\tPASS\t.\tGT";
        
        for (size_t j = 0; j < num_samples_; ++j) {
            float val = genotype_data_[i][j];
            if (std::isnan(val)) {
                vcf_out << "\t./.";
            } else {
                if (val < 0.5) vcf_out << "\t0/0";
                else if (val < 1.5) vcf_out << "\t0/1";
                else vcf_out << "\t1/1";
            }
        }
        vcf_out << "\n";
    }
}

} // namespace cosmic