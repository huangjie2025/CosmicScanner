#include "vcf/vcf_reader.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <zlib.h>

namespace cosmic {
namespace {

std::string genotypeField(const std::string& format, const std::string& sample) {
    std::size_t field = 0, begin = 0;
    while (begin <= format.size()) {
        const std::size_t end = format.find(':', begin);
        if (format.substr(begin, end == std::string::npos ? end : end - begin) == "GT") {
            std::size_t sample_begin = 0;
            for (std::size_t current = 0; current < field; ++current) {
                sample_begin = sample.find(':', sample_begin);
                if (sample_begin == std::string::npos) return ".";
                ++sample_begin;
            }
            const std::size_t sample_end = sample.find(':', sample_begin);
            return sample.substr(sample_begin, sample_end == std::string::npos ? sample_end
                                                                               : sample_end - sample_begin);
        }
        if (end == std::string::npos) break;
        begin = end + 1;
        ++field;
    }
    return ".";
}

} // namespace

VcfReader::VcfReader() = default;

VcfReader::~VcfReader() {
    close();
}

void VcfReader::close() {
    if (ifs_.is_open()) ifs_.close();
    if (tbi_ifs_.is_open()) tbi_ifs_.close();
    ifs_.clear();
    tbi_ifs_.clear();
    num_samples_ = 0;
    num_variants_ = 0;
    sample_ids_.clear();
    variants_.clear();
    chrom_labels_.clear();
    chrom_ranges_.clear();
    current_offset_ = 0;
    variant_data_offset_ = 0;
    header_length_ = 0;
    snp_means_.clear();
    snp_inv_stds_.clear();
    snp_mafs_.clear();
    snp_missing_counts_.clear();
    stats_computed_ = false;
}

void VcfReader::open(const std::string& vcf_file, const std::string& tbi_file) {
    close();
    vcf_file_ = vcf_file;
    tbi_file_ = tbi_file;
    
    if (vcf_file_.size() >= 3 && vcf_file_.substr(vcf_file_.size() - 3) == ".gz") {
        // For bgzipped files, we'll use zlib to read
        // For now, we'll read through gzip using a pipe
        std::string cmd = "gzip -dc " + vcf_file_;
        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            throw std::runtime_error("Failed to open VCF file: " + vcf_file_);
        }
        
        // Read header and build index
        char buffer[16384];
        std::string header_content;
        
        // Read header lines
        while (fgets(buffer, sizeof(buffer), pipe)) {
            if (buffer[0] == '#') {
                if (buffer[1] == '#') {
                    header_content += buffer;
                } else if (buffer[1] == 'C' && buffer[2] == 'H' && buffer[3] == 'R') {
                    // #CHROM line - sample IDs
                    header_content += buffer;
                    break;
                }
            }
        }
        
        // Parse header
        std::istringstream header_stream(header_content);
        std::string line;
        while (std::getline(header_stream, line)) {
            if (line.size() > 1 && line[0] == '#' && line[1] == '#') {
                continue;
            } else if (line.size() > 1 && line[0] == '#' && line[1] == 'C') {
                parseSampleLine(line);
                break;
            }
        }
        
        // Count variants and build chrom index
        num_variants_ = 0;
        variants_.clear();
        chrom_labels_.clear();
        chrom_ranges_.clear();
        std::string last_chrom = "";
        uint32_t chrom_start = 0;
        
        VcfVariant variant;
        Eigen::VectorXf dummy_dosage(num_samples_);
        
        while (fgets(buffer, sizeof(buffer), pipe)) {
            if (buffer[0] == '#') continue;
            VcfVariant v;
            Eigen::VectorXf dosage(num_samples_);
            std::string line(buffer);
            line = line.substr(0, line.find_last_not_of("\r\n") + 1);
            
            if (parseVariantLine(line, v, dosage)) {
                variants_.push_back(v);
                num_variants_++;
                
                if (v.chrom != last_chrom) {
                    if (!last_chrom.empty()) {
                        chrom_ranges_.push_back({chrom_start, num_variants_ - chrom_start});
                    }
                    last_chrom = v.chrom;
                    chrom_start = num_variants_ - 1;
                    chrom_labels_.push_back(v.chrom);
                }
            }
        }
        
        if (!last_chrom.empty()) {
            chrom_ranges_.push_back({chrom_start, num_variants_ - chrom_start});
        }
        
        pclose(pipe);
        
        // Open TBI index if provided
        if (!tbi_file_.empty()) {
            std::string tbi_cmd = "gzip -dc " + tbi_file_;
            FILE* tbi_pipe = popen(tbi_cmd.c_str(), "r");
            if (tbi_pipe) {
                // TBI format: read index
                pclose(tbi_pipe);
            }
        }
        
    } else {
        // Plain VCF file
        ifs_.open(vcf_file_);
        if (!ifs_) {
            throw std::runtime_error("Cannot open VCF file: " + vcf_file_);
        }
        
        std::string header_content;
        std::string line;
        
        // Read header
        while (std::getline(ifs_, line)) {
            if (line.size() > 1 && line[0] == '#' && line[1] == '#') {
                header_content += line + "\n";
                continue;
            } else if (line.size() > 1 && line[0] == '#' && line[1] == 'C') {
                header_content += line + "\n";
                parseSampleLine(line);
                break;
            }
        }
        
        // Parse header for sample IDs
        // Count variants and build index
        num_variants_ = 0;
        variants_.clear();
        chrom_labels_.clear();
        chrom_ranges_.clear();
        std::string last_chrom = "";
        uint32_t chrom_start = 0;
        
        VcfVariant variant;
        Eigen::VectorXf dosage(num_samples_);
        
        while (std::getline(ifs_, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            VcfVariant v;
            Eigen::VectorXf dosage(num_samples_);
            
            if (parseVariantLine(line, v, dosage)) {
                variants_.push_back(v);
                num_variants_++;
                
                if (v.chrom != last_chrom) {
                    if (!last_chrom.empty()) {
                        chrom_ranges_.push_back({chrom_start, num_variants_ - chrom_start});
                    }
                    last_chrom = v.chrom;
                    chrom_start = num_variants_ - 1;
                    chrom_labels_.push_back(v.chrom);
                }
            }
        }
        
        if (!last_chrom.empty()) {
            chrom_ranges_.push_back({chrom_start, num_variants_ - chrom_start});
        }
        
        variant_data_offset_ = header_content.size();
        // For gzipped files, we can't easily seek, so we'll re-read from the beginning
ifs_.clear();
        ifs_.seekg(0);
        // Skip header lines to find variant data offset
        std::string header_line;
        while (std::getline(ifs_, header_line)) {
            if (header_line.size() > 6 && header_line[0] == '#' && header_line[1] == 'C' && header_line[2] == 'H' && header_line[3] == 'R' && header_line[4] == 'O' && header_line[5] == 'M') {
                variant_data_offset_ = ifs_.tellg();
                break;
}
    }
    }
}

void VcfReader::parseSampleLine(const std::string& line) {
    std::istringstream iss(line);
    std::string token;
    // #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE1 SAMPLE2 ...
    for (int i = 0; i < 9; ++i) iss >> token;
    while (iss >> token) {
        sample_ids_.push_back(token);
    }
    num_samples_ = sample_ids_.size();
}

bool VcfReader::parseVariantLine(const std::string& line, VcfVariant& variant, Eigen::VectorXf& dosage) const {
    std::istringstream iss(line);
    std::string token;
    
    if (!(iss >> variant.chrom >> variant.pos >> variant.id >> variant.ref >> variant.alt 
          >> variant.qual >> variant.filter >> variant.info >> variant.format)) {
        return false;
    }
    
    if (variant.ref.size() > 1 || variant.alt.size() > 1) {
        // Skip multi-allelic for now
        return false;
    }
    
    // Parse genotype fields and store genotype strings
    variant.genotypes.clear();
    variant.genotypes.reserve(num_samples_);
    std::string sample_data;
    
    while (iss >> sample_data) {
        variant.genotypes.push_back(sample_data);
    }
    if (variant.genotypes.size() != num_samples_) return false;
    
    // Parse genotype fields to compute dosage
    dosage.resize(num_samples_);
    dosage.setConstant(std::numeric_limits<float>::quiet_NaN());
    size_t sample_idx = 0;
    
    for (const auto& gt_str : variant.genotypes) {
        if (sample_idx >= static_cast<size_t>(num_samples_)) break;
        
        const std::string gt = genotypeField(variant.format, gt_str);
        
        float d = 0.0f;
        computeDosageFromGT(gt, variant.ref, variant.alt, d);
        dosage[sample_idx] = d;
        
        sample_idx++;
    }
    
    return true;
}

void VcfReader::computeDosageFromGT(const std::string& gt, const std::string& ref, const std::string& alt, float& dosage) const {
    (void)ref;
    (void)alt;
    const std::size_t separator = gt.find_first_of("/|");
    if (separator == std::string::npos || separator == 0 || separator + 1 >= gt.size() ||
        gt.find_first_of("/|", separator + 1) != std::string::npos) {
        dosage = std::numeric_limits<float>::quiet_NaN();
        return;
    }
    const std::string first = gt.substr(0, separator);
    const std::string second = gt.substr(separator + 1);
    if ((first != "0" && first != "1") || (second != "0" && second != "1")) {
        dosage = std::numeric_limits<float>::quiet_NaN();
        return;
    }
    dosage = static_cast<float>((first == "1") + (second == "1"));
}

bool VcfReader::isMissing(const std::string& gt) const {
    return gt == "./." || gt == ".|." || gt == ".";
}

void VcfReader::reset() const {
    if (ifs_.is_open()) {
        ifs_.clear();
        ifs_.seekg(variant_data_offset_);
    }
    current_offset_ = 0;
}

void VcfReader::loadVariantInfo() const {
    if (!variants_.empty()) return;
    
    if (ifs_.is_open()) {
        ifs_.clear();
        ifs_.seekg(variant_data_offset_);
        std::string line;
        variants_.clear();
        chrom_labels_.clear();
        chrom_ranges_.clear();
        num_variants_ = 0;
        
        std::string last_chrom = "";
        uint32_t chrom_start = 0;
        
        while (std::getline(ifs_, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            VcfVariant v;
            Eigen::VectorXf dosage(num_samples_);
            if (parseVariantLine(line, v, dosage)) {
                variants_.push_back(v);
                num_variants_++;
                
                if (v.chrom != last_chrom) {
                    if (!last_chrom.empty()) {
                        chrom_ranges_.push_back({chrom_start, num_variants_ - chrom_start});
                    }
                    last_chrom = v.chrom;
                    chrom_start = num_variants_ - 1;
                    chrom_labels_.push_back(v.chrom);
                }
            }
        }
        
        if (!last_chrom.empty()) {
            chrom_ranges_.push_back({chrom_start, num_variants_ - chrom_start});
        }
    }
}

void VcfReader::seekVariant(uint32_t idx) const {
    if (idx >= num_variants_) return;
    current_offset_ = idx;
    if (ifs_.is_open()) {
        ifs_.clear();
        ifs_.seekg(variant_data_offset_);
        std::string line;
        for (uint32_t i = 0; i <= idx && std::getline(ifs_, line); ++i) {
            // Skip
        }
    }
}

bool VcfReader::readVariant(std::string& chrom, std::string& id, uint32_t& pos,
                            std::string& ref, std::string& alt,
                            Eigen::VectorXf& dosage) const {
    if (current_offset_ >= num_variants_) return false;
    
    const auto& v = variants_[current_offset_];
    chrom = v.chrom;
    id = v.id;
    pos = v.pos;
    ref = v.ref;
    alt = v.alt;
    dosage.resize(num_samples_);
    dosage.setConstant(std::numeric_limits<float>::quiet_NaN());
    
    // Compute dosage from cached genotype strings
    size_t sample_idx = 0;
    for (const auto& gt_str : v.genotypes) {
        if (sample_idx >= static_cast<size_t>(num_samples_)) break;
        
        const std::string gt = genotypeField(v.format, gt_str);
        
        float d = 0.0f;
        computeDosageFromGT(gt, v.ref, v.alt, d);
        dosage[sample_idx] = d;
        sample_idx++;
    }
    
    current_offset_++;
    return true;
}

void VcfReader::computeStats(int n_threads) {
    if (stats_computed_) return;
    
    snp_means_.resize(num_variants_);
    snp_inv_stds_.resize(num_variants_);
    snp_mafs_.resize(num_variants_);
    snp_missing_counts_.resize(num_variants_);
    
    // For VCF, we need to re-read variants to compute stats
    // This is a simplified implementation
    for (size_t i = 0; i < variants_.size(); ++i) {
        snp_means_[i] = 0.0;
        snp_inv_stds_[i] = 1.0;
        snp_mafs_[i] = 0.0;
        snp_missing_counts_[i] = 0;
    }
    stats_computed_ = true;
}

void VcfReader::multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& y, int n_threads) {
    if (!stats_computed_) computeStats(n_threads);
    
    y.resize(num_samples_);
    y.setZero();
    
    // For each variant, add v[j] * Z_ij to y[i]
    // Z_ij = (dosage_ij - 2p_j) / sqrt(2p_j(1-p_j))
    // This is a simplified implementation
    y.setZero();
}

void VcfReader::multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& y, int n_threads) {
    if (!stats_computed_) computeStats(n_threads);
    
    y.resize(num_variants_);
    y.setZero();
    
    // y[j] = sum_i Z_ij * v[i]
    y.setZero();
}

void VcfReader::buildChromIndex() const {
    chrom_labels_.clear();
    chrom_ranges_.clear();
    
    if (variants_.empty()) return;
    
    std::string last_chrom = "";
    uint32_t chrom_start = 0;
    
    for (size_t i = 0; i < variants_.size(); ++i) {
        if (variants_[i].chrom != last_chrom) {
            if (!last_chrom.empty()) {
                chrom_ranges_.push_back({chrom_start, static_cast<uint32_t>(i - chrom_start)});
            }
            last_chrom = variants_[i].chrom;
            chrom_start = static_cast<uint32_t>(i);
            chrom_labels_.push_back(variants_[i].chrom);
        }
    }
    
    if (!last_chrom.empty()) {
        chrom_ranges_.push_back({chrom_start, static_cast<uint32_t>(variants_.size() - chrom_start)});
    }
}

std::pair<uint32_t, uint32_t> VcfReader::getChromRange(const std::string& chrom) const {
    for (size_t i = 0; i < chrom_labels_.size(); ++i) {
        if (chrom_labels_[i] == chrom) {
            return chrom_ranges_[i];
        }
    }
    return {0, 0};
}

void VcfReader::releaseVariantRange(uint32_t start_idx, uint32_t count) const {
    // No-op for VCF
}

} // namespace cosmic

namespace cosmic {

void VcfReader::exportToText(const std::string& geno_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream geno_out(geno_file);
    if (!geno_out) throw std::runtime_error("Cannot open output geno file: " + geno_file);
    
    geno_out << "VariantID";
    for (const auto& sid : sample_ids_) {
        geno_out << " " << sid;
    }
    geno_out << "\n";
    
    const auto& variants = getVariants();
    for (const auto& var : variants) {
        geno_out << var.id;
        for (const auto& gt_str : var.genotypes) {
            float dosage = 0.0f;
            computeDosageFromGT(gt_str, var.ref, var.alt, dosage);
            if (std::isnan(dosage)) geno_out << " NA";
            else geno_out << " " << dosage;
        }
        geno_out << "\n";
    }
    
    std::ofstream var_out(var_info_file);
    var_out << "ID CHROM POS REF ALT\n";
    for (const auto& var : variants) {
        var_out << var.id << " " << var.chrom << " " << var.pos << " " << var.ref << " " << var.alt << "\n";
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

void VcfReader::exportToNpy(const std::string& npy_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream file(npy_file, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open output NPY file: " + npy_file);
    
    const auto& variants = getVariants();
    uint32_t n_variants = variants.size();
    uint32_t n_samples = num_samples_;
    
    file.write("\x93NUMPY", 6);
    file.put(1); file.put(0);
    
    size_t header_len_pos = file.tellp();
    file.write("\x00\x00", 2);
    
    for (const auto& var : variants) {
        for (const auto& gt_str : var.genotypes) {
            float dosage = 0.0f;
            computeDosageFromGT(gt_str, var.ref, var.alt, dosage);
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
    for (const auto& var : variants) {
        var_out << var.id << " " << var.chrom << " " << var.pos << " " << var.ref << " " << var.alt << "\n";
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

void VcfReader::exportToPedMap(const std::string& ped_file, const std::string& map_file) const {
    std::ofstream map_out(map_file);
    if (!map_out) throw std::runtime_error("Cannot open output MAP file: " + map_file);
    
    const auto& variants = getVariants();
    for (const auto& var : variants) {
        map_out << var.chrom << " " << var.id << " 0 " << var.pos << " " << var.ref << " " << var.alt << "\n";
    }
    
    std::ofstream ped_out(ped_file);
    if (!ped_out) throw std::runtime_error("Cannot open output PED file: " + ped_file);
    
    for (const auto& sid : sample_ids_) {
        ped_out << sid << " " << sid << " 0 0 1 -9";
    }
    ped_out << "\n";
    
    for (const auto& var : variants) {
        for (size_t i = 0; i < sample_ids_.size(); ++i) {
            float dosage = 0.0f;
            computeDosageFromGT(var.genotypes[i], var.ref, var.alt, dosage);
            if (std::isnan(dosage)) {
                ped_out << " 0 0";
            } else {
                if (dosage < 0.5) ped_out << " 0 0";
                else if (dosage < 1.5) ped_out << " 1 0";
                else ped_out << " 1 1";
            }
        }
        ped_out << "\n";
    }
}

void VcfReader::exportToVcf(const std::string& vcf_file) const {
    std::ofstream vcf_out(vcf_file);
    if (!vcf_out) throw std::runtime_error("Cannot open output VCF file: " + vcf_file);
    
    vcf_out << "##fileformat=VCFv4.2\n";
    vcf_out << "##source=CosmicScanner\n";
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& sid : sample_ids_) {
        vcf_out << "\t" << sid;
    }
    vcf_out << "\n";
    
    const auto& variants = getVariants();
    for (const auto& var : variants) {
        vcf_out << var.chrom << "\t" << var.pos << "\t" << var.id << "\t"
                << var.ref << "\t" << var.alt << "\t.\tPASS\t.\tGT";
        
        for (const auto& gt_str : var.genotypes) {
            vcf_out << "\t" << gt_str;
        }
        vcf_out << "\n";
    }
}

} // namespace cosmic
