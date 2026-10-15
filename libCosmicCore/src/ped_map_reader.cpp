#include "ped_map_reader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <numeric>

namespace cosmic {

void PedMapReader::parseMapFile() {
    std::ifstream file(map_file_);
    if (!file) {
        throw std::runtime_error("Cannot open MAP file: " + map_file_);
    }
    
    std::string line;
    variants_.clear();
    chrom_labels_.clear();
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        MapVariant var;
        ss >> var.chrom >> var.id >> var.cm_pos >> var.bp_pos >> var.allele1 >> var.allele2;
        if (!var.id.empty()) {
            variants_.push_back(var);
            if (std::find(chrom_labels_.begin(), chrom_labels_.end(), var.chrom) == chrom_labels_.end()) {
                chrom_labels_.push_back(var.chrom);
            }
        }
    }
    
    num_variants_ = variants_.size();
}

void PedMapReader::parsePedFile() {
    std::ifstream file(ped_file_);
    if (!file) {
        throw std::runtime_error("Cannot open PED file: " + ped_file_);
    }
    
    std::string line;
    samples_.clear();
    sample_ids_.clear();
    
    std::vector<std::vector<std::string>> all_genotypes;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        PedSample sample;
        ss >> sample.fid >> sample.iid >> sample.father >> sample.mother >> sample.sex >> sample.phenotype;
        
        if (sample.iid.empty()) continue;
        
        sample_ids_.push_back(sample.fid + ":" + sample.iid);
        samples_.push_back(sample);
        
        std::vector<std::string> sample_genotypes;
        std::string allele;
        while (ss >> allele) {
            sample_genotypes.push_back(allele);
        }
        all_genotypes.push_back(sample_genotypes);
    }
    
    num_samples_ = samples_.size();
    
    if (num_variants_ == 0) {
        throw std::runtime_error("MAP file must be parsed before PED file");
    }
    
    // Verify genotype count
    for (size_t i = 0; i < all_genotypes.size(); ++i) {
        if (all_genotypes[i].size() != num_variants_ * 2) {
            std::cerr << "Warning: Sample " << sample_ids_[i] 
                      << " has " << all_genotypes[i].size() 
                      << " alleles, expected " << (num_variants_ * 2) << "\n";
        }
    }
    
    // Encode genotypes: [variant][sample]
    genotypes_.resize(num_variants_, std::vector<uint8_t>(num_samples_));
    for (size_t s = 0; s < num_samples_; ++s) {
        for (size_t v = 0; v < num_variants_; ++v) {
            const std::string& a1 = all_genotypes[s][v * 2];
            const std::string& a2 = all_genotypes[s][v * 2 + 1];
            const std::string& ref = variants_[v].allele1;
            const std::string& alt = variants_[v].allele2;
            genotypes_[v][s] = encodeGenotype(a1, a2, ref, alt);
        }
    }
}

uint8_t PedMapReader::encodeGenotype(const std::string& a1, const std::string& a2, 
                                      const std::string& ref, const std::string& alt) const {
    // Missing alleles
    if (a1 == "0" || a1 == "NA" || a1 == "N" || a1 == "." ||
        a2 == "0" || a2 == "NA" || a2 == "N" || a2 == ".") {
        return 3;
    }
    
    // Count reference alleles
    int ref_count = 0;
    if (a1 == ref) ref_count++;
    if (a2 == ref) ref_count++;
    
    if (ref_count == 2) return 0;      // hom ref
    if (ref_count == 1) return 1;      // het
    if (ref_count == 0) return 2;      // hom alt
    
    return 3; // missing
}

void PedMapReader::decodeGenotype(uint8_t code, const std::string& ref, const std::string& alt,
                                   std::string& a1, std::string& a2) const {
    switch (code) {
        case 0: a1 = ref; a2 = ref; break;
        case 1: a1 = ref; a2 = alt; break;
        case 2: a1 = alt; a2 = alt; break;
        default: a1 = "0"; a2 = "0"; break;
    }
}

void PedMapReader::buildChromIndex() {
    chrom_ranges_.clear();
    if (variants_.empty()) return;
    
    std::string current_chrom = variants_[0].chrom;
    uint32_t start_idx = 0;
    for (size_t i = 1; i <= num_variants_; ++i) {
        if (i == num_variants_ || variants_[i].chrom != current_chrom) {
            chrom_ranges_.push_back({start_idx, i});
            if (i < num_variants_) {
                current_chrom = variants_[i].chrom;
                start_idx = i;
            }
        }
    }
}

void PedMapReader::open(const std::string& ped_file, const std::string& map_file) {
    ped_file_ = ped_file;
    map_file_ = map_file;
    
    parseMapFile();
    parsePedFile();
    buildChromIndex();
    
    // Compute initial stats
    computeStats(1);
}

void PedMapReader::close() {
    samples_.clear();
    variants_.clear();
    sample_ids_.clear();
    chrom_labels_.clear();
    genotypes_.clear();
    snp_means_.clear();
    snp_inv_stds_.clear();
    num_samples_ = 0;
    num_variants_ = 0;
    stats_computed_ = false;
    current_variant_idx_ = 0;
}

bool PedMapReader::readVariant(std::string& chrom, std::string& id, uint32_t& pos,
                               std::string& ref, std::string& alt,
                               Eigen::VectorXf& dosage) const {
    if (current_variant_idx_ >= num_variants_) {
        return false;
    }
    
    size_t idx = current_variant_idx_++;
    const auto& var = variants_[idx];
    
    id = var.id;
    chrom = var.chrom;
    pos = var.bp_pos;
    ref = var.allele1;
    alt = var.allele2;
    
    dosage.resize(num_samples_);
    for (size_t i = 0; i < num_samples_; ++i) {
        uint8_t code = genotypes_[idx][i];
        if (code == 3) {
            dosage[i] = std::numeric_limits<float>::quiet_NaN();
        } else {
            dosage[i] = static_cast<float>(code); // 0, 1, 2
        }
    }
    
    return true;
}

void PedMapReader::computeStats(int n_threads) {
    if (stats_computed_) return;
    
    snp_means_.resize(num_variants_);
    snp_inv_stds_.resize(num_variants_);
    
    #pragma omp parallel for num_threads(n_threads)
    for (size_t i = 0; i < num_variants_; ++i) {
        double sum = 0.0;
        double sum_sq = 0.0;
        int count = 0;
        for (size_t j = 0; j < num_samples_; ++j) {
            uint8_t code = genotypes_[i][j];
            if (code != 3) {
                sum += code;
                sum_sq += static_cast<double>(code) * code;
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

void PedMapReader::multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (!stats_computed_) computeStats(n_threads);
    
    out.resize(num_samples_);
    out.setZero();
    
    #pragma omp parallel for num_threads(n_threads)
    for (size_t i = 0; i < num_samples_; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < num_variants_; ++j) {
            uint8_t code = genotypes_[j][i];
            if (code != 3) {
                sum += snp_inv_stds_[j] * (static_cast<double>(code) - snp_means_[j]) * v[j];
            }
        }
        out[i] = sum;
    }
}

void PedMapReader::multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (!stats_computed_) computeStats(n_threads);
    
    out.resize(num_variants_);
    out.setZero();
    
    #pragma omp parallel for num_threads(n_threads)
    for (size_t j = 0; j < num_variants_; ++j) {
        double sum = 0.0;
        for (size_t i = 0; i < num_samples_; ++i) {
            uint8_t code = genotypes_[j][i];
            if (code != 3) {
                sum += snp_inv_stds_[j] * (static_cast<double>(code) - snp_means_[j]) * v[i];
            }
        }
        out[j] = sum;
    }
}

void PedMapReader::seekVariant(uint32_t idx) const {
    current_variant_idx_ = idx;
}

std::pair<uint32_t, uint32_t> PedMapReader::getChromRange(const std::string& chrom) const {
    for (const auto& range : chrom_ranges_) {
        if (range.first < variants_.size() && variants_[range.first].chrom == chrom) {
            return range;
        }
    }
    return {0, 0};
}

void PedMapReader::exportToPedMap(const std::string& ped_file, const std::string& map_file) const {
    // Export MAP
    std::ofstream map_out(map_file);
    if (!map_out) {
        throw std::runtime_error("Cannot open output MAP file: " + map_file);
    }
    for (const auto& var : variants_) {
        map_out << var.chrom << " " << var.id << " " << var.cm_pos << " "
                << var.bp_pos << " " << var.allele1 << " " << var.allele2 << "\n";
    }
    
    // Export PED
    std::ofstream ped_out(ped_file);
    if (!ped_out) {
        throw std::runtime_error("Cannot open output PED file: " + ped_file);
    }
    for (size_t s = 0; s < num_samples_; ++s) {
        const auto& sample = samples_[s];
        ped_out << sample.fid << " " << sample.iid << " " << sample.father << " "
                << sample.mother << " " << sample.sex << " " << sample.phenotype;
        for (size_t v = 0; v < num_variants_; ++v) {
            std::string a1, a2;
            decodeGenotype(genotypes_[v][s], variants_[v].allele1, variants_[v].allele2, a1, a2);
            ped_out << " " << a1 << " " << a2;
        }
        ped_out << "\n";
    }
}

} // namespace cosmic