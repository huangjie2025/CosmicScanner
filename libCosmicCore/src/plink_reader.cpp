#include "plink_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <omp.h>

namespace cosmic {

PlinkReader::PlinkReader() {
    // PLINK binary encoding (bits 0-1, 2-3, 4-5, 6-7)
    // 00 (0) -> Homozygous 1st allele (0.0)
    // 01 (1) -> Missing (NAN)
    // 10 (2) -> Heterozygous (1.0)
    // 11 (3) -> Homozygous 2nd allele (2.0)
    
    lookup_table[0] = 2.0f; 
    lookup_table[1] = NAN;  
    lookup_table[2] = 1.0f; 
    lookup_table[3] = 0.0f; 
}

PlinkReader::~PlinkReader() {
    bed_file.close();
}

void PlinkReader::load(const std::string& prefix) {
    readFam(prefix + ".fam");
    readBim(prefix + ".bim");
    openBed(prefix + ".bed");
}

void PlinkReader::readFam(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open .fam file: " + path);
    
    fam_data.clear();
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        PlinkFamInfo info;
        ss >> info.fid >> info.iid >> info.father >> info.mother >> info.sex >> info.phenotype;
        fam_data.push_back(info);
    }
    num_samples = fam_data.size();
    if (num_samples == 0) throw std::runtime_error("No samples found in .fam file");
}

void PlinkReader::readBim(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open .bim file: " + path);
    
    bim_data.clear();
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        PlinkBimInfo info;
        ss >> info.chrom >> info.id >> info.cm_pos >> info.bp_pos >> info.alt >> info.ref;
        bim_data.push_back(info);
    }
    num_snps = bim_data.size();
}

void PlinkReader::openBed(const std::string& path) {
    bed_file.open(path);
    const uint8_t* data = static_cast<const uint8_t*>(bed_file.get_data());
    size_t size = bed_file.get_size();
    
    if (size < 3) throw std::runtime_error("Invalid .bed file (too small)");
    
    if (data[0] != 0x6C || data[1] != 0x1B) {
        throw std::runtime_error("Invalid .bed file magic number");
    }
    
    if (data[2] != 0x01) {
        throw std::runtime_error("Only SNP-major mode (0x01) .bed files are supported");
    }
    
    // Calculate bytes per snp
    // ceil(N / 4)
    bytes_per_snp = (num_samples + 3) / 4;
    
    size_t expected_size = 3 + num_snps * bytes_per_snp;
    if (size < expected_size) {
        // Warning or error? Some files might be truncated or have extra data.
        // Let's warn but proceed if possible.
        std::cerr << "Warning: .bed file size (" << size << ") < expected (" << expected_size << "). Data may be truncated." << std::endl;
    }
}

size_t PlinkReader::getSnpGenotypes(size_t snp_idx, Eigen::VectorXf& out_genotypes) const {
    if (snp_idx >= num_snps) return 0;
    out_genotypes.resize(num_samples);
    
    const uint8_t* data = static_cast<const uint8_t*>(bed_file.get_data());
    const uint8_t* snp_data = data + 3 + snp_idx * bytes_per_snp;
    
    size_t valid_count = 0;
    
    for (size_t byte_idx = 0; byte_idx < bytes_per_snp; ++byte_idx) {
        uint8_t byte = snp_data[byte_idx];
        for (int k = 0; k < 4; ++k) {
            size_t sample_idx = byte_idx * 4 + k;
            if (sample_idx >= num_samples) break;
            
            uint8_t code = (byte >> (k * 2)) & 0x03;
            float val = lookup_table[code];
            out_genotypes[sample_idx] = val;
            
            if (!std::isnan(val)) valid_count++;
        }
    }
    return valid_count;
}

void PlinkReader::getSnpBlock(size_t start_snp_idx, size_t end_snp_idx, Eigen::MatrixXf& out_matrix) const {
    if (start_snp_idx >= num_snps) return;
    if (end_snp_idx > num_snps) end_snp_idx = num_snps;
    
    size_t block_size = end_snp_idx - start_snp_idx;
    out_matrix.resize(num_samples, block_size);
    
    const uint8_t* data = static_cast<const uint8_t*>(bed_file.get_data());
    const uint8_t* snp_data = data + 3 + start_snp_idx * bytes_per_snp;
    
    #pragma omp parallel for
    for (long long j = 0; j < (long long)block_size; ++j) {
        const uint8_t* col_data = snp_data + j * bytes_per_snp;
        float* col_ptr = out_matrix.col(j).data();
        
        for (size_t i = 0; i < num_samples; i += 4) {
            size_t byte_idx = i >> 2;
            uint8_t byte = col_data[byte_idx];
            
            if (i + 3 < num_samples) {
                col_ptr[i+0] = lookup_table[(byte >> 0) & 0x03];
                col_ptr[i+1] = lookup_table[(byte >> 2) & 0x03];
                col_ptr[i+2] = lookup_table[(byte >> 4) & 0x03];
                col_ptr[i+3] = lookup_table[(byte >> 6) & 0x03];
            } else {
                for (int k = 0; k < 4 && (i + k) < num_samples; ++k) {
                    col_ptr[i+k] = lookup_table[(byte >> (k * 2)) & 0x03];
                }
            }
        }
    }
}

void PlinkReader::getSnpBlock(size_t start_snp_idx, size_t end_snp_idx, Eigen::MatrixXd& out_matrix) const {
    if (start_snp_idx >= num_snps) return;
    if (end_snp_idx > num_snps) end_snp_idx = num_snps;
    
    size_t block_size = end_snp_idx - start_snp_idx;
    out_matrix.resize(num_samples, block_size);
    
    const uint8_t* data = static_cast<const uint8_t*>(bed_file.get_data());
    const uint8_t* snp_data = data + 3 + start_snp_idx * bytes_per_snp;
    
    #pragma omp parallel for
    for (long long j = 0; j < (long long)block_size; ++j) {
        const uint8_t* col_data = snp_data + j * bytes_per_snp;
        double* col_ptr = out_matrix.col(j).data();
        
        for (size_t i = 0; i < num_samples; i += 4) {
            size_t byte_idx = i >> 2;
            uint8_t byte = col_data[byte_idx];
            
            if (i + 3 < num_samples) {
                col_ptr[i+0] = lookup_table[(byte >> 0) & 0x03];
                col_ptr[i+1] = lookup_table[(byte >> 2) & 0x03];
                col_ptr[i+2] = lookup_table[(byte >> 4) & 0x03];
                col_ptr[i+3] = lookup_table[(byte >> 6) & 0x03];
            } else {
                for (int k = 0; k < 4 && (i + k) < num_samples; ++k) {
                    col_ptr[i+k] = lookup_table[(byte >> (k * 2)) & 0x03];
                }
            }
        }
    }
}

// --- RHE / PCG Support ---

void PlinkReader::computeStats(int n_threads) {
    if (stats_computed) return;
    
    snp_means.resize(num_snps);
    snp_inv_stds.resize(num_snps);
    
    const uint8_t* data = static_cast<const uint8_t*>(bed_file.get_data()) + 3;
    
    #pragma omp parallel for num_threads(n_threads)
    for (long long j = 0; j < (long long)num_snps; ++j) {
        const uint8_t* col_data = data + j * bytes_per_snp;
        
        double sum = 0.0;
        double sum_sq = 0.0;
        int count = 0;
        
        for (size_t i = 0; i < num_samples; ++i) {
            size_t byte_idx = i / 4;
            size_t bit_idx = (i % 4) * 2;
            uint8_t code = (col_data[byte_idx] >> bit_idx) & 0x03;
            
            float val = lookup_table[code];
            if (!std::isnan(val)) {
                sum += val;
                sum_sq += val * val;
                count++;
            }
        }
        
        if (count > 0) {
            double mean = sum / count;
            // var = E[x^2] - (E[x])^2
            // Or sample variance? GRM usually uses 2p(1-p).
            // Let's use 2p(1-p) where p = mean/2.
            double p = mean / 2.0;
            double var = 2.0 * p * (1.0 - p);
            
            // Check for monomorphic
            if (var < 1e-9) {
                snp_means[j] = mean;
                snp_inv_stds[j] = 0.0; // Filtered out effectively
            } else {
                snp_means[j] = mean;
                snp_inv_stds[j] = 1.0 / std::sqrt(var);
            }
        } else {
            snp_means[j] = 0.0;
            snp_inv_stds[j] = 0.0;
        }
    }
    stats_computed = true;
    std::cout << "Computed stats for " << num_snps << " SNPs." << std::endl;
}

void PlinkReader::multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& y, int n_threads) {
    if (!stats_computed) computeStats(n_threads);
    if (v.size() != (long long)num_snps) throw std::runtime_error("Dimension mismatch in multiply_Z_v");
    
    y.setZero(num_samples);
    
    const uint8_t* data = static_cast<const uint8_t*>(bed_file.get_data()) + 3;
    
    // We parallelize over SNPs and reduce to y
    // Each thread needs a local y buffer
    
    #pragma omp parallel num_threads(n_threads)
    {
        Eigen::VectorXd y_local = Eigen::VectorXd::Zero(num_samples);
        
        #pragma omp for
        for (long long j = 0; j < (long long)num_snps; ++j) {
            double inv_std = snp_inv_stds[j];
            if (inv_std == 0.0) continue;
            
            double mean = snp_means[j];
            double weight = v(j) * inv_std;
            
            // Optimization: Precompute values for 0, 1, 2, 3
            // Code 0 (2.0): (2 - mean) * weight
            // Code 1 (NAN): 0 (mean imputed -> 0 in Z)
            // Code 2 (1.0): (1 - mean) * weight
            // Code 3 (0.0): (0 - mean) * weight
            
            double val_map[4];
            val_map[0] = (lookup_table[0] - mean) * weight;
            val_map[1] = 0.0;
            val_map[2] = (lookup_table[2] - mean) * weight;
            val_map[3] = (lookup_table[3] - mean) * weight;
            
            const uint8_t* col_data = data + j * bytes_per_snp;
            
            // Unroll loop over bytes
            for (size_t byte_idx = 0; byte_idx < bytes_per_snp; ++byte_idx) {
                uint8_t byte = col_data[byte_idx];
                
                size_t base_sample = byte_idx * 4;
                
                // Sample 0
                if (base_sample < num_samples) {
                    y_local(base_sample) += val_map[byte & 3];
                }
                // Sample 1
                if (base_sample + 1 < num_samples) {
                    y_local(base_sample + 1) += val_map[(byte >> 2) & 3];
                }
                // Sample 2
                if (base_sample + 2 < num_samples) {
                    y_local(base_sample + 2) += val_map[(byte >> 4) & 3];
                }
                // Sample 3
                if (base_sample + 3 < num_samples) {
                    y_local(base_sample + 3) += val_map[(byte >> 6) & 3];
                }
            }
        }
        
        #pragma omp critical
        {
            y += y_local;
        }
    }
}

void PlinkReader::multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& y, const std::vector<int>& keep_indices, int n_threads) {
    if (!stats_computed) computeStats(n_threads);
    if (v.size() != (long long)num_snps) throw std::runtime_error("Dimension mismatch in multiply_Z_v");
    
    int n_keep = keep_indices.size();
    y.setZero(n_keep);
    
    const uint8_t* data = static_cast<const uint8_t*>(bed_file.get_data()) + 3;
    
    #pragma omp parallel num_threads(n_threads)
    {
        Eigen::VectorXd y_local = Eigen::VectorXd::Zero(n_keep);
        
        #pragma omp for
        for (long long j = 0; j < (long long)num_snps; ++j) {
            double inv_std = snp_inv_stds[j];
            if (inv_std == 0.0) continue;
            
            double mean = snp_means[j];
            double weight = v(j) * inv_std;
            
            double val_map[4];
            val_map[0] = (lookup_table[0] - mean) * weight;
            val_map[1] = 0.0;
            val_map[2] = (lookup_table[2] - mean) * weight;
            val_map[3] = (lookup_table[3] - mean) * weight;
            
            const uint8_t* col_data = data + j * bytes_per_snp;
            
            for (int k = 0; k < n_keep; ++k) {
                int sample_idx = keep_indices[k];
                size_t byte_idx = sample_idx / 4;
                int shift = (sample_idx % 4) * 2;
                uint8_t code = (col_data[byte_idx] >> shift) & 3;
                y_local(k) += val_map[code];
            }
        }
        
        #pragma omp critical
        {
            y += y_local;
        }
    }
}

void PlinkReader::multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& y, int n_threads) {
    if (!stats_computed) computeStats(n_threads);
    if (v.size() != (long long)num_samples) throw std::runtime_error("Dimension mismatch in multiply_Zt_v");
    
    y.resize(num_snps);
    
    const uint8_t* data = static_cast<const uint8_t*>(bed_file.get_data()) + 3;
    
    #pragma omp parallel for num_threads(n_threads)
    for (long long j = 0; j < (long long)num_snps; ++j) {
        double inv_std = snp_inv_stds[j];
        if (inv_std == 0.0) {
            y(j) = 0.0;
            continue;
        }
        
        double mean = snp_means[j];
        
        double sum_centered_v = 0.0;
        
        const uint8_t* col_data = data + j * bytes_per_snp;
        
        // Unroll
        for (size_t byte_idx = 0; byte_idx < bytes_per_snp; ++byte_idx) {
            uint8_t byte = col_data[byte_idx];
            size_t base_sample = byte_idx * 4;
            
            // Helper lambda
            auto process = [&](int offset, int shift) {
                size_t idx = base_sample + offset;
                if (idx < num_samples) {
                    uint8_t code = (byte >> shift) & 3;
                    if (code != 1) { // Not Missing
                        sum_centered_v +=
                            (static_cast<double>(lookup_table[code]) - mean) * v(idx);
                    }
                }
            };
            
            process(0, 0);
            process(1, 2);
            process(2, 4);
            process(3, 6);
        }
        
        y(j) = sum_centered_v * inv_std;
    }
}

void PlinkReader::multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& y, const std::vector<int>& keep_indices, int n_threads) {
    if (!stats_computed) computeStats(n_threads);
    int n_keep = keep_indices.size();
    if (v.size() != n_keep) throw std::runtime_error("Dimension mismatch in multiply_Zt_v (subset)");
    
    y.resize(num_snps);
    const uint8_t* data = static_cast<const uint8_t*>(bed_file.get_data()) + 3;
    
    #pragma omp parallel for num_threads(n_threads)
    for (long long j = 0; j < (long long)num_snps; ++j) {
        double inv_std = snp_inv_stds[j];
        if (inv_std == 0.0) {
            y(j) = 0.0;
            continue;
        }
        
        double mean = snp_means[j];
        double sum_centered_v = 0.0;
        const uint8_t* col_data = data + j * bytes_per_snp;
        
        for (int k = 0; k < n_keep; ++k) {
            int sample_idx = keep_indices[k];
            size_t byte_idx = sample_idx / 4;
            int shift = (sample_idx % 4) * 2;
            uint8_t code = (col_data[byte_idx] >> shift) & 3;
            
            if (code != 1) { // 1 is missing
                sum_centered_v +=
                    (static_cast<double>(lookup_table[code]) - mean) * v(k);
            }
        }
        
        y(j) = sum_centered_v * inv_std;
    }
}

void PlinkReader::readSnp(int snp_idx, Eigen::VectorXd& out) const {
    if (snp_idx < 0 || snp_idx >= num_snps) throw std::runtime_error("Invalid SNP index");
    out.resize(num_samples);
    size_t byte_offset = 3 + (size_t)snp_idx * bytes_per_snp;
    const uint8_t* bed_data = static_cast<const uint8_t*>(bed_file.get_data());
    for (int i = 0; i < num_samples; ++i) {
        int byte_idx = i / 4;
        int bit_idx = (i % 4) * 2;
        unsigned char b = bed_data[byte_offset + byte_idx];
        unsigned char geno = (b >> bit_idx) & 3;
        
        if (geno == 1) { // Missing
            out(i) = std::numeric_limits<double>::quiet_NaN();
        } else if (geno == 0) { // Homozygous Major
            out(i) = 2.0;
        } else if (geno == 2) { // Heterozygous
            out(i) = 1.0;
        } else { // Homozygous Minor
            out(i) = 0.0;
        }
    }
}

} // namespace cosmic

namespace cosmic {

void PlinkReader::exportToText(const std::string& geno_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream geno_out(geno_file);
    if (!geno_out) throw std::runtime_error("Cannot open output geno file: " + geno_file);
    
    // Header: VariantID Sample1 Sample2 ...
    geno_out << "VariantID";
    for (const auto& fam : fam_data) {
        geno_out << " " << fam.fid << ":" << fam.iid;
    }
    geno_out << "\n";
    
    // Write each variant
    for (size_t j = 0; j < num_snps; ++j) {
        geno_out << bim_data[j].id;
        Eigen::VectorXf dosages;
        getSnpGenotypes(j, dosages);
        for (size_t i = 0; i < num_samples; ++i) {
            float val = dosages[i];
            if (std::isnan(val)) {
                geno_out << " NA";
            } else {
                geno_out << " " << val;
            }
        }
        geno_out << "\n";
    }
    
    // Variant info
    std::ofstream var_out(var_info_file);
    if (var_out) {
        var_out << "ID CHROM POS REF ALT\n";
        for (const auto& bim : bim_data) {
            var_out << bim.id << " " << bim.chrom << " " << bim.bp_pos << " " << bim.ref << " " << bim.alt << "\n";
        }
    }
    
    // Sample info
    std::ofstream sample_out(sample_info_file);
    if (sample_out) {
        sample_out << "FID IID\n";
        for (const auto& fam : fam_data) {
            sample_out << fam.fid << " " << fam.iid << "\n";
        }
    }
}

void PlinkReader::exportToNpy(const std::string& npy_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream file(npy_file, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open output NPY file: " + npy_file);
    
    // NPY header
    file.write("\x93NUMPY", 6);
    file.put(1); file.put(0);
    
    size_t header_len_pos = file.tellp();
    file.write("\x00\x00", 2); // placeholder
    
    // Write data first to compute exact size
    size_t data_start = file.tellp();
    
    Eigen::MatrixXf mat(num_snps, num_samples);
    for (size_t j = 0; j < num_snps; ++j) {
        Eigen::VectorXf dosages;
        getSnpGenotypes(j, dosages);
        mat.row(j) = dosages.transpose();
    }
    
    // Write data
    for (size_t j = 0; j < num_snps; ++j) {
        for (size_t i = 0; i < num_samples; ++i) {
            float val = mat(j, i);
            file.write(reinterpret_cast<const char*>(&val), 4);
        }
    }
    
    size_t data_end = file.tellp();
    
    // Now write header
    std::ostringstream header;
    header << "{'descr': '<f4', 'fortran_order': False, 'shape': (" << num_snps << ", " << num_samples << "), }";
    std::string header_str = header.str();
    while ((header_str.size() + 10) % 64 != 0) header_str += " ";
    
    size_t header_len = header_str.size();
    file.seekp(header_len_pos);
    uint16_t header_len_le = static_cast<uint16_t>(header_len);
    file.write(reinterpret_cast<const char*>(&header_len_le), 2);
    file.write(header_str.c_str(), header_len);
    file.seekp(data_end);
    
    // Write var info
    std::ofstream var_out(var_info_file);
    var_out << "ID CHROM POS REF ALT\n";
    for (const auto& bim : bim_data) {
        var_out << bim.id << " " << bim.chrom << " " << bim.bp_pos << " " << bim.ref << " " << bim.alt << "\n";
    }
    
    // Sample info
    std::ofstream sample_out(sample_info_file);
    sample_out << "FID IID\n";
    for (const auto& fam : fam_data) {
        sample_out << fam.fid << " " << fam.iid << "\n";
    }
}

void PlinkReader::exportToPedMap(const std::string& ped_file, const std::string& map_file) const {
    // Export MAP
    std::ofstream map_out(map_file);
    if (!map_out) throw std::runtime_error("Cannot open output MAP file: " + map_file);
    for (const auto& bim : bim_data) {
        map_out << bim.chrom << " " << bim.id << " " << bim.cm_pos << " "
                << bim.bp_pos << " " << bim.ref << " " << bim.alt << "\n";
    }
    
    // Export PED
    std::ofstream ped_out(ped_file);
    if (!ped_out) throw std::runtime_error("Cannot open output PED file: " + ped_file);
    
    for (size_t i = 0; i < num_samples; ++i) {
        const auto& fam = fam_data[i];
        ped_out << fam.fid << " " << fam.iid << " " << fam.father << " "
                << fam.mother << " " << fam.sex << " " << fam.phenotype;
        
        for (size_t j = 0; j < num_snps; ++j) {
            Eigen::VectorXf dosages;
            getSnpGenotypes(j, dosages);
            float val = dosages[i];
            const auto& bim = bim_data[j];
            if (std::isnan(val)) {
                ped_out << " 0 0";
            } else {
                if (val == 0.0) ped_out << " " << bim.ref << " " << bim.ref;
                else if (val == 1.0) ped_out << " " << bim.ref << " " << bim.alt;
                else ped_out << " " << bim.alt << " " << bim.alt;
            }
        }
        ped_out << "\n";
    }
}

void PlinkReader::exportToVcf(const std::string& vcf_file) const {
    std::ofstream vcf_out(vcf_file);
    if (!vcf_out) throw std::runtime_error("Cannot open output VCF file: " + vcf_file);
    
    vcf_out << "##fileformat=VCFv4.2\n";
    vcf_out << "##source=CosmicScanner\n";
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& fam : fam_data) {
        vcf_out << "\t" << fam.fid << ":" << fam.iid;
    }
    vcf_out << "\n";
    
    for (size_t j = 0; j < num_snps; ++j) {
        const auto& bim = bim_data[j];
        Eigen::VectorXf dosages;
        getSnpGenotypes(j, dosages);
        
        vcf_out << bim.chrom << "\t" << bim.bp_pos << "\t" << bim.id << "\t"
                << bim.ref << "\t" << bim.alt << "\t.\tPASS\t.\tGT";
        
        for (size_t i = 0; i < num_samples; ++i) {
            float val = dosages[i];
            if (std::isnan(val)) {
                vcf_out << "\t./.";
            } else {
                if (val == 0.0) vcf_out << "\t0/0";
                else if (val == 1.0) vcf_out << "\t0/1";
                else vcf_out << "\t1/1";
            }
        }
        vcf_out << "\n";
    }
}

} // namespace cosmic
