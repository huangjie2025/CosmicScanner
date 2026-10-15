#include "grm_builder.h"
#include "bgen_dominance_encoding.h"
#include "qc_utils.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace cosmic {

GrmBuilder::GrmBuilder() : options() {}

GrmBuilder::GrmBuilder(const Options& opts) : options(opts) {
}

void GrmBuilder::loadAuxFiles() {
    snp_weights.clear();
    if (!options.snp_weight_file.empty()) {
        std::ifstream fin(options.snp_weight_file);
        if (!fin) throw std::runtime_error("Cannot open snp_weight_file: " + options.snp_weight_file);
        std::string line, snp; float w;
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            if (iss >> snp >> w) {
                snp_weights[snp] = w;
            }
        }
    }

    pop_map.clear();
    num_pops = 0;
    if (!options.pop_class_file.empty()) {
        std::ifstream fin(options.pop_class_file);
        if (!fin) throw std::runtime_error("Cannot open pop_class_file: " + options.pop_class_file);
        std::string line, id, pop;
        std::unordered_map<std::string, int> pop_name_to_id;
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string col1, col2, col3;
            iss >> col1;
            if (iss >> col2) {
                if (iss >> col3) {
                    id = col2; pop = col3; // FID IID POP
                } else {
                    id = col1; pop = col2; // IID POP
                }
                if (pop_name_to_id.find(pop) == pop_name_to_id.end()) {
                    pop_name_to_id[pop] = num_pops++;
                }
                pop_map[id] = pop_name_to_id[pop];
            }
        }
    }
}

void GrmBuilder::compute(const PlinkReader& reader, const std::vector<int>& keep_indices) {
    if (!keep_indices.empty()) {
        throw std::runtime_error("PlinkReader GRM calculation with subsetting not yet implemented. Use PGEN/BGEN or full dataset.");
    }

    if (options.sparse) {
        computeSparse(reader);
        return;
    }

    size_t n = reader.getNumSamples();
    size_t m = reader.getNumSnps();

    if (n == 0 || m == 0) {
        throw std::runtime_error("Empty dataset provided to GrmBuilder");
    }

    // Initialize matrices
    G = Eigen::MatrixXd::Zero(n, n);
    N_snps = Eigen::MatrixXf::Zero(n, n);

    // Copy IDs for saving later (shared with computeSparse so both agree on sample order)
    initSampleIdsFromFam(reader);

    loadAuxFiles();
    ind_pop.assign(n, 0);
    if (num_pops > 0) {
        for (size_t i = 0; i < n; ++i) {
            auto it = pop_map.find(sample_iids[i]);
            if (it != pop_map.end()) ind_pop[i] = it->second;
            else {
                auto it2 = pop_map.find(sample_fids[i] + ":" + sample_iids[i]);
                if (it2 != pop_map.end()) ind_pop[i] = it2->second;
                else ind_pop[i] = 0;
            }
        }
    } else {
        num_pops = 1;
    }
    std::vector<double> pop_n_samples(num_pops, 0.0);
    for (size_t i = 0; i < n; ++i) pop_n_samples[ind_pop[i]] += 1.0;

    total_sum_2pq = 0.0;
    size_t total_valid_snps = 0;
    allele_freqs.assign(m, -1.0f);
    
    if (options.compute_homo_hete) {
        homo_counts.assign(n, 0);
        hete_counts.assign(n, 0);
        missing_counts.assign(n, 0);
    }

    // Get BIM info for chromosome filtering
    const auto& bim = reader.getBimInfo();
    
    // Set threads
#ifdef _OPENMP
    omp_set_num_threads(options.thread_num);
#endif

    std::cout << "Calculating GRM for " << n << " individuals and " << m << " SNPs (Optimized Bit-Parallel)." << std::endl;
    std::cout << "Block size: " << options.block_size << std::endl;
    std::string algo_name;
    switch(options.algorithm) {
        case 2: algo_name = "Zeng"; break;
        case 4: algo_name = "Vitezica"; break;
        case 5: algo_name = "TrueVar (Yang 2010 corrected)"; break;
        case 1: algo_name = "VanRaden"; break;
        default: algo_name = "Yang et al (GCTA default)"; break;
    }
    std::cout << "Algorithm: " << algo_name << std::endl;

    size_t num_blocks = (m + options.block_size - 1) / options.block_size;
    
    // Direct Access to Raw Data
    const uint8_t* raw_data = reader.getData();
    size_t bytes_per_snp = reader.getBytesPerSnp();
    
    for (size_t b = 0; b < num_blocks; ++b) {
        size_t start = b * options.block_size;
        size_t end = std::min(start + options.block_size, m);
        size_t current_block_size = end - start;

        // Check if this block has any relevant SNPs
        if (!options.chrom_filter.empty()) {
            bool block_has_chr = false;
            for(size_t k=start; k<end; ++k) {
                if(bim[k].chrom == options.chrom_filter) {
                    block_has_chr = true;
                    break;
                }
            }
            if (!block_has_chr) continue;
        }

        // 1. First Pass: Compute Stats (Mean, Variance, Filter)
        std::vector<float> snp_means(current_block_size * num_pops);
        std::vector<char> snp_valid(current_block_size, 0);
        std::vector<float> snp_scales(current_block_size * num_pops, 1.0f);
        std::vector<float> snp_2pq(current_block_size, 0.0f);
        std::vector<float> snp_weights_block(current_block_size, 1.0f);
        
        int valid_count_in_block = 0;
        
        #pragma omp parallel for reduction(+:valid_count_in_block)
        for (long long j = 0; j < (long long)current_block_size; ++j) {
            size_t global_idx = start + j;
            
            // Chrom filter
            if (!options.chrom_filter.empty() && bim[global_idx].chrom != options.chrom_filter) continue;
            if (!options.extract_set.empty() && !options.extract_set.count(bim[global_idx].id)) continue;
            if (!options.exclude_set.empty() && options.exclude_set.count(bim[global_idx].id)) continue;
            
            float weight = 1.0f;
            if (!snp_weights.empty()) {
                auto it = snp_weights.find(bim[global_idx].id);
                if (it != snp_weights.end()) weight = it->second;
                else continue;
            }
            if (weight <= 0.0f) continue;
            snp_weights_block[j] = weight;
            
            const uint8_t* snp_ptr = raw_data + 3 + global_idx * bytes_per_snp;
            
            std::vector<long long> pop_sum_x(num_pops, 0);
            std::vector<long long> pop_n_obs(num_pops, 0);
            std::vector<long long> pop_n_het(num_pops, 0);
            std::vector<long long> pop_n_hom1(num_pops, 0);
            long long global_sum_x = 0;
            long long global_n_obs = 0;
            long long n_hom1 = 0; // 00 -> 2
            long long n_het = 0;  // 10 -> 1
            long long n_hom2 = 0; // 11 -> 0
            
            for (size_t i = 0; i < n; i += 4) {
                size_t byte_idx = i >> 2;
                if (byte_idx >= bytes_per_snp) break; 
                uint8_t byte = snp_ptr[byte_idx];
                for (int k = 0; k < 4 && (i + k) < n; ++k) {
                    uint8_t code = (byte >> (k * 2)) & 3;
                    if (code == 1) continue;
                    int p = ind_pop[i + k];
                    pop_n_obs[p]++;
                    global_n_obs++;
                    if (code == 0) { pop_sum_x[p] += 2; pop_n_hom1[p]++; global_sum_x += 2; n_hom1++; } 
                    else if (code == 2) { pop_sum_x[p] += 1; pop_n_het[p]++; global_sum_x += 1; n_het++; } 
                    else { n_hom2++; } 
                }
            }
            
            if (global_n_obs == 0) continue;
            
            double global_freq = (double)global_sum_x / (2.0 * global_n_obs);
            double global_maf = std::min(global_freq, 1.0 - global_freq);
            
            if (global_maf < options.min_maf || std::abs(global_maf) < 1e-8) continue;
            
            if (options.min_hwe > 0.0f) {
                double p_hwe = calculateHWE(n_hom1, n_het, n_hom2);
                if (p_hwe < options.min_hwe) continue;
            }
            
            if (options.algorithm == 4) {
                snp_2pq[j] = (float)(weight * 4.0 * global_freq * global_freq * (1.0 - global_freq) * (1.0 - global_freq));
            } else if (options.algorithm == 5) {
                double p1 = (double)n_het / global_n_obs;
                double p2 = (double)n_hom1 / global_n_obs;
                double mu = p1 + 2.0 * p2;
                double var = p1 + 4.0 * p2 - mu * mu;
                snp_2pq[j] = (float)(weight * var);
            } else {
                snp_2pq[j] = (float)(weight * 2.0 * global_freq * (1.0 - global_freq));
            }
            allele_freqs[global_idx] = (float)global_freq;
            
            for (int p = 0; p < num_pops; ++p) {
                double freq = pop_n_obs[p] > 0 ? (double)pop_sum_x[p] / (2.0 * pop_n_obs[p]) : global_freq;
                snp_means[j * num_pops + p] = (float)(freq * 2.0);
                if (options.algorithm == 0 || options.algorithm == 3) {
                    double pq2 = 2.0 * freq * (1.0 - freq);
                    snp_scales[j * num_pops + p] = (float)((pq2 > 1e-9) ? (1.0 / std::sqrt(pq2)) : 0.0);
                } else if (options.algorithm == 5) {
                    double p1 = pop_n_obs[p] > 0 ? (double)pop_n_het[p] / pop_n_obs[p] : (double)n_het / global_n_obs;
                    double p2 = pop_n_obs[p] > 0 ? (double)pop_n_hom1[p] / pop_n_obs[p] : (double)n_hom1 / global_n_obs;
                    double mu = p1 + 2.0 * p2;
                    double var = p1 + 4.0 * p2 - mu * mu;
                    snp_scales[j * num_pops + p] = (float)((var > 1e-9) ? (1.0 / std::sqrt(var)) : 0.0);
                } else {
                    snp_scales[j * num_pops + p] = 1.0f;
                }
            }
            
            snp_valid[j] = true;
            valid_count_in_block++;
        }
        
        if (valid_count_in_block == 0) continue;

        // 2. Second Pass: Fill Z Matrix
        Eigen::MatrixXd Z(n, valid_count_in_block);
        
        std::vector<int> block_to_valid(current_block_size, -1);
        int cur_valid = 0;
        double block_sum_2pq = 0;
        
        for(size_t j=0; j<current_block_size; ++j) {
            if (snp_valid[j]) {
                block_to_valid[j] = cur_valid++;
                block_sum_2pq += snp_2pq[j];
            }
        }
        
        #pragma omp parallel for
        for (long long j = 0; j < (long long)current_block_size; ++j) {
            if (!snp_valid[j]) continue;
            
            int valid_idx = block_to_valid[j];
            float w = std::sqrt(snp_weights_block[j]);
            
            const uint8_t* snp_ptr = raw_data + 3 + (start + j) * bytes_per_snp;
            
            std::vector<std::array<double, 4>> pop_val_map(num_pops);
            for (int p = 0; p < num_pops; ++p) {
                double mean = snp_means[j * num_pops + p];
                double scale = snp_scales[j * num_pops + p] * w;
                if (options.algorithm == 4) {
                    double freq_p = mean / 2.0;
                    double freq_q = 1.0 - freq_p;
                    pop_val_map[p][0] = (-2.0 * freq_q * freq_q) * scale;
                    pop_val_map[p][1] = 0.0;
                    pop_val_map[p][2] = (2.0 * freq_p * freq_q) * scale;
                    pop_val_map[p][3] = (-2.0 * freq_p * freq_p) * scale;
                } else if (options.algorithm == 2) {
                    double freq_p = mean / 2.0;
                    pop_val_map[p][0] = (1.0 - freq_p) * scale;
                    pop_val_map[p][1] = 0.0;
                    pop_val_map[p][2] = (1.0 - 2.0 * freq_p) * scale;
                    pop_val_map[p][3] = (-freq_p) * scale;
                } else {
                    pop_val_map[p][0] = (2.0 - mean) * scale;
                    pop_val_map[p][1] = 0.0;
                    pop_val_map[p][2] = (1.0 - mean) * scale;
                    pop_val_map[p][3] = (0.0 - mean) * scale;
                }
            }
            
            double* z_col_ptr = Z.col(valid_idx).data();
            
            // Loop unrolling for SIMD auto-vectorization
            for (size_t i = 0; i < n; i += 4) {
                size_t byte_idx = i >> 2;
                uint8_t byte = snp_ptr[byte_idx];
                
                if (i + 3 < n) {
                    z_col_ptr[i+0] = pop_val_map[ind_pop[i+0]][(byte >> 0) & 3];
                    z_col_ptr[i+1] = pop_val_map[ind_pop[i+1]][(byte >> 2) & 3];
                    z_col_ptr[i+2] = pop_val_map[ind_pop[i+2]][(byte >> 4) & 3];
                    z_col_ptr[i+3] = pop_val_map[ind_pop[i+3]][(byte >> 6) & 3];
                    
                    if (options.compute_homo_hete) {
                        for(int k=0; k<4; ++k) {
                            uint8_t g = (byte >> (k * 2)) & 3;
                            if (g == 0 || g == 3) homo_counts[i+k]++;
                            else if (g == 2) hete_counts[i+k]++;
                            else if (g == 1) missing_counts[i+k]++;
                        }
                    }
                } else {
                    for (int k = 0; k < 4 && (i + k) < n; ++k) {
                        uint8_t g = (byte >> (k * 2)) & 3;
                        z_col_ptr[i+k] = pop_val_map[ind_pop[i+k]][g];
                        if (options.compute_homo_hete) {
                            if (g == 0 || g == 3) homo_counts[i+k]++;
                            else if (g == 2) hete_counts[i+k]++;
                            else if (g == 1) missing_counts[i+k]++;
                        }
                    }
                }
            }
        }
        
        // Fast Matrix Multiplication (uses AVX/AVX2/AVX-512 internally via Eigen/OpenBLAS/MKL)
        G.noalias() += Z * Z.transpose();
        
        total_sum_2pq += block_sum_2pq;
        total_valid_snps += valid_count_in_block;
        
        std::cout << "\rProcessed block " << b + 1 << "/" << num_blocks 
                  << " (Valid SNPs: " << total_valid_snps << ")" << std::flush;
    }
    std::cout << std::endl;
    
    if (total_valid_snps == 0 && !options.chrom_filter.empty()) {
        std::cout << "Warning: No valid SNPs found for chromosome " << options.chrom_filter << "\n";
    }

    if (options.chrom_filter.empty()) {
        for(size_t i=0; i<m; ++i) {
            observed_chroms.insert(bim[i].chrom);
        }
    } else {
        observed_chroms.insert(options.chrom_filter);
    }

    // Normalize
    if (!options.pop_class_file.empty() && num_pops > 1) {
        std::vector<double> pop_S(num_pops, 0.0);
        for (int p = 0; p < num_pops; ++p) {
            double sum_diag = 0.0;
            for (size_t i = 0; i < n; ++i) {
                if (ind_pop[i] == p) sum_diag += G(i, i);
            }
            pop_S[p] = (pop_n_samples[p] > 0) ? (sum_diag / pop_n_samples[p]) : 1.0;
            if (pop_S[p] < 1e-9) pop_S[p] = 1.0;
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                double divisor = std::sqrt(pop_S[ind_pop[i]] * pop_S[ind_pop[j]]);
                G(i, j) /= divisor;
                if (i != j) G(j, i) = G(i, j);
            }
        }
    } else {
        if (options.algorithm == 0 || options.algorithm == 3) {
             if (total_valid_snps > 0) G /= (double)total_valid_snps;
        } else if (options.algorithm == 2) {
             double tr = G.trace();
             if (tr > 0) G /= (tr / n);
        } else {
             if (total_sum_2pq > 0) G /= total_sum_2pq;
        }
    }

    // Fill N_snps
    N_snps.fill((float)total_valid_snps);
    total_valid_snps_count = (long long)total_valid_snps;
}

GrmBuilder::SplitGrmResult GrmBuilder::computeSplit(
    const PlinkReader& reader,
    const std::set<std::string>& chrom_labels)
{
    if (chrom_labels.empty()) {
        compute(reader);
        SplitGrmResult res;
        res.global = G;
        res.global_sum_2pq = total_sum_2pq;
        res.global_snp_count = total_valid_snps_count;
        return res;
    }

    size_t n = reader.getNumSamples();
    size_t m = reader.getNumSnps();
    if (n == 0 || m == 0)
        throw std::runtime_error("Empty dataset provided to GrmBuilder::computeSplit");

    // Initialize global GRM
    G = Eigen::MatrixXd::Zero(n, n);
    N_snps = Eigen::MatrixXf::Zero(n, n);

    // Per-chromosome accumulators
    std::vector<std::string> chrom_vec(chrom_labels.begin(), chrom_labels.end());
    std::map<std::string, int> chrom_to_idx;
    size_t n_chrom = chrom_vec.size();
    for (size_t k = 0; k < n_chrom; ++k)
        chrom_to_idx[chrom_vec[k]] = (int)k;

    std::vector<Eigen::MatrixXd> per_chrom_raw(n_chrom, Eigen::MatrixXd::Zero(n, n));
    std::vector<long long> per_chrom_snp_count(n_chrom, 0);
    std::vector<double> per_chrom_sum_2pq(n_chrom, 0.0);

    // Copy IDs
    const auto& fam = reader.getFamInfo();
    sample_fids.resize(n);
    sample_iids.resize(n);
    for (size_t i = 0; i < n; ++i) {
        sample_fids[i] = fam[i].fid;
        sample_iids[i] = fam[i].iid;
    }

    loadAuxFiles();
    ind_pop.assign(n, 0);
    if (num_pops > 0) {
        for (size_t i = 0; i < n; ++i) {
            auto it = pop_map.find(sample_iids[i]);
            if (it != pop_map.end()) ind_pop[i] = it->second;
            else {
                auto it2 = pop_map.find(sample_fids[i] + ":" + sample_iids[i]);
                if (it2 != pop_map.end()) ind_pop[i] = it2->second;
                else ind_pop[i] = 0;
            }
        }
    } else {
        num_pops = 1;
    }
    std::vector<double> pop_n_samples(num_pops, 0.0);
    for (size_t i = 0; i < n; ++i) pop_n_samples[ind_pop[i]] += 1.0;

    total_sum_2pq = 0.0;
    size_t total_valid_snps = 0;
    allele_freqs.assign(m, -1.0f);

    if (options.compute_homo_hete) {
        homo_counts.assign(n, 0);
        hete_counts.assign(n, 0);
        missing_counts.assign(n, 0);
    }

    const auto& bim = reader.getBimInfo();

#ifdef _OPENMP
    omp_set_num_threads(options.thread_num);
#endif

    std::cout << "Calculating global + " << n_chrom << " chromosome GRMs in single pass ("
              << n << " individuals, " << m << " SNPs)." << std::endl;

    size_t num_blocks = (m + options.block_size - 1) / options.block_size;
    const uint8_t* raw_data = reader.getData();
    size_t bytes_per_snp = reader.getBytesPerSnp();

    for (size_t b = 0; b < num_blocks; ++b) {
        size_t start = b * options.block_size;
        size_t end = std::min(start + options.block_size, m);
        size_t current_block_size = end - start;

        // Early skip if no SNP in this block belongs to any target chromosome
        if (!options.chrom_filter.empty()) {
            bool block_has_chr = false;
            for (size_t k = start; k < end; ++k) {
                if (bim[k].chrom == options.chrom_filter) { block_has_chr = true; break; }
            }
            if (!block_has_chr) continue;
        }

        // 1. First Pass: Compute Stats (Mean, Variance, Filter)
        std::vector<float> snp_means(current_block_size * num_pops);
        std::vector<char> snp_valid(current_block_size, 0);
        std::vector<float> snp_scales(current_block_size * num_pops, 1.0f);
        std::vector<float> snp_2pq(current_block_size, 0.0f);
        std::vector<float> snp_weights_block(current_block_size, 1.0f);
        std::vector<int> snp_chrom_idx(current_block_size, -1);

        int valid_count_in_block = 0;

        #pragma omp parallel for reduction(+:valid_count_in_block)
        for (long long j = 0; j < (long long)current_block_size; ++j) {
            size_t global_idx = start + j;

            if (!options.chrom_filter.empty() && bim[global_idx].chrom != options.chrom_filter) continue;
            if (!options.extract_set.empty() && !options.extract_set.count(bim[global_idx].id)) continue;
            if (!options.exclude_set.empty() && options.exclude_set.count(bim[global_idx].id)) continue;

            float weight = 1.0f;
            if (!snp_weights.empty()) {
                auto it = snp_weights.find(bim[global_idx].id);
                if (it != snp_weights.end()) weight = it->second;
                else continue;
            }
            if (weight <= 0.0f) continue;
            snp_weights_block[j] = weight;

            const uint8_t* snp_ptr = raw_data + 3 + global_idx * bytes_per_snp;

            std::vector<long long> pop_sum_x(num_pops, 0);
            std::vector<long long> pop_n_obs(num_pops, 0);
            std::vector<long long> pop_n_het(num_pops, 0);
            std::vector<long long> pop_n_hom1(num_pops, 0);
            long long global_sum_x = 0;
            long long global_n_obs = 0;
            long long n_hom1 = 0;
            long long n_het = 0;
            long long n_hom2 = 0;

            for (size_t i = 0; i < n; i += 4) {
                size_t byte_idx = i >> 2;
                if (byte_idx >= bytes_per_snp) break;
                uint8_t byte = snp_ptr[byte_idx];
                for (int k = 0; k < 4 && (i + k) < n; ++k) {
                    uint8_t code = (byte >> (k * 2)) & 3;
                    if (code == 1) continue;
                    int p = ind_pop[i + k];
                    pop_n_obs[p]++;
                    global_n_obs++;
                    if (code == 0) { pop_sum_x[p] += 2; pop_n_hom1[p]++; global_sum_x += 2; n_hom1++; }
                    else if (code == 2) { pop_sum_x[p] += 1; pop_n_het[p]++; global_sum_x += 1; n_het++; }
                    else { n_hom2++; }
                }
            }

            if (global_n_obs == 0) continue;

            double global_freq = (double)global_sum_x / (2.0 * global_n_obs);
            double global_maf = std::min(global_freq, 1.0 - global_freq);

            if (global_maf < options.min_maf || std::abs(global_maf) < 1e-8) continue;

            if (options.min_hwe > 0.0f) {
                double p_hwe = calculateHWE(n_hom1, n_het, n_hom2);
                if (p_hwe < options.min_hwe) continue;
            }

            if (options.algorithm == 4) {
                snp_2pq[j] = (float)(weight * 4.0 * global_freq * global_freq * (1.0 - global_freq) * (1.0 - global_freq));
            } else if (options.algorithm == 5) {
                double p1 = (double)n_het / global_n_obs;
                double p2 = (double)n_hom1 / global_n_obs;
                double mu = p1 + 2.0 * p2;
                double var = p1 + 4.0 * p2 - mu * mu;
                snp_2pq[j] = (float)(weight * var);
            } else {
                snp_2pq[j] = (float)(weight * 2.0 * global_freq * (1.0 - global_freq));
            }
            allele_freqs[global_idx] = (float)global_freq;

            for (int p = 0; p < num_pops; ++p) {
                double freq = pop_n_obs[p] > 0 ? (double)pop_sum_x[p] / (2.0 * pop_n_obs[p]) : global_freq;
                snp_means[j * num_pops + p] = (float)(freq * 2.0);
                if (options.algorithm == 0 || options.algorithm == 3) {
                    double pq2 = 2.0 * freq * (1.0 - freq);
                    snp_scales[j * num_pops + p] = (float)((pq2 > 1e-9) ? (1.0 / std::sqrt(pq2)) : 0.0);
                } else if (options.algorithm == 5) {
                    double p1 = pop_n_obs[p] > 0 ? (double)pop_n_het[p] / pop_n_obs[p] : (double)n_het / global_n_obs;
                    double p2 = pop_n_obs[p] > 0 ? (double)pop_n_hom1[p] / pop_n_obs[p] : (double)n_hom1 / global_n_obs;
                    double mu = p1 + 2.0 * p2;
                    double var = p1 + 4.0 * p2 - mu * mu;
                    snp_scales[j * num_pops + p] = (float)((var > 1e-9) ? (1.0 / std::sqrt(var)) : 0.0);
                } else {
                    snp_scales[j * num_pops + p] = 1.0f;
                }
            }

            snp_valid[j] = true;
            // Record which chromosome this valid SNP belongs to
            auto it = chrom_to_idx.find(bim[global_idx].chrom);
            if (it != chrom_to_idx.end())
                snp_chrom_idx[j] = it->second;
            valid_count_in_block++;
        }

        if (valid_count_in_block == 0) continue;

        // 2. Second Pass: Fill Z Matrix
        Eigen::MatrixXd Z(n, valid_count_in_block);

        std::vector<int> block_to_valid(current_block_size, -1);
        int cur_valid = 0;
        double block_sum_2pq = 0;

        for (size_t j = 0; j < current_block_size; ++j) {
            if (snp_valid[j]) {
                block_to_valid[j] = cur_valid++;
                block_sum_2pq += snp_2pq[j];
            }
        }

        #pragma omp parallel for
        for (long long j = 0; j < (long long)current_block_size; ++j) {
            if (!snp_valid[j]) continue;

            int valid_idx = block_to_valid[j];
            float w = std::sqrt(snp_weights_block[j]);

            const uint8_t* snp_ptr = raw_data + 3 + (start + j) * bytes_per_snp;

            std::vector<std::array<double, 4>> pop_val_map(num_pops);
            for (int p = 0; p < num_pops; ++p) {
                double mean = snp_means[j * num_pops + p];
                double scale = snp_scales[j * num_pops + p] * w;
                if (options.algorithm == 4) {
                    double freq_p = mean / 2.0;
                    double freq_q = 1.0 - freq_p;
                    pop_val_map[p][0] = (-2.0 * freq_q * freq_q) * scale;
                    pop_val_map[p][1] = 0.0;
                    pop_val_map[p][2] = (2.0 * freq_p * freq_q) * scale;
                    pop_val_map[p][3] = (-2.0 * freq_p * freq_p) * scale;
                } else if (options.algorithm == 2) {
                    double freq_p = mean / 2.0;
                    pop_val_map[p][0] = (1.0 - freq_p) * scale;
                    pop_val_map[p][1] = 0.0;
                    pop_val_map[p][2] = (1.0 - 2.0 * freq_p) * scale;
                    pop_val_map[p][3] = (-freq_p) * scale;
                } else {
                    pop_val_map[p][0] = (2.0 - mean) * scale;
                    pop_val_map[p][1] = 0.0;
                    pop_val_map[p][2] = (1.0 - mean) * scale;
                    pop_val_map[p][3] = (0.0 - mean) * scale;
                }
            }

            double* z_col_ptr = Z.col(valid_idx).data();

            for (size_t i = 0; i < n; i += 4) {
                size_t byte_idx = i >> 2;
                uint8_t byte = snp_ptr[byte_idx];

                if (i + 3 < n) {
                    z_col_ptr[i+0] = pop_val_map[ind_pop[i+0]][(byte >> 0) & 3];
                    z_col_ptr[i+1] = pop_val_map[ind_pop[i+1]][(byte >> 2) & 3];
                    z_col_ptr[i+2] = pop_val_map[ind_pop[i+2]][(byte >> 4) & 3];
                    z_col_ptr[i+3] = pop_val_map[ind_pop[i+3]][(byte >> 6) & 3];

                    if (options.compute_homo_hete) {
                        for (int k = 0; k < 4; ++k) {
                            uint8_t g = (byte >> (k * 2)) & 3;
                            if (g == 0 || g == 3) homo_counts[i+k]++;
                            else if (g == 2) hete_counts[i+k]++;
                            else if (g == 1) missing_counts[i+k]++;
                        }
                    }
                } else {
                    for (int k = 0; k < 4 && (i + k) < n; ++k) {
                        uint8_t g = (byte >> (k * 2)) & 3;
                        z_col_ptr[i+k] = pop_val_map[ind_pop[i+k]][g];
                        if (options.compute_homo_hete) {
                            if (g == 0 || g == 3) homo_counts[i+k]++;
                            else if (g == 2) hete_counts[i+k]++;
                            else if (g == 1) missing_counts[i+k]++;
                        }
                    }
                }
            }
        }

        // Accumulate global GRM
        G.noalias() += Z * Z.transpose();

        // For each chromosome, build Z_sub and accumulate into per-chromosome raw G
        // Build per-chromosome column index lists
        std::map<int, std::vector<int>> chrom_col_lists;
        int col = 0;
        for (size_t j = 0; j < current_block_size; ++j) {
            if (snp_valid[j]) {
                int cidx = snp_chrom_idx[j];
                if (cidx >= 0)
                    chrom_col_lists[cidx].push_back(col);
                col++;
            }
        }

        for (auto& kv : chrom_col_lists) {
            int cidx = kv.first;
            auto& cols = kv.second;
            if (cols.empty()) continue;

            Eigen::MatrixXd Z_sub(n, cols.size());
            for (size_t c = 0; c < cols.size(); ++c)
                Z_sub.col(c) = Z.col(cols[c]);

            per_chrom_raw[cidx].noalias() += Z_sub * Z_sub.transpose();
        }

        // Update per-chromosome stats
        for (size_t j = 0; j < current_block_size; ++j) {
            if (snp_valid[j]) {
                int cidx = snp_chrom_idx[j];
                if (cidx >= 0) {
                    per_chrom_snp_count[cidx]++;
                    per_chrom_sum_2pq[cidx] += snp_2pq[j];
                }
            }
        }

        total_sum_2pq += block_sum_2pq;
        total_valid_snps += valid_count_in_block;

        std::cout << "\rProcessed block " << b + 1 << "/" << num_blocks
                  << " (Valid SNPs: " << total_valid_snps << ")" << std::flush;
    }
    std::cout << std::endl;

    for (size_t i = 0; i < m; ++i) {
        observed_chroms.insert(bim[i].chrom);
    }

    // Normalize global GRM
    if (!options.pop_class_file.empty() && num_pops > 1) {
        std::vector<double> pop_S(num_pops, 0.0);
        for (int p = 0; p < num_pops; ++p) {
            double sum_diag = 0.0;
            for (size_t i = 0; i < n; ++i)
                if (ind_pop[i] == p) sum_diag += G(i, i);
            pop_S[p] = (pop_n_samples[p] > 0) ? (sum_diag / pop_n_samples[p]) : 1.0;
            if (pop_S[p] < 1e-9) pop_S[p] = 1.0;
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                double divisor = std::sqrt(pop_S[ind_pop[i]] * pop_S[ind_pop[j]]);
                G(i, j) /= divisor;
                if (i != j) G(j, i) = G(i, j);
            }
        }
    } else {
        if (options.algorithm == 0 || options.algorithm == 3) {
            if (total_valid_snps > 0) G /= (double)total_valid_snps;
        } else if (options.algorithm == 2) {
            double tr = G.trace();
            if (tr > 0) G /= (tr / n);
        } else {
            if (total_sum_2pq > 0) G /= total_sum_2pq;
        }
    }

    N_snps.fill((float)total_valid_snps);
    total_valid_snps_count = (long long)total_valid_snps;

    // Assemble result
    SplitGrmResult res;
    res.global = G;
    res.global_sum_2pq = total_sum_2pq;
    res.global_snp_count = total_valid_snps_count;
    res.chrom_labels = chrom_vec;
    res.per_chrom_raw.resize(n_chrom);
    res.per_chrom_snp_count.resize(n_chrom);
    res.per_chrom_sum_2pq.resize(n_chrom);
    for (size_t k = 0; k < n_chrom; ++k) {
        res.per_chrom_raw[k] = std::move(per_chrom_raw[k]);
        res.per_chrom_snp_count[k] = per_chrom_snp_count[k];
        res.per_chrom_sum_2pq[k] = per_chrom_sum_2pq[k];
    }

    return res;
}

void GrmBuilder::compute(BgenReader& reader, const std::vector<int>& keep_indices) {
    size_t n_total = reader.getNumSamples();
    size_t n = keep_indices.empty() ? n_total : keep_indices.size();
    size_t m = reader.getNumVariants();
    
    if (n == 0) throw std::runtime_error("BGEN reader empty");
    if (options.algorithm == 4 && reader.getCompressionType() == 0) {
        throw std::runtime_error("Probabilistic BGEN dominance GRM does not support uncompressed Layout-2 input");
    }
    if (options.algorithm == 4 && options.min_hwe > 0.0f) {
        throw std::runtime_error("Probabilistic BGEN dominance GRM does not support exact HWE filtering");
    }
    if (options.algorithm == 4 && !options.pop_class_file.empty()) {
        throw std::runtime_error("Probabilistic BGEN dominance GRM does not support population-specific frequencies");
    }

    G = Eigen::MatrixXd::Zero(n, n);
    N_snps = Eigen::MatrixXf::Zero(n, n);
    
    const auto& ids = reader.getSampleIds();
    sample_fids.resize(n);
    sample_iids.resize(n);
    for(size_t i=0; i<n; ++i) {
        size_t idx = keep_indices.empty() ? i : keep_indices[i];
        sample_fids[i] = ids[idx];
        sample_iids[i] = ids[idx];
    }
    
    loadAuxFiles();
    ind_pop.assign(n, 0);
    if (num_pops > 0) {
        for (size_t i = 0; i < n; ++i) {
            auto it = pop_map.find(sample_iids[i]);
            if (it != pop_map.end()) ind_pop[i] = it->second;
            else {
                auto it2 = pop_map.find(sample_fids[i] + ":" + sample_iids[i]);
                if (it2 != pop_map.end()) ind_pop[i] = it2->second;
                else ind_pop[i] = 0;
            }
        }
    } else {
        num_pops = 1;
    }
    std::vector<double> pop_n_samples(num_pops, 0.0);
    for (size_t i = 0; i < n; ++i) pop_n_samples[ind_pop[i]] += 1.0;

    total_sum_2pq = 0.0;
    size_t total_valid_snps = 0;
    allele_freqs.assign(m, -1.0f);
    
    if (options.compute_homo_hete) {
        homo_counts.assign(n, 0);
        hete_counts.assign(n, 0);
        missing_counts.assign(n, 0);
    }
    
    reader.reset();
    
    std::cout << "Calculating GRM for " << n << " individuals and " << m << " variants (BGEN)." << std::endl;
    std::cout << "Algorithm: " << ((options.algorithm == 0 || options.algorithm == 3) ? "Yang et al (GCTA default)" : (options.algorithm == 2 ? "Zeng" : (options.algorithm == 4 ? "Vitezica" : "VanRaden"))) << std::endl;
    
    size_t block_size = options.block_size;
    std::vector<Eigen::VectorXf> Z_cols;
    Z_cols.reserve(block_size);
    
    size_t processed_snps = 0;
    double sum_2pq_accum = 0.0;
    
    while (processed_snps < m) {
        Z_cols.clear();
        double block_2pq = 0.0;
        
        for (size_t k = 0; k < block_size && processed_snps < m; ++k) {
            std::string chrom, rsid, a1, a2;
            uint32_t pos;
            Eigen::VectorXf dosage;
            BgenGenotypeProbs probs;

            const bool read_ok = options.algorithm == 4
                ? reader.readVariantProbs(chrom, rsid, pos, a1, a2, probs)
                : reader.readVariant(chrom, rsid, pos, a1, a2, dosage);
            if (!read_ok) {
                processed_snps = m; 
                break;
            }
            if (options.algorithm == 4) dosage = probs.dosage;
            processed_snps++;
            
            if (options.chrom_filter.empty()) {
                observed_chroms.insert(chrom);
            }

            if (!options.extract_set.empty() && options.extract_set.find(rsid) == options.extract_set.end()) continue;
            if (!options.exclude_set.empty() && options.exclude_set.find(rsid) != options.exclude_set.end()) continue;

            float weight = 1.0f;
            if (!snp_weights.empty()) {
                auto it = snp_weights.find(rsid);
                if (it != snp_weights.end()) weight = it->second;
                else continue;
            }
            if (weight <= 0.0f) continue;
            
            if (!options.chrom_filter.empty() && chrom != options.chrom_filter) continue;
            
            std::vector<double> pop_sum(num_pops, 0.0);
            std::vector<int> pop_count(num_pops, 0);
            std::vector<int> pop_het(num_pops, 0);
            std::vector<int> pop_hom1(num_pops, 0);
            double global_sum = 0.0;
            int global_count = 0;
            int n_het = 0;
            int n_hom1 = 0;
            bgen::DominanceQcMetrics dominance_qc;
            if (options.algorithm == 4) {
                dominance_qc = bgen::dominanceQcMetrics(
                    probs.p0, probs.p1, probs.p2, probs.missing, keep_indices);
                if (dominance_qc.observed == 0 ||
                    dominance_qc.missing_rate > options.bgen_dom_max_missing) continue;
            }
            
            if (keep_indices.empty()) {
                for(int i=0; i<n; ++i) {
                    if (!std::isnan(dosage[i])) {
                        double val = dosage[i];
                        global_sum += val;
                        global_count++;
                        int p = ind_pop[i];
                        pop_sum[p] += val;
                        pop_count[p]++;
                        if (val > 1.5) { pop_hom1[p]++; n_hom1++; }
                        else if (val >= 0.5) { pop_het[p]++; n_het++; }
                    }
                }
            } else {
                for(int i=0; i<n; ++i) {
                    float val = dosage[keep_indices[i]];
                    if (!std::isnan(val)) {
                        global_sum += val;
                        global_count++;
                        int p = ind_pop[i];
                        pop_sum[p] += val;
                        pop_count[p]++;
                        if (val > 1.5f) { pop_hom1[p]++; n_hom1++; }
                        else if (val >= 0.5f) { pop_het[p]++; n_het++; }
                    }
                }
            }
            
            if (global_count == 0) continue;
            
            double global_freq = options.algorithm == 4
                ? dominance_qc.frequency
                : (global_sum / global_count) / 2.0;
            double global_maf = std::min(global_freq, 1.0 - global_freq);
            
            if (global_maf < options.min_maf || std::abs(global_maf) < 1e-8) continue;
            if (options.algorithm == 4 && dominance_qc.info_d < options.bgen_dom_min_info) continue;
            
            std::vector<double> pop_mean(num_pops);
            std::vector<double> pop_scale(num_pops, 1.0);
            for (int p = 0; p < num_pops; ++p) {
                double p_freq = pop_count[p] > 0 ? (pop_sum[p] / pop_count[p]) / 2.0 : global_freq;
                pop_mean[p] = 2.0 * p_freq;
                if (options.algorithm == 0 || options.algorithm == 3) {
                    double pq2 = 2.0 * p_freq * (1.0 - p_freq);
                    pop_scale[p] = (pq2 > 1e-9) ? (1.0 / std::sqrt(pq2)) : 0.0;
                } else if (options.algorithm == 5) {
                    double p1 = pop_count[p] > 0 ? (double)pop_het[p] / pop_count[p] : (double)n_het / global_count;
                    double p2 = pop_count[p] > 0 ? (double)pop_hom1[p] / pop_count[p] : (double)n_hom1 / global_count;
                    double mu = p1 + 2.0 * p2;
                    double var = p1 + 4.0 * p2 - mu * mu;
                    pop_scale[p] = (var > 1e-9) ? (1.0 / std::sqrt(var)) : 0.0;
                }
            }
            
            Eigen::VectorXf z_col(n);
            float w = std::sqrt(weight);
            
            for(int i=0; i<n; ++i) {
                float val = keep_indices.empty() ? dosage[i] : dosage[keep_indices[i]];
                if (std::isnan(val)) {
                    z_col[i] = 0.0f;
                    if (options.compute_homo_hete) missing_counts[i]++;
                } else {
                    int p = ind_pop[i];
                    double mean = pop_mean[p];
                    double scale = pop_scale[p] * w;

                    if (options.algorithm == 4) {
                        const int source = keep_indices.empty() ? i : keep_indices[i];
                        z_col[i] = probs.missing[source] != 0.0f
                            ? 0.0f
                            : static_cast<float>(bgen::posteriorMeanVitezica(
                                  probs.p0[source], probs.p1[source], probs.p2[source],
                                  global_freq) * w);
                    } else if (options.algorithm == 2) {
                        double freq_p = mean / 2.0;
                        if (val > 1.5f) z_col[i] = (float)((1.0 - freq_p) * scale);
                        else if (val < 0.5f) z_col[i] = (float)(-freq_p * scale);
                        else z_col[i] = (float)((1.0 - 2.0 * freq_p) * scale);
                    } else {
                        z_col[i] = (float)((val - mean) * scale);
                    }
                    
                    if (options.compute_homo_hete) {
                        if (val < 0.5f || val > 1.5f) homo_counts[i]++;
                        else hete_counts[i]++;
                    }
                }
            }
            
            Z_cols.push_back(z_col);
            if (options.algorithm == 4) {
                block_2pq += weight * 4.0 * global_freq * global_freq * (1.0 - global_freq) * (1.0 - global_freq);
            } else if (options.algorithm == 5) {
                double p1 = (double)n_het / global_count;
                double p2 = (double)n_hom1 / global_count;
                double mu = p1 + 2.0 * p2;
                double var = p1 + 4.0 * p2 - mu * mu;
                block_2pq += weight * var;
            } else {
                block_2pq += weight * 2.0 * global_freq * (1.0 - global_freq);
            }
            allele_freqs[processed_snps - 1] = (float)global_freq;
        }
        
        if (!Z_cols.empty()) {
            size_t valid_count = Z_cols.size();
            Eigen::MatrixXd Z(n, valid_count);
            
            #pragma omp parallel for
            for(int k=0; k<(int)valid_count; ++k) {
                Z.col(k) = Z_cols[k].cast<double>();
            }
            
            G.noalias() += Z * Z.transpose();
            
            total_valid_snps += valid_count;
            sum_2pq_accum += block_2pq;
        }
        
        std::cout << "\rProcessed " << processed_snps << " variants (Valid: " << total_valid_snps << ")" << std::flush;
    }
    std::cout << std::endl;
    
    total_sum_2pq = sum_2pq_accum;
    
    if (!options.pop_class_file.empty() && num_pops > 1) {
        std::vector<double> pop_S(num_pops, 0.0);
        for (int p = 0; p < num_pops; ++p) {
            double sum_diag = 0.0;
            for (size_t i = 0; i < n; ++i) {
                if (ind_pop[i] == p) sum_diag += G(i, i);
            }
            pop_S[p] = (pop_n_samples[p] > 0) ? (sum_diag / pop_n_samples[p]) : 1.0;
            if (pop_S[p] < 1e-9) pop_S[p] = 1.0;
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                double divisor = std::sqrt(pop_S[ind_pop[i]] * pop_S[ind_pop[j]]);
                G(i, j) /= divisor;
                if (i != j) G(j, i) = G(i, j);
            }
        }
    } else {
        if (options.algorithm == 0 || options.algorithm == 3) {
             if (total_valid_snps > 0) G /= (double)total_valid_snps;
        } else if (options.algorithm == 2) {
             double tr = G.trace();
             if (tr > 0) G /= (tr / n);
        } else {
             if (total_sum_2pq > 0) G /= total_sum_2pq;
        }
    }
    
    N_snps.fill((float)total_valid_snps);
    total_valid_snps_count = (long long)total_valid_snps;
}

GrmBuilder::SplitGrmResult GrmBuilder::computeSplit(
    BgenReader& reader,
    const std::set<std::string>& chrom_labels,
    const std::vector<int>& keep_indices)
{
    if (chrom_labels.empty()) {
        compute(reader, keep_indices);
        SplitGrmResult res;
        res.global = G; res.global_sum_2pq = total_sum_2pq; res.global_snp_count = total_valid_snps_count;
        return res;
    }

    size_t n_total = reader.getNumSamples();
    size_t n = keep_indices.empty() ? n_total : keep_indices.size();
    size_t m = reader.getNumVariants();
    if (n == 0) throw std::runtime_error("BGEN reader empty");
    if (options.algorithm == 4 && reader.getCompressionType() == 0) {
        throw std::runtime_error("Probabilistic BGEN dominance GRM does not support uncompressed Layout-2 input");
    }
    if (options.algorithm == 4 && options.min_hwe > 0.0f) {
        throw std::runtime_error("Probabilistic BGEN dominance GRM does not support exact HWE filtering");
    }
    if (options.algorithm == 4 && !options.pop_class_file.empty()) {
        throw std::runtime_error("Probabilistic BGEN dominance GRM does not support population-specific frequencies");
    }

    G = Eigen::MatrixXd::Zero(n, n);
    N_snps = Eigen::MatrixXf::Zero(n, n);

    std::vector<std::string> chrom_vec(chrom_labels.begin(), chrom_labels.end());
    std::map<std::string, int> chrom_to_idx;
    size_t n_chrom = chrom_vec.size();
    for (size_t k = 0; k < n_chrom; ++k) chrom_to_idx[chrom_vec[k]] = (int)k;

    std::vector<Eigen::MatrixXd> per_chrom_raw(n_chrom, Eigen::MatrixXd::Zero(n, n));
    std::vector<long long> per_chrom_snp_count(n_chrom, 0);
    std::vector<double> per_chrom_sum_2pq(n_chrom, 0.0);

    const auto& ids = reader.getSampleIds();
    sample_fids.resize(n); sample_iids.resize(n);
    for (size_t i = 0; i < n; ++i) {
        size_t idx = keep_indices.empty() ? i : keep_indices[i];
        sample_fids[i] = ids[idx]; sample_iids[i] = ids[idx];
    }

    loadAuxFiles();
    ind_pop.assign(n, 0);
    if (num_pops > 0) {
        for (size_t i = 0; i < n; ++i) {
            auto it = pop_map.find(sample_iids[i]);
            if (it != pop_map.end()) ind_pop[i] = it->second;
            else {
                auto it2 = pop_map.find(sample_fids[i] + ":" + sample_iids[i]);
                if (it2 != pop_map.end()) ind_pop[i] = it2->second;
                else ind_pop[i] = 0;
            }
        }
    } else { num_pops = 1; }
    std::vector<double> pop_n_samples(num_pops, 0.0);
    for (size_t i = 0; i < n; ++i) pop_n_samples[ind_pop[i]] += 1.0;

    total_sum_2pq = 0.0;
    size_t total_valid_snps = 0;
    allele_freqs.assign(m, -1.0f);
    if (options.compute_homo_hete) {
        homo_counts.assign(n, 0); hete_counts.assign(n, 0); missing_counts.assign(n, 0);
    }

    reader.reset();
    std::cout << "Calculating global + " << n_chrom << " chromosome GRMs in single pass ("
              << n << " individuals, " << m << " BGEN variants)." << std::endl;

    size_t block_size = options.block_size;
    std::vector<Eigen::VectorXf> Z_cols;
    Z_cols.reserve(block_size);
    std::vector<int> col_chrom;

    size_t processed_snps = 0;
    double sum_2pq_accum = 0.0;
    long long bgen_qc_candidates = 0;
    long long bgen_qc_rejected_missing = 0;
    long long bgen_qc_missing_warnings = 0;
    long long bgen_qc_rejected_maf = 0;
    long long bgen_qc_rejected_info = 0;
    long long bgen_qc_info_warnings = 0;
    long long bgen_qc_entropy_warnings = 0;

    while (processed_snps < m) {
        Z_cols.clear();
        col_chrom.clear();
        double block_2pq = 0.0;

        for (size_t k = 0; k < block_size && processed_snps < m; ++k) {
            std::string chrom, rsid, a1, a2;
            uint32_t pos;
            Eigen::VectorXf dosage;
            BgenGenotypeProbs probs;

            const bool read_ok = options.algorithm == 4
                ? reader.readVariantProbs(chrom, rsid, pos, a1, a2, probs)
                : reader.readVariant(chrom, rsid, pos, a1, a2, dosage);
            if (!read_ok) {
                processed_snps = m; break;
            }
            if (options.algorithm == 4) dosage = probs.dosage;
            processed_snps++;

            if (options.chrom_filter.empty())
                observed_chroms.insert(chrom);

            if (!options.extract_set.empty() && options.extract_set.find(rsid) == options.extract_set.end()) continue;
            if (!options.exclude_set.empty() && options.exclude_set.find(rsid) != options.exclude_set.end()) continue;

            float weight = 1.0f;
            if (!snp_weights.empty()) {
                auto it = snp_weights.find(rsid);
                if (it != snp_weights.end()) weight = it->second;
                else continue;
            }
            if (weight <= 0.0f) continue;
            if (!options.chrom_filter.empty() && chrom != options.chrom_filter) continue;

            int cidx = -1;
            auto cit = chrom_to_idx.find(chrom);
            if (cit != chrom_to_idx.end()) cidx = cit->second;

            std::vector<double> pop_sum(num_pops, 0.0);
            std::vector<int> pop_count(num_pops, 0);
            std::vector<int> pop_het(num_pops, 0);
            std::vector<int> pop_hom1(num_pops, 0);
            double global_sum = 0.0;
            int global_count = 0;
            int n_het = 0, n_hom1 = 0;
            bgen::DominanceQcMetrics dominance_qc;
            if (options.algorithm == 4) {
                ++bgen_qc_candidates;
                dominance_qc = bgen::dominanceQcMetrics(
                    probs.p0, probs.p1, probs.p2, probs.missing, keep_indices);
                if (dominance_qc.observed == 0 ||
                    dominance_qc.missing_rate > options.bgen_dom_max_missing) {
                    ++bgen_qc_rejected_missing;
                    continue;
                }
                if (dominance_qc.missing_rate > bgen::kDominanceMissingWarning) {
                    ++bgen_qc_missing_warnings;
                }
            }

            if (keep_indices.empty()) {
                for (int i = 0; i < n; ++i) {
                    if (!std::isnan(dosage[i])) {
                        double val = dosage[i]; global_sum += val; global_count++;
                        int p = ind_pop[i]; pop_sum[p] += val; pop_count[p]++;
                        if (val > 1.5) { pop_hom1[p]++; n_hom1++; }
                        else if (val >= 0.5) { pop_het[p]++; n_het++; }
                    }
                }
            } else {
                for (int i = 0; i < n; ++i) {
                    float val = dosage[keep_indices[i]];
                    if (!std::isnan(val)) {
                        global_sum += val; global_count++;
                        int p = ind_pop[i]; pop_sum[p] += val; pop_count[p]++;
                        if (val > 1.5f) { pop_hom1[p]++; n_hom1++; }
                        else if (val >= 0.5f) { pop_het[p]++; n_het++; }
                    }
                }
            }

            if (global_count == 0) continue;

            double global_freq = options.algorithm == 4
                ? dominance_qc.frequency
                : (global_sum / global_count) / 2.0;
            double global_maf = std::min(global_freq, 1.0 - global_freq);
            if (global_maf < options.min_maf || std::abs(global_maf) < 1e-8) {
                if (options.algorithm == 4) ++bgen_qc_rejected_maf;
                continue;
            }
            if (options.algorithm == 4 && dominance_qc.info_d < options.bgen_dom_min_info) {
                ++bgen_qc_rejected_info;
                continue;
            }
            if (options.algorithm == 4 &&
                dominance_qc.info_d < bgen::kDominanceInfoWarning) {
                ++bgen_qc_info_warnings;
            }
            if (options.algorithm == 4 &&
                dominance_qc.mean_entropy > options.bgen_dom_entropy_warning) {
                ++bgen_qc_entropy_warnings;
            }

            std::vector<double> pop_mean(num_pops);
            std::vector<double> pop_scale(num_pops, 1.0);
            for (int p = 0; p < num_pops; ++p) {
                double p_freq = pop_count[p] > 0 ? (pop_sum[p] / pop_count[p]) / 2.0 : global_freq;
                pop_mean[p] = 2.0 * p_freq;
                if (options.algorithm == 0 || options.algorithm == 3) {
                    double pq2 = 2.0 * p_freq * (1.0 - p_freq);
                    pop_scale[p] = (pq2 > 1e-9) ? (1.0 / std::sqrt(pq2)) : 0.0;
                } else if (options.algorithm == 5) {
                    double p1 = pop_count[p] > 0 ? (double)pop_het[p] / pop_count[p] : (double)n_het / global_count;
                    double p2 = pop_count[p] > 0 ? (double)pop_hom1[p] / pop_count[p] : (double)n_hom1 / global_count;
                    double mu = p1 + 2.0 * p2;
                    double var = p1 + 4.0 * p2 - mu * mu;
                    pop_scale[p] = (var > 1e-9) ? (1.0 / std::sqrt(var)) : 0.0;
                }
            }

            Eigen::VectorXf z_col(n);
            float w = std::sqrt(weight);
            for (int i = 0; i < n; ++i) {
                float val = keep_indices.empty() ? dosage[i] : dosage[keep_indices[i]];
                if (std::isnan(val)) {
                    z_col[i] = 0.0f;
                    if (options.compute_homo_hete) missing_counts[i]++;
                } else {
                    int p = ind_pop[i];
                    double mean = pop_mean[p];
                    double scale = pop_scale[p] * w;
                    if (options.algorithm == 4) {
                        const int source = keep_indices.empty() ? i : keep_indices[i];
                        z_col[i] = probs.missing[source] != 0.0f
                            ? 0.0f
                            : static_cast<float>(bgen::posteriorMeanVitezica(
                                  probs.p0[source], probs.p1[source], probs.p2[source],
                                  global_freq) * w);
                    } else if (options.algorithm == 2) {
                        double freq_p = mean / 2.0;
                        if (val > 1.5f) z_col[i] = (float)((1.0 - freq_p) * scale);
                        else if (val < 0.5f) z_col[i] = (float)(-freq_p * scale);
                        else z_col[i] = (float)((1.0 - 2.0 * freq_p) * scale);
                    } else {
                        z_col[i] = (float)((val - mean) * scale);
                    }
                    if (options.compute_homo_hete) {
                        if (val < 0.5f || val > 1.5f) homo_counts[i]++;
                        else hete_counts[i]++;
                    }
                }
            }

            Z_cols.push_back(z_col);
            col_chrom.push_back(cidx);

            if (options.algorithm == 4)
                block_2pq += weight * 4.0 * global_freq * global_freq * (1.0 - global_freq) * (1.0 - global_freq);
            else if (options.algorithm == 5) {
                double p1 = (double)n_het / global_count;
                double p2 = (double)n_hom1 / global_count;
                double mu = p1 + 2.0 * p2;
                double var = p1 + 4.0 * p2 - mu * mu;
                block_2pq += weight * var;
            } else
                block_2pq += weight * 2.0 * global_freq * (1.0 - global_freq);

            if (cidx >= 0) {
                per_chrom_snp_count[cidx]++;
                per_chrom_sum_2pq[cidx] += (options.algorithm == 4)
                    ? weight * 4.0 * global_freq * global_freq * (1.0 - global_freq) * (1.0 - global_freq)
                    : (options.algorithm == 5
                        ? weight * ((double)n_het / global_count + 4.0 * (double)n_hom1 / global_count
                            - std::pow((double)n_het / global_count + 2.0 * (double)n_hom1 / global_count, 2))
                        : weight * 2.0 * global_freq * (1.0 - global_freq));
            }
            allele_freqs[processed_snps - 1] = (float)global_freq;
        }

        if (!Z_cols.empty()) {
            size_t valid_count = Z_cols.size();
            Eigen::MatrixXd Z(n, valid_count);
            #pragma omp parallel for
            for (int k = 0; k < (int)valid_count; ++k)
                Z.col(k) = Z_cols[k].cast<double>();

            G.noalias() += Z * Z.transpose();

            std::map<int, std::vector<int>> chrom_col_lists;
            for (int c = 0; c < (int)valid_count; ++c) {
                int cidx = col_chrom[c];
                if (cidx >= 0) chrom_col_lists[cidx].push_back(c);
            }
            for (auto& kv : chrom_col_lists) {
                int cidx = kv.first;
                auto& cols = kv.second;
                Eigen::MatrixXd Z_sub(n, cols.size());
                for (size_t c = 0; c < cols.size(); ++c)
                    Z_sub.col(c) = Z.col(cols[c]);
                per_chrom_raw[cidx].noalias() += Z_sub * Z_sub.transpose();
            }

            total_valid_snps += valid_count;
            sum_2pq_accum += block_2pq;
        }

        std::cout << "\rProcessed " << processed_snps << " variants (Valid: " << total_valid_snps << ")" << std::flush;
    }
    std::cout << std::endl;

    total_sum_2pq = sum_2pq_accum;

    if (!options.pop_class_file.empty() && num_pops > 1) {
        std::vector<double> pop_S(num_pops, 0.0);
        for (int p = 0; p < num_pops; ++p) {
            double sum_diag = 0.0;
            for (size_t i = 0; i < n; ++i)
                if (ind_pop[i] == p) sum_diag += G(i, i);
            pop_S[p] = (pop_n_samples[p] > 0) ? (sum_diag / pop_n_samples[p]) : 1.0;
            if (pop_S[p] < 1e-9) pop_S[p] = 1.0;
        }
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j <= i; ++j) {
                double divisor = std::sqrt(pop_S[ind_pop[i]] * pop_S[ind_pop[j]]);
                G(i, j) /= divisor;
                if (i != j) G(j, i) = G(i, j);
            }
    } else {
        if (options.algorithm == 0 || options.algorithm == 3) {
            if (total_valid_snps > 0) G /= (double)total_valid_snps;
        } else if (options.algorithm == 2) {
            double tr = G.trace();
            if (tr > 0) G /= (tr / n);
        } else {
            if (total_sum_2pq > 0) G /= total_sum_2pq;
        }
    }

    N_snps.fill((float)total_valid_snps);
    total_valid_snps_count = (long long)total_valid_snps;

    SplitGrmResult res;
    res.global = G;
    res.global_sum_2pq = total_sum_2pq;
    res.global_snp_count = total_valid_snps_count;
    res.bgen_qc_candidates = bgen_qc_candidates;
    res.bgen_qc_rejected_missing = bgen_qc_rejected_missing;
    res.bgen_qc_missing_warnings = bgen_qc_missing_warnings;
    res.bgen_qc_rejected_maf = bgen_qc_rejected_maf;
    res.bgen_qc_rejected_info = bgen_qc_rejected_info;
    res.bgen_qc_info_warnings = bgen_qc_info_warnings;
    res.bgen_qc_entropy_warnings = bgen_qc_entropy_warnings;
    res.chrom_labels = chrom_vec;
    res.per_chrom_raw.resize(n_chrom);
    res.per_chrom_snp_count.resize(n_chrom);
    res.per_chrom_sum_2pq.resize(n_chrom);
    for (size_t k = 0; k < n_chrom; ++k) {
        res.per_chrom_raw[k] = std::move(per_chrom_raw[k]);
        res.per_chrom_snp_count[k] = per_chrom_snp_count[k];
        res.per_chrom_sum_2pq[k] = per_chrom_sum_2pq[k];
    }
    return res;
}

void GrmBuilder::compute(PgenReader& reader, const std::vector<int>& keep_indices) {
    size_t n_total = reader.getNumSamples();
    size_t n = keep_indices.empty() ? n_total : keep_indices.size();
    size_t m = reader.getNumVariants();
    
    if (n == 0) throw std::runtime_error("PGEN reader empty");

    G = Eigen::MatrixXd::Zero(n, n);
    N_snps = Eigen::MatrixXf::Zero(n, n);
    
    const auto& ids = reader.getSampleIds();
    sample_fids.resize(n);
    sample_iids.resize(n);
    for(size_t i=0; i<n; ++i) {
        size_t idx = keep_indices.empty() ? i : keep_indices[i];
        sample_fids[i] = ids[idx];
        sample_iids[i] = ids[idx];
    }
    
    loadAuxFiles();
    ind_pop.assign(n, 0);
    if (num_pops > 0) {
        for (size_t i = 0; i < n; ++i) {
            auto it = pop_map.find(sample_iids[i]);
            if (it != pop_map.end()) ind_pop[i] = it->second;
            else {
                auto it2 = pop_map.find(sample_fids[i] + ":" + sample_iids[i]);
                if (it2 != pop_map.end()) ind_pop[i] = it2->second;
                else ind_pop[i] = 0;
            }
        }
    } else {
        num_pops = 1;
    }
    std::vector<double> pop_n_samples(num_pops, 0.0);
    for (size_t i = 0; i < n; ++i) pop_n_samples[ind_pop[i]] += 1.0;

    total_sum_2pq = 0.0;
    size_t total_valid_snps = 0;
    allele_freqs.assign(m, -1.0f);
    
    if (options.compute_homo_hete) {
        homo_counts.assign(n, 0);
        hete_counts.assign(n, 0);
        missing_counts.assign(n, 0);
    }
    
    reader.reset();
    
    std::cout << "Calculating GRM for " << n << " individuals and " << m << " variants (PGEN)." << std::endl;
    std::cout << "Algorithm: " << ((options.algorithm == 0 || options.algorithm == 3) ? "Yang et al (GCTA default)" : (options.algorithm == 2 ? "Zeng" : (options.algorithm == 4 ? "Vitezica" : "VanRaden"))) << std::endl;
    
    size_t block_size = options.block_size;
    std::vector<Eigen::VectorXf> Z_cols;
    Z_cols.reserve(block_size);
    
    size_t processed_snps = 0;
    double sum_2pq_accum = 0.0;
    
    while (processed_snps < m) {
        Z_cols.clear();
        double block_2pq = 0.0;
        
        for (size_t k = 0; k < block_size && processed_snps < m; ++k) {
            std::string chrom, rsid, a1, a2;
            uint32_t pos;
            Eigen::VectorXf dosage;
            
            if (!reader.readVariant(chrom, rsid, pos, a1, a2, dosage)) {
                processed_snps = m; 
                break;
            }
            processed_snps++;
            
            float weight = 1.0f;
            if (!snp_weights.empty()) {
                auto it = snp_weights.find(rsid);
                if (it != snp_weights.end()) weight = it->second;
                else continue;
            }
            if (weight <= 0.0f) continue;
            
            if (!options.chrom_filter.empty() && chrom != options.chrom_filter) continue;
            
            std::vector<double> pop_sum(num_pops, 0.0);
            std::vector<int> pop_count(num_pops, 0);
            std::vector<int> pop_het(num_pops, 0);
            std::vector<int> pop_hom1(num_pops, 0);
            double global_sum = 0.0;
            int global_count = 0;
            int n_het = 0;
            int n_hom1 = 0;
            
            if (keep_indices.empty()) {
                for(int i=0; i<n; ++i) {
                    if (!std::isnan(dosage[i])) {
                        double val = dosage[i];
                        global_sum += val;
                        global_count++;
                        int p = ind_pop[i];
                        pop_sum[p] += val;
                        pop_count[p]++;
                        if (val > 1.5f) { pop_hom1[p]++; n_hom1++; }
                        else if (val >= 0.5f) { pop_het[p]++; n_het++; }
                    }
                }
            } else {
                for(int i=0; i<n; ++i) {
                    float val = dosage[keep_indices[i]];
                    if (!std::isnan(val)) {
                        global_sum += val;
                        global_count++;
                        int p = ind_pop[i];
                        pop_sum[p] += val;
                        pop_count[p]++;
                        if (val > 1.5f) { pop_hom1[p]++; n_hom1++; }
                        else if (val >= 0.5f) { pop_het[p]++; n_het++; }
                    }
                }
            }
            
            if (global_count == 0) continue;
            
            double global_freq = (global_sum / global_count) / 2.0;
            double global_maf = std::min(global_freq, 1.0 - global_freq);
            
            if (global_maf < options.min_maf || std::abs(global_maf) < 1e-8) continue;

            // Exact HWE on decoded hardcalls, matching the PLINK BED GRM paths.
            if (options.min_hwe > 0.0f) {
                const long long n_hom2 = (long long)global_count - n_het - n_hom1;
                if (calculateHWE(n_hom1, n_het, n_hom2) < options.min_hwe) continue;
            }
            
            std::vector<double> pop_mean(num_pops);
            std::vector<double> pop_scale(num_pops, 1.0);
            for (int p = 0; p < num_pops; ++p) {
                double p_freq = pop_count[p] > 0 ? (pop_sum[p] / pop_count[p]) / 2.0 : global_freq;
                pop_mean[p] = 2.0 * p_freq;
                if (options.algorithm == 0 || options.algorithm == 3) {
                    double pq2 = 2.0 * p_freq * (1.0 - p_freq);
                    pop_scale[p] = (pq2 > 1e-9) ? (1.0 / std::sqrt(pq2)) : 0.0;
                } else if (options.algorithm == 5) {
                    double p1 = pop_count[p] > 0 ? (double)pop_het[p] / pop_count[p] : (double)n_het / global_count;
                    double p2 = pop_count[p] > 0 ? (double)pop_hom1[p] / pop_count[p] : (double)n_hom1 / global_count;
                    double mu = p1 + 2.0 * p2;
                    double var = p1 + 4.0 * p2 - mu * mu;
                    pop_scale[p] = (var > 1e-9) ? (1.0 / std::sqrt(var)) : 0.0;
                }
            }
            
            Eigen::VectorXf z_col(n);
            float w = std::sqrt(weight);
            
            for(int i=0; i<n; ++i) {
                float val = keep_indices.empty() ? dosage[i] : dosage[keep_indices[i]];
                if (std::isnan(val)) {
                    z_col[i] = 0.0f;
                    if (options.compute_homo_hete) missing_counts[i]++;
                } else {
                    int p = ind_pop[i];
                    double mean = pop_mean[p];
                    double scale = pop_scale[p] * w;
                    
                    if (options.algorithm == 4) {
                        double freq_p = mean / 2.0;
                        double freq_q = 1.0 - freq_p;
                        if (val > 1.5f) z_col[i] = (float)(-2.0 * freq_q * freq_q * scale);
                        else if (val < 0.5f) z_col[i] = (float)(-2.0 * freq_p * freq_p * scale);
                        else z_col[i] = (float)(2.0 * freq_p * freq_q * scale);
                    } else if (options.algorithm == 2) {
                        double freq_p = mean / 2.0;
                        if (val > 1.5f) z_col[i] = (float)((1.0 - freq_p) * scale);
                        else if (val < 0.5f) z_col[i] = (float)(-freq_p * scale);
                        else z_col[i] = (float)((1.0 - 2.0 * freq_p) * scale);
                    } else {
                        z_col[i] = (float)((val - mean) * scale);
                    }
                    
                    if (options.compute_homo_hete) {
                        if (val < 0.5f || val > 1.5f) homo_counts[i]++;
                        else hete_counts[i]++;
                    }
                }
            }
            
            Z_cols.push_back(z_col);
            if (options.algorithm == 4) {
                block_2pq += weight * 4.0 * global_freq * global_freq * (1.0 - global_freq) * (1.0 - global_freq);
            } else if (options.algorithm == 5) {
                double p1 = (double)n_het / global_count;
                double p2 = (double)n_hom1 / global_count;
                double mu = p1 + 2.0 * p2;
                double var = p1 + 4.0 * p2 - mu * mu;
                block_2pq += weight * var;
            } else {
                block_2pq += weight * 2.0 * global_freq * (1.0 - global_freq);
            }
            allele_freqs[processed_snps - 1] = (float)global_freq;
        }
        
        if (!Z_cols.empty()) {
            size_t valid_count = Z_cols.size();
            Eigen::MatrixXd Z(n, valid_count);
            
            #pragma omp parallel for
            for(int k=0; k<(int)valid_count; ++k) {
                Z.col(k) = Z_cols[k].cast<double>();
            }
            
            G.noalias() += Z * Z.transpose();
            
            total_valid_snps += valid_count;
            sum_2pq_accum += block_2pq;
        }
        
        std::cout << "\rProcessed " << processed_snps << " variants (Valid: " << total_valid_snps << ")" << std::flush;
    }
    std::cout << std::endl;
    
    total_sum_2pq = sum_2pq_accum;
    
    if (!options.pop_class_file.empty() && num_pops > 1) {
        std::vector<double> pop_S(num_pops, 0.0);
        for (int p = 0; p < num_pops; ++p) {
            double sum_diag = 0.0;
            for (size_t i = 0; i < n; ++i) {
                if (ind_pop[i] == p) sum_diag += G(i, i);
            }
            pop_S[p] = (pop_n_samples[p] > 0) ? (sum_diag / pop_n_samples[p]) : 1.0;
            if (pop_S[p] < 1e-9) pop_S[p] = 1.0;
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                double divisor = std::sqrt(pop_S[ind_pop[i]] * pop_S[ind_pop[j]]);
                G(i, j) /= divisor;
                if (i != j) G(j, i) = G(i, j);
            }
        }
    } else {
        if (options.algorithm == 0 || options.algorithm == 3) {
             if (total_valid_snps > 0) G /= (double)total_valid_snps;
        } else if (options.algorithm == 2) {
             double tr = G.trace();
             if (tr > 0) G /= (tr / n);
        } else {
             if (total_sum_2pq > 0) G /= total_sum_2pq;
        }
    }
    
    N_snps.fill((float)total_valid_snps);
    total_valid_snps_count = (long long)total_valid_snps;
}

GrmBuilder::SplitGrmResult GrmBuilder::computeSplit(
    PgenReader& reader,
    const std::set<std::string>& chrom_labels,
    const std::vector<int>& keep_indices)
{
    if (chrom_labels.empty()) {
        compute(reader, keep_indices);
        SplitGrmResult res;
        res.global = G; res.global_sum_2pq = total_sum_2pq; res.global_snp_count = total_valid_snps_count;
        return res;
    }

    size_t n_total = reader.getNumSamples();
    size_t n = keep_indices.empty() ? n_total : keep_indices.size();
    size_t m = reader.getNumVariants();
    if (n == 0) throw std::runtime_error("PGEN reader empty");

    G = Eigen::MatrixXd::Zero(n, n);
    N_snps = Eigen::MatrixXf::Zero(n, n);

    std::vector<std::string> chrom_vec(chrom_labels.begin(), chrom_labels.end());
    std::map<std::string, int> chrom_to_idx;
    size_t n_chrom = chrom_vec.size();
    for (size_t k = 0; k < n_chrom; ++k) chrom_to_idx[chrom_vec[k]] = (int)k;

    std::vector<Eigen::MatrixXd> per_chrom_raw(n_chrom, Eigen::MatrixXd::Zero(n, n));
    std::vector<long long> per_chrom_snp_count(n_chrom, 0);
    std::vector<double> per_chrom_sum_2pq(n_chrom, 0.0);

    const auto& ids = reader.getSampleIds();
    sample_fids.resize(n); sample_iids.resize(n);
    for (size_t i = 0; i < n; ++i) {
        size_t idx = keep_indices.empty() ? i : keep_indices[i];
        sample_fids[i] = ids[idx]; sample_iids[i] = ids[idx];
    }

    loadAuxFiles();
    ind_pop.assign(n, 0);
    if (num_pops > 0) {
        for (size_t i = 0; i < n; ++i) {
            auto it = pop_map.find(sample_iids[i]);
            if (it != pop_map.end()) ind_pop[i] = it->second;
            else {
                auto it2 = pop_map.find(sample_fids[i] + ":" + sample_iids[i]);
                if (it2 != pop_map.end()) ind_pop[i] = it2->second;
                else ind_pop[i] = 0;
            }
        }
    } else { num_pops = 1; }
    std::vector<double> pop_n_samples(num_pops, 0.0);
    for (size_t i = 0; i < n; ++i) pop_n_samples[ind_pop[i]] += 1.0;

    total_sum_2pq = 0.0;
    size_t total_valid_snps = 0;
    allele_freqs.assign(m, -1.0f);
    if (options.compute_homo_hete) {
        homo_counts.assign(n, 0); hete_counts.assign(n, 0); missing_counts.assign(n, 0);
    }

    reader.reset();
    std::cout << "Calculating global + " << n_chrom << " chromosome GRMs in single pass ("
              << n << " individuals, " << m << " PGEN variants)." << std::endl;

    size_t block_size = options.block_size;
    std::vector<Eigen::VectorXf> Z_cols;
    Z_cols.reserve(block_size);
    std::vector<int> col_chrom;

    size_t processed_snps = 0;
    double sum_2pq_accum = 0.0;

    while (processed_snps < m) {
        Z_cols.clear();
        col_chrom.clear();
        double block_2pq = 0.0;

        for (size_t k = 0; k < block_size && processed_snps < m; ++k) {
            std::string chrom, rsid, a1, a2;
            uint32_t pos;
            Eigen::VectorXf dosage;

            if (!reader.readVariant(chrom, rsid, pos, a1, a2, dosage)) {
                processed_snps = m; break;
            }
            processed_snps++;

            if (options.chrom_filter.empty())

            // Keep GRM SNP membership aligned with the scan filters so that
            // --extract/--exclude constrain the global GRM, the per-chromosome
            // split GRMs, and therefore the LOCO denominators.
            if (!options.extract_set.empty() && !options.extract_set.count(rsid)) continue;
            if (!options.exclude_set.empty() && options.exclude_set.count(rsid)) continue;
                observed_chroms.insert(chrom);

            float weight = 1.0f;
            if (!snp_weights.empty()) {
                auto it = snp_weights.find(rsid);
                if (it != snp_weights.end()) weight = it->second;
                else continue;
            }
            if (weight <= 0.0f) continue;
            if (!options.chrom_filter.empty() && chrom != options.chrom_filter) continue;

            int cidx = -1;
            auto cit = chrom_to_idx.find(chrom);
            if (cit != chrom_to_idx.end()) cidx = cit->second;

            std::vector<double> pop_sum(num_pops, 0.0);
            std::vector<int> pop_count(num_pops, 0);
            std::vector<int> pop_het(num_pops, 0);
            std::vector<int> pop_hom1(num_pops, 0);
            double global_sum = 0.0;
            int global_count = 0;
            int n_het = 0, n_hom1 = 0;

            if (keep_indices.empty()) {
                for (int i = 0; i < n; ++i) {
                    if (!std::isnan(dosage[i])) {
                        double val = dosage[i]; global_sum += val; global_count++;
                        int p = ind_pop[i]; pop_sum[p] += val; pop_count[p]++;
                        if (val > 1.5) { pop_hom1[p]++; n_hom1++; }
                        else if (val >= 0.5) { pop_het[p]++; n_het++; }
                    }
                }
            } else {
                for (int i = 0; i < n; ++i) {
                    float val = dosage[keep_indices[i]];
                    if (!std::isnan(val)) {
                        global_sum += val; global_count++;
                        int p = ind_pop[i]; pop_sum[p] += val; pop_count[p]++;
                        if (val > 1.5f) { pop_hom1[p]++; n_hom1++; }
                        else if (val >= 0.5f) { pop_het[p]++; n_het++; }
                    }
                }
            }

            if (global_count == 0) continue;

            double global_freq = (global_sum / global_count) / 2.0;
            double global_maf = std::min(global_freq, 1.0 - global_freq);
            if (global_maf < options.min_maf || std::abs(global_maf) < 1e-8) continue;

            // Exact HWE on decoded hardcalls, matching the PLINK BED split-GRM path
            // so K_D / K_A SNP membership is identical between BED and PGEN.
            if (options.min_hwe > 0.0f) {
                const long long n_hom2 = (long long)global_count - n_het - n_hom1;
                if (calculateHWE(n_hom1, n_het, n_hom2) < options.min_hwe) continue;
            }

            std::vector<double> pop_mean(num_pops);
            std::vector<double> pop_scale(num_pops, 1.0);
            for (int p = 0; p < num_pops; ++p) {
                double p_freq = pop_count[p] > 0 ? (pop_sum[p] / pop_count[p]) / 2.0 : global_freq;
                pop_mean[p] = 2.0 * p_freq;
                if (options.algorithm == 0 || options.algorithm == 3) {
                    double pq2 = 2.0 * p_freq * (1.0 - p_freq);
                    pop_scale[p] = (pq2 > 1e-9) ? (1.0 / std::sqrt(pq2)) : 0.0;
                } else if (options.algorithm == 5) {
                    double p1 = pop_count[p] > 0 ? (double)pop_het[p] / pop_count[p] : (double)n_het / global_count;
                    double p2 = pop_count[p] > 0 ? (double)pop_hom1[p] / pop_count[p] : (double)n_hom1 / global_count;
                    double mu = p1 + 2.0 * p2;
                    double var = p1 + 4.0 * p2 - mu * mu;
                    pop_scale[p] = (var > 1e-9) ? (1.0 / std::sqrt(var)) : 0.0;
                }
            }

            Eigen::VectorXf z_col(n);
            float w = std::sqrt(weight);
            for (int i = 0; i < n; ++i) {
                float val = keep_indices.empty() ? dosage[i] : dosage[keep_indices[i]];
                if (std::isnan(val)) {
                    z_col[i] = 0.0f;
                    if (options.compute_homo_hete) missing_counts[i]++;
                } else {
                    int p = ind_pop[i];
                    double mean = pop_mean[p];
                    double scale = pop_scale[p] * w;
                    if (options.algorithm == 4) {
                        double freq_p = mean / 2.0;
                        double freq_q = 1.0 - freq_p;
                        if (val > 1.5f) z_col[i] = (float)(-2.0 * freq_q * freq_q * scale);
                        else if (val < 0.5f) z_col[i] = (float)(-2.0 * freq_p * freq_p * scale);
                        else z_col[i] = (float)(2.0 * freq_p * freq_q * scale);
                    } else if (options.algorithm == 2) {
                        double freq_p = mean / 2.0;
                        if (val > 1.5f) z_col[i] = (float)((1.0 - freq_p) * scale);
                        else if (val < 0.5f) z_col[i] = (float)(-freq_p * scale);
                        else z_col[i] = (float)((1.0 - 2.0 * freq_p) * scale);
                    } else {
                        z_col[i] = (float)((val - mean) * scale);
                    }
                    if (options.compute_homo_hete) {
                        if (val < 0.5f || val > 1.5f) homo_counts[i]++;
                        else hete_counts[i]++;
                    }
                }
            }

            Z_cols.push_back(z_col);
            col_chrom.push_back(cidx);

            if (options.algorithm == 4)
                block_2pq += weight * 4.0 * global_freq * global_freq * (1.0 - global_freq) * (1.0 - global_freq);
            else if (options.algorithm == 5) {
                double p1 = (double)n_het / global_count;
                double p2 = (double)n_hom1 / global_count;
                double mu = p1 + 2.0 * p2;
                double var = p1 + 4.0 * p2 - mu * mu;
                block_2pq += weight * var;
            } else
                block_2pq += weight * 2.0 * global_freq * (1.0 - global_freq);

            if (cidx >= 0) {
                per_chrom_snp_count[cidx]++;
                per_chrom_sum_2pq[cidx] += (options.algorithm == 4)
                    ? weight * 4.0 * global_freq * global_freq * (1.0 - global_freq) * (1.0 - global_freq)
                    : (options.algorithm == 5
                        ? weight * ((double)n_het / global_count + 4.0 * (double)n_hom1 / global_count
                            - std::pow((double)n_het / global_count + 2.0 * (double)n_hom1 / global_count, 2))
                        : weight * 2.0 * global_freq * (1.0 - global_freq));
            }
            allele_freqs[processed_snps - 1] = (float)global_freq;
        }

        if (!Z_cols.empty()) {
            size_t valid_count = Z_cols.size();
            Eigen::MatrixXd Z(n, valid_count);
            #pragma omp parallel for
            for (int k = 0; k < (int)valid_count; ++k)
                Z.col(k) = Z_cols[k].cast<double>();

            G.noalias() += Z * Z.transpose();

            std::map<int, std::vector<int>> chrom_col_lists;
            for (int c = 0; c < (int)valid_count; ++c) {
                int cidx = col_chrom[c];
                if (cidx >= 0) chrom_col_lists[cidx].push_back(c);
            }
            for (auto& kv : chrom_col_lists) {
                int cidx = kv.first;
                auto& cols = kv.second;
                Eigen::MatrixXd Z_sub(n, cols.size());
                for (size_t c = 0; c < cols.size(); ++c)
                    Z_sub.col(c) = Z.col(cols[c]);
                per_chrom_raw[cidx].noalias() += Z_sub * Z_sub.transpose();
            }

            total_valid_snps += valid_count;
            sum_2pq_accum += block_2pq;
        }

        std::cout << "\rProcessed " << processed_snps << " variants (Valid: " << total_valid_snps << ")" << std::flush;
    }
    std::cout << std::endl;

    total_sum_2pq = sum_2pq_accum;

    if (!options.pop_class_file.empty() && num_pops > 1) {
        std::vector<double> pop_S(num_pops, 0.0);
        for (int p = 0; p < num_pops; ++p) {
            double sum_diag = 0.0;
            for (size_t i = 0; i < n; ++i)
                if (ind_pop[i] == p) sum_diag += G(i, i);
            pop_S[p] = (pop_n_samples[p] > 0) ? (sum_diag / pop_n_samples[p]) : 1.0;
            if (pop_S[p] < 1e-9) pop_S[p] = 1.0;
        }
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j <= i; ++j) {
                double divisor = std::sqrt(pop_S[ind_pop[i]] * pop_S[ind_pop[j]]);
                G(i, j) /= divisor;
                if (i != j) G(j, i) = G(i, j);
            }
    } else {
        if (options.algorithm == 0 || options.algorithm == 3) {
            if (total_valid_snps > 0) G /= (double)total_valid_snps;
        } else if (options.algorithm == 2) {
            double tr = G.trace();
            if (tr > 0) G /= (tr / n);
        } else {
            if (total_sum_2pq > 0) G /= total_sum_2pq;
        }
    }

    N_snps.fill((float)total_valid_snps);
    total_valid_snps_count = (long long)total_valid_snps;

    SplitGrmResult res;
    res.global = G;
    res.global_sum_2pq = total_sum_2pq;
    res.global_snp_count = total_valid_snps_count;
    res.chrom_labels = chrom_vec;
    res.per_chrom_raw.resize(n_chrom);
    res.per_chrom_snp_count.resize(n_chrom);
    res.per_chrom_sum_2pq.resize(n_chrom);
    for (size_t k = 0; k < n_chrom; ++k) {
        res.per_chrom_raw[k] = std::move(per_chrom_raw[k]);
        res.per_chrom_snp_count[k] = per_chrom_snp_count[k];
        res.per_chrom_sum_2pq[k] = per_chrom_sum_2pq[k];
    }
    return res;
}

void GrmBuilder::save(const std::string& out_prefix) {
    std::string bin_file = out_prefix + ".grm.bin";
    std::string n_file = out_prefix + ".grm.N.bin";
    std::string id_file = out_prefix + ".grm.id";

    // 1. Write .grm.id
    std::ofstream ofs_id(id_file);
    if (!ofs_id) throw std::runtime_error("Cannot create " + id_file);
    for (size_t i = 0; i < sample_fids.size(); ++i) {
        ofs_id << sample_fids[i] << "\t" << sample_iids[i] << "\n";
    }
    ofs_id.close();

    // 2. Write .grm.bin (Lower triangle, float)
    // 3. Write .grm.N.bin (Lower triangle, float)
    std::ofstream ofs_bin(bin_file, std::ios::binary);
    std::ofstream ofs_n(n_file, std::ios::binary);

    if (!ofs_bin || !ofs_n) throw std::runtime_error("Cannot create binary GRM files");

    int n = G.rows();
    
    const size_t buffer_size = 1024 * 1024; // 1M floats
    std::vector<float> g_buffer;
    std::vector<float> n_buffer;
    g_buffer.reserve(buffer_size);
    n_buffer.reserve(buffer_size);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            g_buffer.push_back(G(i, j));
            n_buffer.push_back(N_snps(i, j));

            if (g_buffer.size() >= buffer_size) {
                ofs_bin.write(reinterpret_cast<const char*>(g_buffer.data()), g_buffer.size() * sizeof(float));
                ofs_n.write(reinterpret_cast<const char*>(n_buffer.data()), n_buffer.size() * sizeof(float));
                g_buffer.clear();
                n_buffer.clear();
            }
        }
    }

    // Flush remaining
    if (!g_buffer.empty()) {
        ofs_bin.write(reinterpret_cast<const char*>(g_buffer.data()), g_buffer.size() * sizeof(float));
        ofs_n.write(reinterpret_cast<const char*>(n_buffer.data()), n_buffer.size() * sizeof(float));
    }

    std::cout << "GRM saved to " << out_prefix << ".grm.bin / .id / .N.bin" << std::endl;
}

void GrmBuilder::initSampleIdsFromFam(const PlinkReader& reader) {
    const auto& fam = reader.getFamInfo();
    const size_t n = reader.getNumSamples();
    if (fam.size() != n) {
        throw std::runtime_error("PLINK sample metadata count (" + std::to_string(fam.size()) +
                                 ") does not match genotype sample count (" + std::to_string(n) + ")");
    }
    sample_fids.resize(n);
    sample_iids.resize(n);
    for (size_t i = 0; i < n; ++i) {
        sample_fids[i] = fam[i].fid;
        sample_iids[i] = fam[i].iid;
        if (sample_fids[i].empty() || sample_iids[i].empty()) {
            throw std::runtime_error("PLINK sample " + std::to_string(i) +
                                     " has an empty FID or IID; cannot write sample metadata");
        }
    }
}

void GrmBuilder::computeSparse(const PlinkReader& reader) {
    size_t n = reader.getNumSamples();
    size_t m = reader.getNumSnps();
    if (n == 0 || m == 0) throw std::runtime_error("Empty dataset");

    if (options.sparse_out_prefix.empty()) {
        throw std::runtime_error("Output prefix must be set for sparse GRM computation");
    }

    // Sample identity must be filled before any sparse output is written: the .grm.id
    // file is emitted at the end of the strip loop and is indexed by sample position.
    initSampleIdsFromFam(reader);

    // 1. Pass 1: Compute sum 2p(1-p)
    std::cout << "Pass 1: Computing allele frequencies and total sum 2p(1-p)..." << std::endl;
    total_sum_2pq = 0.0;
    
    const auto& bim = reader.getBimInfo();
    size_t num_blocks = (m + options.block_size - 1) / options.block_size;
    
    #ifdef _OPENMP
    omp_set_num_threads(options.thread_num);
    #endif

    for (size_t b = 0; b < num_blocks; ++b) {
        size_t start = b * options.block_size;
        size_t end = std::min(start + options.block_size, m);
        size_t current_block_size = end - start;

        // Chrom filter
        if (!options.chrom_filter.empty()) {
            bool block_has_chr = false;
            for(size_t k=start; k<end; ++k) {
                if(bim[k].chrom == options.chrom_filter) {
                    block_has_chr = true;
                    break;
                }
            }
            if (!block_has_chr) continue;
        }

        Eigen::MatrixXf raw_block;
        reader.getSnpBlock(start, end, raw_block);

        double block_2pq = 0;
        
        #pragma omp parallel for reduction(+:block_2pq)
        for (long long j = 0; j < (long long)current_block_size; ++j) {
            size_t global_idx = start + j;
            if (!options.chrom_filter.empty() && bim[global_idx].chrom != options.chrom_filter) continue;

            double sum = 0.0;
            int count = 0;
            long long obs_hom1 = 0, obs_hets = 0, obs_hom2 = 0;

            for (size_t i = 0; i < n; ++i) {
                float val = raw_block(i, j);
                if (!std::isnan(val)) {
                    sum += val;
                    count++;
                    if (val > 1.5f) obs_hom1++;
                    else if (val < 0.5f) obs_hom2++;
                    else obs_hets++;
                }
            }

            if (count == 0) continue;

            double freq = (sum / count) / 2.0;
            double maf = std::min(freq, 1.0 - freq);

            if (maf < options.min_maf || std::abs(maf) < 1e-8) continue;

            if (options.min_hwe > 0.0f) {
                double p_hwe = calculateHWE(obs_hom1, obs_hets, obs_hom2);
                if (p_hwe < options.min_hwe) continue;
            }

            block_2pq += 2.0 * freq * (1.0 - freq);
        }
        total_sum_2pq += block_2pq;
        std::cout << "\rProcessed block " << b + 1 << "/" << num_blocks << std::flush;
    }
    std::cout << std::endl;
    std::cout << "Total sum 2p(1-p): " << total_sum_2pq << std::endl;

    if (total_sum_2pq <= 0) {
        throw std::runtime_error("No valid SNPs found or sum 2p(1-p) is zero");
    }

    // 2. Pass 2: Strip-wise Compute
    std::cout << "Pass 2: Computing Sparse GRM (cutoff " << options.sparse_threshold << ")..." << std::endl;
    
    std::string out_file = options.sparse_out_prefix + ".grm.sp";
    std::ofstream ofs(out_file, std::ios::binary);
    if (!ofs) throw std::runtime_error("Cannot create output file " + out_file);

    // Header: Magic "SPGRM", Version 1, N, threshold
    const char magic[] = "SPGRM";
    int32_t version = 1;
    int32_t n_i32 = (int32_t)n;
    ofs.write(magic, 5);
    ofs.write((char*)&version, sizeof(int32_t));
    ofs.write((char*)&n_i32, sizeof(int32_t));
    ofs.write((char*)&options.sparse_threshold, sizeof(float));

    // Iterate strips of samples
    // Choose strip size based on available RAM.
    // Assume 16GB RAM. G_strip (2000 x 500k) is 4GB (float).
    // Let's use 1000 rows per strip (2GB).
    size_t strip_size = 1000; 
    
    for (size_t row_start = 0; row_start < n; row_start += strip_size) {
        size_t row_end = std::min(row_start + strip_size, n);
        size_t current_strip_rows = row_end - row_start;
        
        // Accumulator for this strip: (current_strip_rows x n)
        // We initialize with zeros
        Eigen::MatrixXf G_strip = Eigen::MatrixXf::Zero(current_strip_rows, n);

        // Iterate all SNPs again
        for (size_t b = 0; b < num_blocks; ++b) {
            size_t start = b * options.block_size;
            size_t end = std::min(start + options.block_size, m);
            size_t current_block_size = end - start;

            // Chrom filter check
            if (!options.chrom_filter.empty()) {
                bool block_has_chr = false;
                for(size_t k=start; k<end; ++k) {
                    if(bim[k].chrom == options.chrom_filter) {
                        block_has_chr = true;
                        break;
                    }
                }
                if (!block_has_chr) continue;
            }

            Eigen::MatrixXf raw_block;
            reader.getSnpBlock(start, end, raw_block); // (n x current_block_size)
            
            // Convert to centered Z for valid SNPs
            // Similar logic as Pass 1 but we need to store Z
            std::vector<int> valid_cols;
            std::vector<float> p_vec;
            
            for (size_t j = 0; j < current_block_size; ++j) {
                size_t global_idx = start + j;
                if (!options.chrom_filter.empty() && bim[global_idx].chrom != options.chrom_filter) continue;

                double sum = 0;
                int count = 0;
                // ... same stats ...
                // Optimization: store stats from Pass 1? Too much memory for 1M SNPs?
                // Just recompute. It's fast.
                for (size_t i = 0; i < n; ++i) {
                    float val = raw_block(i, j);
                    if (!std::isnan(val)) { sum += val; count++; }
                }
                if (count == 0) continue;
                double freq = (sum / count) / 2.0;
                double maf = std::min(freq, 1.0 - freq);
                if (maf < options.min_maf || std::abs(maf) < 1e-8) continue;
                // HWE check skipped for speed or re-check?
                // Let's assume passed.
                
                valid_cols.push_back(j);
                p_vec.push_back((float)freq);
            }

            if (valid_cols.empty()) continue;

            size_t valid_count = valid_cols.size();
            Eigen::MatrixXf Z(n, valid_count);
            
            #pragma omp parallel for collapse(2)
            for (long long j_idx = 0; j_idx < (long long)valid_count; ++j_idx) {
                for (long long i = 0; i < (long long)n; ++i) {
                    float val = raw_block(i, valid_cols[j_idx]);
                    float p = p_vec[j_idx];
                    if (std::isnan(val)) Z(i, j_idx) = 0.0f;
                    else Z(i, j_idx) = val - 2.0f * p;
                }
            }

            // Update G_strip
            // G_strip += Z[row_start:row_end, :] * Z.transpose()
            // Z_sub is (current_strip_rows x valid_count)
            Eigen::MatrixXf Z_sub = Z.block(row_start, 0, current_strip_rows, valid_count);
            
            // Parallel matrix multiplication is handled by Eigen/BLAS
            G_strip.noalias() += Z_sub * Z.transpose();
        }

        // Normalize
        G_strip /= (float)total_sum_2pq;

        // Threshold and write
        // We only write lower triangle (col <= row)
        // Global row index: row_start + r
        // Column index: c
        // Condition: c <= row_start + r
        
        std::vector<int32_t> rows_buf;
        std::vector<int32_t> cols_buf;
        std::vector<float> vals_buf;
        rows_buf.reserve(10000);
        cols_buf.reserve(10000);
        vals_buf.reserve(10000);

        for(size_t r = 0; r < current_strip_rows; ++r) {
            size_t global_row = row_start + r;
            for(size_t c = 0; c <= global_row; ++c) {
                float val = G_strip(r, c);
                if (std::abs(val) >= options.sparse_threshold) {
                    rows_buf.push_back((int32_t)global_row);
                    cols_buf.push_back((int32_t)c);
                    vals_buf.push_back(val);
                }
            }
        }
        
        if (!rows_buf.empty()) {
             size_t count = rows_buf.size();
             // Write block of triplets
             // Or write individually? Let's write array of triplets.
             // Format: row, col, val repeated.
             for(size_t k=0; k<count; ++k) {
                 ofs.write((char*)&rows_buf[k], sizeof(int32_t));
                 ofs.write((char*)&cols_buf[k], sizeof(int32_t));
                 ofs.write((char*)&vals_buf[k], sizeof(float));
             }
        }

        std::cout << "\rProcessed strip rows " << row_end << "/" << n << std::flush;
    }
    std::cout << std::endl;
    
    ofs.close();
    
    // Save IDs
    if (sample_fids.size() != n || sample_iids.size() != n) {
        throw std::runtime_error("Sparse GRM sample metadata size (" + std::to_string(sample_fids.size()) +
                                 ") does not match the computed sample count (" + std::to_string(n) + ")");
    }
    std::ofstream ofs_id(options.sparse_out_prefix + ".grm.id");
    if (!ofs_id) throw std::runtime_error("Cannot create output file " + options.sparse_out_prefix + ".grm.id");
    for(size_t i=0; i<n; ++i) ofs_id << sample_fids[i] << "\t" << sample_iids[i] << "\n";
    ofs_id.close();
    if (!ofs_id) throw std::runtime_error("Failed while writing " + options.sparse_out_prefix + ".grm.id");

    std::cout << "Sparse GRM saved to " << out_file << std::endl;
}

void GrmBuilder::computeSparse(BgenReader& reader) {
    throw std::runtime_error("Sparse GRM for BGEN not yet implemented");
}

void GrmBuilder::computeSparse(PgenReader& reader) {
    throw std::runtime_error("Sparse GRM for PGEN not yet implemented (use PLINK format or wait for update)");
}

void GrmBuilder::compute(const VcfReader& reader, const std::vector<int>& keep_indices) {
    size_t n = reader.getNumSamples();
    size_t m = reader.getNumVariants();
    if (n == 0 || m == 0) throw std::runtime_error("Empty dataset provided to GrmBuilder");

    if (options.sparse_out_prefix.empty() && options.sparse) {
        throw std::runtime_error("Output prefix must be set for sparse GRM computation");
    }

    // Apply keep_indices filter
    std::vector<int> sample_indices;
    if (keep_indices.empty()) {
        sample_indices.resize(n);
        std::iota(sample_indices.begin(), sample_indices.end(), 0);
    } else {
        sample_indices = keep_indices;
    }
    size_t n_keep = sample_indices.size();

    // Allocate GRM matrix
    G.resize(n_keep, n_keep);
    G.setZero();
    
    // Initialize N_snps matrix (required for save())
    N_snps = Eigen::MatrixXf::Zero(n_keep, n_keep);
    
    total_sum_2pq = 0.0;
    total_valid_snps_count = 0;

    // Load all variants into memory first (VCF is typically small enough)
    std::vector<VcfVariant> variants;
    std::vector<Eigen::VectorXd> dosages;
    variants.reserve(reader.getNumVariants());
    
    // Read all variants
    VcfReader& mutable_reader = const_cast<VcfReader&>(reader);
    mutable_reader.reset();
    
    std::string chrom, id, ref, alt;
    uint32_t pos = 0;
    Eigen::VectorXf dosage_f;
    
    while (reader.readVariant(chrom, id, pos, ref, alt, dosage_f)) {
        // Apply MAF/HWE filters
        double sum = 0.0;
        int obs = 0;
        Eigen::VectorXd dosage(n);
        for (int i = 0; i < dosage_f.size(); ++i) {
            float val = dosage_f[i];
            if (std::isfinite(val) && val >= 0.0 && val <= 2.0) {
                sum += val;
                ++obs;
                dosage[i] = val;
            } else {
                dosage[i] = std::numeric_limits<double>::quiet_NaN();
            }
        }
        if (obs < 3) continue;
        double freq = sum / (2.0 * obs);
        double maf = std::min(freq, 1.0 - freq);
        if (maf < options.min_maf) continue;
        
        // Apply extract/exclude filters
        if (!options.extract_set.empty() && !options.extract_set.count(id)) continue;
        if (!options.exclude_set.empty() && options.exclude_set.count(id)) continue;
        
        VcfVariant v;
        v.chrom = chrom;
        v.id = id;
        v.pos = pos;
        v.ref = ref;
        v.alt = alt;
        variants.push_back(v);
        dosages.push_back(std::move(dosage));
        total_valid_snps_count++;
    }
    
    if (variants.empty()) throw std::runtime_error("No valid variants after filtering");
    
    // Compute GRM
    G.setZero();
    size_t valid_m = variants.size();
    
    // Process in blocks to manage memory
    size_t block_size = options.block_size;
    for (size_t block_start = 0; block_start < valid_m; block_start += block_size) {
        size_t block_end = std::min(block_start + block_size, valid_m);
        size_t block_m = block_end - block_start;
        
        // Build genotype matrix for this block
        Eigen::MatrixXd Z(n, block_m);
        double block_sum_2pq = 0.0;
        int block_valid = 0;
        
        for (size_t j = 0; j < block_m; ++j) {
            const auto& dosage = dosages[block_start + j];
            const auto& variant = variants[block_start + j];
            
            // Compute allele frequency from observed dosages
            double sum = 0.0;
            int obs = 0;
            for (size_t i = 0; i < n; ++i) {
                if (std::isfinite(dosage[i]) && dosage[i] >= 0.0 && dosage[i] <= 2.0) {
                    sum += dosage[i];
                    ++obs;
                }
            }
            if (obs < 3) continue;
            double freq = sum / (2.0 * obs);
            double maf = std::min(freq, 1.0 - freq);
            if (maf < options.min_maf) continue;
            
            // Compute scale and Z based on algorithm (matching PLINK implementation)
            double scale = 1.0;
            double pq2 = 2.0 * freq * (1.0 - freq);
            if (options.algorithm == 0 || options.algorithm == 3) {
                // Yang et al (GCTA default): Z = (x - 2p) / sqrt(2pq)
                scale = (pq2 > 1e-9) ? (1.0 / std::sqrt(pq2)) : 0.0;
            } else if (options.algorithm == 2) {
                // Zeng: Z = (x - 2p) * w
                scale = 1.0;
            } else if (options.algorithm == 4) {
                // Vitezica: Z = (x - 2p) / sqrt(2pq)
                scale = (pq2 > 1e-9) ? (1.0 / std::sqrt(pq2)) : 0.0;
            } else if (options.algorithm == 5) {
                // TrueVar: use observed variance
                double p1 = 0.0, p2 = 0.0;
                for (size_t i = 0; i < n; ++i) {
                    if (std::isfinite(dosage[i])) {
                        if (dosage[i] == 1.0) p1 += 1.0;
                        else if (dosage[i] == 2.0) p2 += 1.0;
                    }
                }
                p1 /= obs; p2 /= obs;
                double mu = p1 + 2.0 * p2;
                double var = p1 + 4.0 * p2 - mu * mu;
                scale = (var > 1e-9) ? (1.0 / std::sqrt(var)) : 0.0;
            } else {
                // VanRaden (algorithm 1): PLINK uses Z = 2(1-p) - x
                // This is different from standard VanRaden (x - 2p)
                scale = 1.0;
            }
            
            double mean = 2.0 * freq;
            
            for (size_t i = 0; i < n; ++i) {
                if (!std::isfinite(dosage[i])) {
                    // Missing data: PLINK sets Z=0 for all algorithms
                    Z(i, block_valid) = 0.0;
                } else if (options.algorithm == 1) {
                    Z(i, block_valid) = (dosage[i] - mean) * scale;
                } else {
                    // Standard: Z = (x - 2p) / scale_denom
                    Z(i, block_valid) = (dosage[i] - mean) * scale;
                }
            }
            
            block_sum_2pq += pq2;
            block_valid++;
        }
        
        if (block_valid > 0) {
            // Resize Z to only valid columns
            Eigen::MatrixXd Z_valid = Z.leftCols(block_valid);
            
            // Z * Z' contribution to GRM
            G.noalias() += Z_valid * Z_valid.transpose();
            total_sum_2pq += block_sum_2pq;
        }
        
        if (options.sparse) {
            // Write sparse block
        }
    }
    
    // Normalize
    if (total_sum_2pq > 0) {
        G /= total_sum_2pq;
    }
    
    // Fill N_snps (required for save())
    if (total_valid_snps_count > 0) {
        N_snps.fill((float)total_valid_snps_count);
    }
    
    // Extract keep_indices
    if (!keep_indices.empty() && keep_indices.size() != n) {
        Eigen::MatrixXd G_sub(keep_indices.size(), keep_indices.size());
        Eigen::MatrixXf N_sub(keep_indices.size(), keep_indices.size());
        for (size_t i = 0; i < keep_indices.size(); ++i) {
            for (size_t j = 0; j < keep_indices.size(); ++j) {
                G_sub(i, j) = G(keep_indices[i], keep_indices[j]);
                N_sub(i, j) = N_snps(keep_indices[i], keep_indices[j]);
            }
        }
        G = std::move(G_sub);
        N_snps = std::move(N_sub);
    }
    
    // Store sample IDs
    sample_fids.resize(n_keep);
    sample_iids.resize(n_keep);
    const auto& sample_ids = reader.getSampleIds();
    for (size_t i = 0; i < n_keep; ++i) {
        size_t idx = keep_indices.empty() ? i : keep_indices[i];
        size_t sep = sample_ids[idx].find('_');
        if (sep == std::string::npos) sep = sample_ids[idx].find(':');
        if (sep != std::string::npos) {
            sample_fids[i] = sample_ids[idx].substr(0, sep);
            sample_iids[i] = sample_ids[idx].substr(sep + 1);
        } else {
            sample_fids[i] = sample_ids[idx];
            sample_iids[i] = sample_ids[idx];
        }
    }
}

void GrmBuilder::computeSparse(const VcfReader& reader) {
    throw std::runtime_error("Sparse GRM for VCF not yet implemented");
}

GrmBuilder::SplitGrmResult GrmBuilder::computeSplit(const VcfReader& reader,
                                                    const std::set<std::string>& chrom_labels,
                                                    const std::vector<int>& keep_indices) {
    if (chrom_labels.empty()) {
        compute(reader, keep_indices);
        SplitGrmResult res;
        res.global = G;
        res.global_sum_2pq = total_sum_2pq;
        res.global_snp_count = total_valid_snps_count;
        return res;
    }

    size_t n_total = reader.getNumSamples();
    size_t m = reader.getNumVariants();
    if (m == 0) {
        throw std::runtime_error("Empty VCF dataset provided to GrmBuilder::computeSplit");
    }

    const std::vector<int>& effective_keep = keep_indices.empty() ? 
        [&]() {
            std::vector<int> all(n_total);
            std::iota(all.begin(), all.end(), 0);
            static thread_local std::vector<int> cached;
            cached = std::move(all);
            return cached;
        }() : keep_indices;
    size_t n = effective_keep.size();

    // Initialize global GRM
    G = Eigen::MatrixXd::Zero(n, n);
    N_snps = Eigen::MatrixXf::Zero(n, n);

    // Per-chromosome accumulators
    std::vector<std::string> chrom_vec(chrom_labels.begin(), chrom_labels.end());
    std::map<std::string, int> chrom_to_idx;
    size_t n_chrom = chrom_vec.size();
    for (size_t k = 0; k < n_chrom; ++k)
        chrom_to_idx[chrom_vec[k]] = (int)k;

    std::vector<Eigen::MatrixXd> per_chrom_raw(n_chrom, Eigen::MatrixXd::Zero(n, n));
    std::vector<long long> per_chrom_snp_count(n_chrom, 0);
    std::vector<double> per_chrom_sum_2pq(n_chrom, 0.0);

    // Copy IDs
    const auto& sample_ids = reader.getSampleIds();
    sample_fids.resize(n);
    sample_iids.resize(n);
    for (size_t i = 0; i < n; ++i) {
        size_t idx = effective_keep[i];
        const std::string& id = sample_ids[idx];
        size_t sep = id.find(':');
        if (sep != std::string::npos) {
            sample_fids[i] = id.substr(0, sep);
            sample_iids[i] = id.substr(sep + 1);
        } else {
            sample_fids[i] = id;
            sample_iids[i] = id;
        }
    }

    loadAuxFiles();
    ind_pop.assign(n, 0);
    if (num_pops > 0) {
        for (size_t i = 0; i < n; ++i) {
            auto it = pop_map.find(sample_iids[i]);
            if (it != pop_map.end()) ind_pop[i] = it->second;
            else {
                auto it2 = pop_map.find(sample_fids[i] + ":" + sample_iids[i]);
                if (it2 != pop_map.end()) ind_pop[i] = it2->second;
                else ind_pop[i] = 0;
            }
        }
    } else {
        num_pops = 1;
    }
    std::vector<double> pop_n_samples(num_pops, 0.0);
    for (size_t i = 0; i < n; ++i) pop_n_samples[ind_pop[i]] += 1.0;

    total_sum_2pq = 0.0;
    size_t total_valid_snps = 0;
    allele_freqs.assign(m, -1.0f);

    if (options.compute_homo_hete) {
        homo_counts.assign(n, 0);
        hete_counts.assign(n, 0);
        missing_counts.assign(n, 0);
    }

    // Ensure variant info is loaded
    reader.loadVariantInfo();
    const auto& variants = reader.getVariants();

    std::cout << "Calculating global + " << n_chrom << " chromosome GRMs from VCF ("
              << n << " individuals, " << m << " variants)." << std::endl;

    // Process variants in blocks for efficiency
    size_t block_size = options.block_size;
    size_t num_blocks = (m + block_size - 1) / block_size;

    for (size_t b = 0; b < num_blocks; ++b) {
        size_t start = b * block_size;
        size_t end = std::min(start + block_size, m);
        size_t current_block_size = end - start;

        // First Pass: Compute Stats (Mean, Variance, Filter)
        std::vector<float> snp_means(current_block_size * num_pops);
        std::vector<char> snp_valid(current_block_size, 0);
        std::vector<float> snp_scales(current_block_size * num_pops, 1.0f);
        std::vector<float> snp_2pq(current_block_size, 0.0f);
        std::vector<float> snp_weights_block(current_block_size, 1.0f);
        std::vector<int> snp_chrom_idx(current_block_size, -1);

        int valid_count_in_block = 0;

        #pragma omp parallel for reduction(+:valid_count_in_block)
        for (long long j = 0; j < (long long)current_block_size; ++j) {
            size_t global_idx = start + j;
            if (global_idx >= variants.size()) continue;

            const auto& variant = variants[global_idx];

            // Chrom filter
            if (!options.chrom_filter.empty() && variant.chrom != options.chrom_filter) continue;

            float weight = 1.0f;
            if (!snp_weights.empty()) {
                auto it = snp_weights.find(variant.id);
                if (it != snp_weights.end()) weight = it->second;
                else continue;
            }
            if (weight <= 0.0f) continue;
            snp_weights_block[j] = weight;

            // Parse genotypes and compute stats
            std::vector<long long> pop_sum_x(num_pops, 0);
            std::vector<long long> pop_n_obs(num_pops, 0);
            std::vector<long long> pop_n_het(num_pops, 0);
            std::vector<long long> pop_n_hom1(num_pops, 0);
            long long global_sum_x = 0;
            long long global_n_obs = 0;
            long long n_hom1 = 0;
            long long n_het = 0;
            long long n_hom2 = 0;

            for (size_t si = 0; si < n; ++si) {
                size_t sample_idx = effective_keep[si];
                if (sample_idx >= variant.genotypes.size()) continue;
                const std::string& gt = variant.genotypes[sample_idx];
                if (gt.empty() || gt[0] == '.') continue;

                // Parse GT field: 0/1, 1|0, 0/0, 1/1, etc.
                int alt_count = 0;
                if (gt.size() >= 3) {
                    char a1 = gt[0];
                    char a2 = gt[2];
                    if (a1 == '1') alt_count++;
                    if (a2 == '1') alt_count++;
                } else if (gt.size() == 1) {
                    if (gt[0] == '1') alt_count = 2;
                    else if (gt[0] == '0') alt_count = 0;
                    else continue;
                } else continue;

                int p = ind_pop[si];
                pop_n_obs[p]++;
                global_n_obs++;
                if (alt_count == 2) { pop_sum_x[p] += 2; pop_n_hom1[p]++; global_sum_x += 2; n_hom1++; }
                else if (alt_count == 1) { pop_sum_x[p] += 1; pop_n_het[p]++; global_sum_x += 1; n_het++; }
                else { n_hom2++; }
            }

            if (global_n_obs == 0) continue;

            double global_freq = (double)global_sum_x / (2.0 * global_n_obs);
            double global_maf = std::min(global_freq, 1.0 - global_freq);

            if (global_maf < options.min_maf || std::abs(global_maf) < 1e-8) continue;

            if (options.min_hwe > 0.0f) {
                double p_hwe = calculateHWE(n_hom1, n_het, n_hom2);
                if (p_hwe < options.min_hwe) continue;
            }

            if (options.algorithm == 4) {
                snp_2pq[j] = (float)(weight * 4.0 * global_freq * global_freq * (1.0 - global_freq) * (1.0 - global_freq));
            } else if (options.algorithm == 5) {
                double p1 = (double)n_het / global_n_obs;
                double p2 = (double)n_hom1 / global_n_obs;
                double mu = p1 + 2.0 * p2;
                double var = p1 + 4.0 * p2 - mu * mu;
                snp_2pq[j] = (float)(weight * var);
            } else {
                snp_2pq[j] = (float)(weight * 2.0 * global_freq * (1.0 - global_freq));
            }
            allele_freqs[global_idx] = (float)global_freq;

            for (int p = 0; p < num_pops; ++p) {
                double freq = pop_n_obs[p] > 0 ? (double)pop_sum_x[p] / (2.0 * pop_n_obs[p]) : global_freq;
                snp_means[j * num_pops + p] = (float)(freq * 2.0);
                if (options.algorithm == 0 || options.algorithm == 3) {
                    double pq2 = 2.0 * freq * (1.0 - freq);
                    snp_scales[j * num_pops + p] = (float)((pq2 > 1e-9) ? (1.0 / std::sqrt(pq2)) : 0.0);
                } else if (options.algorithm == 5) {
                    double p1 = pop_n_obs[p] > 0 ? (double)pop_n_het[p] / pop_n_obs[p] : (double)n_het / global_n_obs;
                    double p2 = pop_n_obs[p] > 0 ? (double)pop_n_hom1[p] / pop_n_obs[p] : (double)n_hom1 / global_n_obs;
                    double mu = p1 + 2.0 * p2;
                    double var = p1 + 4.0 * p2 - mu * mu;
                    snp_scales[j * num_pops + p] = (float)((var > 1e-9) ? (1.0 / std::sqrt(var)) : 0.0);
                } else {
                    snp_scales[j * num_pops + p] = 1.0f;
                }
            }

            snp_valid[j] = true;
            auto it = chrom_to_idx.find(variant.chrom);
            if (it != chrom_to_idx.end())
                snp_chrom_idx[j] = it->second;
            valid_count_in_block++;
        }

        if (valid_count_in_block == 0) continue;

        // Second Pass: Fill Z Matrix
        Eigen::MatrixXd Z(n, valid_count_in_block);

        std::vector<int> block_to_valid(current_block_size, -1);
        int cur_valid = 0;
        double block_sum_2pq = 0.0;

        for (size_t j = 0; j < current_block_size; ++j) {
            if (snp_valid[j]) {
                block_to_valid[j] = cur_valid++;
                block_sum_2pq += snp_2pq[j];
            }
        }

        #pragma omp parallel for
        for (long long j = 0; j < (long long)current_block_size; ++j) {
            if (!snp_valid[j]) continue;

            int valid_idx = block_to_valid[j];
            float w = std::sqrt(snp_weights_block[j]);
            const auto& variant = variants[start + j];

            std::vector<std::array<double, 4>> pop_val_map(num_pops);
            for (int p = 0; p < num_pops; ++p) {
                double mean = snp_means[j * num_pops + p];
                double scale = snp_scales[j * num_pops + p] * w;
                if (options.algorithm == 4) {
                    double freq_p = mean / 2.0;
                    double freq_q = 1.0 - freq_p;
                    pop_val_map[p][0] = (-2.0 * freq_q * freq_q) * scale;
                    pop_val_map[p][1] = 0.0;
                    pop_val_map[p][2] = (2.0 * freq_p * freq_q) * scale;
                    pop_val_map[p][3] = (-2.0 * freq_p * freq_p) * scale;
                } else if (options.algorithm == 2) {
                    double freq_p = mean / 2.0;
                    pop_val_map[p][0] = (1.0 - freq_p) * scale;
                    pop_val_map[p][1] = 0.0;
                    pop_val_map[p][2] = (1.0 - 2.0 * freq_p) * scale;
                    pop_val_map[p][3] = (-freq_p) * scale;
                } else if (options.algorithm == 1) {
                    pop_val_map[p][0] = (0.0 - mean) * scale;
                    pop_val_map[p][1] = (1.0 - mean) * scale;
                    pop_val_map[p][2] = (2.0 - mean) * scale;
                    pop_val_map[p][3] = 0.0;
                } else {
                    pop_val_map[p][0] = (2.0 - mean) * scale;
                    pop_val_map[p][1] = 0.0;
                    pop_val_map[p][2] = (1.0 - mean) * scale;
                    pop_val_map[p][3] = (0.0 - mean) * scale;
                }
            }

            double* z_col_ptr = Z.col(valid_idx).data();

            for (size_t si = 0; si < n; ++si) {
                size_t sample_idx = effective_keep[si];
                double val = 0.0;
                if (sample_idx < variant.genotypes.size()) {
                    const std::string& gt = variant.genotypes[sample_idx];
                    if (gt.empty() || gt[0] == '.') {
                        z_col_ptr[si] = 0.0;
                        continue;
                    }
                    int alt_count = 0;
                    if (gt.size() >= 3) {
                        char a1 = gt[0];
                        char a2 = gt[2];
                        if (a1 == '1') alt_count++;
                        if (a2 == '1') alt_count++;
                    } else if (gt.size() == 1) {
                        if (gt[0] == '1') alt_count = 2;
                        else if (gt[0] == '0') alt_count = 0;
                        else alt_count = 1;
                    }
                    val = pop_val_map[ind_pop[si]][alt_count];
                } else {
                    val = 0.0;
                }
                z_col_ptr[si] = val;
            }
        }

        // Accumulate global GRM
        G.noalias() += Z * Z.transpose();

        // For each chromosome, build Z_sub and accumulate into per-chromosome raw G
        std::map<int, std::vector<int>> chrom_col_lists;
        int col = 0;
        for (size_t j = 0; j < current_block_size; ++j) {
            if (snp_valid[j]) {
                int cidx = snp_chrom_idx[j];
                if (cidx >= 0)
                    chrom_col_lists[cidx].push_back(col);
                col++;
            }
        }

        for (auto& kv : chrom_col_lists) {
            int cidx = kv.first;
            auto& cols = kv.second;
            if (cols.empty()) continue;

            Eigen::MatrixXd Z_sub(n, cols.size());
            for (size_t c = 0; c < cols.size(); ++c)
                Z_sub.col(c) = Z.col(cols[c]);

            per_chrom_raw[cidx].noalias() += Z_sub * Z_sub.transpose();
        }

        // Update per-chromosome stats
        for (size_t j = 0; j < current_block_size; ++j) {
            if (snp_valid[j]) {
                int cidx = snp_chrom_idx[j];
                if (cidx >= 0) {
                    per_chrom_snp_count[cidx]++;
                    per_chrom_sum_2pq[cidx] += snp_2pq[j];
                }
            }
        }

        total_sum_2pq += block_sum_2pq;
        total_valid_snps += valid_count_in_block;

        std::cout << "\rProcessed block " << b + 1 << "/" << num_blocks
                  << " (Valid SNPs: " << total_valid_snps << ")" << std::flush;
    }
    std::cout << std::endl;

    // Normalize global GRM
    if (!options.pop_class_file.empty() && num_pops > 1) {
        std::vector<double> pop_S(num_pops, 0.0);
        for (int p = 0; p < num_pops; ++p) {
            double sum_diag = 0.0;
            for (size_t i = 0; i < n; ++i)
                if (ind_pop[i] == p) sum_diag += G(i, i);
            pop_S[p] = (pop_n_samples[p] > 0) ? (sum_diag / pop_n_samples[p]) : 1.0;
            if (pop_S[p] < 1e-9) pop_S[p] = 1.0;
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                double divisor = std::sqrt(pop_S[ind_pop[i]] * pop_S[ind_pop[j]]);
                G(i, j) /= divisor;
                if (i != j) G(j, i) = G(i, j);
            }
        }
    } else {
        if (options.algorithm == 0 || options.algorithm == 3) {
            if (total_valid_snps > 0) G /= (double)total_valid_snps;
        } else if (options.algorithm == 2) {
            double tr = G.trace();
            if (tr > 0) G /= (tr / n);
        } else {
            if (total_sum_2pq > 0) G /= total_sum_2pq;
        }
    }

    N_snps.fill((float)total_valid_snps);
    total_valid_snps_count = (long long)total_valid_snps;

    // Assemble result
    SplitGrmResult res;
    res.global = G;
    res.global_sum_2pq = total_sum_2pq;
    res.global_snp_count = total_valid_snps_count;
    res.chrom_labels = chrom_vec;
    res.per_chrom_raw.resize(n_chrom);
    res.per_chrom_snp_count.resize(n_chrom);
    res.per_chrom_sum_2pq.resize(n_chrom);
    for (size_t k = 0; k < n_chrom; ++k) {
        res.per_chrom_raw[k] = std::move(per_chrom_raw[k]);
        res.per_chrom_snp_count[k] = per_chrom_snp_count[k];
        res.per_chrom_sum_2pq[k] = per_chrom_sum_2pq[k];
    }

    return res;
}

GrmBuilder::SplitGrmResult GrmBuilder::computeSplit(
    SequentialGenotypeMatrix& reader,
    const std::set<std::string>& requested_chroms,
    const std::vector<int>& keep_indices) {
    if (options.algorithm != 1) {
        throw std::runtime_error("HapMap/text GRM currently supports only the default VanRaden algorithm");
    }
    if (!options.snp_weight_file.empty() || !options.pop_class_file.empty()) {
        throw std::runtime_error("SNP weights and population-specific GRM are not supported for HapMap/text input");
    }

    const size_t n_total = reader.getNumSamples();
    const size_t n = keep_indices.empty() ? n_total : keep_indices.size();
    if (n == 0) throw std::runtime_error("HapMap/text genotype reader is empty");
    for (int index : keep_indices) {
        if (index < 0 || static_cast<size_t>(index) >= n_total) {
            throw std::runtime_error("HapMap/text GRM sample index is out of range");
        }
    }

    std::set<std::string> chrom_set = requested_chroms;
    if (chrom_set.empty()) {
        reader.reset();
        std::string chrom, id, ref, alt;
        uint32_t pos = 0;
        Eigen::VectorXf dosage;
        while (reader.readVariant(chrom, id, pos, ref, alt, dosage)) chrom_set.insert(chrom);
    }
    std::vector<std::string> chrom_vec(chrom_set.begin(), chrom_set.end());
    std::unordered_map<std::string, size_t> chrom_index;
    for (size_t i = 0; i < chrom_vec.size(); ++i) chrom_index[chrom_vec[i]] = i;

    Eigen::MatrixXd global_raw = Eigen::MatrixXd::Zero(n, n);
    std::vector<Eigen::MatrixXd> per_chrom_raw(
        chrom_vec.size(), Eigen::MatrixXd::Zero(n, n));
    std::vector<long long> per_chrom_count(chrom_vec.size(), 0);
    std::vector<double> per_chrom_sum(chrom_vec.size(), 0.0);
    allele_freqs.assign(reader.getNumVariants(), -1.0f);
    observed_chroms.clear();

    const auto& ids = reader.getSampleIds();
    sample_fids.resize(n);
    sample_iids.resize(n);
    for (size_t i = 0; i < n; ++i) {
        const size_t source = keep_indices.empty() ? i : static_cast<size_t>(keep_indices[i]);
        const size_t separator = ids[source].find(':');
        sample_fids[i] = separator == std::string::npos ? ids[source] : ids[source].substr(0, separator);
        sample_iids[i] = separator == std::string::npos ? ids[source] : ids[source].substr(separator + 1);
    }

    total_sum_2pq = 0.0;
    total_valid_snps_count = 0;
    reader.reset();
    std::string chrom, id, ref, alt;
    uint32_t pos = 0;
    Eigen::VectorXf dosage;
    size_t variant_index = 0;
    while (reader.readVariant(chrom, id, pos, ref, alt, dosage)) {
        const size_t current_index = variant_index++;
        observed_chroms.insert(chrom);
        if (!options.chrom_filter.empty() && chrom != options.chrom_filter) continue;
        if (!options.extract_set.empty() && !options.extract_set.count(id)) continue;
        if (options.exclude_set.count(id)) continue;
        auto chrom_found = chrom_index.find(chrom);
        if (chrom_found == chrom_index.end()) continue;
        if (dosage.size() != static_cast<int>(n_total)) {
            throw std::runtime_error("HapMap/text variant dosage length does not match sample count");
        }

        double sum = 0.0;
        int count = 0;
        long long hom0 = 0, het = 0, hom2 = 0;
        for (size_t i = 0; i < n; ++i) {
            const size_t source = keep_indices.empty() ? i : static_cast<size_t>(keep_indices[i]);
            if (std::isfinite(dosage[source])) {
                const double value = dosage[source];
                sum += value;
                ++count;
                const long long hardcall = std::llround(value);
                if (hardcall == 0) ++hom0;
                else if (hardcall == 1) ++het;
                else ++hom2;
            }
        }
        if (count == 0) continue;
        const double frequency = sum / (2.0 * count);
        const double maf = std::min(frequency, 1.0 - frequency);
        const double sum_2pq = 2.0 * frequency * (1.0 - frequency);
        if (maf < options.min_maf || sum_2pq <= 1e-12) continue;
        if (options.min_hwe > 0.0f && calculateHWE(hom0, het, hom2) < options.min_hwe) continue;
        allele_freqs[current_index] = static_cast<float>(frequency);

        Eigen::VectorXd centered(n);
        const double mean = 2.0 * frequency;
        for (size_t i = 0; i < n; ++i) {
            const size_t source = keep_indices.empty() ? i : static_cast<size_t>(keep_indices[i]);
            centered[i] = std::isfinite(dosage[source]) ? dosage[source] - mean : 0.0;
        }
        const Eigen::MatrixXd contribution = centered * centered.transpose();
        global_raw.noalias() += contribution;
        const size_t ci = chrom_found->second;
        per_chrom_raw[ci].noalias() += contribution;
        ++per_chrom_count[ci];
        per_chrom_sum[ci] += sum_2pq;
        ++total_valid_snps_count;
        total_sum_2pq += sum_2pq;
    }
    if (total_valid_snps_count == 0 || total_sum_2pq <= 0.0) {
        throw std::runtime_error("No valid HapMap/text variants remained for GRM construction");
    }

    G = global_raw / total_sum_2pq;
    N_snps = Eigen::MatrixXf::Constant(n, n, static_cast<float>(total_valid_snps_count));
    SplitGrmResult result;
    result.global = G;
    result.chrom_labels = std::move(chrom_vec);
    result.per_chrom_raw = std::move(per_chrom_raw);
    result.per_chrom_snp_count = std::move(per_chrom_count);
    result.per_chrom_sum_2pq = std::move(per_chrom_sum);
    result.global_sum_2pq = total_sum_2pq;
    result.global_snp_count = total_valid_snps_count;
    return result;
}

void GrmBuilder::compute(SequentialGenotypeMatrix& reader,
                         const std::vector<int>& keep_indices) {
    (void)computeSplit(reader, {}, keep_indices);
}

} // namespace cosmic
