#include "mvlmm_engine.h"
#include <iostream>
#include <cmath>
#include <omp.h>

namespace cosmic {

static double gamma_q(double a, double x) {
    if (x <= 0.0) return 1.0;
    if (x < a + 1.0) {
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int i = 0; i < 100; ++i) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::abs(del) < std::abs(sum) * 1e-7) {
                return 1.0 - sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
            }
        }
    } else {
        double b = x + 1.0 - a;
        double c = 1.0 / 1e-30;
        double d = 1.0 / b;
        double h = d;
        for (int i = 1; i <= 100; ++i) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (std::abs(d) < 1e-30) d = 1e-30;
            c = b + an / c;
            if (std::abs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (std::abs(del - 1.0) < 1e-7) {
                return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
            }
        }
    }
    return 0.0;
}

static double chisq_prob(double stat, int df) {
    if (stat <= 0) return 1.0;
    return gamma_q(df / 2.0, stat / 2.0);
}

MvLmmEngine::MvLmmEngine(const Options& opts) : options(opts) {}

void MvLmmEngine::prepareFromVce(const MvEigenVce& vce) {
    Vg = vce.getVg();
    Ve = vce.getVe();
    D = vce.getEigenvalues();
    W_proj = vce.getProjectedW();
    Uty = vce.getUty();
    
    n_samples = D.size();
    n_traits = Vg.rows();
    
    // In covariate projection mode, UtX might be empty
    UtX = Eigen::MatrixXd::Zero(n_samples, 0); 
    n_covars = 0;
}

// Sparse Null Model Solver using PCG and MC-REML
void MvLmmEngine::processBlock(const Eigen::MatrixXd& S_block, 
                               const std::vector<double>& means, 
                               size_t global_start_idx, 
                               std::vector<MvGwasResult>& results,
                               const std::vector<struct PlinkBimInfo>& bim,
                               const Eigen::MatrixXd& H_XX_inv,
                               const Eigen::VectorXd& RHS_X,
                               const std::vector<Eigen::MatrixXd>& Sigma_inv) {
    size_t current_block_size = S_block.cols();
    
    // S_block is already rotated by W_proj_t in scan()
    const Eigen::MatrixXd& UtS_block = S_block;
    
    int pd = n_covars * n_traits;
    
    // Initialize Accumulators for all SNPs in block
    std::vector<Eigen::MatrixXd> H_XS(current_block_size, Eigen::MatrixXd::Zero(pd, n_traits));
    std::vector<Eigen::MatrixXd> H_SS(current_block_size, Eigen::MatrixXd::Zero(n_traits, n_traits));
    std::vector<Eigen::VectorXd> RHS_S(current_block_size, Eigen::VectorXd::Zero(n_traits));
    
    // Optimized Loop: Sample-Major
    // Pre-calculate X_Si for sample i if possible, but Si is d x d.
    // Inner loop over k (SNPs) should benefit from cache locality of Si, xi, yi.
    
    #pragma omp parallel
    {
        // Thread-local accumulators
        std::vector<Eigen::MatrixXd> H_XS_loc(current_block_size, Eigen::MatrixXd::Zero(pd, n_traits));
        std::vector<Eigen::MatrixXd> H_SS_loc(current_block_size, Eigen::MatrixXd::Zero(n_traits, n_traits));
        std::vector<Eigen::VectorXd> RHS_S_loc(current_block_size, Eigen::VectorXd::Zero(n_traits));

        #pragma omp for schedule(static)
        for (int i = 0; i < n_samples; ++i) {
            const Eigen::MatrixXd& Si = Sigma_inv[i];
            Eigen::VectorXd yi = Uty.row(i).transpose(); // d x 1
            
            // Pre-compute common terms for this sample
            Eigen::VectorXd Si_yi = Si * yi; // d x 1
            
            for (size_t k = 0; k < current_block_size; ++k) {
                double si = UtS_block(i, k);
                
                // H_SS += si^2 * Si
                H_SS_loc[k] += (si * si) * Si;
                
                // RHS_S += si * Si * yi
                RHS_S_loc[k] += si * Si_yi;
                
                // H_XS
                // Block j: xi(j) * si * Si
                if (n_covars > 0) {
                    Eigen::RowVectorXd xi = UtX.row(i); // 1 x p
                    for (int j = 0; j < n_covars; ++j) {
                        H_XS_loc[k].block(j*n_traits, 0, n_traits, n_traits) += (xi(j) * si) * Si;
                    }
                }
            }
        }
        
        // Reduction
        #pragma omp critical
        {
            for (size_t k = 0; k < current_block_size; ++k) {
                if (n_covars > 0) {
                    H_XS[k] += H_XS_loc[k];
                }
                H_SS[k] += H_SS_loc[k];
                RHS_S[k] += RHS_S_loc[k];
            }
        }
    }
    
    // Finalize per SNP (independent)
    #pragma omp parallel for
    for (long long k = 0; k < (long long)current_block_size; ++k) {
        size_t global_idx = global_start_idx + k;
        
        // Schur Complement
        Eigen::MatrixXd Schur;
        Eigen::VectorXd adjusted_RHS;
        if (n_covars > 0) {
            Eigen::MatrixXd H_XS_T = H_XS[k].transpose();
            Schur = H_SS[k] - H_XS_T * H_XX_inv * H_XS[k];
            adjusted_RHS = RHS_S[k] - H_XS_T * H_XX_inv * RHS_X;
        } else {
            Schur = H_SS[k];
            adjusted_RHS = RHS_S[k];
        }
        
        Eigen::MatrixXd Var_beta;
        Eigen::LLT<Eigen::MatrixXd> lltSchur(Schur);
        if(lltSchur.info() == Eigen::Success) {
            Var_beta = lltSchur.solve(Eigen::MatrixXd::Identity(Schur.rows(), Schur.cols()));
        } else {
            Var_beta = Schur.inverse();
        }
        Eigen::VectorXd beta = Var_beta * adjusted_RHS;
        
        double wald_stat = beta.transpose() * Schur * beta;
        double p_val = chisq_prob(wald_stat, n_traits);
        
        MvGwasResult res;
        if (!bim.empty()) {
            // bim is the block_bim passed from caller, so use local index k
            res.chrom = bim[k].chrom;
            res.snp_id = bim[k].id;
            res.pos = bim[k].bp_pos;
            res.a1 = bim[k].alt;
            res.a2 = bim[k].ref;
        } 
        // For PGEN, we might fill these outside or pass them in differently.
        // Currently scan(PgenReader) handles metadata separately.
        // But to unify, we can pass metadata vector.
        
        res.freq = means[k] / 2.0;
        res.beta = beta;
        res.se = Var_beta; 
        res.wald_stat = wald_stat;
        res.p_value = p_val;
        
        results[global_idx] = res;
    }
}


// --- Sparse GRM Support ---

void MvLmmEngine::scan(PgenReader& reader, 
                       std::vector<MvGwasResult>& results, 
                       const std::vector<int>& keep_indices,
                       const std::string& chrom_filter) {
    throw std::runtime_error("MvLmmEngine::scan(PgenReader) not implemented");
}



/*
void MvLmmEngine::scan_internal(...) {
    results.resize(m);
    
    if (use_sparse_grm) {
        // ...
    }
    // ...
}
*/

    // Orphaned code removed


void MvLmmEngine::scan(PlinkReader& reader, 
                       std::vector<MvGwasResult>& results, 
                       const std::vector<int>& keep_indices,
                       const std::string& chrom_filter) {
    size_t m = reader.getNumSnps();
    size_t n = keep_indices.size();
    
    int n_original = W_proj.rows();
    if (n != (size_t)n_original) {
        std::cerr << "MvLmmEngine::scan: n = " << n << ", n_original = " << n_original << std::endl;
        throw std::runtime_error("MvLmmEngine::scan: Sample size mismatch");
    }
    
    results.resize(m);
    const auto& bim = reader.getBimInfo();
    
    std::vector<Eigen::MatrixXd> Sigma_inv_local;

    // Pre-compute Sigma_inv
    Sigma_inv_local.resize(n_samples);
    #pragma omp parallel for
    for (int i = 0; i < n_samples; ++i) {
        Eigen::MatrixXd Sigma = D(i) * Vg + Ve;
        Sigma_inv_local[i] = Sigma.llt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
    }

    // Pre-compute H_XX and H_XX_inv
    int pd = n_covars * n_traits;
    Eigen::MatrixXd H_XX = Eigen::MatrixXd::Zero(pd, pd);
    Eigen::VectorXd RHS_X = Eigen::VectorXd::Zero(pd);
    
    for (int i = 0; i < n_samples; ++i) {
        Eigen::MatrixXd Si = Sigma_inv_local[i];
        Eigen::RowVectorXd xi = UtX.row(i);
        Eigen::VectorXd yi = Uty.row(i).transpose();
        
        for (int j = 0; j < n_covars; ++j) {
            for (int k = 0; k < n_covars; ++k) {
                H_XX.block(j*n_traits, k*n_traits, n_traits, n_traits) += xi(j) * xi(k) * Si;
            }
            RHS_X.segment(j*n_traits, n_traits) += xi(j) * (Si * yi);
        }
    }
    
    Eigen::MatrixXd H_XX_inv;
    if (pd > 0) {
        H_XX_inv = H_XX.llt().solve(Eigen::MatrixXd::Identity(pd, pd));
    } else {
        H_XX_inv = Eigen::MatrixXd::Zero(0, 0);
    }
    
    // Scan Loop
    bool is_projected = (W_proj.cols() < W_proj.rows());
    std::cout << "Starting Multivariate Scan (PLINK" << (is_projected ? ", Projected" : "") << ") on " << m << " SNPs..." << std::endl;
    
    size_t processed = 0;
    size_t block_size = 2048; // Optimized block size
    
    // U transpose for rotation
    Eigen::MatrixXd W_proj_t = W_proj.transpose(); // W_proj is n_original x n_samples

    while (processed < m) {
        size_t current_block_size = std::min(block_size, m - processed);
        size_t start_idx = processed;
        size_t end_idx = processed + current_block_size;
        
        // 1. Read Block (Raw Doubles)
        Eigen::MatrixXd S_block_raw;
        reader.getSnpBlock(start_idx, end_idx, S_block_raw);
        
        // 2. Impute and Center
        std::vector<double> means(current_block_size);
        std::vector<struct PlinkBimInfo> block_bim(current_block_size);
        
        #pragma omp parallel for
        for (long long k = 0; k < (long long)current_block_size; ++k) {
            size_t snp_idx = start_idx + k;
            block_bim[k] = bim[snp_idx];
            
            double sum = 0;
            int count = 0;
            for(int i=0; i<n_original; ++i) {
                double val = S_block_raw(i, k);
                if (!std::isnan(val)) {
                    sum += val;
                    count++;
                }
            }
            
            double mean = (count > 0) ? sum / count : 0.0;
            means[k] = mean;
            
            for(int i=0; i<n_original; ++i) {
                if (std::isnan(S_block_raw(i, k))) S_block_raw(i, k) = 0.0;
                else S_block_raw(i, k) -= mean;
            }
        }
        
        // 3. Rotate Block: S_rot = W_proj' * S_block
        // This is the heavy lifting: Matrix Multiplication
        Eigen::MatrixXd S_rot = W_proj_t * S_block_raw;
        
        // 4. Process Block (using rotated genotypes)
        // Apply Chromosome Filter here if needed
        bool block_relevant = false;
        if (chrom_filter.empty()) {
            block_relevant = true;
        } else {
            for(const auto& info : block_bim) {
                if (info.chrom == chrom_filter) {
                    block_relevant = true;
                    break;
                }
            }
        }
        
        if (block_relevant) {
            processBlock(S_rot, means, processed, results, block_bim, H_XX_inv, RHS_X, Sigma_inv_local);
        }
        
        processed += current_block_size;
        std::cout << "\rProcessed " << processed << "/" << m << " SNPs..." << std::flush;
    }
    std::cout << std::endl;
}

} // namespace cosmic
