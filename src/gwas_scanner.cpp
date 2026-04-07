#include "gwas_scanner.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <omp.h>
#include <stdexcept>
#include <unordered_map>
#include "design.h"

namespace cosmic {

GwasScanner::GwasScanner(const Options& opts) : options(opts) {}

void GwasScanner::writeHeader(const std::string& filename) {
    std::ofstream ofs(filename);
    if (options.run_skat || options.run_burden) {
        ofs << "SET_ID\tN_SNPS\tSTAT\tP\n";
    } else if (options.out_ssf) {
        ofs << "chromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tsample_size\tbeta\tstandard_error\tp_value\tvariant_id\n";
    } else {
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tLOG10P\n";
    }
    ofs.close();
}

void GwasScanner::appendResults(const std::string& filename, const std::vector<GwasResult>& results) {
    std::ofstream ofs(filename, std::ios::app);
    for(const auto& res : results) {
        if (std::isnan(res.p_value)) continue;
        if (options.run_skat || options.run_burden) {
            ofs << res.snp_id << "\t" << res.n_obs << "\t" << res.beta << "\t" << res.p_value << "\n"; // using snp_id as set_id, beta as stat
        } else if (options.out_ssf) {
            ofs << res.chrom << "\t" << res.pos << "\t"
                << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
                << res.p_value << "\t" << res.snp_id << "\n";
        } else {
            ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
                << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
                << res.p_value << "\t" << res.log10_p << "\n";
        }
    }
}

void GwasScanner::runLocoScan(
    const Eigen::VectorXd& y,
    const Eigen::MatrixXf& X,
    PlinkReader& plink_reader,
    BgenReader& bgen_reader,
    PgenReader& pgen_reader,
    bool is_bgen,
    bool is_pgen,
    const DataFilter::KeepInfo& keep
) {
    std::cout << "Calculating Global GRM (Pass 1)...\n";
    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options.min_maf;
    grm_opts.min_hwe = (float)options.min_hwe;
    grm_opts.algorithm = options.grm_alg;
    #ifdef _OPENMP
    grm_opts.thread_num = options.n_threads;
    #endif
    
    GrmBuilder global_builder(grm_opts);
    if (is_bgen) {
            bgen_reader.reset();
            global_builder.compute(bgen_reader);
    } else if (is_pgen) {
            pgen_reader.reset();
            global_builder.compute(pgen_reader);
    } else {
            global_builder.compute(plink_reader);
    }
    
    const Eigen::MatrixXd& G_global_norm = global_builder.getGrm();
    double S_global = 0.0;
    
    if (options.grm_alg == 0) {
        S_global = (double)global_builder.getTotalSnps();
    } else {
        S_global = global_builder.getSum2pq();
    }
    
    std::cout << "Global GRM computed. Scaling Factor = " << S_global << std::endl;
    auto chroms_set = global_builder.getObservedChroms();
    
    // Prepare Output File
    std::string out_file = options.out_prefix + ".assoc.txt";
    writeHeader(out_file);

    // Subset y and X for analysis
    Eigen::VectorXd y_sub = DataFilter::subset(y, keep);
    Eigen::MatrixXf X_sub = DataFilter::subsetRows(X, keep);
    
    // Apply Observation Weights (e.g., from Repeatability Model) if provided
    bool has_weights = (options.obs_weights.size() == keep.n_keep);
    if (has_weights) {
        y_sub.array() *= options.obs_weights.array();
        for(int j=0; j<X_sub.cols(); ++j) {
            X_sub.col(j).array() *= options.obs_weights.cast<float>().array();
        }
    }
    
    // LMM Options
    LmmEngine::Options lmm_opts;
    lmm_opts.fast_mode = options.fast_lmm;
    lmm_opts.min_maf = options.min_maf;
    lmm_opts.min_hwe = options.min_hwe;
    lmm_opts.use_covariate_projection = options.use_covariate_projection;
    lmm_opts.use_spa = options.use_spa;
    if (has_weights) lmm_opts.obs_weights = options.obs_weights;
    // Enable GCTA-style residuals if using Yang et al. algorithm (GCTA mode)
    if (options.grm_alg == 0) {
        lmm_opts.use_residuals = true;
    }
    #ifdef _OPENMP
    lmm_opts.thread_num = options.n_threads;
    #endif

    if (chroms_set.size() <= 1 || options.no_loco) {
        std::cout << (options.no_loco ? "LOCO disabled via flag. " : "Only 1 chromosome found. ");
        std::cout << "Using standard LMM (no LOCO)...\n";
        
        // Subset Global GRM
        Eigen::MatrixXd G_sub = DataFilter::subsetSquare(G_global_norm, keep);
        
        if (has_weights) {
            G_sub.array().colwise() *= options.obs_weights.array();
            G_sub.array().rowwise() *= options.obs_weights.transpose().array();
        }
        
        EigenVce::Options vce_opts;
        vce_opts.use_covariate_projection = options.use_covariate_projection;
        EigenVce vce(vce_opts);
        vce.prepare(y_sub, X_sub, G_sub);
        vce.runNullModel();
        
        std::cout << "  Null Model VCs: V(G)=" << vce.getSigma2g() << ", V(e)=" << vce.getSigma2e() << std::endl;
        
        LmmEngine engine(lmm_opts);
        engine.prepareFromVce(vce);
        
        std::vector<GwasResult> results;
        if (is_bgen) {
            bgen_reader.reset();
            engine.scan(bgen_reader, results, keep.keep_indices);
        } else if (is_pgen) {
            pgen_reader.reset();
            engine.scan(pgen_reader, results, keep.keep_indices);
        } else {
            engine.scan(plink_reader, results, keep.keep_indices);
        }
        
        appendResults(out_file, results);
    } else {
        std::cout << "Running LOCO Scan on " << chroms_set.size() << " chromosomes...\n";
        
        // Sort chroms
        std::vector<std::string> sorted_chroms(chroms_set.begin(), chroms_set.end());
        std::sort(sorted_chroms.begin(), sorted_chroms.end(), [](const std::string& a, const std::string& b) {
            try {
                int ia = std::stoi(a);
                int ib = std::stoi(b);
                return ia < ib;
            } catch (...) {
                return a < b;
            }
        });

        for(const std::string& chrom : sorted_chroms) {
            std::cout << "Processing Chromosome " << chrom << "...\n";
            
            // Calculate Chr GRM Contribution
            GrmBuilder::Options chr_opts = grm_opts;
            chr_opts.chrom_filter = chrom;
            GrmBuilder chr_builder(chr_opts);
            
            if (is_bgen) { bgen_reader.reset(); chr_builder.compute(bgen_reader); }
            else if (is_pgen) { pgen_reader.reset(); chr_builder.compute(pgen_reader); }
            else { chr_builder.compute(plink_reader); }
            
            double S_chr = 0.0;
            if (options.grm_alg == 0) {
                 S_chr = (double)chr_builder.getTotalSnps();
            } else {
                 S_chr = chr_builder.getSum2pq();
            }
            
            const auto& G_chr_norm = chr_builder.getGrm();

            // Calculate LOCO GRM Subset
            int n_sub = keep.n_keep;
            Eigen::MatrixXd G_loco_sub(n_sub, n_sub);
            
            double S_loco = S_global - S_chr;
            
            if (S_chr > 0 && S_loco > 1e-6) {
                double w1 = (S_global / S_loco);
                double w2 = (S_chr / S_loco);
                
                #pragma omp parallel for collapse(2)
                for(int i=0; i<n_sub; ++i) {
                    for(int j=0; j<n_sub; ++j) {
                        int gi = keep.keep_indices[i];
                        int gj = keep.keep_indices[j];
                        G_loco_sub(i, j) = G_global_norm(gi, gj) * w1 - G_chr_norm(gi, gj) * w2;
                    }
                }
            } else {
                    G_loco_sub = DataFilter::subsetSquare(G_global_norm, keep);
            }
            
            if (has_weights) {
                G_loco_sub.array().colwise() *= options.obs_weights.array();
                G_loco_sub.array().rowwise() *= options.obs_weights.transpose().array();
            }
            
            EigenVce::Options vce_opts;
            vce_opts.use_covariate_projection = options.use_covariate_projection;
            EigenVce vce(vce_opts);
            vce.prepare(y_sub, X_sub, G_loco_sub);
            vce.runNullModel();
            
            std::cout << "  Null Model VCs: V(G)=" << vce.getSigma2g() << ", V(e)=" << vce.getSigma2e() << std::endl;
            
            // Run LMM for this chromosome
            LmmEngine engine(lmm_opts);
            engine.prepareFromVce(vce);
            
            // Scan
            std::vector<GwasResult> chrom_results;
            if (is_bgen) { bgen_reader.reset(); engine.scan(bgen_reader, chrom_results, keep.keep_indices, chrom); }
            else if (is_pgen) { pgen_reader.reset(); engine.scan(pgen_reader, chrom_results, keep.keep_indices, chrom); }
            else { engine.scan(plink_reader, chrom_results, keep.keep_indices, chrom); }
            
            appendResults(out_file, chrom_results);
        }
    }
    std::cout << "Results saved to " << out_file << std::endl;
}

void GwasScanner::runBiobankScan(
    const Eigen::VectorXd& y,
    const Eigen::MatrixXf& X,
    PlinkReader& plink_reader,
    BgenReader& bgen_reader,
    PgenReader& pgen_reader,
    bool is_bgen,
    bool is_pgen,
    const DataFilter::KeepInfo& keep
) {
    std::cout << "\n=========================================\n";
    std::cout << "  Biobank Mode (REGENIE-style Ridge)     \n";
    std::cout << "=========================================\n";
    
    int n = y.size();
    int m = is_bgen ? bgen_reader.getNumVariants() : (is_pgen ? pgen_reader.getNumVariants() : plink_reader.getNumSnps());
    
    // Step 1: Run Block-wise Ridge Regression
    RegenieStep1::Options step1_opts;
    step1_opts.block_size = options.regenie_block_size;
    step1_opts.use_loco = !options.no_loco;
    step1_opts.n_threads = options.n_threads;
    
    RegenieStep1 regenie_step1(step1_opts);
    
    std::cout << "  Loading Genotypes in blocks for Step 1..." << std::endl;
    
    if (is_bgen) bgen_reader.reset();
    else if (is_pgen) pgen_reader.reset();
    else plink_reader.reset();
    
    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage_f;
    
    int valid_snps = 0;
    int current_block_snps = 0;
    Eigen::MatrixXd X_block(n, step1_opts.block_size);
    int current_block_chrom = 0;
    
    while(true) {
        bool ok = false;
        if (is_bgen) {
            ok = bgen_reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f);
        } else if (is_pgen) {
            ok = pgen_reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f);
        } else {
            if (valid_snps >= m) break;
            ok = true;
            const auto& bim = plink_reader.getBimInfo();
            chrom = bim[valid_snps].chrom;
            Eigen::VectorXd d_d(plink_reader.getNumSamples());
            plink_reader.readSnp(valid_snps, d_d);
            dosage_f = d_d.cast<float>();
        }
        
        if (!ok || valid_snps >= m) break;
        
        int chr_int = 0;
        try { chr_int = std::stoi(chrom); } catch (...) { chr_int = 0; }
        
        if (current_block_snps == 0) {
            current_block_chrom = chr_int;
        }
        
        // Fill NaNs with mean and insert into block
        float sum = 0.0f; int count = 0;
        for (int i = 0; i < n; ++i) {
            float val = dosage_f[keep.keep_indices[i]];
            if (!std::isnan(val)) { sum += val; count++; }
        }
        float mean = count > 0 ? sum / count : 0.0f;
        for (int i = 0; i < n; ++i) {
            float val = dosage_f[keep.keep_indices[i]];
            X_block(i, current_block_snps) = std::isnan(val) ? mean : val;
        }
        
        current_block_snps++;
        valid_snps++;
        
        // Add block if full or chromosome changes
        if (current_block_snps == step1_opts.block_size || (chr_int != current_block_chrom && current_block_snps > 0)) {
            regenie_step1.add_block(X_block.leftCols(current_block_snps), y, current_block_chrom);
            current_block_snps = 0;
        }
    }
    
    if (current_block_snps > 0) {
        regenie_step1.add_block(X_block.leftCols(current_block_snps), y, current_block_chrom);
    }
    
    // Run Level 1 stacking
    regenie_step1.fit_level1_stacking(y);
    
    std::cout << "\n[Step 2] Running Single-SNP Association Tests...\n";
    std::string out_file = options.out_prefix + ".biobank.assoc.txt";
    writeHeader(out_file);
    
    LmmEngine::Options lmm_opts;
    lmm_opts.fast_mode = true; // Always fast mode in biobank
    lmm_opts.min_maf = options.min_maf;
    lmm_opts.min_hwe = options.min_hwe;
    lmm_opts.use_covariate_projection = options.use_covariate_projection;
    lmm_opts.use_spa = options.use_spa;
    lmm_opts.thread_num = options.n_threads;
    
    LmmEngine engine(lmm_opts);
    
    if (options.no_loco) {
        Eigen::VectorXd loco_pred = regenie_step1.get_loco_prediction(0); // 0 means all chromosomes included
        Eigen::MatrixXf X_new(n, X.cols() + 1);
        X_new.leftCols(X.cols()) = X;
        X_new.col(X.cols()) = loco_pred.cast<float>();
        
        Eigen::MatrixXd G_dummy = Eigen::MatrixXd::Identity(n, n);
        
        EigenVce::Options ev_opts;
        ev_opts.use_covariate_projection = options.use_covariate_projection;
        EigenVce vce(ev_opts);
        vce.prepare(y, X_new, G_dummy);
        vce.runNullModel();
        engine.prepareFromVce(vce);
        
        std::vector<GwasResult> results;
        if (is_bgen) engine.scan(bgen_reader, results, keep.keep_indices);
        else if (is_pgen) engine.scan(pgen_reader, results, keep.keep_indices);
        else engine.scan(plink_reader, results, keep.keep_indices);
        
        appendResults(out_file, results);
    } else {
        // LOCO scanning
        for (int chr = 1; chr <= 22; ++chr) {
            Eigen::VectorXd loco_pred = regenie_step1.get_loco_prediction(chr);
            if (loco_pred.norm() < 1e-8) continue; // Skip empty chromosomes
            
            std::cout << "  Scanning Chromosome " << chr << " using LOCO predictions...\n";
            Eigen::MatrixXf X_new(n, X.cols() + 1);
            X_new.leftCols(X.cols()) = X;
            X_new.col(X.cols()) = loco_pred.cast<float>();
            
            Eigen::MatrixXd G_dummy = Eigen::MatrixXd::Identity(n, n);
            
            EigenVce::Options ev_opts;
            ev_opts.use_covariate_projection = options.use_covariate_projection;
            EigenVce vce(ev_opts);
            vce.prepare(y, X_new, G_dummy);
            vce.runNullModel();
            engine.prepareFromVce(vce);
            
            std::vector<GwasResult> results;
            if (is_bgen) engine.scan(bgen_reader, results, keep.keep_indices);
            else if (is_pgen) engine.scan(pgen_reader, results, keep.keep_indices);
            else engine.scan(plink_reader, results, keep.keep_indices);
            
            // Filter results for current chromosome
            std::vector<GwasResult> chr_results;
            for (const auto& res : results) {
                if (res.chrom == std::to_string(chr) || res.chrom == "0" + std::to_string(chr)) {
                    chr_results.push_back(res);
                }
            }
            appendResults(out_file, chr_results);
        }
    }
    
    std::cout << "  Biobank GWAS Scan complete. Results saved to: " << out_file << std::endl;
}

void GwasScanner::runRrmScan(
    const std::vector<PhenoReader::PhenoRecord>& records,
    int rrm_order,
    const Eigen::MatrixXf& X_base,
    PlinkReader& plink_reader,
    BgenReader& bgen_reader,
    PgenReader& pgen_reader,
    bool is_bgen,
    bool is_pgen,
    const DataFilter::KeepInfo& keep
) {
    int n_ind = keep.keep_indices.size();
    int n_obs = records.size();
    
    std::cout << "\n=== Starting Random Regression Model (RRM) GWAS Scan ===\n";
    std::cout << "Observations: " << n_obs << ", Individuals: " << n_ind << ", Legendre Order: " << rrm_order << "\n";
    
    // 1. Build Time Vector and ID Map
    Eigen::VectorXd time_vec(n_obs);
    Eigen::VectorXd y_vec(n_obs);
    std::vector<int> id_map(n_obs, -1);
    
    // Create a map from FID:IID to kept index
    std::unordered_map<std::string, int> ind_idx_map;
    const auto* fam_ptr = (!is_bgen && !is_pgen) ? &plink_reader.getFamInfo() : nullptr;
    
    for (int i = 0; i < n_ind; ++i) {
        int orig_idx = keep.keep_indices[i];
        std::string fid, iid;
        if (is_pgen) {
            fid = pgen_reader.getSampleIds()[orig_idx];
            iid = fid;
        } else if (is_bgen) {
            fid = bgen_reader.getSampleIds()[orig_idx];
            iid = fid;
        } else if (fam_ptr) {
            fid = (*fam_ptr)[orig_idx].fid;
            iid = (*fam_ptr)[orig_idx].iid;
        }
        ind_idx_map[fid + ":" + iid] = i;
    }
    
    int valid_obs = 0;
    std::vector<int> valid_obs_indices;
    double tmin = std::numeric_limits<double>::max();
    double tmax = std::numeric_limits<double>::lowest();
    
    for (int i = 0; i < n_obs; ++i) {
        std::string key = records[i].fid + ":" + records[i].iid;
        auto it = ind_idx_map.find(key);
        if (it != ind_idx_map.end() && !std::isnan(records[i].value)) {
            time_vec(valid_obs) = records[i].time;
            y_vec(valid_obs) = records[i].value;
            id_map[valid_obs] = it->second;
            
            if (records[i].time < tmin) tmin = records[i].time;
            if (records[i].time > tmax) tmax = records[i].time;
            
            valid_obs_indices.push_back(i);
            valid_obs++;
        }
    }
    
    time_vec.conservativeResize(valid_obs);
    y_vec.conservativeResize(valid_obs);
    id_map.resize(valid_obs);
    
    std::cout << "Matched " << valid_obs << " valid observations to " << n_ind << " genotyped individuals.\n";
    std::cout << "Time range: [" << tmin << ", " << tmax << "]\n";
    
    if (valid_obs == 0) {
        std::cerr << "Error: No valid observations after matching.\n";
        return;
    }
    
    // 2. Build Legendre Polynomials Matrix P
    Eigen::MatrixXd P = buildLegendreMatrix(time_vec, tmin, tmax, rrm_order);
    
    // Note: The rest of the RRM GWAS logic (VCE, MME construction, Score Test) 
    // is highly complex and requires dense matrix operations on large dimensions.
    // For now, we output the data processing success.
    std::cout << "Successfully built Legendre design matrix of size " << P.rows() << " x " << P.cols() << ".\n";
    
    // Future work:
    // 1. Compute Global GRM (K)
    // 2. Setup VCESolver with multiple RandomComponents (one for each Legendre order)
    // 3. Solve for Variance Components
    // 4. Construct MME LHS and invert it to get C_aa
    // 5. Run Score Test for each SNP using C_aa and BLUPs (a_hat)
    
    std::cout << "RRM exact LMM scan is under development. Data processing complete.\n";
}

} // namespace cosmic
