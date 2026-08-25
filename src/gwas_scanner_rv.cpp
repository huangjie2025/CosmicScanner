#ifdef COSMICSCANNER_ENABLE_EXPERIMENTAL

#include "gwas_scanner.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <unordered_set>

namespace cosmic {

void GwasScanner::runRareVariantScan(
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
    std::cout << "  Rare Variant Association Test Module   \n";
    std::cout << "=========================================\n";

    if (options.set_file.empty()) {
        throw std::runtime_error("Rare variant test requires a --set-file to define SNP sets (e.g., genes).");
    }

    // 1. Load Set File
    std::vector<VariantSet> sets = SetFileReader::read(options.set_file);
    if (sets.empty()) {
        throw std::runtime_error("No valid sets found in " + options.set_file);
    }

    std::unordered_set<std::string> target_snps;
    for (const auto& set : sets) {
        target_snps.insert(set.snp_ids.begin(), set.snp_ids.end());
    }

    // 2. Prepare Null Model (LMM)
    std::cout << "\n[Step 1] Fitting Null LMM Model...\n";
    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options.min_maf;
    grm_opts.min_hwe = (float)options.min_hwe;
    grm_opts.algorithm = options.grm_alg;
    #ifdef _OPENMP
    grm_opts.thread_num = options.n_threads;
    #endif
    
    GrmBuilder global_builder(grm_opts);
    if (is_bgen) { bgen_reader.reset(); global_builder.compute(bgen_reader); }
    else if (is_pgen) { pgen_reader.reset(); global_builder.compute(pgen_reader); }
    else { global_builder.compute(plink_reader); }
    
    const Eigen::MatrixXd& G_global_norm = global_builder.getGrm();
    
    Eigen::VectorXd y_sub = DataFilter::subset(y, keep);
    Eigen::MatrixXf X_sub = DataFilter::subsetRows(X, keep);
    Eigen::MatrixXd G_sub = DataFilter::subsetSquare(G_global_norm, keep);

    EigenVce::Options vce_opts;
    vce_opts.use_covariate_projection = options.use_covariate_projection;
    EigenVce vce(vce_opts);
    vce.prepare(y_sub, X_sub, G_sub);
    vce.runNullModel();

    std::cout << "  Null Model VCs: V(G)=" << vce.getSigma2g() << ", V(e)=" << vce.getSigma2e() << "\n";

    // 3. Extract Null Model Projection Matrix P and Scores
    // P = V^{-1} - V^{-1}X(X'V^{-1}X)^{-1}X'V^{-1}
    int n = y_sub.size();
    Eigen::MatrixXd V = vce.getSigma2g() * G_sub + vce.getSigma2e() * Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd V_inv = V.inverse();
    Eigen::MatrixXd X_d = X_sub.cast<double>();
    Eigen::MatrixXd XtVinvX = X_d.transpose() * V_inv * X_d;
    Eigen::MatrixXd XtVinvX_inv = XtVinvX.inverse();
    Eigen::MatrixXd P_matrix = V_inv - V_inv * X_d * XtVinvX_inv * X_d.transpose() * V_inv;
    
    Eigen::VectorXd P_y = P_matrix * y_sub;

    // 4. Setup Rare Variant Test
    RareVariantTest::Method rv_method = options.run_skat ? RareVariantTest::Method::SKAT : RareVariantTest::Method::BURDEN;
    RareVariantTest rv_test(rv_method);
    
    std::string method_name = options.run_skat ? "SKAT" : "Burden Test";
    std::cout << "\n[Step 2] Running " << method_name << " on " << sets.size() << " sets...\n";
    
    std::string out_file = options.out_prefix + ".rare.txt";
    writeHeader(out_file);
    
    // Create map from SNP ID to Genotype vector for fast lookup
    std::cout << "  Caching genotype data for target SNPs...\n";
    std::unordered_map<std::string, Eigen::VectorXd> snp_genotypes;
    std::unordered_map<std::string, double> snp_mafs;
    
    if (is_bgen) bgen_reader.reset();
    else if (is_pgen) pgen_reader.reset();
    
    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage_f;
    
    size_t m_snps = plink_reader.getNumSnps();
    size_t current_snp = 0;
    const auto& bim = plink_reader.getBimInfo();
    
    while(true) {
        bool ok = false;
        if (is_bgen) {
            ok = bgen_reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f);
        } else if (is_pgen) {
            ok = pgen_reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f);
        } else {
            if (current_snp >= m_snps) break;
            ok = true;
            chrom = bim[current_snp].chrom;
            rsid = bim[current_snp].id;
            pos = bim[current_snp].bp_pos;
            ref = bim[current_snp].ref;
            alt = bim[current_snp].alt;
            
            Eigen::VectorXd d_d(plink_reader.getNumSamples());
            plink_reader.readSnp(current_snp, d_d);
            dosage_f = d_d.cast<float>();
            current_snp++;
        }
        
        if (!ok) break;
        
        if (target_snps.find(rsid) == target_snps.end()) continue;
        
        Eigen::VectorXd g(n);
        double sum_g = 0.0;
        int valid_n = 0;
        for(int i=0; i<n; ++i) {
            float val = dosage_f[keep.keep_indices[i]];
            g(i) = val;
            if (!std::isnan(val)) {
                sum_g += val;
                valid_n++;
            }
        }
        
        if (valid_n == 0) continue;
        double mean_g = sum_g / valid_n;
        double freq = mean_g / 2.0;
        double maf = std::min(freq, 1.0 - freq);
        
        // Impute missing
        for(int i=0; i<n; ++i) {
            if (std::isnan(g(i))) g(i) = mean_g;
        }
        
        snp_genotypes[rsid] = g;
        snp_mafs[rsid] = maf;
    }
    
    std::cout << "  Found " << snp_genotypes.size() << " valid target SNPs in genotype data.\n";
    
    // 5. Run tests for each set
    std::vector<GwasResult> all_results;
    int processed = 0;
    
    for (const auto& set : sets) {
        std::vector<Eigen::VectorXd> valid_gs;
        std::vector<double> valid_mafs;
        
        for (const auto& snp : set.snp_ids) {
            if (snp_genotypes.count(snp)) {
                valid_gs.push_back(snp_genotypes[snp]);
                valid_mafs.push_back(snp_mafs[snp]);
            }
        }
        
        int k = valid_gs.size();
        if (k == 0) {
            // No valid SNPs found for this set
            GwasResult res;
            res.snp_id = set.set_id;
            res.n_obs = 0;
            res.beta = 0.0;
            res.p_value = 1.0;
            all_results.push_back(res);
            continue;
        }
        
        // Build X_set (N x K)
        Eigen::MatrixXd X_set(n, k);
        Eigen::VectorXd MAF(k);
        for (int j = 0; j < k; ++j) {
            X_set.col(j) = valid_gs[j];
            MAF(j) = valid_mafs[j];
        }
        
        // Marginal Scores U = X_set' * P_y
        Eigen::VectorXd score_vec = X_set.transpose() * P_y;
        
        // Run SKAT or Burden
        RareVariantTest::Result rv_res = rv_test.runTest(set.set_id, X_set, score_vec, P_matrix, MAF);
        
        GwasResult out_res;
        out_res.snp_id = rv_res.set_id; // Repurpose snp_id for set_id
        out_res.n_obs = rv_res.n_snps;  // Repurpose n_obs for valid SNP count
        out_res.beta = rv_res.stat;     // Repurpose beta for test statistic
        out_res.p_value = rv_res.p_value;
        all_results.push_back(out_res);
        
        processed++;
        if (processed % 100 == 0) {
            std::cout << "\r  Processed " << processed << " / " << sets.size() << " sets..." << std::flush;
        }
    }
    
    std::cout << "\n[Step 3] Saving Results...\n";
    appendResults(out_file, all_results);
    std::cout << "  Rare Variant Results saved to: " << out_file << std::endl;
}

} // namespace cosmic

#endif // COSMICSCANNER_ENABLE_EXPERIMENTAL