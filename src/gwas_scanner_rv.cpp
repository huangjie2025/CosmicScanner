#include "gwas_scanner.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <unordered_set>
#include <limits>

namespace cosmic {

void GwasScanner::runRareVariantScan(
    const Eigen::VectorXd& y,
    const Eigen::MatrixXf& X,
    PlinkReader& plink_reader,
    BgenReader& bgen_reader,
    PgenReader& pgen_reader,
    VcfReader& vcf_reader,
    bool is_bgen,
    bool is_pgen,
    bool is_vcf,
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
    else if (is_vcf) { vcf_reader.reset(); global_builder.compute(vcf_reader); }
    else { global_builder.compute(plink_reader); }
    
    const Eigen::MatrixXd& G_global_norm = global_builder.getGrm();
    
    Eigen::VectorXd y_sub = DataFilter::subset(y, keep);
    Eigen::MatrixXf X_sub = DataFilter::subsetRows(X, keep);
    Eigen::MatrixXd G_sub = DataFilter::subsetSquare(G_global_norm, keep);

    // Reject invalid/singular designs before VCE and the legacy inverse-based
    // projection. Valid-input arithmetic and all set-test statistics stay unchanged.
    if (!y_sub.allFinite() || !X_sub.allFinite() || !G_sub.allFinite() ||
        y_sub.size() <= X_sub.cols())
        throw std::runtime_error("Rare-variant null requires finite inputs and residual degrees of freedom");
    const Eigen::MatrixXd X_d = X_sub.cast<double>();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> design_rank(X_d);
    // Accumulation error grows with the observation count, including Eigen-only
    // builds. The default min(rows, cols) threshold can accept duplicate columns.
    design_rank.setThreshold(std::numeric_limits<double>::epsilon() *
                             static_cast<double>(std::max(X_d.rows(), X_d.cols())));
    if (design_rank.rank() != X_d.cols())
        throw std::runtime_error("Rare-variant covariate design is rank deficient");

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
    if (!V.allFinite() || Eigen::LLT<Eigen::MatrixXd>(V).info() != Eigen::Success)
        throw std::runtime_error("Rare-variant null covariance is not finite positive definite");
    Eigen::MatrixXd V_inv = V.inverse();
    Eigen::MatrixXd XtVinvX = X_d.transpose() * V_inv * X_d;
    if (!XtVinvX.allFinite() || (X_d.cols() > 0 &&
        Eigen::LLT<Eigen::MatrixXd>(XtVinvX).info() != Eigen::Success))
        throw std::runtime_error("Rare-variant covariate information is not positive definite");
    Eigen::MatrixXd XtVinvX_inv = XtVinvX.inverse();
    Eigen::MatrixXd P_matrix = V_inv - V_inv * X_d * XtVinvX_inv * X_d.transpose() * V_inv;
    
    Eigen::VectorXd P_y = P_matrix * y_sub;
    if (!P_matrix.allFinite() || !P_y.allFinite())
        throw std::runtime_error("Rare-variant null projection is nonfinite");

    // 4. Setup Rare Variant Test
    RareVariantTest::Method rv_method = options.run_skat ? RareVariantTest::Method::SKAT : RareVariantTest::Method::BURDEN;
    const bool extension = !options.rare_method.empty();
    if (options.rare_method == "skat-o") rv_method = RareVariantTest::Method::SKAT_O;
    if (options.rare_method == "acat-v") rv_method = RareVariantTest::Method::ACAT_V;
    RareVariantTest rv_test(rv_method);
    
    std::string method_name = options.run_skat ? "SKAT" : "BURDEN";
    if (extension) method_name = options.rare_method == "skat-o"
        ? "SKAT_O_CANONICAL" : "ACAT_V";
    std::ofstream rv_meta;
    if (extension) {
        rv_meta.open(options.out_prefix + ".rare.meta.tsv");
        if (!rv_meta) throw std::runtime_error("Cannot write rare-variant metadata");
        rv_meta << "SET_ID\tMETHOD\tN_SAMPLES\tN_VARIANTS\tSTAT\tP\tRHO_GRID\tBEST_RHO\tMIN_P\tP_METHOD\tFALLBACK_USED\n";
        rv_meta.precision(17);
    }
    std::cout << "\n[Step 2] Running " << method_name << " on " << sets.size() << " sets...\n";
    
    std::string out_file = options.out_prefix + ".rare.txt";
    
    // Create map from SNP ID to Genotype vector for fast lookup
    std::cout << "  Caching genotype data for target SNPs...\n";
    std::unordered_map<std::string, Eigen::VectorXd> snp_genotypes;
    std::unordered_map<std::string, double> snp_mafs;
    std::unordered_map<std::string, std::string> snp_chroms;
    std::unordered_map<std::string, long long> snp_positions;
    
    if (is_bgen) bgen_reader.reset();
    else if (is_pgen) pgen_reader.reset();
    else if (is_vcf) vcf_reader.reset();
    
    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage_f;
    
    size_t m_snps = is_vcf ? vcf_reader.getNumVariants() : plink_reader.getNumSnps();
    size_t current_snp = 0;
    const auto& bim = plink_reader.getBimInfo();
    
    while(true) {
        bool ok = false;
        if (is_bgen) {
            ok = bgen_reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f);
        } else if (is_pgen) {
            ok = pgen_reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f);
        } else if (is_vcf) {
            ok = vcf_reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f);
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
        if (extension && ((!extract_ids_.empty() && !extract_ids_.count(rsid)) || exclude_ids_.count(rsid))) continue;

        snp_chroms[rsid] = chrom;
        snp_positions[rsid] = static_cast<long long>(pos);
        
        Eigen::VectorXd g(n);
        double sum_g = 0.0;
        int valid_n = 0;
        for(int i=0; i<n; ++i) {
            float val = dosage_f[keep.keep_indices[i]];
            g(i) = val;
            if (std::isfinite(val) && val >= 0.0 && val <= 2.0) {
                sum_g += val;
                valid_n++;
            } else {
                g(i) = std::numeric_limits<double>::quiet_NaN();
            }
        }
        
        if (valid_n == 0) continue;
        double mean_g = sum_g / valid_n;
        double freq = mean_g / 2.0;
        double maf = std::min(freq, 1.0 - freq);
        const double mac = freq <= 0.5 ? sum_g : 2.0 * valid_n - sum_g;
        if (maf <= 0.0 || maf > options.rare_max_maf || mac < options.rare_min_mac) continue;

        if (freq > 0.5) {
            for (int i = 0; i < n; ++i) if (std::isfinite(g(i))) g(i) = 2.0 - g(i);
            mean_g = 2.0 - mean_g;
        }
        
        // Impute missing
        for(int i=0; i<n; ++i) {
            if (!std::isfinite(g(i))) g(i) = mean_g;
        }
        
        snp_genotypes[rsid] = g;
        snp_mafs[rsid] = maf;
    }
    
    std::cout << "  Found " << snp_genotypes.size() << " valid target SNPs in genotype data.\n";
    
    // 5. Run tests for each set
    struct RareOutputRow {
        std::string chrom;
        long long start = 0;
        long long end = 0;
        std::string set_id;
        std::string method;
        double statistic = 0.0;
        double p_value = std::numeric_limits<double>::quiet_NaN();
    };
    std::vector<RareOutputRow> output_rows;
    int processed = 0;
    
    for (const auto& set : sets) {
        std::vector<Eigen::VectorXd> valid_gs;
        std::vector<double> valid_mafs;
        std::string set_chrom;
        long long set_start = std::numeric_limits<long long>::max();
        long long set_end = 0;
        bool multiple_chromosomes = false;
        
        for (const auto& snp : set.snp_ids) {
            if (extension && std::count(set.snp_ids.begin(), set.snp_ids.end(), snp) > 1)
                throw std::runtime_error("Duplicate SNP in variant set: " + set.set_id);
            if (snp_genotypes.count(snp)) {
                valid_gs.push_back(snp_genotypes[snp]);
                valid_mafs.push_back(snp_mafs[snp]);
                const auto chrom_it = snp_chroms.find(snp);
                const auto pos_it = snp_positions.find(snp);
                if (chrom_it != snp_chroms.end() && pos_it != snp_positions.end()) {
                    if (set_chrom.empty()) set_chrom = chrom_it->second;
                    else if (set_chrom != chrom_it->second) multiple_chromosomes = true;
                    set_start = std::min(set_start, pos_it->second);
                    set_end = std::max(set_end, pos_it->second);
                }
            }
        }

        if (multiple_chromosomes) set_chrom = "MULTI";
        if (set_chrom.empty()) set_chrom = "NA";
        if (set_start == std::numeric_limits<long long>::max()) set_start = 0;
        
        int k = valid_gs.size();
        if (k == 0) {
            if (extension) rv_meta << set.set_id << '\t' << method_name << '\t' << n << "\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
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
        if (!std::isfinite(rv_res.p_value) || rv_res.p_value < 0.0 || rv_res.p_value > 1.0 ||
            std::isnan(rv_res.stat) ||
            (rv_method != RareVariantTest::Method::ACAT_V && !std::isfinite(rv_res.stat)))
            throw std::runtime_error("Invalid rare-variant inference for set: " + set.set_id);
        
        RareOutputRow output_row;
        output_row.chrom = set_chrom;
        output_row.start = set_start;
        output_row.end = set_end;
        output_row.set_id = rv_res.set_id;
        output_row.method = extension ? rv_res.method_name : method_name;
        output_row.statistic = rv_res.stat;
        output_row.p_value = rv_res.p_value;
        if (extension) {
            std::string rho_grid;
            for (size_t r = 0; r < rv_res.rho_grid.size(); ++r) {
                if (r) rho_grid += ',';
                rho_grid += std::to_string(rv_res.rho_grid[r]);
            }
            rv_meta << set.set_id << '\t' << method_name << '\t' << n << '\t' << k << '\t'
                    << rv_res.stat << '\t' << rv_res.p_value << '\t' << (rho_grid.empty() ? "NA" : rho_grid)
                    << '\t' << rv_res.selected_rho << '\t' << rv_res.min_p_value << '\t'
                    << (rv_res.p_method.empty() ? "NA" : rv_res.p_method) << '\t'
                    << (rv_res.fallback_used ? 1 : 0) << '\n';
        }
        output_rows.push_back(std::move(output_row));
        
        processed++;
        if (processed % 100 == 0) {
            std::cout << "\r  Processed " << processed << " / " << sets.size() << " sets..." << std::flush;
        }
    }
    
    std::cout << "\n[Step 3] Saving Results...\n";
    std::ofstream output(out_file);
    if (!output) throw std::runtime_error("Cannot write rare-variant results: " + out_file);
    output << "CHR\tSTART\tEND\tVARIANT_SET_ID\tMETHOD\tSTATISTIC\tP_VALUE\n";
    output << std::setprecision(17);
    for (const auto& row : output_rows) {
        if (!std::isfinite(row.p_value)) continue;
        output << row.chrom << '\t' << row.start << '\t' << row.end << '\t'
               << row.set_id << '\t' << row.method << '\t' << row.statistic << '\t'
               << row.p_value << '\n';
    }
    output.close();
    if (!output) throw std::runtime_error("Failed while writing rare-variant results: " + out_file);
    if (extension) {
        rv_meta.close();
        if (!rv_meta) throw std::runtime_error("Failed while writing rare-variant metadata");
    }
    std::cout << "  Rare Variant Results saved to: " << out_file << std::endl;
}

} // namespace cosmic
