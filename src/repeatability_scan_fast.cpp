#include "repeatability_model.h"
#include "repeatability_helpers.h"

#include "data_filter.h"
#include "eigen_vce.h"
#include "grm_builder.h"
#include "lmm_engine.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>

namespace cosmic {

void RepeatabilityModel::scanFast() {
    auto scan_total_start = std::chrono::high_resolution_clock::now();
    std::cout << "\nRunning Repeatability GWAS Scan (fast aggregated mean model)...\n";

    struct SubjectData { double sum_y = 0.0, count = 0.0; std::vector<double> sum_X; };
    std::vector<SubjectData> subj_data(n_samples_);
    std::vector<bool> subj_present(n_samples_, false);
    for (int i = 0; i < n_recs_kept_; ++i) {
        int gen_idx = gen_records_[i].aid - 1;
        if (gen_idx < 0 || gen_idx >= n_samples_) continue;
        subj_present[gen_idx] = true; subj_data[gen_idx].sum_y += gen_records_[i].y;
        subj_data[gen_idx].count += 1.0;
        if (subj_data[gen_idx].sum_X.empty()) subj_data[gen_idx].sum_X.resize(fd_.p - 1, 0.0);
        for (const auto& kv : fd_.rows[i])
            if (kv.first > 0 && kv.first - 1 < (int)subj_data[gen_idx].sum_X.size())
                subj_data[gen_idx].sum_X[kv.first - 1] += kv.second;
    }

    std::vector<int> agg_keep_indices; std::vector<double> y_agg_vec, n_counts;
    std::vector<std::vector<double>> X_agg_vec;
    for (int i = 0; i < n_samples_; ++i) {
        if (!subj_present[i]) continue;
        agg_keep_indices.push_back(i); double n_i = subj_data[i].count;
        n_counts.push_back(n_i); y_agg_vec.push_back(subj_data[i].sum_y / n_i);
        std::vector<double> x_row; x_row.push_back(1.0);
        for (double val : subj_data[i].sum_X) x_row.push_back(val / n_i);
        X_agg_vec.push_back(std::move(x_row));
    }

    int n_agg = (int)agg_keep_indices.size();
    std::cout << "Aggregated to " << n_agg << " unique subjects for fast repeatability scan.\n";

    Eigen::VectorXd y_agg(n_agg); Eigen::MatrixXf X_agg(n_agg, fd_.p); Eigen::VectorXd obs_weights(n_agg);
    double n_mean = 0.0, n_var = 0.0, n_min = 1e30, n_max = 0.0;
    for (int i = 0; i < n_agg; ++i) {
        y_agg(i) = y_agg_vec[i];
        for (int j = 0; j < fd_.p; ++j) X_agg(i, j) = (float)X_agg_vec[i][j];
        double cnt = n_counts[i]; n_mean += cnt; n_min = std::min(n_min, cnt); n_max = std::max(n_max, cnt);
        double var_resid = v_pe_ + v_e_ / cnt;
        obs_weights(i) = (var_resid > 1e-12) ? (1.0 / std::sqrt(var_resid)) : 1.0;
    }
    if (n_agg > 0) { n_mean /= n_agg; for (double c : n_counts) n_var += (c - n_mean) * (c - n_mean); n_var /= n_agg; }
    double mean_w = obs_weights.mean(); if (mean_w > 0.0) obs_weights /= mean_w;

    std::cout << "Repeated Measures Balance Check:\n  Mean records/subject: " << n_mean
              << "\n  Min/Max records: " << n_min << " / " << n_max
              << "\n  Variance of records: " << n_var << "\n";
    if (n_var > 0.1) std::cout << "WARNING: Fast repeatability mode is approximate for unbalanced data.\n";

    DataFilter::KeepInfo agg_keep; agg_keep.keep_indices = agg_keep_indices; agg_keep.n_keep = n_agg;

    std::string agg_pheno_file = options_.out_prefix + ".agg.pheno";
    { std::ofstream ofs(agg_pheno_file); ofs << "FID\tIID\tPHENO_MEAN\tN_RECS\n";
      for (size_t i = 0; i < agg_keep_indices.size(); ++i) {
          int idx = agg_keep_indices[i]; std::string fid = "0", iid = "0";
          if (is_pgen_) { iid = pgen_reader_->getSampleIds()[idx]; fid = iid; }
          else if (is_bgen_) { iid = bgen_reader_->getSampleIds()[idx]; fid = iid; }
          else { fid = plink_reader_->getFamInfo()[idx].fid; iid = plink_reader_->getFamInfo()[idx].iid; }
          ofs << fid << "\t" << iid << "\t" << y_agg(i) << "\t" << n_counts[i] << "\n";
      }
    }

    std::vector<std::pair<std::string, std::string>> fast_summary_fields = {
        {"mean_records", std::to_string(n_mean)},
        {"min_records", std::to_string(n_min)},
        {"max_records", std::to_string(n_max)}
    };
    writeNullSummary(options_.out_prefix + ".repeatability.null.txt", "fast_repeatability_mean_model", fast_summary_fields);

    Eigen::VectorXd y_w = y_agg; Eigen::MatrixXf X_w = X_agg;
    for (int i = 0; i < n_agg; ++i) { double w = obs_weights(i); y_w(i) *= w; X_w.row(i) *= (float)w; }

    LmmEngine::Options lmm_fast_opts;
    lmm_fast_opts.fast_mode = true; lmm_fast_opts.min_maf = options_.min_maf;
    lmm_fast_opts.min_hwe = options_.min_hwe; lmm_fast_opts.use_covariate_projection = options_.use_covariate_projection;
    lmm_fast_opts.use_spa = options_.use_spa; lmm_fast_opts.thread_num = options_.n_threads;

    std::string assoc_file = options_.out_prefix + ".assoc.txt";
    writeAssocHeader(assoc_file, options_.out_ssf, options_.eb_lmm);
    double eb_effective_tests = 0.0;
    long long eb_tested_markers = 0;

    GrmBuilder::Options grm_opts; grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe; grm_opts.algorithm = options_.grm_alg;
    #ifdef _OPENMP
    grm_opts.thread_num = omp_get_max_threads();
    #endif

    GrmBuilder global_builder(grm_opts);
    if (is_bgen_) { bgen_reader_->reset(); global_builder.compute(*bgen_reader_); }
    else if (is_pgen_) { pgen_reader_->reset(); global_builder.compute(*pgen_reader_); }
    else { global_builder.compute(*plink_reader_); }
    auto chroms_set_fast = global_builder.getObservedChroms();

    if (!options_.no_loco && chroms_set_fast.size() > 1) {
        std::cout << "Fast repeatability GWAS: LOCO scan...\n";
        double S_global_fast = (uses_per_snp_scaling(options_.grm_alg)) ? (double)global_builder.getTotalSnps() : global_builder.getSum2pq();
        const auto& G_global = global_builder.getGrm();
        auto sorted_chroms = sortChroms(chroms_set_fast);

        for (int chr_idx = 0; chr_idx < (int)sorted_chroms.size(); ++chr_idx) {
            const auto& chrom = sorted_chroms[chr_idx];
            std::cout << "  Fast LOCO chromosome " << chrom << "...\n";

            GrmBuilder::Options chr_opts = grm_opts;
            chr_opts.chrom_filter = chrom;
            GrmBuilder chr_builder(chr_opts);
            if (is_bgen_) { bgen_reader_->reset(); chr_builder.compute(*bgen_reader_); }
            else if (is_pgen_) { pgen_reader_->reset(); chr_builder.compute(*pgen_reader_); }
            else { chr_builder.compute(*plink_reader_); }

            double S_chr = (uses_per_snp_scaling(options_.grm_alg)) ? (double)chr_builder.getTotalSnps() : chr_builder.getSum2pq();
            const auto& G_chr_norm = chr_builder.getGrm();
            double S_loco = S_global_fast - S_chr;
            int n_sub = agg_keep.n_keep;
            Eigen::MatrixXd G_loco_sub(n_sub, n_sub);
            if (S_chr > 0.0 && S_loco > 1e-6) {
                double w1 = S_global_fast / S_loco;
                double w2 = S_chr / S_loco;
                #pragma omp parallel for collapse(2)
                for (int i = 0; i < n_sub; ++i) for (int j = 0; j < n_sub; ++j) {
                    int gi = agg_keep.keep_indices[i], gj = agg_keep.keep_indices[j];
                    G_loco_sub(i, j) = G_global(gi, gj) * w1 - G_chr_norm(gi, gj) * w2;
                }
            } else {
                G_loco_sub = DataFilter::subsetSquare(G_global, agg_keep);
            }
            G_loco_sub.array().colwise() *= obs_weights.array();
            G_loco_sub.array().rowwise() *= obs_weights.transpose().array();

            EigenVce::Options vce_fast_opts; vce_fast_opts.use_covariate_projection = options_.use_covariate_projection;
            vce_fast_opts.null_solver = options_.lmm_null_solver;
            vce_fast_opts.max_iter = options_.reml_maxit;
            vce_fast_opts.tol = options_.reml_tol;
            EigenVce vce_fast(vce_fast_opts); vce_fast.prepare(y_w, X_w, G_loco_sub); vce_fast.runNullModel();
            LmmEngine engine_fast(lmm_fast_opts); engine_fast.prepareFromVce(vce_fast);

            std::vector<GwasResult> chrom_results;
            if (is_bgen_) { bgen_reader_->reset(); if (options_.eb_lmm) engine_fast.scanEB(*bgen_reader_, chrom_results, agg_keep.keep_indices, chrom); else engine_fast.scan(*bgen_reader_, chrom_results, agg_keep.keep_indices, chrom); }
            else if (is_pgen_) { pgen_reader_->reset(); if (options_.eb_lmm) engine_fast.scanEB(*pgen_reader_, chrom_results, agg_keep.keep_indices, chrom); else engine_fast.scan(*pgen_reader_, chrom_results, agg_keep.keep_indices, chrom); }
            else { if (options_.eb_lmm) engine_fast.scanEB(*plink_reader_, chrom_results, agg_keep.keep_indices, chrom); else engine_fast.scan(*plink_reader_, chrom_results, agg_keep.keep_indices, chrom); }

            appendAssocResults(assoc_file, chrom_results, options_.out_ssf, options_.eb_lmm);
            if (options_.eb_lmm) accumulateEbResults(chrom_results, eb_effective_tests, eb_tested_markers);
        }
    } else {
        std::cout << "Fast repeatability GWAS: global scan on " << n_agg << " aggregated subjects.\n";
        Eigen::MatrixXd G_fast_sub = DataFilter::subsetSquare(G_subject_, agg_keep);
        G_fast_sub.array().colwise() *= obs_weights.array();
        G_fast_sub.array().rowwise() *= obs_weights.transpose().array();

        EigenVce::Options vce_fast_opts; vce_fast_opts.use_covariate_projection = options_.use_covariate_projection;
        vce_fast_opts.null_solver = options_.lmm_null_solver;
        vce_fast_opts.max_iter = options_.reml_maxit;
        vce_fast_opts.tol = options_.reml_tol;
        EigenVce vce_fast(vce_fast_opts); vce_fast.prepare(y_w, X_w, G_fast_sub); vce_fast.runNullModel();
        LmmEngine engine_fast(lmm_fast_opts); engine_fast.prepareFromVce(vce_fast);

        std::vector<GwasResult> results_fast;
        if (is_bgen_) { bgen_reader_->reset(); if (options_.eb_lmm) engine_fast.scanEB(*bgen_reader_, results_fast, agg_keep.keep_indices); else engine_fast.scan(*bgen_reader_, results_fast, agg_keep.keep_indices); }
        else if (is_pgen_) { pgen_reader_->reset(); if (options_.eb_lmm) engine_fast.scanEB(*pgen_reader_, results_fast, agg_keep.keep_indices); else engine_fast.scan(*pgen_reader_, results_fast, agg_keep.keep_indices); }
        else { if (options_.eb_lmm) engine_fast.scanEB(*plink_reader_, results_fast, agg_keep.keep_indices); else engine_fast.scan(*plink_reader_, results_fast, agg_keep.keep_indices); }
        appendAssocResults(assoc_file, results_fast, options_.out_ssf, options_.eb_lmm);
        if (options_.eb_lmm) accumulateEbResults(results_fast, eb_effective_tests, eb_tested_markers);
    }
    auto scan_total_end = std::chrono::high_resolution_clock::now();
    scan_snp_time_sec_ = std::chrono::duration<double>(scan_total_end - scan_total_start).count();
    fast_summary_fields.push_back({"time_fast_total_sec", std::to_string(scan_snp_time_sec_)});
    writeNullSummary(options_.out_prefix + ".repeatability.null.txt", "fast_repeatability_mean_model", fast_summary_fields);
    if (options_.eb_lmm) {
        std::string eb_summary = options_.out_prefix + ".repeatability.eb.summary.txt";
        writeRepeatabilityEbSummary(eb_summary, eb_effective_tests, eb_tested_markers);
        std::cout << "Repeatability EB summary saved to " << eb_summary << std::endl;
    }
    std::cout << "Fast repeatability GWAS results saved to " << assoc_file << std::endl;
}

} // namespace cosmic
