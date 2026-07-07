#include "gwas_scanner.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <numeric>
#include <omp.h>
#include <stdexcept>
#include <unordered_map>
#include <chrono>
#include <iomanip>
#include <unordered_set>
#include <iterator>
#include "design.h"
#include "data_filter.h"
#include "lmm_engine.h"
#include "mvlmm_engine.h"

namespace cosmic {

// Forward declaration
std::unordered_set<std::string> load_snp_id_set(const std::string& path);

GwasScanner::GwasScanner(const Options& opts) : options(opts) {
    if (!opts.extract_file.empty()) {
        extract_ids_ = load_snp_id_set(opts.extract_file);
        std::cout << "  --extract: " << extract_ids_.size() << " variant IDs loaded." << std::endl;
    }
    if (!opts.exclude_file.empty()) {
        exclude_ids_ = load_snp_id_set(opts.exclude_file);
        std::cout << "  --exclude: " << exclude_ids_.size() << " variant IDs loaded." << std::endl;
    }
}

double elapsed_seconds(const std::chrono::steady_clock::time_point& start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

void log_stage_time(const std::string& label, const std::chrono::steady_clock::time_point& start) {
    std::cout << "  [time] " << label << ": " << std::fixed << std::setprecision(2)
              << elapsed_seconds(start) << "s" << std::defaultfloat << std::endl;
}

std::unordered_set<std::string> load_snp_id_set(const std::string& path) {
    std::unordered_set<std::string> ids;
    if (path.empty()) return ids;
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open SNP list: " + path);
    std::string id;
    while (ifs >> id) ids.insert(id);
    if (ids.empty()) throw std::runtime_error("SNP list is empty: " + path);
    return ids;
}

std::vector<std::string> sort_chrom_labels(const std::set<std::string>& chroms_set) {
    std::vector<std::string> chroms(chroms_set.begin(), chroms_set.end());
    std::sort(chroms.begin(), chroms.end(), [](const std::string& a, const std::string& b) {
        try {
            return std::stoi(a) < std::stoi(b);
        } catch (...) {
            return a < b;
        }
    });
    return chroms;
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
    if (!is_bgen && !is_pgen && options.grm_alg != 2) {
        auto total_start = std::chrono::steady_clock::now();
        GrmBuilder::Options grm_opts;
        grm_opts.min_maf = (float)options.min_maf;
        grm_opts.min_hwe = (float)options.min_hwe;
        grm_opts.algorithm = options.grm_alg;
        grm_opts.extract_set = extract_ids_;
        grm_opts.exclude_set = exclude_ids_;
        #ifdef _OPENMP
        grm_opts.thread_num = options.n_threads;
        #endif

        std::set<std::string> chrom_labels;
        for (const auto& snp : plink_reader.getBimInfo()) chrom_labels.insert(snp.chrom);

        const bool use_loco = chrom_labels.size() > 1 && !options.no_loco;

        Eigen::MatrixXd G_global;
        double S_global = 0.0;
        GrmBuilder::SplitGrmResult split;
        if (use_loco) {
            std::cout << "Calculating Global + per-chromosome GRMs in one PLINK pass...\n";
            auto grm_start = std::chrono::steady_clock::now();
            GrmBuilder split_builder(grm_opts);
            split = split_builder.computeSplit(plink_reader, chrom_labels);
            S_global = (options.grm_alg == 0 || options.grm_alg == 3)
                ? (double)split.global_snp_count
                : split.global_sum_2pq;
            G_global = split.global;
            log_stage_time("PLINK split GRM build", grm_start);
        } else {
            std::cout << "Calculating Global GRM in one PLINK pass...\n";
            auto grm_start = std::chrono::steady_clock::now();
            GrmBuilder global_builder(grm_opts);
            global_builder.compute(plink_reader);
            S_global = (options.grm_alg == 0 || options.grm_alg == 3)
                ? (double)global_builder.getTotalSnps()
                : global_builder.getSum2pq();
            G_global = global_builder.getGrm();
            log_stage_time("PLINK global GRM build", grm_start);
        }
        std::cout << "Global GRM computed. Scaling Factor = " << S_global << std::endl;

        std::string out_file = options.out_prefix + ".assoc.txt";
        writeHeader(out_file);

        Eigen::VectorXd y_sub = DataFilter::subset(y, keep);
        Eigen::MatrixXf X_sub = DataFilter::subsetRows(X, keep);
        bool has_weights = (options.obs_weights.size() == keep.n_keep);
        if (has_weights) {
            y_sub.array() *= options.obs_weights.array();
            for (int j = 0; j < X_sub.cols(); ++j) {
                X_sub.col(j).array() *= options.obs_weights.cast<float>().array();
            }
        }

        LmmEngine::Options lmm_opts;
        lmm_opts.fast_mode = options.fast_lmm;
        lmm_opts.gemma_refit = options.lmm_gemma_refit;
        lmm_opts.min_maf = options.min_maf;
        lmm_opts.min_hwe = options.min_hwe;
        lmm_opts.use_covariate_projection = options.use_covariate_projection;
        lmm_opts.use_spa = options.use_spa;
        lmm_opts.extract_set = extract_ids_;
        lmm_opts.exclude_set = exclude_ids_;
        if (has_weights) lmm_opts.obs_weights = options.obs_weights;
        if (options.grm_alg == 0) lmm_opts.use_residuals = true;
        #ifdef _OPENMP
        lmm_opts.thread_num = options.n_threads;
        #endif

        auto run_scan_with_grm = [&](const Eigen::MatrixXd& G_sub, const std::string& chrom_filter) {
            auto null_start = std::chrono::steady_clock::now();
            Eigen::MatrixXd G_work = G_sub;
            if (has_weights) {
                G_work.array().colwise() *= options.obs_weights.array();
                G_work.array().rowwise() *= options.obs_weights.transpose().array();
            }

            EigenVce::Options vce_opts;
            vce_opts.use_covariate_projection = options.use_covariate_projection;
            vce_opts.null_solver = options.lmm_null_solver;
            vce_opts.max_iter = options.reml_maxit;
            vce_opts.tol = options.reml_tol;
            EigenVce vce(vce_opts);
            vce.prepare(y_sub, X_sub, G_work);
            vce.runNullModel();
            std::cout << "  Null Model VCs: V(G)=" << vce.getSigma2g()
                      << ", V(e)=" << vce.getSigma2e() << std::endl;
            log_stage_time(chrom_filter.empty() ? "null model" : ("null model chr " + chrom_filter), null_start);

            auto scan_start = std::chrono::steady_clock::now();
            LmmEngine engine(lmm_opts);
            engine.prepareFromVce(vce);
            std::vector<GwasResult> results;
            if (options.eb_lmm) engine.scanEB(plink_reader, results, keep.keep_indices, chrom_filter);
            else engine.scan(plink_reader, results, keep.keep_indices, chrom_filter);
            appendResults(out_file, results);
            log_stage_time(chrom_filter.empty() ? "SNP scan" : ("SNP scan chr " + chrom_filter), scan_start);
        };

        if (!use_loco) {
            std::cout << (options.no_loco ? "LOCO disabled via flag. " : "Only 1 chromosome found. ");
            std::cout << "Using standard LMM (no LOCO)...\n";
            run_scan_with_grm(DataFilter::subsetSquare(G_global, keep), "");
        } else {
            std::cout << "Running LOCO Scan on " << split.chrom_labels.size()
                      << " chromosomes using cached split GRMs...\n";
            const int n_sub = keep.n_keep;
            Eigen::MatrixXd G_global_sub = DataFilter::subsetSquare(G_global, keep);
            for (size_t c = 0; c < split.chrom_labels.size(); ++c) {
                const std::string& chrom = split.chrom_labels[c];
                std::cout << "Processing Chromosome " << chrom << "...\n";
                double S_chr = (options.grm_alg == 0 || options.grm_alg == 3)
                    ? (double)split.per_chrom_snp_count[c]
                    : split.per_chrom_sum_2pq[c];
                double S_loco = S_global - S_chr;
                Eigen::MatrixXd G_loco_sub = G_global_sub;
                if (S_chr > 0.0 && S_loco > 1e-6) {
                    Eigen::MatrixXd raw_chr_sub(n_sub, n_sub);
                    #pragma omp parallel for collapse(2)
                    for (int i = 0; i < n_sub; ++i) {
                        for (int j = 0; j < n_sub; ++j) {
                            raw_chr_sub(i, j) = split.per_chrom_raw[c](keep.keep_indices[i], keep.keep_indices[j]);
                        }
                    }
                    G_loco_sub = G_global_sub * (S_global / S_loco) - raw_chr_sub / S_loco;
                }
                run_scan_with_grm(G_loco_sub, chrom);
            }
        }

        std::cout << "Results saved to " << out_file << std::endl;
        if (options.eb_lmm) {
            std::string summary_file = options.out_prefix + ".eb.summary.txt";
            writeEbSummary(summary_file);
            std::cout << "EB summary saved to " << summary_file << std::endl;
        }
        log_stage_time("PLINK optimized GWAS total", total_start);
        return;
    }

    std::cout << "Calculating Global GRM (Pass 1)...\n";
    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options.min_maf;
    grm_opts.min_hwe = (float)options.min_hwe;
    grm_opts.algorithm = options.grm_alg;
    grm_opts.extract_set = extract_ids_;
    grm_opts.exclude_set = exclude_ids_;
    #ifdef _OPENMP
    grm_opts.thread_num = options.n_threads;
    #endif

    auto chroms_set = [&]() -> std::set<std::string> {
        if (is_bgen) {
            bgen_reader.loadVariantInfo();
            std::set<std::string> s;
            for (const auto& v : bgen_reader.getVariants()) s.insert(v.chrom);
            return s;
        } else if (is_pgen) {
            std::set<std::string> s;
            for (const auto& v : pgen_reader.getVariants()) s.insert(v.chrom);
            return s;
        } else {
            std::set<std::string> s;
            for (const auto& snp : plink_reader.getBimInfo()) s.insert(snp.chrom);
            return s;
        }
    }();

    Eigen::MatrixXd G_global_norm;
    double S_global = 0.0;
    GrmBuilder::SplitGrmResult split_result;
    const bool use_loco = chroms_set.size() > 1 && !options.no_loco;

    if (is_bgen && use_loco) {
        GrmBuilder split_builder(grm_opts);
        bgen_reader.reset();
        split_result = split_builder.computeSplit(bgen_reader, chroms_set);
        S_global = (options.grm_alg == 0 || options.grm_alg == 3)
            ? (double)split_result.global_snp_count
            : split_result.global_sum_2pq;
        G_global_norm = split_result.global;
    } else if (is_pgen && use_loco) {
        GrmBuilder split_builder(grm_opts);
        pgen_reader.reset();
        split_result = split_builder.computeSplit(pgen_reader, chroms_set);
        S_global = (options.grm_alg == 0 || options.grm_alg == 3)
            ? (double)split_result.global_snp_count
            : split_result.global_sum_2pq;
        G_global_norm = split_result.global;
    } else {
        GrmBuilder global_builder(grm_opts);
        if (is_bgen) { bgen_reader.reset(); global_builder.compute(bgen_reader); }
        else if (is_pgen) { pgen_reader.reset(); global_builder.compute(pgen_reader); }
        else { global_builder.compute(plink_reader); }
        G_global_norm = global_builder.getGrm();
        S_global = (options.grm_alg == 0) ? (double)global_builder.getTotalSnps() : global_builder.getSum2pq();
    }

    std::cout << "Global GRM computed. Scaling Factor = " << S_global << std::endl;

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
    lmm_opts.gemma_refit = options.lmm_gemma_refit;
    lmm_opts.min_maf = options.min_maf;
    lmm_opts.min_hwe = options.min_hwe;
    lmm_opts.use_covariate_projection = options.use_covariate_projection;
    lmm_opts.use_spa = options.use_spa;
    lmm_opts.extract_set = extract_ids_;
    lmm_opts.exclude_set = exclude_ids_;
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
        vce_opts.null_solver = options.lmm_null_solver;
        vce_opts.max_iter = options.reml_maxit;
        vce_opts.tol = options.reml_tol;
        EigenVce vce(vce_opts);
        vce.prepare(y_sub, X_sub, G_sub);
        vce.runNullModel();

        std::cout << "  Null Model VCs: V(G)=" << vce.getSigma2g() << ", V(e)=" << vce.getSigma2e() << std::endl;

        LmmEngine engine(lmm_opts);
        engine.prepareFromVce(vce);

        std::vector<GwasResult> results;
        if (is_bgen) {
            bgen_reader.reset();
            if (options.eb_lmm) engine.scanEB(bgen_reader, results, keep.keep_indices);
            else engine.scan(bgen_reader, results, keep.keep_indices);
        } else if (is_pgen) {
            pgen_reader.reset();
            if (options.eb_lmm) engine.scanEB(pgen_reader, results, keep.keep_indices);
            else engine.scan(pgen_reader, results, keep.keep_indices);
        } else {
            if (options.eb_lmm) engine.scanEB(plink_reader, results, keep.keep_indices);
            else engine.scan(plink_reader, results, keep.keep_indices);
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

            double S_chr = 0.0;
            const Eigen::MatrixXd* G_chr_norm_ptr = nullptr;
            Eigen::MatrixXd G_chr_single;

            if (is_bgen && use_loco && !split_result.chrom_labels.empty()) {
                // Use cached split GRM result
                for (size_t ci = 0; ci < split_result.chrom_labels.size(); ++ci) {
                    if (split_result.chrom_labels[ci] == chrom) {
                        G_chr_norm_ptr = &split_result.per_chrom_raw[ci];
                        S_chr = (options.grm_alg == 0 || options.grm_alg == 3)
                            ? (double)split_result.per_chrom_snp_count[ci]
                            : split_result.per_chrom_sum_2pq[ci];
                        break;
                    }
                }
            } else if (is_pgen && use_loco && !split_result.chrom_labels.empty()) {
                for (size_t ci = 0; ci < split_result.chrom_labels.size(); ++ci) {
                    if (split_result.chrom_labels[ci] == chrom) {
                        G_chr_norm_ptr = &split_result.per_chrom_raw[ci];
                        S_chr = (options.grm_alg == 0 || options.grm_alg == 3)
                            ? (double)split_result.per_chrom_snp_count[ci]
                            : split_result.per_chrom_sum_2pq[ci];
                        break;
                    }
                }
            }

            if (!G_chr_norm_ptr) {
                // Fallback: compute per-chromosome GRM (single pass)
                GrmBuilder::Options chr_opts = grm_opts;
                chr_opts.chrom_filter = chrom;
                GrmBuilder chr_builder(chr_opts);
                if (is_bgen) { bgen_reader.reset(); chr_builder.compute(bgen_reader); }
                else if (is_pgen) { pgen_reader.reset(); chr_builder.compute(pgen_reader); }
                else { chr_builder.compute(plink_reader); }
                S_chr = (options.grm_alg == 0) ? (double)chr_builder.getTotalSnps() : chr_builder.getSum2pq();
                G_chr_single = chr_builder.getGrm();
                G_chr_norm_ptr = &G_chr_single;
            }

            const auto& G_chr_norm = *G_chr_norm_ptr;

            // Calculate LOCO GRM Subset
            int n_sub = keep.n_keep;
            Eigen::MatrixXd G_loco_sub(n_sub, n_sub);

            double S_loco = S_global - S_chr;

            if (S_chr > 0 && S_loco > 1e-6) {
                bool is_split_chr = (G_chr_norm_ptr &&
                    (is_bgen || is_pgen) && use_loco && !split_result.chrom_labels.empty());

                if (is_split_chr) {
                    // per_chrom_raw is raw (un-normalized), use PLINK-style formula
                    #pragma omp parallel for collapse(2)
                    for (int i = 0; i < n_sub; ++i) {
                        for (int j = 0; j < n_sub; ++j) {
                            int gi = keep.keep_indices[i];
                            int gj = keep.keep_indices[j];
                            G_loco_sub(i, j) = G_global_norm(gi, gj) * (S_global / S_loco)
                                              - G_chr_norm(gi, gj) / S_loco;
                        }
                    }
                } else {
                    // Normalized per-chromosome GRM
                    double w1 = (S_global / S_loco);
                    double w2 = (S_chr / S_loco);
                    #pragma omp parallel for collapse(2)
                    for (int i = 0; i < n_sub; ++i) {
                        for (int j = 0; j < n_sub; ++j) {
                            int gi = keep.keep_indices[i];
                            int gj = keep.keep_indices[j];
                            G_loco_sub(i, j) = G_global_norm(gi, gj) * w1 - G_chr_norm(gi, gj) * w2;
                        }
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
            vce_opts.null_solver = options.lmm_null_solver;
            vce_opts.max_iter = options.reml_maxit;
            vce_opts.tol = options.reml_tol;
            EigenVce vce(vce_opts);
            vce.prepare(y_sub, X_sub, G_loco_sub);
            vce.runNullModel();

            std::cout << "  Null Model VCs: V(G)=" << vce.getSigma2g() << ", V(e)=" << vce.getSigma2e() << std::endl;

            // Run LMM for this chromosome
            LmmEngine engine(lmm_opts);
            engine.prepareFromVce(vce);

            // Scan
            std::vector<GwasResult> chrom_results;
            if (is_bgen) { bgen_reader.reset(); if (options.eb_lmm) engine.scanEB(bgen_reader, chrom_results, keep.keep_indices, chrom); else engine.scan(bgen_reader, chrom_results, keep.keep_indices, chrom); }
            else if (is_pgen) { pgen_reader.reset(); if (options.eb_lmm) engine.scanEB(pgen_reader, chrom_results, keep.keep_indices, chrom); else engine.scan(pgen_reader, chrom_results, keep.keep_indices, chrom); }
            else { if (options.eb_lmm) engine.scanEB(plink_reader, chrom_results, keep.keep_indices, chrom); else engine.scan(plink_reader, chrom_results, keep.keep_indices, chrom); }

            appendResults(out_file, chrom_results);
        }
    }
    std::cout << "Results saved to " << out_file << std::endl;
    if (options.eb_lmm) {
        std::string summary_file = options.out_prefix + ".eb.summary.txt";
        writeEbSummary(summary_file);
        std::cout << "EB summary saved to " << summary_file << std::endl;
    }
}


} // namespace cosmic
