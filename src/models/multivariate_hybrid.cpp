#include "models/multivariate_model.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <omp.h>
#include "grm_builder.h"
#include "mvlmm_engine.h"

namespace cosmic {
namespace {

std::ofstream mv_open_output_file_hybrid(const std::string& path, std::ios_base::openmode mode = std::ios::out) {
    std::ofstream ofs(path, mode);
    if (!ofs) throw std::runtime_error("Failed to open output file: " + path);
    return ofs;
}

void mv_validate_output_parent(const std::string& out_prefix) {
    std::filesystem::path parent = std::filesystem::path(out_prefix).parent_path();
    if (!parent.empty() && !std::filesystem::exists(parent)) {
        throw std::runtime_error("Output directory does not exist: " + parent.string());
    }
}

void mv_write_vbeta_header_hybrid(std::ostream& os, int n_traits) {
    for (int i = 0; i < n_traits; ++i) for (int j = i; j < n_traits; ++j) os << "\tVBETA_" << (i + 1) << "_" << (j + 1);
}

std::vector<std::string> mv_split_tsv_line_hybrid(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, '\t')) fields.push_back(token);
    return fields;
}

std::unordered_map<std::string, double> mv_read_score_prepass_hybrid(const std::string& path) {
    std::unordered_map<std::string, double> scores;
    std::ifstream ifs(path);
    if (!ifs) return scores;
    std::string line;
    if (!std::getline(ifs, line)) return scores;
    auto header = mv_split_tsv_line_hybrid(line);
    int snp_col = -1, score_col = -1;
    for (int i = 0; i < (int)header.size(); ++i) {
        if (header[i] == "SNP") snp_col = i;
        if (header[i] == "P_SCORE") score_col = i;
    }
    if (snp_col < 0 || score_col < 0) return scores;
    while (std::getline(ifs, line)) {
        auto fields = mv_split_tsv_line_hybrid(line);
        if ((int)fields.size() <= std::max(snp_col, score_col)) continue;
        try { double p = std::stod(fields[score_col]); if (std::isfinite(p)) scores[fields[snp_col]] = p; } catch (...) {}
    }
    return scores;
}

} // namespace

void MultivariateModel::runHybrid() {
    mv_validate_output_parent(options_.out_prefix);

    if (is_bgen_ || is_pgen_) {
        throw std::runtime_error("--mv-method hybrid currently supports PLINK --bfile input only.");
    }
    std::cout << "Multivariate GWAS method: "
              << (options_.no_loco ? "hybrid EMmvGWAS scan + GEMMA-full top-SNP refinement"
                                   : "LOCO-aware hybrid EMmvGWAS scan + GEMMA-full top-SNP refinement")
              << "\n";
    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe;
    grm_opts.algorithm = options_.grm_alg;
    grm_opts.extract_set = options_.extract_set;
    grm_opts.exclude_set = options_.exclude_set;
    #ifdef _OPENMP
    grm_opts.thread_num = options_.n_threads;
    #endif

    std::cout << "Hybrid stage 1: " << (options_.no_loco ? "EMmvGWAS full scan" : "LOCO-aware EMmvGWAS full scan") << "...\n";
    GrmBuilder global_builder(grm_opts);
    global_builder.compute(*plink_reader_);
    const Eigen::MatrixXd& G_global = global_builder.getGrm();
    double S_global = global_builder.getSum2pq();
    auto chroms_set = global_builder.getObservedChroms();

    Eigen::MatrixXd G_global_sub(n_keep_, n_keep_);
    for (int i = 0; i < n_keep_; ++i) {
        for (int j = 0; j < n_keep_; ++j) {
            G_global_sub(i, j) = G_global(keep_.keep_indices[i], keep_.keep_indices[j]);
        }
    }

    MvEigenVce::Options vce_opts;
    vce_opts.use_projection = options_.use_covariate_projection;
    vce_opts.max_iter = options_.reml_maxit;
    vce_opts.tol = options_.reml_tol;
    vce_opts.solver = options_.use_aireml ? MvEigenVce::VceSolver::AIREML : MvEigenVce::VceSolver::EM;
    vce_opts.em_min_iter = options_.em_min_iter;
    vce_opts.ai_max_iter = options_.ai_max_iter;
    vce_opts.ai_tol = options_.ai_tol;
    vce_opts.use_emmvgwas_lambda = options_.use_emmvgwas_lambda;
    MvLmmEngine::Options mv_opts;
    mv_opts.use_projection = options_.use_covariate_projection;
    mv_opts.em_mvgwas_test = true;
    mv_opts.use_general_sigma_inverse = options_.use_emmvgwas_lambda;
    mv_opts.center_genotypes = !options_.use_emmvgwas_lambda;
    mv_opts.use_dosage_minus_one = options_.use_emmvgwas_lambda;
    mv_opts.scan_chrom = options_.scan_chrom;
    mv_opts.scan_from_bp = options_.scan_from_bp;
    mv_opts.scan_to_bp = options_.scan_to_bp;
    mv_opts.extract_set = options_.extract_set;
    mv_opts.exclude_set = options_.exclude_set;
    #ifdef _OPENMP
    mv_opts.thread_num = options_.n_threads;
    #endif
    std::string out_file = options_.out_prefix + ".assoc.mv.txt";
    {
        auto ofs = mv_open_output_file_hybrid(out_file);
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tWALD_STAT\tP";
        for (int t = 0; t < Y_.cols(); ++t) ofs << "\tBETA_" << (t + 1) << "\tSE_" << (t + 1);
        if (options_.output_vbeta) mv_write_vbeta_header_hybrid(ofs, Y_.cols());
        ofs << "\n";
    }

    std::vector<MvGwasResult> results(plink_reader_->getNumSnps());
    std::map<std::string, Eigen::MatrixXd> loco_refit_grms;
    if (options_.no_loco || chroms_set.size() <= 1) {
        MvEigenVce vce(vce_opts);
        vce.prepare(Y_, X_, G_global_sub);
        vce.runNullModel();

        MvLmmEngine engine(mv_opts);
        engine.prepareFromVce(vce);
        engine.scan(*plink_reader_, results, keep_.keep_indices);
        writeResults(out_file, results);
    } else {
        std::vector<std::string> sorted_chroms(chroms_set.begin(), chroms_set.end());
        std::sort(sorted_chroms.begin(), sorted_chroms.end(), [](const std::string& a, const std::string& b) {
            try { return std::stoi(a) < std::stoi(b); } catch (...) { return a < b; }
        });

        for (const std::string& chrom : sorted_chroms) {
            if (!options_.scan_chrom.empty() && chrom != options_.scan_chrom) continue;
            std::cout << "Hybrid LOCO screening chromosome " << chrom << "...\n";

            GrmBuilder::Options chr_opts = grm_opts;
            chr_opts.chrom_filter = chrom;
            GrmBuilder chr_builder(chr_opts);
            chr_builder.compute(*plink_reader_);

            double S_chr = chr_builder.getSum2pq();
            double S_loco = S_global - S_chr;
            const auto& G_chr = chr_builder.getGrm();
            Eigen::MatrixXd G_loco_sub = G_global_sub;
            if (S_chr > 0.0 && S_loco > 1e-6) {
                double w1 = S_global / S_loco;
                double w2 = S_chr / S_loco;
                #pragma omp parallel for collapse(2)
                for (int i = 0; i < n_keep_; ++i) {
                    for (int j = 0; j < n_keep_; ++j) {
                        int gi = keep_.keep_indices[i];
                        int gj = keep_.keep_indices[j];
                        G_loco_sub(i, j) = G_global(gi, gj) * w1 - G_chr(gi, gj) * w2;
                    }
                }
            }
            loco_refit_grms[chrom] = G_loco_sub;

            MvEigenVce vce(vce_opts);
            vce.prepare(Y_, X_, G_loco_sub);
            vce.runNullModel();

            MvLmmEngine engine(mv_opts);
            engine.prepareFromVce(vce);
            std::vector<MvGwasResult> chrom_results;
            engine.scan(*plink_reader_, chrom_results, keep_.keep_indices, chrom);
            writeResults(out_file, chrom_results);
            for (size_t i = 0; i < chrom_results.size() && i < results.size(); ++i) {
                if (!chrom_results[i].snp_id.empty() && std::isfinite(chrom_results[i].p_value)) {
                    results[i] = std::move(chrom_results[i]);
                }
            }
        }
    }

    std::vector<std::pair<double, size_t>> valid_by_p;
    valid_by_p.reserve(results.size());
    for (size_t i = 0; i < results.size(); ++i) {
        if (std::isfinite(results[i].p_value) && !results[i].snp_id.empty()) {
            valid_by_p.emplace_back(results[i].p_value, i);
        }
    }
    std::sort(valid_by_p.begin(), valid_by_p.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    std::unordered_map<std::string, double> prepass_score_p;
    size_t score_prepass_candidates = 0;
    if (options_.refine_score_prepass > 0 && options_.refine_discrepancy_top > 0 && !valid_by_p.empty()) {
        score_prepass_candidates = std::min(valid_by_p.size(), (size_t)options_.refine_score_prepass);
        refine_mask_.assign(plink_reader_->getNumSnps(), 0);
        for (size_t j = 0; j < score_prepass_candidates; ++j) {
            refine_mask_[valid_by_p[j].second] = 1;
        }

        bool old_full = options_.use_gemma_full;
        bool old_refit = options_.use_gemma_refit;
        bool old_no_loco = options_.no_loco;
        double old_refit_p = options_.refit_p;
        std::string old_out_prefix = options_.out_prefix;
        options_.use_gemma_full = true;
        options_.use_gemma_refit = true;
        options_.no_loco = true;
        options_.refit_p = 0.0;
        options_.out_prefix = old_out_prefix + ".score_prepass";

        std::cout << "Hybrid Stage 1.5 Score-only prepass on top " << score_prepass_candidates << " SNPs...\n";
        if (!loco_refit_grms.empty()) {
            bool first_refit = true;
            for (auto& entry : loco_refit_grms) {
                refit_grm_override_ = &entry.second;
                refit_chrom_filter_ = entry.first;
                refit_append_mode_ = !first_refit;
                runGemmaRefit();
                first_refit = false;
            }
            refit_grm_override_ = nullptr;
            refit_chrom_filter_.clear();
            refit_append_mode_ = false;
        } else {
            runGemmaRefit();
        }

        prepass_score_p = mv_read_score_prepass_hybrid(options_.out_prefix + ".refined.assoc.mv.txt");
        options_.out_prefix = old_out_prefix;
        options_.use_gemma_full = old_full;
        options_.use_gemma_refit = old_refit;
        options_.no_loco = old_no_loco;
        options_.refit_p = old_refit_p;
        refine_mask_.clear();
        std::cout << "Hybrid Score-only prepass produced " << prepass_score_p.size() << " score p-values.\n";
    }

    std::vector<char> selected(results.size(), 0), selected_by_p(results.size(), 0), selected_by_top(results.size(), 0), selected_by_discrepancy(results.size(), 0);
    size_t candidates_by_p = 0;
    size_t candidates_by_top = 0;
    bool p_threshold_enabled = options_.refine_p < 1.0;
    for (const auto& item : valid_by_p) {
        if (p_threshold_enabled && item.first <= options_.refine_p) {
            selected[item.second] = 1;
            selected_by_p[item.second] = 1;
            ++candidates_by_p;
        }
    }
    if (options_.refine_top > 0) {
        size_t top_n = std::min(valid_by_p.size(), (size_t)options_.refine_top);
        for (size_t j = 0; j < top_n; ++j) {
            selected[valid_by_p[j].second] = 1;
            selected_by_top[valid_by_p[j].second] = 1;
        }
        candidates_by_top = top_n;
    }

    auto neg_log10_p = [](double p) {
        if (!std::isfinite(p) || p <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        return -std::log10(std::max(p, 1e-300));
    };
    std::vector<std::pair<double, size_t>> valid_by_discrepancy;
    valid_by_discrepancy.reserve(results.size());
    for (size_t i = 0; i < results.size(); ++i) {
        if (results[i].snp_id.empty()) continue;
        double p_wald = std::isfinite(results[i].p_wald) ? results[i].p_wald : results[i].p_value;
        double p_score = results[i].p_score;
        auto prepass_it = prepass_score_p.find(results[i].snp_id);
        if (prepass_it != prepass_score_p.end()) p_score = prepass_it->second;
        double lw = neg_log10_p(p_wald);
        double ls = neg_log10_p(p_score);
        if (std::isfinite(lw) && std::isfinite(ls)) {
            double discrepancy = std::abs(ls - lw);
            if (discrepancy >= options_.refine_discrepancy) {
                valid_by_discrepancy.emplace_back(discrepancy, i);
            }
        }
    }
    std::sort(valid_by_discrepancy.begin(), valid_by_discrepancy.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    size_t candidates_by_discrepancy = 0;
    if (options_.refine_discrepancy_top > 0) {
        size_t disc_n = std::min(valid_by_discrepancy.size(), (size_t)options_.refine_discrepancy_top);
        for (size_t j = 0; j < disc_n; ++j) {
            selected[valid_by_discrepancy[j].second] = 1;
            selected_by_discrepancy[valid_by_discrepancy[j].second] = 1;
        }
        candidates_by_discrepancy = disc_n;
    }

    std::vector<std::pair<double, size_t>> candidates;
    candidates.reserve(valid_by_p.size());
    for (const auto& item : valid_by_p) {
        if (selected[item.second]) candidates.push_back(item);
    }
    size_t candidates_before_top = candidates_by_p;

    refine_mask_.assign(plink_reader_->getNumSnps(), 0);
    std::string snplist = options_.out_prefix + ".refine.snplist";
    std::string summary_file = options_.out_prefix + ".hybrid.summary.tsv";
    {
        auto ofs = mv_open_output_file_hybrid(snplist);
        for (const auto& item : candidates) {
            refine_mask_[item.second] = 1;
            ofs << results[item.second].snp_id << "\t" << item.first << "\n";
        }
    }
    {
        auto ofs = mv_open_output_file_hybrid(summary_file);
        ofs << "stage1_mode\t" << ((options_.no_loco || chroms_set.size() <= 1) ? "global" : "loco") << "\n";
        ofs << "screening_method\tEMmvGWAS semi-exact Wald\n";
        ofs << "refinement_method\t" << (loco_refit_grms.empty() ? "GEMMA-full Wald/LRT/Score"
                                                                 : "per-chrom LOCO GEMMA-full Wald/LRT/Score") << "\n";
        ofs << "refine_p\t" << options_.refine_p << "\n";
        ofs << "refine_p_threshold_enabled\t" << (p_threshold_enabled ? "yes" : "no") << "\n";
        ofs << "refine_top\t" << options_.refine_top << "\n";
        ofs << "refine_discrepancy\t" << options_.refine_discrepancy << "\n";
        ofs << "refine_discrepancy_top\t" << options_.refine_discrepancy_top << "\n";
        ofs << "refine_score_prepass\t" << options_.refine_score_prepass << "\n";
        ofs << "score_prepass_candidates\t" << score_prepass_candidates << "\n";
        ofs << "score_prepass_scores\t" << prepass_score_p.size() << "\n";
        ofs << "discrepancy_available\t" << (!prepass_score_p.empty() || !valid_by_discrepancy.empty() ? "yes" : "no") << "\n";
        ofs << "candidates_by_p\t" << candidates_by_p << "\n";
        ofs << "candidates_by_top\t" << candidates_by_top << "\n";
        ofs << "candidates_by_discrepancy\t" << candidates_by_discrepancy << "\n";
        ofs << "candidates_before_top_filter\t" << candidates_before_top << "\n";
        ofs << "candidates_selected\t" << candidates.size() << "\n";
        ofs << "snplist\t" << snplist << "\n";
    }
    std::cout << "Hybrid selected " << candidates.size() << " SNPs for GEMMA-full refinement: " << snplist << "\n";
    if (candidates.empty()) {
        std::cout << "Hybrid refinement skipped because no SNP passed --mv-refine-p.\n";
        return;
    }

    bool old_full = options_.use_gemma_full;
    bool old_refit = options_.use_gemma_refit;
    bool old_no_loco = options_.no_loco;
    options_.use_gemma_full = true;
    options_.use_gemma_refit = true;
    options_.no_loco = true;

    if (!loco_refit_grms.empty()) {
        bool first_refit = true;
        for (auto& entry : loco_refit_grms) {
            const std::string& chrom = entry.first;
            refit_grm_override_ = &entry.second;
            refit_chrom_filter_ = chrom;
            refit_append_mode_ = !first_refit;
            std::cout << "Hybrid Stage 2 LOCO GEMMA-full refinement chromosome " << chrom << "...\n";
            runGemmaRefit();
            first_refit = false;
        }
        refit_grm_override_ = nullptr;
        refit_chrom_filter_.clear();
        refit_append_mode_ = false;
    } else {
        runGemmaRefit();
    }

    options_.use_gemma_full = old_full;
    options_.use_gemma_refit = old_refit;
    options_.no_loco = old_no_loco;
    refine_mask_.clear();
}


} // namespace cosmic
