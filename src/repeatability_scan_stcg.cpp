#include "repeatability_model.h"
#include "repeatability_helpers.h"

#include "data_filter.h"
#include "grm_builder.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

namespace cosmic {

void RepeatabilityModel::scanSTCG() {
    std::cout << "\nRunning Repeatability GWAS Scan (STCG PCG-based)...\n";

    int n_obs = n_recs_kept_;
    std::vector<int> keep_subject_indices = obs_to_sample_;
    std::sort(keep_subject_indices.begin(), keep_subject_indices.end());
    keep_subject_indices.erase(std::unique(keep_subject_indices.begin(), keep_subject_indices.end()), keep_subject_indices.end());
    int n_subjects_keep = (int)keep_subject_indices.size();

    std::vector<int> sample_to_subject_rank(n_samples_, -1);
    for (int i = 0; i < n_subjects_keep; ++i) sample_to_subject_rank[keep_subject_indices[i]] = i;
    std::vector<int> obs_to_subject(n_obs);
    for (int i = 0; i < n_obs; ++i) obs_to_subject[i] = sample_to_subject_rank[obs_to_sample_[i]];

    std::cout << "  Records: " << n_obs << ", Subjects: " << n_subjects_keep << "\n";
    std::cout << "  V(G)=" << v_g_ << " V(PE)=" << v_pe_ << " V(E)=" << v_e_ << "\n";

    std::vector<std::vector<int>> subj_to_rec(n_subjects_keep);
    for (int i = 0; i < n_obs; ++i) subj_to_rec[obs_to_subject[i]].push_back(i);

    std::vector<int> n_obs_per_subj(n_subjects_keep, 0);
    for (int i = 0; i < n_obs; ++i) n_obs_per_subj[obs_to_subject[i]]++;

    auto scan_stcg_with_grm = [&](const Eigen::MatrixXd& G_loco,
                                  const std::string& chrom_filter,
                                  std::ofstream& ofs) {
        int q = n_subjects_keep;
        int p = X_obs_.cols();
        double vg = v_g_, vpe = v_pe_, ve = v_e_;

        Eigen::VectorXd Zty(q);
        Eigen::MatrixXd ZtX(q, p);
        Zty.setZero(); ZtX.setZero();
        for (int i = 0; i < n_obs; ++i) {
            int s = obs_to_subject[i];
            Zty[s] += y_obs_[i];
            for (int j = 0; j < p; ++j) ZtX(s, j) += X_obs_(i, j);
        }

        Eigen::VectorXd D_diag(q);
        for (int s = 0; s < q; ++s) D_diag[s] = n_obs_per_subj[s];

        Eigen::MatrixXd K = vg * G_loco + vpe * Eigen::MatrixXd::Identity(q, q);
        Eigen::LDLT<Eigen::MatrixXd> K_ldlt(K);
        Eigen::MatrixXd Kinv_Zty = K_ldlt.solve(Zty);
        Eigen::MatrixXd Kinv_ZtX = K_ldlt.solve(ZtX);

        Eigen::MatrixXd M = K_ldlt.solve(Eigen::MatrixXd::Identity(q, q));
        for (int s = 0; s < q; ++s) M(s, s) += D_diag[s] / ve;
        Eigen::LDLT<Eigen::MatrixXd> M_ldlt(M);

        auto woodbury_Vinv_vec = [&](const Eigen::VectorXd& v) -> Eigen::VectorXd {
            Eigen::VectorXd Ztv(q);
            Ztv.setZero();
            for (int i = 0; i < n_obs; ++i) Ztv[obs_to_subject[i]] += v[i];
            Eigen::VectorXd Minv_Ztv = M_ldlt.solve(Ztv);
            Eigen::VectorXd out(n_obs);
            for (int i = 0; i < n_obs; ++i) {
                out[i] = v[i] / ve - Minv_Ztv[obs_to_subject[i]] / (ve * ve);
            }
            return out;
        };

        Eigen::VectorXd Vinv_y = woodbury_Vinv_vec(y_obs_);
        Eigen::MatrixXd Vinv_X(n_obs, p);
        for (int j = 0; j < p; ++j) {
            Eigen::VectorXd X_col = X_obs_.col(j).cast<double>();
            Vinv_X.col(j) = woodbury_Vinv_vec(X_col);
        }
        Eigen::MatrixXd XtVinvX = X_obs_.cast<double>().transpose() * Vinv_X;
        Eigen::LDLT<Eigen::MatrixXd> x_solver(XtVinvX);
        Eigen::VectorXd beta;
        if (x_solver.info() == Eigen::Success) {
            Eigen::VectorXd Xty = X_obs_.cast<double>().transpose() * Vinv_y;
            beta = x_solver.solve(Xty);
        } else {
            beta = Eigen::VectorXd::Zero(p);
        }
        Eigen::VectorXd resid = y_obs_ - X_obs_.cast<double>() * beta;
        Eigen::VectorXd P_y = woodbury_Vinv_vec(resid);

        if (chrom_filter.empty()) {
            std::cout << "  Null model (Woodbury) completed (global). Beta[0]=" << (beta.size() > 0 ? beta[0] : 0.0) << "\n";
        } else {
            std::cout << "  Null model (Woodbury) completed (LOCO chr " << chrom_filter << ").\n";
        }

        auto process_snp = [&](const std::string& chrom, const std::string& rsid, uint32_t pos,
                               const std::string& ref, const std::string& alt,
                               const Eigen::VectorXf& dosage) {
            Eigen::VectorXd g(n_obs);
            double sum_g = 0.0;
            int valid_n = 0;
            for (int i = 0; i < n_obs; ++i) {
                float val = dosage[obs_to_sample_[i]];
                g[i] = val;
                if (!std::isnan(val)) { sum_g += val; valid_n++; }
            }
            if (valid_n < 10) return;
            double maf = sum_g / (2.0 * valid_n);
            if (maf > 0.5) maf = 1.0 - maf;
            if (maf < options_.min_maf) return;
            double g_mean = sum_g / valid_n;
            for (int i = 0; i < n_obs; ++i) {
                g[i] = std::isnan(g[i]) ? 0.0 : g[i] - g_mean;
            }

            double score = g.dot(P_y);
            Eigen::VectorXd Vinv_g = woodbury_Vinv_vec(g);
            Eigen::VectorXd XtVinvX_XtVinv_g = x_solver.solve(X_obs_.cast<double>().transpose() * Vinv_g);
            Eigen::VectorXd P_g = Vinv_g - Vinv_X * XtVinvX_XtVinv_g;
            double score_var = g.dot(P_g);
            if (score_var < 1e-15) return;
            double z = score / std::sqrt(score_var);
            double pval = std::erfc(std::abs(z) / std::sqrt(2.0));
            double snp_beta = score / score_var;
            double se = 1.0 / std::sqrt(score_var);
            double log10p = (pval > 0.0 && pval <= 1.0) ? -std::log10(pval) : (pval > 1.0 ? 0.0 : 300.0);

            ofs << chrom << "\t" << rsid << "\t" << pos << "\t"
                << alt << "\t" << ref << "\t" << maf << "\t" << n_subjects_keep << "\t"
                << snp_beta << "\t" << se << "\t" << pval << "\t" << log10p << "\n";
        };

        if (is_bgen_) {
            bgen_reader_->reset();
            size_t m = bgen_reader_->getNumVariants();
            for (size_t j = 0; j < m; ++j) {
                std::string chrom, rsid, ref, alt; uint32_t pos; Eigen::VectorXf dosage;
                if (!bgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage)) break;
                if (!chrom_filter.empty() && chrom != chrom_filter) continue;
                if (options_.scan_from_bp >= 0 && static_cast<long long>(pos) < options_.scan_from_bp) continue;
                if (options_.scan_to_bp >= 0 && static_cast<long long>(pos) > options_.scan_to_bp) continue;
                if (!options_.extract_set.empty() && options_.extract_set.find(rsid) == options_.extract_set.end()) continue;
                if (!options_.exclude_set.empty() && options_.exclude_set.find(rsid) != options_.exclude_set.end()) continue;
                process_snp(chrom, rsid, pos, ref, alt, dosage);
            }
        } else if (is_pgen_) {
            pgen_reader_->reset();
            size_t m = pgen_reader_->getNumVariants();
            for (size_t j = 0; j < m; ++j) {
                std::string chrom, rsid, ref, alt; uint32_t pos; Eigen::VectorXf dosage;
                if (!pgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage)) break;
                if (!chrom_filter.empty() && chrom != chrom_filter) continue;
                if (options_.scan_from_bp >= 0 && static_cast<long long>(pos) < options_.scan_from_bp) continue;
                if (options_.scan_to_bp >= 0 && static_cast<long long>(pos) > options_.scan_to_bp) continue;
                if (!options_.extract_set.empty() && options_.extract_set.find(rsid) == options_.extract_set.end()) continue;
                if (!options_.exclude_set.empty() && options_.exclude_set.find(rsid) != options_.exclude_set.end()) continue;
                process_snp(chrom, rsid, pos, ref, alt, dosage);
            }
        } else {
            const auto& bim = plink_reader_->getBimInfo();
            int m = plink_reader_->getNumSnps();
            for (int j = 0; j < m; ++j) {
                if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;
                if (options_.scan_from_bp >= 0 && bim[j].bp_pos < options_.scan_from_bp) continue;
                if (options_.scan_to_bp >= 0 && bim[j].bp_pos > options_.scan_to_bp) continue;
                if (!options_.extract_set.empty() && options_.extract_set.find(bim[j].id) == options_.extract_set.end()) continue;
                if (!options_.exclude_set.empty() && options_.exclude_set.find(bim[j].id) != options_.exclude_set.end()) continue;
                Eigen::VectorXf dosage;
                plink_reader_->getSnpGenotypes(j, dosage);
                process_snp(bim[j].chrom, bim[j].id, bim[j].bp_pos, bim[j].ref, bim[j].alt, dosage);
            }
        }
    };

    GrmBuilder::Options grm_opts; grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe; grm_opts.algorithm = options_.grm_alg;
    grm_opts.thread_num = options_.n_threads;
    grm_opts.extract_set = options_.extract_set;
    grm_opts.exclude_set = options_.exclude_set;

    bool loco_enabled = !options_.no_loco;
    std::set<std::string> observed_chroms;

    if (is_bgen_) {
        const auto& variants = bgen_reader_->getVariants();
        for (const auto& v : variants) observed_chroms.insert(v.chrom);
    } else if (is_pgen_) {
        const auto& variants = pgen_reader_->getVariants();
        for (const auto& v : variants) observed_chroms.insert(v.chrom);
    } else {
        const auto& bim = plink_reader_->getBimInfo();
        for (const auto& b : bim) observed_chroms.insert(b.chrom);
    }
    std::vector<std::string> chroms = sortChroms(observed_chroms);
    if (chroms.size() <= 1) loco_enabled = false;

    std::string out_file = options_.out_prefix + ".stcg.assoc.txt";
    std::string null_file = options_.out_prefix + ".stcg.null.summary.txt";
    std::ofstream ofs(out_file);
    ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tLOG10P\n";

    if (!loco_enabled) {
        std::cout << "  Running global GRM scan (LOCO disabled or single chromosome).\n";
        scan_stcg_with_grm(G_subject_, options_.scan_chrom, ofs);
        writeNullSummary(null_file, "STCG",
            {{"V_G", std::to_string(v_g_)}, {"V_PE", std::to_string(v_pe_)},
             {"V_E", std::to_string(v_e_)}, {"n_records", std::to_string(n_obs)},
             {"n_subjects", std::to_string(n_subjects_keep)}, {"loco", "false"}});
    } else {
        std::cout << "  Running LOCO scan on " << chroms.size() << " chromosomes.\n";

        double S_global = 0.0;
        bool use_single_pass = options_.grm_cache_prefix.empty() && options_.grm_alg != 2;
        GrmBuilder::SplitGrmResult split_res;
        std::vector<Eigen::MatrixXf> per_chrom_raw;
        std::vector<double> chr_scaling;
        DataFilter::KeepInfo ki; ki.keep_indices = keep_subject_indices; ki.n_keep = keep_subject_indices.size();

        if (use_single_pass) {
            GrmBuilder split_builder(grm_opts);
            if (is_bgen_) {
                bgen_reader_->reset();
                split_res = split_builder.computeSplit(*bgen_reader_, observed_chroms, keep_subject_indices);
            } else if (is_pgen_) {
                pgen_reader_->reset();
                split_res = split_builder.computeSplit(*pgen_reader_, observed_chroms, keep_subject_indices);
            } else {
                split_res = split_builder.computeSplit(*plink_reader_, observed_chroms);
            }
            S_global = (uses_per_snp_scaling(options_.grm_alg))
                ? (double)split_res.global_snp_count : split_res.global_sum_2pq;

            per_chrom_raw.assign(chroms.size(), Eigen::MatrixXf());
            chr_scaling.assign(chroms.size(), 0.0);
            for (size_t i = 0; i < chroms.size(); ++i) {
                auto it = std::find(split_res.chrom_labels.begin(), split_res.chrom_labels.end(), chroms[i]);
                if (it != split_res.chrom_labels.end()) {
                    size_t k = std::distance(split_res.chrom_labels.begin(), it);
                    per_chrom_raw[i] = (is_bgen_ || is_pgen_)
                        ? split_res.per_chrom_raw[k].cast<float>()
                        : DataFilter::subsetSquare(split_res.per_chrom_raw[k], ki).cast<float>();
                    chr_scaling[i] = (uses_per_snp_scaling(options_.grm_alg))
                        ? (double)split_res.per_chrom_snp_count[k] : split_res.per_chrom_sum_2pq[k];
                }
            }
        } else {
            S_global = (uses_per_snp_scaling(options_.grm_alg))
                ? (double)(is_bgen_ ? bgen_reader_->getNumVariants() : is_pgen_ ? pgen_reader_->getNumVariants() : plink_reader_->getNumSnps())
                : 1.0;

            per_chrom_raw.assign(chroms.size(), Eigen::MatrixXf());
            chr_scaling.assign(chroms.size(), 0.0);
            for (size_t i = 0; i < chroms.size(); ++i) {
                GrmBuilder chr_builder(grm_opts);
                GrmBuilder::Options lo = grm_opts; lo.chrom_filter = chroms[i];
                chr_builder = GrmBuilder(lo);
                if (is_bgen_) { bgen_reader_->reset(); chr_builder.compute(*bgen_reader_); }
                else if (is_pgen_) { pgen_reader_->reset(); chr_builder.compute(*pgen_reader_); }
                else { chr_builder.compute(*plink_reader_); }
                per_chrom_raw[i] = DataFilter::subsetSquare(chr_builder.getGrm(), ki).cast<float>();
                chr_scaling[i] = (uses_per_snp_scaling(options_.grm_alg))
                    ? (double)chr_builder.getTotalSnps() : chr_builder.getSum2pq();
            }
        }

        Eigen::MatrixXd G_global = (use_single_pass)
            ? ((is_bgen_ || is_pgen_)
                ? split_res.global
                : DataFilter::subsetSquare(split_res.global, ki))
            : G_subject_;

        for (size_t i = 0; i < chroms.size(); ++i) {
            const auto& chrom = chroms[i];
            if (!options_.scan_chrom.empty() && chrom != options_.scan_chrom) continue;
            std::cout << "  STCG LOCO chromosome " << chrom << "...\n";

            double S_loco = S_global - chr_scaling[i];
            Eigen::MatrixXd G_loco = G_global;
            if (chr_scaling[i] > 0.0 && S_loco > 1e-6) {
                G_loco = (G_global * S_global - per_chrom_raw[i].cast<double>()) / S_loco;
            }
            scan_stcg_with_grm(G_loco, chrom, ofs);
        }

        writeNullSummary(null_file, "STCG_LOCO",
            {{"V_G", std::to_string(v_g_)}, {"V_PE", std::to_string(v_pe_)},
             {"V_E", std::to_string(v_e_)}, {"n_records", std::to_string(n_obs)},
             {"n_subjects", std::to_string(n_subjects_keep)}, {"loco", "true"},
             {"n_chromosomes", std::to_string(chroms.size())}});
    }

    ofs.close();
    std::cout << "STCG repeatability GWAS results saved to " << out_file << "\n";
}

} // namespace cosmic
