#include "repeatability_model.h"
#include "repeatability_helpers.h"

#include "matrix_io.h"
#include "qc_utils.h"
#include "vce.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace cosmic {

namespace {

double wald_pvalue_1df(double z_score) {
    return std::erfc(std::abs(z_score) / std::sqrt(2.0));
}

} // namespace

void RepeatabilityModel::scanRefitValidation() {
    if (is_bgen_ || is_pgen_) throw std::runtime_error("Repeatability refit validation currently supports PLINK input only.");
    const bool primary_exact_scan = options_.scan_tier == "exact";
    std::cout << (primary_exact_scan
        ? "\nRunning Repeatability GWAS exact per-SNP full VC refit...\n"
        : "\nRunning Repeatability GWAS validation (per-SNP refit; small-data gold standard)...\n");

    const auto& bim = plink_reader_->getBimInfo();
    int n_variants = (int)bim.size();
    int n_obs = n_recs_kept_;
    int n_subjects = n_samples_;
    int variants_to_run = (options_.refit_limit > 0) ? std::min(n_variants, options_.refit_limit) : n_variants;

    std::string out_file = options_.out_prefix + (primary_exact_scan ? ".assoc.txt" : ".refit.assoc.txt");
    std::string meta_file = options_.out_prefix + ".refit.meta.txt";
    writeAssocHeader(out_file, options_.out_ssf);

    { std::ofstream meta(meta_file); meta << "Item\tValue\nmethod\t"
        << (primary_exact_scan ? "repeatability_per_snp_full_vc_refit" : "repeatability_per_snp_refit_validation")
        << "\nvariant_reader\tplink\n"
        << "variants_total\t" << n_variants << "\nvariants_attempted\t" << variants_to_run
        << "\nn_records\t" << n_obs << "\nn_subjects\t" << n_subjects
        << "\ninit_Vg\t" << v_g_ << "\ninit_Vpe\t" << v_pe_ << "\ninit_Ve\t" << v_e_ << "\n"; }

    Eigen::VectorXf g_raw; int written = 0, skipped_no_info = 0, skipped_solver = 0, skipped_nonpositive = 0;
    for (int j = 0; j < variants_to_run; ++j) {
        if (!options_.extract_set.empty() && options_.extract_set.find(bim[j].id) == options_.extract_set.end()) continue;
        if (!options_.exclude_set.empty() && options_.exclude_set.find(bim[j].id) != options_.exclude_set.end()) continue;
        plink_reader_->getSnpGenotypes(j, g_raw);
        if (g_raw.size() == 0) { skipped_no_info++; continue; }
        double sum_g = 0.0; int valid_n = 0;
        for (int s = 0; s < n_subjects; ++s) { float val = g_raw[s]; if (!std::isnan(val)) { sum_g += val; valid_n++; } }
        if (valid_n == 0) { skipped_no_info++; continue; }
        double mean_dosage = sum_g / valid_n; double freq = sum_g / (2.0 * valid_n);
        double maf = std::min(freq, 1.0 - freq);
        if (maf < options_.min_maf || freq <= 1e-8 || freq >= 1.0 - 1e-8) { skipped_no_info++; continue; }
        if (options_.min_hwe > 0.0) {
            long long n_hom1 = 0, n_het = 0, n_hom2 = 0;
            for (int s = 0; s < n_subjects; ++s) {
                if (std::isnan(g_raw[s])) continue;
                int gt = (int)std::llround(g_raw[s]);
                if (gt == 0) ++n_hom2;
                else if (gt == 1) ++n_het;
                else if (gt == 2) ++n_hom1;
            }
            if (calculateHWE(n_hom1, n_het, n_hom2) < options_.min_hwe) { skipped_no_info++; continue; }
        }

        FixedDesignG fd_snp = fd_; const int snp_col = fd_.p; fd_snp.p = fd_.p + 1; fd_snp.rows.resize(fd_.rows.size());
        for (size_t i = 0; i < fd_.rows.size(); ++i) {
            fd_snp.rows[i] = fd_.rows[i]; int sample_idx = obs_to_sample_[i];
            double dosage = mean_dosage;
            if (sample_idx >= 0 && sample_idx < g_raw.size() && !std::isnan(g_raw[sample_idx])) dosage = g_raw[sample_idx];
            fd_snp.rows[i].push_back({snp_col, dosage});
        }

        VCEConfig vce_cfg; vce_cfg.verbose = false; vce_cfg.print_report = false;
        vce_cfg.max_iter = options_.reml_maxit; vce_cfg.tol = options_.reml_tol; vce_cfg.vce_mode = "em";
        vce_cfg.set_algorithm_from_string(vce_cfg.vce_mode);
        AI_REML refit_solver(gen_records_, fd_snp, components_, y_obs_, vce_cfg);
        refit_solver.setInitialVars(v_e_, {v_g_, v_pe_}); refit_solver.solve();
        if (!refit_solver.isConverged()) { skipped_solver++; continue; }
        const auto& vars = refit_solver.getVarsU();
        if (vars.size() < 2) { skipped_solver++; continue; }
        double vg = vars[0], vpe = vars[1], ve = refit_solver.getVarE();
        if (vg <= 0.0 || vpe < 0.0 || ve <= 0.0) { skipped_nonpositive++; continue; }

        Eigen::MatrixXd X_full(n_obs, fd_snp.p); X_full.leftCols(fd_.p) = X_obs_.cast<double>();
        for (int i = 0; i < n_obs; ++i) {
            int sample_idx = obs_to_sample_[i]; double dosage = mean_dosage;
            if (sample_idx >= 0 && sample_idx < g_raw.size() && !std::isnan(g_raw[sample_idx])) dosage = g_raw[sample_idx];
            X_full(i, snp_col) = dosage;
        }

        Eigen::MatrixXd V(n_obs, n_obs);
        for (int i = 0; i < n_obs; ++i) { int si = obs_to_sample_[i]; for (int k = 0; k <= i; ++k) { int sk = obs_to_sample_[k];
            double value = vg * G_subject_(si, sk); if (si == sk) value += vpe; if (i == k) value += ve; V(i, k) = value; V(k, i) = value; } }

        Eigen::LDLT<Eigen::MatrixXd> v_solver;
        if (!factorize_with_jitter(V, v_solver)) { skipped_solver++; continue; }
        Eigen::MatrixXd VinvX = v_solver.solve(X_full); Eigen::VectorXd Vinvy = v_solver.solve(y_obs_);
        Eigen::MatrixXd XtVinvX = X_full.transpose() * VinvX;
        Eigen::LDLT<Eigen::MatrixXd> x_solver;
        if (!factorize_with_jitter(XtVinvX, x_solver)) { skipped_solver++; continue; }
        Eigen::MatrixXd XtVinvX_inv = x_solver.solve(Eigen::MatrixXd::Identity(XtVinvX.rows(), XtVinvX.cols()));
        Eigen::VectorXd beta = XtVinvX_inv * (X_full.transpose() * Vinvy);
        double beta_j = beta(snp_col), var_beta = XtVinvX_inv(snp_col, snp_col);
        if (!(var_beta > 0.0)) { skipped_solver++; continue; }
        double se_j = std::sqrt(var_beta), z = beta_j / se_j, p = wald_pvalue_1df(z);

        GwasResult res; res.chrom = bim[j].chrom; res.snp_id = bim[j].id; res.pos = bim[j].bp_pos;
        res.a1 = bim[j].alt; res.a2 = bim[j].ref; res.freq = freq; res.n_obs = valid_n;
        res.beta = beta_j; res.se = se_j; res.p_value = p; res.log10_p = (p > 0.0) ? -std::log10(p) : 300.0;
        appendAssocResult(out_file, res, options_.out_ssf); written++;
    }

    { std::ofstream meta(meta_file, std::ios::app); meta << "variants_written\t" << written << "\nskipped_no_info\t" << skipped_no_info
        << "\nskipped_solver\t" << skipped_solver << "\nskipped_nonpositive\t" << skipped_nonpositive << "\n"; }
    std::cout << "Repeatability refit validation results saved to " << out_file << "\n";
}

} // namespace cosmic
