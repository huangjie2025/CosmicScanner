// Additive-Dominance (a+d+A+D) GWAS scan orchestration.
//
// New --gwas-model additive_dominance path. Builds the additive (VanRaden,
// algorithm 1) and dominance (Vitezica, algorithm 4) GRMs over the same SNP
// set, fits an EigenVceAD null per chromosome under LOCO, and scans every SNP
// with the joint additive + dominance design implemented by LmmEngineAD.
// Requires PLINK BED hardcall input in this phase.
#include "gwas_scanner.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>
#include "data_filter.h"
#include "grm_builder.h"
#include "eigen_vce_ad.h"
#include "gwas_output_schema.h"
#include "validation.h"

namespace cosmic {

namespace {

double ad_elapsed_seconds(const std::chrono::steady_clock::time_point& start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

void ad_log_stage_time(const std::string& label,
                       const std::chrono::steady_clock::time_point& start) {
    std::cout << "  [time] " << label << ": " << std::fixed << std::setprecision(2)
              << ad_elapsed_seconds(start) << "s" << std::defaultfloat << std::endl;
}

std::vector<std::string> ad_chrom_labels(const PlinkReader& reader) {
    std::set<std::string> chroms;
    for (const auto& snp : reader.getBimInfo()) chroms.insert(snp.chrom);
    std::vector<std::string> out(chroms.begin(), chroms.end());
    std::sort(out.begin(), out.end(), [](const std::string& a, const std::string& b) {
        try {
            return std::stoi(a) < std::stoi(b);
        } catch (...) {
            return a < b;
        }
    });
    return out;
}

std::vector<std::string> ad_chrom_labels(const PgenReader& reader) {
    std::set<std::string> chroms;
    for (const auto& snp : reader.getVariants()) chroms.insert(snp.chrom);
    std::vector<std::string> out(chroms.begin(), chroms.end());
    std::sort(out.begin(), out.end(), [](const std::string& a, const std::string& b) {
        try {
            return std::stoi(a) < std::stoi(b);
        } catch (...) {
            return a < b;
        }
    });
    return out;
}

std::vector<std::string> ad_chrom_labels(const PlinkReader& bed, const PgenReader& pgen,
                                         bool is_pgen) {
    return is_pgen ? ad_chrom_labels(pgen) : ad_chrom_labels(bed);
}

} // namespace

void GwasScanner::runAdLocoScan(const Eigen::VectorXd& y,
                                const Eigen::MatrixXf& X,
                                PlinkReader& plink_reader,
                                PgenReader& pgen_reader,
                                bool is_pgen,
                                const DataFilter::KeepInfo& keep) {
    if (options.gwas_model_md != "additive_dominance") {
        throw std::runtime_error("runAdLocoScan requires --gwas-model additive_dominance");
    }
    if (!options.use_covariate_projection) {
        throw std::runtime_error("additive_dominance requires covariate projection; "
                                 "--no-covar-proj is not supported for this model");
    }
    if (!options.extract_file.empty() || !options.exclude_file.empty()) {
        // extract/exclude sets flow through to both GRM passes and the scan.
    }
    if (options.obs_weights.size() == keep.n_keep) {
        throw std::runtime_error("--gwas-model additive_dominance does not yet support observation weights");
    }

    {
        // Hardcall-only contract: this reader decodes the hardcall track only,
        // so a PGEN whose records also carry phase, dosage, or genotype
        // probabilities must be rejected rather than silently reduced.
        if (is_pgen) {
            const PgenRecordCapabilities caps = pgen_reader.getRecordCapabilities();
            if (!caps.hardcall_only) {
                throw std::runtime_error(
                    "PGEN dosage records are not supported for dominance/AD; provide "
                    "hardcall PGEN or BED. Detected " + caps.describe() + ". This model is "
                    "defined on exact AA/AB/BB hardcalls, and the PGEN reader does not "
                    "consume dosage, phase, or genotype-probability payloads, so the file "
                    "cannot be reduced to hardcalls without changing the statistical "
                    "definition.");
            }
        }
        std::ofstream meta(options.out_prefix + ".scan.meta.txt");
        meta << "Item\tValue\n"
             << "analysis_model\tunivariate_lmm_ad\n"
             << "requested_scan_tier\tfixed_null\n"
             << "resolved_scan_tier\tfixed_null\n"
             << "snp_vc_strategy\tfixed_null\n"
             << "GWAS_MODEL\t" << options.gwas_model_md << "\n"
             << "GENOTYPE_FORMAT\t" << (is_pgen ? "PLINK_PGEN_HARDCALL" : "PLINK_BED_HARDCALL") << "\n"
             << "PGEN_MODE\t" << (is_pgen ? "HARDCALL" : "NOT_APPLICABLE") << "\n"
             << "PGEN_HARDCALL_ONLY\t" << (is_pgen ? "true" : "NOT_APPLICABLE") << "\n"
             << "KA_GRM_ALGORITHM\t1_VANRADEN\n"
             << "KD_GRM_ALGORITHM\t4_VITEZICA_DOMINANCE\n"
             << "AD_NULL_SOLVER\tEIGENVCEAD_COORDINATE_DESCENT\n"
             << "HWE_APPLIED\t" << (options.min_hwe > 0.0 ? "true" : "false") << "\n"
             << "HWE_METHOD\t" << (options.min_hwe > 0.0 ? "EXACT_HARDCALL" : "NONE") << "\n"
             << "DOMINANCE_ENCODING\tCENTERED_HETEROZYGOTE_INDICATOR\n"
             << "JOINT_TEST\tP_ADJUSTED_TWO_COLUMN_SCORE_WITH_2DF_WALD\n";
    }

    const auto start_total = std::chrono::steady_clock::now();
    const std::vector<std::string> chrom_labels = ad_chrom_labels(plink_reader, pgen_reader, is_pgen);
    const long long n_total = is_pgen
        ? static_cast<long long>(pgen_reader.getNumSamples())
        : static_cast<long long>(plink_reader.getNumSamples());
    if (n_total == 0) {
        throw std::runtime_error(is_pgen
            ? "No samples loaded from PGEN/PSAM"
            : "No samples loaded from BED/bim/fam");
    }

    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options.min_maf;
    grm_opts.min_hwe = (float)options.min_hwe;
    grm_opts.extract_set = extract_ids_;
    grm_opts.exclude_set = exclude_ids_;
#ifdef _OPENMP
    grm_opts.thread_num = options.n_threads;
#endif

    const std::set<std::string> label_set(chrom_labels.begin(), chrom_labels.end());

    // Pass A: additive VanRaden GRM (+ per-chromosome raw blocks for LOCO)
    auto start_grm = std::chrono::steady_clock::now();
    GrmBuilder::Options grm_opts_a = grm_opts;
    grm_opts_a.algorithm = 1;
    GrmBuilder::SplitGrmResult split_a;
    std::cout << "Building additive GRM (algorithm 1, VanRaden"
              << (is_pgen ? ", PGEN hardcalls" : "") << ")...\n";
    {
        GrmBuilder builder(grm_opts_a);
        split_a = is_pgen ? builder.computeSplit(pgen_reader, label_set)
                          : builder.computeSplit(plink_reader, label_set);
    }
    ad_log_stage_time("additive GRM build", start_grm);
    if (split_a.global_sum_2pq <= 0.0 || split_a.global_snp_count == 0) {
        throw std::runtime_error("additive GRM constructed no valid SNPs");
    }
    const double S_A_global = split_a.global_sum_2pq;

    // Pass D: dominance Vitezica GRM (+ per-chromosome raw blocks for LOCO)
    start_grm = std::chrono::steady_clock::now();
    GrmBuilder::Options grm_opts_d = grm_opts;
    grm_opts_d.algorithm = 4;
    GrmBuilder::SplitGrmResult split_d;
    std::cout << "Building dominance GRM (algorithm 4, Vitezica"
              << (is_pgen ? ", PGEN hardcalls" : "") << ")...\n";
    {
        GrmBuilder builder(grm_opts_d);
        split_d = is_pgen ? builder.computeSplit(pgen_reader, label_set)
                          : builder.computeSplit(plink_reader, label_set);
    }
    ad_log_stage_time("dominance GRM build", start_grm);
    if (split_d.global_sum_2pq <= 0.0 || split_d.global_snp_count == 0) {
        throw std::runtime_error("dominance GRM constructed no valid SNPs");
    }
    const double S_D_global = split_d.global_sum_2pq;
    if (split_d.chrom_labels != split_a.chrom_labels ||
        split_d.per_chrom_snp_count != split_a.per_chrom_snp_count) {
        throw std::runtime_error("additive/dominance GRM chromosome partitions disagree");
    }

    std::cout << "Global GRMs computed: S_A(g)= " << S_A_global
              << ", S_A(snps)= " << split_a.global_snp_count
              << "; S_D(g)= " << S_D_global
              << ", S_D(snps)= " << split_d.global_snp_count << std::endl;

    const Eigen::MatrixXd GA_global = DataFilter::subsetSquare(split_a.global, keep);
    const Eigen::MatrixXd GD_global = DataFilter::subsetSquare(split_d.global, keep);
    const Eigen::VectorXd y_sub = DataFilter::subset(y, keep);
    const Eigen::MatrixXf X_sub = DataFilter::subsetRows(X, keep);

    const bool use_loco = chrom_labels.size() > 1 && !options.no_loco;
    const std::string out_file = options.out_prefix + ".assoc.ad.txt";

    EigenVceAD::Options vce_opts;
    vce_opts.use_covariate_projection = options.use_covariate_projection;
    vce_opts.max_em_iter = options.reml_maxit;
    vce_opts.tol = options.reml_tol;
    vce_opts.verbose = false;

    LmmEngineAD::Options engine_opts;
    engine_opts.min_maf = options.min_maf;
    engine_opts.min_hwe = options.min_hwe;
    engine_opts.extract_set = extract_ids_;
    engine_opts.exclude_set = exclude_ids_;
    engine_opts.scan_from_bp = options.scan_from_bp;
    engine_opts.scan_to_bp = options.scan_to_bp;
#ifdef _OPENMP
    engine_opts.thread_num = options.n_threads;
#endif

    std::vector<bool> ad_null_converged;
    std::vector<int> ad_null_iterations;
    std::vector<std::string> ad_null_labels;
    std::vector<double> ad_sigma_a2;
    std::vector<double> ad_sigma_d2;
    std::vector<double> ad_sigma_e2;

    const auto run_ad_null_and_scan = [&](const Eigen::MatrixXd& K_A,
                                          const Eigen::MatrixXd& K_D,
                                          const std::string& chrom_filter) {
        auto null_start = std::chrono::steady_clock::now();
        EigenVceAD vce(vce_opts);
        vce.prepare(y_sub, X_sub, K_A, K_D);
        if (!vce.runNullModel()) {
            throw std::runtime_error("EigenVceAD null model failed to converge");
        }
        const double s2a = vce.getSigma2a();
        const double s2d = vce.getSigma2d();
        const double s2e = vce.getSigma2e();
        std::cout << "  [AD null" << (chrom_filter.empty() ? "" : (" chr " + chrom_filter))
                  << "] sigma_A^2=" << s2a
                  << ", sigma_D^2=" << s2d
                  << ", sigma_e^2=" << s2e
                  << ", logL=" << vce.getLogLikelihood() << std::endl;
        const bool ad_conv = vce.converged();
        const int ad_conv_iter = vce.lastIterations();
        ad_null_converged.push_back(ad_conv);
        ad_null_iterations.push_back(ad_conv_iter);
        ad_null_labels.push_back(chrom_filter.empty() ? "global" : chrom_filter);
        ad_sigma_a2.push_back(s2a);
        ad_sigma_d2.push_back(s2d);
        ad_sigma_e2.push_back(s2e);
        if (!ad_conv) {
            const std::string tag = chrom_filter.empty() ? "" : (" chr " + chrom_filter);
            std::cerr << "[warn] AD null model" << tag << " did NOT converge"
                      << " (maxit=" << vce_opts.max_em_iter << ", completed " << ad_conv_iter
                      << " iteration(s)); variance components are initial/final values "
                      << "and may be unreliable (see meta AD_NULL_*)\n";
        }
        ValidationReport runtime_validation;
        validate_variance_components(runtime_validation,
                                     {s2a, s2d, s2e},
                                     "EigenVceAD null model",
                                     ValidationStage::Runtime);
        runtime_validation.throw_if_fail();
        validation_report_.merge(runtime_validation);
        ad_log_stage_time(chrom_filter.empty() ? "AD null model" : ("AD null model chr " + chrom_filter),
                          null_start);

        auto scan_start = std::chrono::steady_clock::now();
        LmmEngineAD engine(engine_opts);
        engine.prepareFromVce(vce);
        std::vector<AdResult> results;
        if (is_pgen) engine.scan(pgen_reader, results, keep.keep_indices, chrom_filter);
        else engine.scan(plink_reader, results, keep.keep_indices, chrom_filter);
        ad_log_stage_time(chrom_filter.empty() ? "AD SNP scan" : ("AD SNP scan chr " + chrom_filter),
                          scan_start);
        appendAdResults(results);
    };

    if (!use_loco) {
        std::cout << (options.no_loco ? "LOCO disabled via flag. " : "Only 1 chromosome found. ")
                  << "Using AD LMM on global GRMs (no LOCO)...\n";
        run_ad_null_and_scan(GA_global, GD_global, options.scan_chr);
    } else {
        std::cout << "Running AD LOCO scan on " << split_a.chrom_labels.size()
                  << " chromosomes...\n";
        const int n_sub = (int)keep.n_keep;
        for (size_t c = 0; c < split_a.chrom_labels.size(); ++c) {
            const std::string& chrom = split_a.chrom_labels[c];
            if (!options.scan_chr.empty() && chrom != options.scan_chr) continue;

            Eigen::MatrixXd K_A_loco = GA_global;
            Eigen::MatrixXd K_D_loco = GD_global;
            const double S_A_chr = split_a.per_chrom_sum_2pq[c];
            const double S_D_chr = split_d.per_chrom_sum_2pq[c];
            const double S_A_loco = S_A_global - S_A_chr;
            const double S_D_loco = S_D_global - S_D_chr;

            if (S_A_chr > 0.0 && S_D_chr > 0.0 && S_A_loco > 1e-6 && S_D_loco > 1e-6) {
                Eigen::MatrixXd rawA_sub(n_sub, n_sub);
                Eigen::MatrixXd rawD_sub(n_sub, n_sub);
                for (int i = 0; i < n_sub; ++i) {
                    for (int j = 0; j < n_sub; ++j) {
                        rawA_sub(i, j) = split_a.per_chrom_raw[c](keep.keep_indices[i], keep.keep_indices[j]);
                        rawD_sub(i, j) = split_d.per_chrom_raw[c](keep.keep_indices[i], keep.keep_indices[j]);
                    }
                }
                K_A_loco = GA_global * (S_A_global / S_A_loco) - rawA_sub / S_A_loco;
                K_D_loco = GD_global * (S_D_global / S_D_loco) - rawD_sub / S_D_loco;
            }
            std::cout << "Processing Chromosome " << chrom << " (S_A_loco=" << S_A_loco
                      << ", S_D_loco=" << S_D_loco << ")...\n";
            run_ad_null_and_scan(K_A_loco, K_D_loco, chrom);
        }
    }

    if (!ad_null_converged.empty()) {
        std::ofstream meta_ad(options.out_prefix + ".scan.meta.txt", std::ios::app);
        if (meta_ad.is_open()) {
            const bool all_conv = std::all_of(ad_null_converged.begin(), ad_null_converged.end(),
                                              [](bool c) { return c; });
            const int max_iter = *std::max_element(ad_null_iterations.begin(), ad_null_iterations.end());
            meta_ad << "AD_NULL_MODELS\t" << ad_null_converged.size() << "\n"
                    << "AD_NULL_CONVERGED_ALL\t" << (all_conv ? "true" : "false") << "\n"
                    << "AD_NULL_MAX_ITERATIONS\t" << max_iter << "\n"
                    << "AD_NULL_MAX_ITERATIONS_REACHED\t"
                    << (!all_conv ? "true" : "false") << "\n";
            for (size_t i = 0; i < ad_null_labels.size(); ++i) {
                meta_ad << "AD_NULL_" << ad_null_labels[i] << "_CONVERGED\t"
                        << (ad_null_converged[i] ? "true" : "false") << "\n"
                        << "AD_NULL_" << ad_null_labels[i] << "_ITERATIONS\t"
                        << ad_null_iterations[i] << "\n";
            }
        }
    }
    if (ad_null_converged.empty()) {
        throw std::runtime_error("No chromosome matched the requested additive-dominance scan filters");
    }
    const bool all_converged = std::all_of(ad_null_converged.begin(), ad_null_converged.end(),
                                           [](bool value) { return value; });
    const std::string metadata_file = options.out_prefix + ".scan.meta.txt";
    gwas_output::appendCommonMetadata(metadata_file, "additive_dominance", keep.n_keep, use_loco,
                                      options.use_covariate_projection, all_converged);
    gwas_output::appendVarianceMetadata(metadata_file, "SIGMA_A2", ad_null_labels, ad_sigma_a2);
    gwas_output::appendVarianceMetadata(metadata_file, "SIGMA_D2", ad_null_labels, ad_sigma_d2);
    gwas_output::appendVarianceMetadata(metadata_file, "SIGMA_E2", ad_null_labels, ad_sigma_e2);

    flushAdResults(out_file);
    std::cout << "AD results saved to " << out_file << std::endl;
    std::cout << "  [AD summary] tested=" << ad_tested_markers_
              << ", rank-deficient/excluded=" << ad_rank_deficient_markers_ << std::endl;
    ad_log_stage_time("AD GWAS total", start_total);
}

void GwasScanner::appendAdResults(const std::vector<AdResult>& results) {
    ValidationReport report;
    long long tested = 0;
    long long rank_def = 0;
    for (const auto& res : results) {
        if (std::isnan(res.p_joint)) continue;
        ++tested;
        validate_probability_value(report, res.p_add, "AD_INVALID_P_ADD", "AD association output");
        validate_probability_value(report, res.p_dom, "AD_INVALID_P_DOM", "AD association output");
        validate_probability_value(report, res.p_joint, "AD_INVALID_P_JOINT", "AD association output");
        validate_finite_scalar(report, res.beta_add, "AD_NONFINITE_BETA_ADD", "AD association output");
        validate_finite_scalar(report, res.beta_dom, "AD_NONFINITE_BETA_DOM", "AD association output");
        validate_nonnegative_finite_scalar(report, res.se_add, "AD_INVALID_SE_ADD", "AD association output");
        validate_nonnegative_finite_scalar(report, res.se_dom, "AD_INVALID_SE_DOM", "AD association output");
        if (report.has_fail()) break;
        if (res.rank_deficient) ++rank_def;
    }
    if (!report.has_fail()) {
        report.pass(ValidationStage::PostRun, "AD_ASSOC_RESULTS_VALID",
                    "AD association output",
                    std::to_string(tested) + " tested marker results passed numerical sanity checks");
    }
    std::cout << report.summary(false);
    report.throw_if_fail();
    validation_report_.merge(report);
    ad_tested_markers_ += tested;
    ad_rank_deficient_markers_ += rank_def;
    ad_results_.insert(ad_results_.end(), results.begin(), results.end());
}

void GwasScanner::flushAdResults(const std::string& filename) {
    if (ad_results_.empty()) return;
    std::ofstream ofs(filename);
    if (!ofs) throw std::runtime_error("Cannot open AD output file: " + filename);
    ofs << std::setprecision(9);
    ofs << gwas_output::kUnifiedHeader;
    for (const auto& res : ad_results_) {
        if (std::isnan(res.p_joint)) continue;
        gwas_output::writeAdRow(ofs, res);
    }
    ad_results_.clear();
}

} // namespace cosmic
