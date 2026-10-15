#include "gwas_scanner.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>

#include "eigen_vce_d.h"
#include "bgen_dominance_encoding.h"
#include "gwas_output_schema.h"
#include "validation.h"

namespace cosmic {

namespace {

std::vector<std::string> dominanceChromLabels(const PlinkReader& reader) {
    std::set<std::string> labels;
    for (const auto& snp : reader.getBimInfo()) labels.insert(snp.chrom);
    std::vector<std::string> result(labels.begin(), labels.end());
    std::sort(result.begin(), result.end(), [](const std::string& a, const std::string& b) {
        try {
            return std::stoi(a) < std::stoi(b);
        } catch (...) {
            return a < b;
        }
    });
    return result;
}

std::vector<std::string> dominanceChromLabels(const PgenReader& reader) {
    std::set<std::string> labels;
    for (const auto& snp : reader.getVariants()) labels.insert(snp.chrom);
    std::vector<std::string> result(labels.begin(), labels.end());
    std::sort(result.begin(), result.end(), [](const std::string& a, const std::string& b) {
        try {
            return std::stoi(a) < std::stoi(b);
        } catch (...) {
            return a < b;
        }
    });
    return result;
}

std::vector<std::string> dominanceChromLabels(BgenReader& reader) {
    reader.loadVariantInfo();
    std::set<std::string> labels;
    for (const auto& variant : reader.getVariants()) labels.insert(variant.chrom);
    std::vector<std::string> result(labels.begin(), labels.end());
    std::sort(result.begin(), result.end(), [](const std::string& a, const std::string& b) {
        try {
            return std::stoi(a) < std::stoi(b);
        } catch (...) {
            return a < b;
        }
    });
    return result;
}

double elapsedSeconds(const std::chrono::steady_clock::time_point& start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

} // namespace

void GwasScanner::runDominanceLocoScan(const Eigen::VectorXd& y,
                                       const Eigen::MatrixXf& X,
                                       PlinkReader& plink_reader,
                                       BgenReader& bgen_reader,
                                       PgenReader& pgen_reader,
                                       bool is_bgen,
                                       bool is_pgen,
                                       const DataFilter::KeepInfo& keep) {
    if (options.gwas_model_md != "dominance") {
        throw std::runtime_error("runDominanceLocoScan requires --gwas-model dominance");
    }
    if (!options.use_covariate_projection) {
        throw std::runtime_error("dominance requires covariate projection; --no-covar-proj is not supported");
    }
    if (options.obs_weights.size() == keep.n_keep) {
        throw std::runtime_error("--gwas-model dominance does not support observation weights");
    }
    if (is_bgen && bgen_reader.getCompressionType() == 0) {
        throw std::runtime_error(
            "BGEN dominance does not support uncompressed Layout-2 input; use zlib or zstd compression");
    }
    if (is_bgen && options.min_hwe > 0.0) {
        throw std::runtime_error(
            "BGEN probability dominance does not support exact --hwe filtering; use --hwe 0");
    }

    const auto total_start = std::chrono::steady_clock::now();

    // Hardcall-only contract: this reader decodes the hardcall track only, so a
    // PGEN whose records also carry phase, dosage, or genotype probabilities
    // must be rejected rather than silently reduced to its hardcalls.
    if (is_pgen) {
        const PgenRecordCapabilities caps = pgen_reader.getRecordCapabilities();
        if (!caps.hardcall_only) {
            throw std::runtime_error(
                "PGEN dosage records are not supported for dominance/AD; provide hardcall "
                "PGEN or BED. Detected " + caps.describe() + ". This model is defined on "
                "exact AA/AB/BB hardcalls, and the PGEN reader does not consume dosage, "
                "phase, or genotype-probability payloads, so the file cannot be reduced "
                "to hardcalls without changing the statistical definition.");
        }
    }

    const std::vector<std::string> chrom_labels = is_bgen
        ? dominanceChromLabels(bgen_reader)
        : (is_pgen ? dominanceChromLabels(pgen_reader) : dominanceChromLabels(plink_reader));
    const long long n_total = is_bgen
        ? static_cast<long long>(bgen_reader.getNumSamples())
        : (is_pgen ? static_cast<long long>(pgen_reader.getNumSamples())
                   : static_cast<long long>(plink_reader.getNumSamples()));
    if (n_total == 0) {
        throw std::runtime_error(is_bgen ? "No samples loaded from BGEN/sample"
            : (is_pgen ? "No samples loaded from PGEN/PSAM"
                       : "No samples loaded from BED/bim/fam"));
    }

    GrmBuilder::Options grm_opts;
    grm_opts.algorithm = 4;
    grm_opts.min_maf = static_cast<float>(options.min_maf);
    grm_opts.min_hwe = static_cast<float>(options.min_hwe);
    grm_opts.extract_set = extract_ids_;
    grm_opts.exclude_set = exclude_ids_;
    grm_opts.thread_num = options.n_threads;
    grm_opts.bgen_dom_max_missing = options.bgen_dom_max_missing;
    grm_opts.bgen_dom_min_info = options.bgen_dom_min_info;
    grm_opts.bgen_dom_entropy_warning = options.bgen_dom_entropy_warning;

    const std::set<std::string> label_set(chrom_labels.begin(), chrom_labels.end());
    std::cout << "Building dominance GRM (algorithm 4, Vitezica"
              << (is_bgen ? ", BGEN posterior means"
                           : (is_pgen ? ", PGEN hardcalls" : "")) << ")...\n";
    GrmBuilder builder(grm_opts);
    GrmBuilder::SplitGrmResult split = is_bgen
        ? builder.computeSplit(bgen_reader, label_set, keep.keep_indices)
        : (is_pgen ? builder.computeSplit(pgen_reader, label_set)
                   : builder.computeSplit(plink_reader, label_set));
    if (split.global_sum_2pq <= 0.0 || split.global_snp_count == 0) {
        throw std::runtime_error("dominance GRM constructed no valid SNPs");
    }
    if (is_bgen && (split.bgen_qc_missing_warnings > 0 ||
                   split.bgen_qc_info_warnings > 0 || split.bgen_qc_entropy_warnings > 0)) {
        std::cerr << "[warn] BGEN dominance QC " << bgen::kDominanceQcPolicyVersion
                  << ": missingness warnings=" << split.bgen_qc_missing_warnings
                  << ", INFO_D warnings=" << split.bgen_qc_info_warnings
                  << ", entropy warnings=" << split.bgen_qc_entropy_warnings
                  << "; thresholds and counts recorded in .scan.meta.txt\n";
    }

    const Eigen::MatrixXd K_D_global = is_bgen
        ? split.global
        : DataFilter::subsetSquare(split.global, keep);
    const Eigen::VectorXd y_sub = DataFilter::subset(y, keep);
    const Eigen::MatrixXf X_sub = DataFilter::subsetRows(X, keep);
    const double S_D_global = split.global_sum_2pq;
    const bool use_loco = split.chrom_labels.size() > 1 && !options.no_loco;

    {
        std::ofstream meta(options.out_prefix + ".scan.meta.txt");
        if (!meta) throw std::runtime_error("Cannot open dominance scan metadata output");
        meta << "Item\tValue\n"
             << "analysis_model\tunivariate_lmm_d\n"
             << "requested_scan_tier\tfixed_null\n"
             << "resolved_scan_tier\tfixed_null\n"
             << "snp_vc_strategy\tfixed_null\n"
             << "GWAS_MODEL\tdominance\n"
              << "GENOTYPE_FORMAT\t" << (is_bgen ? "BGEN_PROBABILITY"
                  : (is_pgen ? "PLINK_PGEN_HARDCALL" : "PLINK_BED_HARDCALL")) << "\n"
              << "PGEN_MODE\t" << (is_pgen ? "HARDCALL" : "NOT_APPLICABLE") << "\n"
              << "PGEN_HARDCALL_ONLY\t" << (is_pgen ? "true" : "NOT_APPLICABLE") << "\n"
              << "KD_GRM_ALGORITHM\t4_VITEZICA_DOMINANCE\n"
              << "KD_ENCODING\t" << (is_bgen ? "POSTERIOR_MEAN_VITEZICA_OUTER_PRODUCT"
                                                    : "HARDCALL_VITEZICA") << "\n"
              << "KD_NORMALIZATION\tSUM_4P2Q2\n"
              << "KD_DIAGONAL_UNCERTAINTY_CORRECTION\tfalse\n"
              << "D_NULL_SOLVER\tEIGENVCE_D_PROFILE_REML\n"
              << "HWE_APPLIED\t" << (options.min_hwe > 0.0 ? "true" : "false") << "\n"
              << "DOMINANCE_ENCODING\t" << (is_bgen
                    ? "posterior_mean_vitezica"
                    : "CENTERED_HETEROZYGOTE_INDICATOR") << "\n"
              << "COUNTED_ALLELE_CONVENTION\t" << (is_bgen
                    ? "FIRST_STORED_BGEN_ALLELE_REPORTED_AS_A1"
                    : "A1") << "\n"
              << "TEST\tONE_DF_FIXED_NULL_WALD\n"
              << "COVARIATE_COLUMNS\t" << X_sub.cols() << "\n"
              << "KD_MARKERS_GLOBAL\t" << split.global_snp_count << "\n"
              << std::setprecision(17)
              << "KD_DENOMINATOR_GLOBAL\t" << split.global_sum_2pq << "\n";
        if (is_bgen) {
            meta << "GENOTYPE_BACKEND\tBGEN\n"
                 << "PROBABILITY_MODE\ttrue\n"
                 << "DOMINANCE_MARKER_ENCODING\tposterior_heterozygosity_centered_sample_mean\n"
                 << "QC_POLICY_VERSION\t" << bgen::kDominanceQcPolicyVersion << "\n"
                 << "QC_SAMPLE_SET\tANALYSIS_SAMPLES\n"
                 << "FREQUENCY_SAMPLE_SET\tANALYSIS_SAMPLES\n"
                 << "INFO_D_DEFINITION\tVAR_OBSERVED_POSTERIOR_MEAN_VITEZICA_DIV_4P2Q2\n"
                 << "PROBABILITY_ENTROPY_DEFINITION\tMEAN_OBSERVED_NEG_SUM_P_LOG_P_DIV_LOG3\n"
                 << "BGEN_DOM_MAX_MISSING\t" << options.bgen_dom_max_missing << "\n"
                 << "BGEN_DOM_MISSING_WARNING_ABOVE\t" << bgen::kDominanceMissingWarning << "\n"
                 << "BGEN_DOM_MISSING_ACTION\tfilter_above_max_warn_above_warning\n"
                 << "BGEN_DOM_MIN_INFO_D\t" << options.bgen_dom_min_info << "\n"
                 << "BGEN_DOM_INFO_D_WARNING_BELOW\t" << bgen::kDominanceInfoWarning << "\n"
                 << "BGEN_DOM_INFO_D_ACTION\tfilter_below_min_warn_below_warning\n"
                 << "BGEN_DOM_ENTROPY_WARNING\t" << options.bgen_dom_entropy_warning << "\n"
                 << "BGEN_DOM_ENTROPY_ACTION\twarning_only\n"
                 << "BGEN_DOM_HWE_ACTION\tretain_not_computed\n"
                 << "BGEN_DOM_QC_CANDIDATES\t" << split.bgen_qc_candidates << "\n"
                 << "BGEN_DOM_QC_REJECTED_MISSING\t" << split.bgen_qc_rejected_missing << "\n"
                 << "BGEN_DOM_QC_MISSING_WARNINGS\t" << split.bgen_qc_missing_warnings << "\n"
                 << "BGEN_DOM_QC_REJECTED_MAF\t" << split.bgen_qc_rejected_maf << "\n"
                 << "BGEN_DOM_QC_REJECTED_INFO_D\t" << split.bgen_qc_rejected_info << "\n"
                 << "BGEN_DOM_QC_INFO_D_WARNINGS\t" << split.bgen_qc_info_warnings << "\n"
                 << "BGEN_DOM_QC_ENTROPY_WARNINGS\t" << split.bgen_qc_entropy_warnings << "\n";
        }
        for (size_t c = 0; c < split.chrom_labels.size(); ++c) {
            meta << "KD_MARKERS_" << split.chrom_labels[c] << "\t"
                 << split.per_chrom_snp_count[c] << "\n"
                 << "KD_DENOMINATOR_" << split.chrom_labels[c] << "\t"
                 << split.per_chrom_sum_2pq[c] << "\n";
        }
    }

    EigenVceD::Options vce_opts;
    vce_opts.use_covariate_projection = true;
    vce_opts.max_em_iter = options.reml_maxit;
    vce_opts.tol = options.reml_tol;
    vce_opts.verbose = false;

    LmmEngineD::Options engine_opts;
    engine_opts.min_maf = options.min_maf;
    engine_opts.min_hwe = options.min_hwe;
    engine_opts.extract_set = extract_ids_;
    engine_opts.exclude_set = exclude_ids_;
    engine_opts.scan_from_bp = options.scan_from_bp;
    engine_opts.scan_to_bp = options.scan_to_bp;
    engine_opts.thread_num = options.n_threads;
    engine_opts.bgen_dom_max_missing = options.bgen_dom_max_missing;
    engine_opts.bgen_dom_min_info = options.bgen_dom_min_info;
    engine_opts.bgen_dom_entropy_warning = options.bgen_dom_entropy_warning;

    std::vector<std::string> null_labels;
    std::vector<bool> null_converged;
    std::vector<int> null_iterations;
    std::vector<double> null_sigma_d2;
    std::vector<double> null_sigma_e2;

    const auto run_null_and_scan = [&](const Eigen::MatrixXd& K_D,
                                       const std::string& chrom_filter) {
        EigenVceD vce(vce_opts);
        vce.prepare(y_sub, X_sub, K_D);
        if (!vce.runNullModel()) {
            throw std::runtime_error("EigenVceD null model failed");
        }
        ValidationReport report;
        validate_variance_components(report, {vce.getSigma2d(), vce.getSigma2e()},
                                     "EigenVceD null model", ValidationStage::Runtime);
        report.throw_if_fail();
        validation_report_.merge(report);

        null_labels.push_back(chrom_filter.empty() ? "global" : chrom_filter);
        null_converged.push_back(vce.converged());
        null_iterations.push_back(vce.lastIterations());
        null_sigma_d2.push_back(vce.getSigma2d());
        null_sigma_e2.push_back(vce.getSigma2e());
        if (!vce.converged()) {
            std::cerr << "[warn] dominance null model did not converge for "
                      << null_labels.back() << "\n";
        }

        std::cout << "  [D null " << null_labels.back() << "] sigma_D^2="
                  << vce.getSigma2d() << ", sigma_e^2=" << vce.getSigma2e()
                  << ", logL=" << vce.getLogLikelihood() << "\n";
        LmmEngineD engine(engine_opts);
        engine.prepareFromVce(vce);
        std::vector<DomResult> results;
        if (is_bgen) engine.scan(bgen_reader, results, keep.keep_indices, chrom_filter);
        else if (is_pgen) engine.scan(pgen_reader, results, keep.keep_indices, chrom_filter);
        else engine.scan(plink_reader, results, keep.keep_indices, chrom_filter);
        appendDomResults(results);
    };

    if (!use_loco) {
        run_null_and_scan(K_D_global, options.scan_chr);
    } else {
        for (size_t c = 0; c < split.chrom_labels.size(); ++c) {
            const std::string& chrom = split.chrom_labels[c];
            if (!options.scan_chr.empty() && chrom != options.scan_chr) continue;
            const double S_D_loco = S_D_global - split.per_chrom_sum_2pq[c];
            if (S_D_loco <= 1e-6) {
                throw std::runtime_error("dominance LOCO GRM has no valid SNPs outside chromosome " + chrom);
            }
            Eigen::MatrixXd raw_sub = is_bgen
                ? split.per_chrom_raw[c]
                : DataFilter::subsetSquare(split.per_chrom_raw[c], keep);
            Eigen::MatrixXd K_D_loco =
                K_D_global * (S_D_global / S_D_loco) - raw_sub / S_D_loco;
            run_null_and_scan(K_D_loco, chrom);
        }
    }
    if (null_converged.empty()) {
        throw std::runtime_error("No chromosome matched the requested dominance scan filters");
    }

    {
        std::ofstream meta(options.out_prefix + ".scan.meta.txt", std::ios::app);
        const bool all_converged = std::all_of(null_converged.begin(), null_converged.end(),
                                               [](bool value) { return value; });
        meta << "D_NULL_MODELS\t" << null_converged.size() << "\n"
             << "D_NULL_CONVERGED_ALL\t" << (all_converged ? "true" : "false") << "\n";
        for (size_t i = 0; i < null_labels.size(); ++i) {
            meta << "D_NULL_" << null_labels[i] << "_CONVERGED\t"
                 << (null_converged[i] ? "true" : "false") << "\n"
                 << "D_NULL_" << null_labels[i] << "_ITERATIONS\t"
                 << null_iterations[i] << "\n";
        }
    }
    const bool all_converged = std::all_of(null_converged.begin(), null_converged.end(),
                                           [](bool value) { return value; });
    const std::string metadata_file = options.out_prefix + ".scan.meta.txt";
    gwas_output::appendCommonMetadata(metadata_file, "dominance", keep.n_keep, use_loco,
                                      options.use_covariate_projection, all_converged);
    gwas_output::appendVarianceMetadata(metadata_file, "SIGMA_D2", null_labels, null_sigma_d2);
    gwas_output::appendVarianceMetadata(metadata_file, "SIGMA_E2", null_labels, null_sigma_e2);

    if (dom_tested_markers_ == 0) {
        throw std::runtime_error("Dominance scan produced no testable markers");
    }

    const std::string output = options.out_prefix + ".assoc.dom.txt";
    flushDomResults(output);
    std::cout << "Dominance results saved to " << output << "\n"
              << "  [D summary] tested=" << dom_tested_markers_
              << ", rank-deficient/excluded=" << dom_rank_deficient_markers_ << "\n"
              << "  [time] dominance GWAS total: " << std::fixed << std::setprecision(2)
              << elapsedSeconds(total_start) << "s\n";
}

void GwasScanner::appendDomResults(const std::vector<DomResult>& results) {
    ValidationReport report;
    long long tested = 0;
    long long rank_deficient = 0;
    for (const auto& result : results) {
        if (result.rank_deficient && !result.snp_id.empty()) ++rank_deficient;
        if (std::isnan(result.p_dom)) continue;
        ++tested;
        validate_probability_value(report, result.p_dom, "D_INVALID_P_DOM", "dominance association output");
        validate_finite_scalar(report, result.beta_dom, "D_NONFINITE_BETA_DOM", "dominance association output");
        validate_nonnegative_finite_scalar(report, result.se_dom, "D_INVALID_SE_DOM", "dominance association output");
        if (report.has_fail()) break;
    }
    report.throw_if_fail();
    validation_report_.merge(report);
    dom_tested_markers_ += tested;
    dom_rank_deficient_markers_ += rank_deficient;
    dom_results_.insert(dom_results_.end(), results.begin(), results.end());
}

void GwasScanner::flushDomResults(const std::string& filename) {
    std::ofstream output(filename);
    if (!output) throw std::runtime_error("Cannot open dominance output file: " + filename);
    output << std::setprecision(9)
           << gwas_output::kUnifiedHeader;
    for (const auto& result : dom_results_) {
        if (std::isnan(result.p_dom)) continue;
        gwas_output::writeDominanceRow(output, result);
    }
    dom_results_.clear();
}

} // namespace cosmic
