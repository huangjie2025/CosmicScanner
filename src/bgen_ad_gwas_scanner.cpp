#include "gwas_scanner.h"
#include "bgen_ad_integration.h"
#include "bgen_ad_scan_engine.h"
#include "bgen_dominance_encoding.h"
#include "gwas_output_schema.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>

namespace cosmic {
std::vector<std::string> sort_chrom_labels(const std::set<std::string>&);

void GwasScanner::runBgenAdResearchScan(const Eigen::VectorXd& y, const Eigen::MatrixXf& X,
                                       BgenReader& reader, const DataFilter::KeepInfo& keep,
                                       const BgenAdResearchOptions& research,
                                       BgenAdIntegrationDiagnostics* diagnostics) {
    const std::string prefix = options.out_prefix;
    const std::string association = prefix + ".assoc.ad.txt";
    const std::string metadata = prefix + ".scan.meta.txt";
    if (std::filesystem::exists(association) || std::filesystem::exists(metadata)) {
        throw std::runtime_error("Internal BGEN AD integration requires a fresh output prefix");
    }
    std::ofstream meta(metadata);
    if (!meta) throw std::runtime_error("Cannot open internal BGEN AD metadata");
    meta << std::setprecision(17) << "Item\tValue\nMODEL\tadditive_dominance\n"
         << "GENOTYPE_BACKEND\tBGEN\nPROBABILITY_MODE\ttrue\n"
         << "INTERNAL_RESEARCH\ttrue\nPUBLIC_CLI_SUPPORTED\tfalse\n"
         << "QC_POLICY\tBGEN_DOM_V1_SHARED_AD_RESEARCH\n"
         << "QC_THRESHOLD_SOURCE\tE2.7_BGEN_DOM_V1\n"
         << "INFO_A_ACTION\tdiagnostic_only_no_new_cutoff\n"
         << "ENTROPY_ACTION\twarning_only\nHWE_ACTION\tretain_not_computed\n"
         << "AD_JOINT_TEST\tTWO_DF_FIXED_NULL_WALD\n"
         << "AD_CONDITIONAL_TESTS\tONE_DF_A_GIVEN_D_AND_D_GIVEN_A_WALD\n"
         << "KA_ENCODING\tPOSTERIOR_MEAN_ADDITIVE_OUTER_PRODUCT\nKA_NORMALIZATION\tSUM_2PQ\n"
         << "KD_ENCODING\tPOSTERIOR_MEAN_VITEZICA_OUTER_PRODUCT\nKD_NORMALIZATION\tSUM_4P2Q2\n"
         << "AD_MARKER_ENCODING\tPOSTERIOR_DOSAGE_AND_HETEROZYGOSITY_CENTERED_OBSERVED_MEAN\n"
         << "COUNTED_ALLELE_CONVENTION\tFIRST_STORED_BGEN_ALLELE_REPORTED_AS_A1\n"
         << "FREQ_SEMANTICS\tA1_allele_frequency_(not_MAF)\n"
         << "FREQUENCY_SAMPLE_SET\tANALYSIS_SAMPLES\nGRM_SAMPLE_SET\tANALYSIS_SAMPLES\n"
         << "POSTERIOR_VARIANCE_DIAGONAL_CORRECTION\tfalse\n"
         << "MISSING_POLICY\tobserved_mean_marker_imputation_zero_after_centering;zero_GRM_column\n"
         << "MAX_MISSING\t" << options.bgen_dom_max_missing
         << "\nMIN_INFO_D\t" << options.bgen_dom_min_info
         << "\nENTROPY_WARNING_ABOVE\t" << options.bgen_dom_entropy_warning
         << "\nMIN_MAF\t" << options.min_maf << "\n";
    std::string phase = "INPUT_GUARDS";
    bool null_summary_written = false;
    try {
        if (options.gwas_model_md != "additive_dominance" || options.snp_model != "additive")
            throw std::invalid_argument("Internal BGEN AD entry requires the additive_dominance model");
        if (reader.getCompressionType() == 0) throw std::invalid_argument("BGEN AD does not support uncompressed Layout-2 input");
        if (options.min_hwe != 0.0) throw std::invalid_argument("BGEN AD does not support exact HWE filtering");
        if (!std::isfinite(options.bgen_dom_entropy_warning) || options.bgen_dom_entropy_warning < 0.0 ||
            options.bgen_dom_entropy_warning > 1.0) throw std::invalid_argument("BGEN AD entropy warning threshold must be in [0,1]");
        if (!research.population_file.empty()) throw std::invalid_argument("BGEN AD does not support population-specific frequencies");
        if (!research.snp_weight_file.empty()) throw std::invalid_argument("BGEN AD does not support SNP weights");
        if (!options.use_covariate_projection || options.obs_weights.size() != 0 || options.use_spa ||
            options.biobank_mode || options.eb_lmm || options.run_skat || options.run_burden ||
            options.fast_lmm || options.fastlmm || options.lmm_low_rank > 0 ||
            options.lmm_gemma_refit || options.lmm_gemma_mode != 0 || options.profile_scale ||
            options.out_gz || options.out_ssf || options.manhattan_tsv || options.qq_tsv ||
            (!options.scan_tier.empty() && options.scan_tier != "approx")) {
            throw std::invalid_argument("Unsupported internal BGEN AD mode; standard fixed-null quantitative LMM only");
        }
        if (y.size() != reader.getNumSamples() || X.rows() != y.size() ||
            keep.keep_indices.size() != keep.n_keep || keep.n_keep == 0) {
            throw std::invalid_argument("BGEN AD analysis sample dimensions");
        }
        std::set<int> unique;
        for (int index : keep.keep_indices) {
            if (index < 0 || index >= y.size() || !unique.insert(index).second)
                throw std::invalid_argument("BGEN AD duplicate/out-of-range analysis sample");
        }
        const Eigen::VectorXd ys = DataFilter::subset(y, keep);
        const Eigen::MatrixXf xs = DataFilter::subsetRows(X, keep);
        if (!ys.allFinite() || !xs.allFinite()) throw std::invalid_argument("BGEN AD analysis samples contain non-finite phenotype/covariate");
        ad_results_.clear(); ad_tested_markers_ = 0; ad_rank_deficient_markers_ = 0;
        if (diagnostics) *diagnostics = BgenAdIntegrationDiagnostics{};
        meta << "N\t" << keep.n_keep
             << "\nN_SEMANTICS\tanalysis_samples_after_phenotype_and_covariate_filtering;per-SNP_missing_genotypes_mean-imputed\n"
             << "EXTRACT_FILE\t" << options.extract_file << "\nEXCLUDE_FILE\t" << options.exclude_file << '\n';

        BgenAdKernelOptions kernel_options;
        kernel_options.min_maf = options.min_maf;
        kernel_options.max_missing = options.bgen_dom_max_missing;
        kernel_options.min_info_d = options.bgen_dom_min_info;
        kernel_options.extract = extract_ids_; kernel_options.exclude = exclude_ids_;
        phase = "BUILD_SHARED_KERNELS";
        auto kernels = BgenAdNullBackend::build(reader, kernel_options, keep.keep_indices);
        const auto global = kernels.global();
        meta << "KA_MARKERS_GLOBAL\t" << global.markers << "\nKD_MARKERS_GLOBAL\t" << global.markers
             << "\nSA_GLOBAL\t" << global.denominator_a << "\nSD_GLOBAL\t" << global.denominator_d
             << "\nQC_CANDIDATES\t" << kernels.candidates
             << "\nQC_REJECTED_MISSING\t" << kernels.rejected_missing
             << "\nQC_REJECTED_MAF\t" << kernels.rejected_maf
             << "\nQC_REJECTED_INFO_D\t" << kernels.rejected_info << '\n';
        std::ofstream membership(prefix + ".ad.kernel.markers.txt");
        for (const auto& id : kernels.accepted_variant_ids) membership << id << '\n';
        std::unordered_set<std::string> accepted(kernels.accepted_variant_ids.begin(), kernels.accepted_variant_ids.end());
        // Diagnostic-only pass, using the same selected samples/probability API.
        long long entropy_warnings = 0, info_warnings = 0, missing_warnings = 0;
        std::ofstream qc(prefix + ".ad.qc.tsv");
        qc << std::setprecision(17) << "SNP\tCHR\tFREQ\tOBSERVED\tMISSING\tINFO_A\tINFO_D\tENTROPY\tGRM_ACCEPTED\n";
        reader.reset();
        for (uint32_t j = 0; j < reader.getNumVariants(); ++j) {
            std::string chromosome, id, a1, a2; uint32_t position; BgenGenotypeProbs probabilities;
            if (!reader.readVariantProbs(chromosome, id, position, a1, a2, probabilities)) throw std::runtime_error("Truncated BGEN QC stream");
            if ((!extract_ids_.empty() && !extract_ids_.count(id)) || exclude_ids_.count(id)) continue;
            const auto marker = BgenAdScanEngine::encode(probabilities, keep.keep_indices);
            const auto probability_qc = bgen::dominanceQcMetrics(probabilities.p0, probabilities.p1,
                probabilities.p2, probabilities.missing, keep.keep_indices);
            const bool retained = accepted.count(id) != 0;
            qc << id << '\t' << chromosome << '\t' << marker.frequency << '\t' << marker.observed
               << '\t' << marker.missing_rate << '\t' << marker.info_a << '\t' << marker.info_d
               << '\t' << probability_qc.mean_entropy << '\t' << (retained ? "true" : "false") << '\n';
            if (retained) {
                entropy_warnings += probability_qc.mean_entropy > options.bgen_dom_entropy_warning;
                info_warnings += marker.info_d < bgen::kDominanceInfoWarning;
                missing_warnings += marker.missing_rate > bgen::kDominanceMissingWarning;
            }
        }
        meta << "QC_ENTROPY_WARNINGS\t" << entropy_warnings << "\nQC_INFO_D_WARNINGS\t" << info_warnings
             << "\nQC_MISSING_WARNINGS\t" << missing_warnings << '\n';
        if (entropy_warnings || info_warnings || missing_warnings)
            std::cerr << "[warn] internal BGEN AD QC: entropy=" << entropy_warnings << ", INFO_D=" << info_warnings
                      << ", missingness=" << missing_warnings << " (warning-only details in .ad.qc.tsv)\n";
        std::set<std::string> chromosomes;
        for (const auto& entry : kernels.chromosome_raw) {
            chromosomes.insert(entry.first);
            meta << "KA_MARKERS_" << entry.first << '\t' << entry.second.markers
                 << "\nKD_MARKERS_" << entry.first << '\t' << entry.second.markers
                 << "\nSA_" << entry.first << '\t' << entry.second.denominator_a
                 << "\nSD_" << entry.first << '\t' << entry.second.denominator_d << '\n';
        }
        const bool loco = !options.no_loco;
        meta << "LOCO\t" << (loco ? "true" : "false") << "\nCOVARIATE_PROJECTION\ttrue\n";
        std::vector<std::string> labels;
        if (loco) {
            for (const auto& chr : sort_chrom_labels(chromosomes))
                if (options.scan_chr.empty() || chr == options.scan_chr) labels.push_back(chr);
        } else labels.push_back("global");
        if (labels.empty()) throw std::invalid_argument("No BGEN AD chromosome matched scan filter");
        std::ofstream certificates(prefix + ".ad.null.tsv"), candidates(prefix + ".ad.null.candidates.tsv");
        certificates << std::setprecision(17) << "NULL\tSIGMA_A2\tSIGMA_D2\tSIGMA_E2\tLOGL\tCONVERGED\tITERATIONS\tLIKELIHOOD_DELTA\tGRADIENT_NORM\tKKT_RESIDUAL\tBOUNDARY_STATUS\tINFORMATION_RANK\tD_LOCAL_SE\n";
        candidates << std::setprecision(17) << "NULL\tMODEL\tSTART\tSIGMA_A2\tSIGMA_D2\tSIGMA_E2\tLOGL\tFACE_CONVERGED\tKKT_VALID\n";
        std::vector<double> a_values, d_values, e_values;
        std::vector<std::string> boundaries;
        long long tested_total = 0, rank_total = 0, weak_count = 0;
        for (const auto& label : labels) {
            phase = "FIT_NULL_" + label;
            const auto pair = loco ? kernels.loco(label) : kernels.global();
            EigenVceAD::Options vce_options;
            vce_options.verbose = false; vce_options.max_em_iter = options.reml_maxit; vce_options.tol = options.reml_tol;
            const auto fit = BgenAdNullBackend::fit(pair, ys, xs, vce_options);
            certificates << label << '\t' << fit.sigma_a2 << '\t' << fit.sigma_d2 << '\t' << fit.sigma_e2
                << '\t' << fit.log_likelihood << '\t' << (fit.converged ? "true" : "false") << '\t' << fit.iterations
                << '\t' << fit.likelihood_delta << '\t' << fit.gradient_norm << '\t' << fit.kkt_residual
                << '\t' << fit.boundary_status << '\t' << fit.information_rank << '\t' << fit.sigma_d2_local_se << '\n';
            certificates.flush();
            for (const auto& candidate : fit.candidates)
                candidates << label << '\t' << candidate.model << '\t' << candidate.start << '\t' << candidate.sigma_a2
                    << '\t' << candidate.sigma_d2 << '\t' << candidate.sigma_e2 << '\t' << candidate.log_likelihood
                    << '\t' << (candidate.face_converged ? "true" : "false") << '\t' << (candidate.kkt_valid ? "true" : "false") << '\n';
            candidates.flush();
            meta << "NULL_" << label << "_CONVERGED\t" << (fit.converged ? "true" : "false")
                 << "\nNULL_" << label << "_GRADIENT_NORM\t" << fit.gradient_norm
                 << "\nNULL_" << label << "_KKT_RESIDUAL\t" << fit.kkt_residual
                 << "\nBOUNDARY_STATUS_" << label << '\t' << fit.boundary_status
                 << "\nWEAK_INFORMATION_DIAGNOSTICS_" << label << "\trank=" << fit.information_rank
                 << ";D_local_SE=" << fit.sigma_d2_local_se << ";D_point=" << fit.sigma_d2 << '\n';
            const bool weak = fit.information_rank < 3 || !std::isfinite(fit.sigma_d2_local_se) ||
                              (fit.sigma_d2 > 0 && fit.sigma_d2_local_se > fit.sigma_d2);
            weak_count += weak;
            if (weak) std::cerr << "[warn] BGEN AD " << label << " variance precision risk: rank=" << fit.information_rank
                << ", D_local_SE=" << fit.sigma_d2_local_se << ", D_point=" << fit.sigma_d2 << "; estimates not clipped\n";
            BgenAdScanEngine::Options scan_options;
            scan_options.min_maf = options.min_maf; scan_options.max_missing = options.bgen_dom_max_missing;
            scan_options.min_info_d = options.bgen_dom_min_info;
            scan_options.extract = extract_ids_; scan_options.exclude = exclude_ids_;
            scan_options.chromosome = loco ? label : options.scan_chr;
            scan_options.from_bp = options.scan_from_bp; scan_options.to_bp = options.scan_to_bp;
            BgenAdScanEngine engine(scan_options);
            // This call certifies and rotates the already fitted optimum; no
            // second null fit and no warm-start rotation reuse.
            engine.prepareCertified(kernels, ys, xs, fit, loco ? label : "");
            phase = "SCAN_" + label;
            const auto before = engine.nullFitCount();
            std::vector<AdResult> results; engine.scan(reader, results);
            if (engine.nullFitCount() != before) throw std::runtime_error("BGEN AD per-SNP null refit detected");
            for (const auto& result : results) if (!result.rank_deficient &&
                (!std::isfinite(result.log10_p_dom) || !std::isfinite(result.log10_p_joint)))
                throw std::runtime_error("BGEN AD non-finite LOG10P");
            appendAdResults(results);
            const auto stats = engine.statistics();
            tested_total += stats.tested; rank_total += stats.rank_deficient;
            a_values.push_back(fit.sigma_a2); d_values.push_back(fit.sigma_d2); e_values.push_back(fit.sigma_e2);
            boundaries.push_back(fit.boundary_status);
            meta << "NULL_BACKGROUND_MARKERS_" << label << '\t' << pair.markers
                 << "\nNULL_SA_" << label << '\t' << pair.denominator_a
                 << "\nNULL_SD_" << label << '\t' << pair.denominator_d
                 << "\nSCAN_TESTED_" << label << '\t' << stats.tested
                 << "\nSCAN_RANK_DEFICIENT_" << label << '\t' << stats.rank_deficient << '\n';
            if (diagnostics) diagnostics->nulls.push_back({label, fit, pair.markers, pair.denominator_a, pair.denominator_d,
                                                         stats.tested, stats.rank_deficient});
        }
        if (tested_total == 0) throw std::runtime_error("BGEN AD produced no testable markers");
        meta << "AD_NULL_MODELS\t" << labels.size() << "\nNULL_CONVERGED\ttrue\nVCE_CONVERGED\ttrue\n"
             << "BOUNDARY_STATUS\t" << (labels.size() == 1 ? boundaries.front() : "PER_NULL_MODEL")
             << "\nWEAK_INFORMATION_DIAGNOSTICS\tSEE_PER_NULL_RANK_AND_LOCAL_SE\nWEAK_INFORMATION_NULLS\t" << weak_count
             << "\nPER_SNP_NULL_REFITS\t0\nTESTED_MARKERS\t" << tested_total
             << "\nRANK_DEFICIENT_MARKERS\t" << rank_total << '\n';
        meta.flush();
        null_summary_written = true;
        meta.close();
        gwas_output::appendVarianceMetadata(metadata, "SIGMA_A2", labels, a_values);
        gwas_output::appendVarianceMetadata(metadata, "SIGMA_D2", labels, d_values);
        gwas_output::appendVarianceMetadata(metadata, "SIGMA_E2", labels, e_values);
        meta.open(metadata, std::ios::app);
        if (!meta) throw std::runtime_error("Cannot append BGEN AD final metadata");
        phase = "WRITE_UNIFIED_AD";
        flushAdResults(association);
        meta << "RUN_STATUS\tSUCCESS\n"; meta.flush();
        if (diagnostics) diagnostics->kernels = std::move(kernels);
    } catch (const std::exception& error) {
        ad_results_.clear();
        if (!meta.is_open()) { meta.clear(); meta.open(metadata, std::ios::app); }
        if (!null_summary_written) meta << "NULL_CONVERGED\tfalse\n";
        meta << "RUN_STATUS\tFAILED\nERROR_PHASE\t" << phase << "\nERROR_MESSAGE\t" << error.what()
             << "\nFAILED_NULL_INFERENCE_USABLE\tfalse\n";
        meta.flush();
        throw;
    }
}
}
