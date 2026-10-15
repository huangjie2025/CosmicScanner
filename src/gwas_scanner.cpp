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
#include "matrix_free_vce.h"
#include "kinship_operator.h"
#include "design.h"
#include "data_filter.h"
#include "lmm_engine.h"
#include "gwas_output_schema.h"
#include "mvlmm_engine.h"
#include "validation.h"
#include "association/sumstats/sumstats.h"

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
    VcfReader& vcf_reader,
    bool is_bgen,
    bool is_pgen,
    bool is_vcf,
    const DataFilter::KeepInfo& keep,
    SequentialGenotypeMatrix* generic_reader
) {
    const bool is_generic = generic_reader != nullptr;
    std::vector<bool> additive_null_converged;
    const std::string resolved_tier = options.scan_tier.empty() ? "approx" : options.scan_tier;
    const std::string vc_strategy = resolved_tier == "exact" ? "per_snp_full_vc_refit"
                                  : resolved_tier == "semi" ? "fixed_ratio_profile_scale"
                                                             : "fixed_null";
    {
        std::ofstream meta(options.out_prefix + ".scan.meta.txt");
        meta << "Item\tValue\n"
             << "analysis_model\tunivariate_lmm\n"
             << "requested_scan_tier\t" << (options.scan_tier.empty() ? "legacy_default" : options.scan_tier) << "\n"
             << "resolved_scan_tier\t" << resolved_tier << "\n"
             << "snp_vc_strategy\t" << vc_strategy << "\n"
             << "SNP_MODEL\t" << options.snp_model << "\n"
             << "HWE_APPLIED\t" << (options.min_hwe > 0.0 && !is_bgen ? "true" : "false") << "\n"
             << "HWE_METHOD\t" << (options.min_hwe > 0.0 && !is_bgen ? "EXACT_HARDCALL" : "NONE") << "\n"
             << "HWE_REASON_IF_SKIPPED\t" << (options.min_hwe <= 0.0 ? "FILTER_DISABLED" : is_bgen ? "CONTINUOUS_DOSAGE_WITHOUT_PROBABILITIES" : "NONE") << "\n"
             << "matrix_approximation\t" << (options.fastlmm ? "fastlmm_auto_trace" : options.lmm_low_rank > 0 ? "randomized_low_rank" : "full_eigen") << "\n"
             << "requested_rank\t" << (options.fastlmm ? options.fastlmm_max_rank : options.lmm_low_rank) << "\n"
             << "trace_threshold\t" << (options.fastlmm ? options.fastlmm_trace_threshold : 0.0) << "\n"
             << "RANK_SELECTION\t" << (options.fastlmm ? "AUTO_TRACE" : options.lmm_low_rank > 0 ? "FIXED" : "FULL") << "\n";
    }
    if (options.lmm_low_rank > 0 || options.fastlmm) {
        std::ofstream meta(options.out_prefix + ".low_rank.meta.tsv");
        meta << "SCOPE\tALGORITHM\tREQUESTED_RANK\tEFFECTIVE_RANK\tSKETCH_SIZE"
             << "\tRESIDUAL_EIGENVALUE\tRELATIVE_FROBENIUS_ERROR\tERROR_TOLERANCE"
             << "\tDECOMPOSITION_SECONDS\tESTIMATED_PEAK_BYTES"
              << "\tFULL_EIGEN_ESTIMATED_PEAK_BYTES\tSELECTED_RANK\tTOTAL_RANK"
              << "\tTRACE_FRACTION\tRANK_SELECTION\tFALLBACK\n";
    }
std::cout << "Single-trait scan tier: " << resolved_tier << " (" << vc_strategy << ")\n";
    // Matrix-free path: no N x N GRM is ever built.  All GRM algebra goes
    // through a PlinkKinshipOperator that already restricts samples and
    // variants, so LOCO is handled by excluding one chromosome's SNPs.
    const bool matrix_free_scan = options.fastlmm && options.fastlmm_matrix_free;
    if (matrix_free_scan) {
        if (is_bgen || is_pgen || is_vcf || is_generic) {
            throw std::runtime_error(
                "--fastlmm-matrix-free currently requires PLINK BED input");
        }
        if (options.grm_alg == 2) {
            throw std::runtime_error(
                "--make-grm-alg 2 (external GRM) is not compatible with --fastlmm-matrix-free");
        }
    }
    // Build GRM for PLINK BED and VCF formats
    if ((!is_bgen && !is_pgen && !is_vcf) || is_vcf || is_generic) {
        if (options.grm_alg != 2) {
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
            // Get chromosome labels from the appropriate reader
            if (is_generic) {
                generic_reader->reset();
                std::string chrom, id, ref, alt;
                uint32_t pos = 0;
                Eigen::VectorXf dosage;
                while (generic_reader->readVariant(chrom, id, pos, ref, alt, dosage)) chrom_labels.insert(chrom);
            } else if (is_vcf) {
                for (const auto& chrom : vcf_reader.getChromLabels()) chrom_labels.insert(chrom);
            } else if (is_pgen) {
                for (const auto& snp : pgen_reader.getVariants()) chrom_labels.insert(snp.chrom);
            } else if (is_bgen) {
                // BGEN doesn't have a simple bim-like interface, use PLINK reader as fallback
                for (const auto& snp : plink_reader.getBimInfo()) chrom_labels.insert(snp.chrom);
            } else {
                for (const auto& snp : plink_reader.getBimInfo()) chrom_labels.insert(snp.chrom);
            }

            const bool use_loco = chrom_labels.size() > 1 && !options.no_loco;

            Eigen::MatrixXd G_global;
            double S_global = 0.0;
            GrmBuilder::SplitGrmResult split;
            if (matrix_free_scan) {
                // No dense GRM is built; the matrix-free operator streams
                // Z-multiplies per LOCO chromosome instead.
                S_global = (double)plink_reader.getNumSnps();
            } else if (use_loco) {
                std::cout << "Calculating Global + per-chromosome GRMs...\n";
                auto grm_start = std::chrono::steady_clock::now();
                GrmBuilder split_builder(grm_opts);
                if (is_generic) {
                    split = split_builder.computeSplit(*generic_reader, chrom_labels);
                } else if (is_vcf) {
                    split = split_builder.computeSplit(vcf_reader, chrom_labels);
                } else if (is_pgen) {
                    split = split_builder.computeSplit(pgen_reader, chrom_labels);
                } else if (is_bgen) {
                    split = split_builder.computeSplit(bgen_reader, chrom_labels);
                } else {
                    split = split_builder.computeSplit(plink_reader, chrom_labels);
                }
                S_global = (options.grm_alg == 0 || options.grm_alg == 3)
                    ? (double)split.global_snp_count
                    : split.global_sum_2pq;
                G_global = split.global;
                log_stage_time("split GRM build", grm_start);
            } else {
                std::cout << "Calculating Global GRM...\n";
                auto grm_start = std::chrono::steady_clock::now();
                GrmBuilder global_builder(grm_opts);
                if (is_generic) {
                    global_builder.compute(*generic_reader);
                } else if (is_vcf) {
                    global_builder.compute(vcf_reader);
                } else if (is_pgen) {
                    global_builder.compute(pgen_reader);
                } else if (is_bgen) {
                    global_builder.compute(bgen_reader);
                } else {
                    global_builder.compute(plink_reader);
                }
                S_global = (options.grm_alg == 0 || options.grm_alg == 3)
                    ? (double)global_builder.getTotalSnps()
                    : global_builder.getSum2pq();
                G_global = global_builder.getGrm();
                log_stage_time("global GRM build", grm_start);
            }

        std::cout << (matrix_free_scan ? "Matrix-free mode: dense GRM not built (M="
                                       : "Global GRM computed. Scaling Factor = ")
                  << S_global << std::endl;

        std::string out_file = options.out_prefix + ".assoc.txt";
        // Header + data written in one pass by flushAssocResults

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
        lmm_opts.profile_scale = options.profile_scale;
        lmm_opts.gemma_refit = options.lmm_gemma_refit;
        lmm_opts.min_maf = options.min_maf;
        lmm_opts.min_hwe = options.min_hwe;
        lmm_opts.snp_model = options.snp_model;
        lmm_opts.use_covariate_projection = options.use_covariate_projection;
        lmm_opts.use_spa = options.use_spa;
        lmm_opts.extract_set = extract_ids_;
        lmm_opts.exclude_set = exclude_ids_;
        lmm_opts.scan_from_bp = options.scan_from_bp;
        lmm_opts.scan_to_bp = options.scan_to_bp;
        if (has_weights) lmm_opts.obs_weights = options.obs_weights;
        if (options.grm_alg == 0) lmm_opts.use_residuals = true;
        #ifdef _OPENMP
        lmm_opts.thread_num = options.n_threads;
        #endif

        const bool loco_active = use_loco;

        // Variant selection shared by the dense and matrix-free paths.
        std::vector<int> grm_variant_indices;   // empty => all SNPs
        std::vector<std::string> grm_variant_chrom;  // chromosome label per entry
        if (matrix_free_scan) {
            plink_reader.computeStats(1);
            const auto& bim = plink_reader.getBimInfo();
            const auto& inv_stds = plink_reader.getInvStds();
            const auto& means = plink_reader.getMeans();
            grm_variant_indices.reserve(bim.size());
            grm_variant_chrom.reserve(bim.size());
            for (std::size_t j = 0; j < bim.size(); ++j) {
                if (inv_stds[j] == 0.0) continue;                 // monomorphic / all-missing
                const double freq = means[j] / 2.0;
                const double maf = std::min(freq, 1.0 - freq);
                if (maf < options.min_maf || std::abs(maf) < 1e-8) continue;
                if (!extract_ids_.empty() && !extract_ids_.count(bim[j].id)) continue;
                if (!exclude_ids_.empty() && exclude_ids_.count(bim[j].id)) continue;
                grm_variant_indices.push_back(static_cast<int>(j));
                grm_variant_chrom.push_back(bim[j].chrom);
            }
            if (grm_variant_indices.empty()) {
                throw std::runtime_error("Matrix-free GRM selection produced no variants");
            }
        }

        auto run_scan_with_grm = [&](const Eigen::MatrixXd& G_sub, const std::string& chrom_filter) {
            auto null_start = std::chrono::steady_clock::now();
            Eigen::MatrixXd G_work;
            if (!matrix_free_scan) {
                G_work = G_sub;
                if (has_weights) {
                    G_work.array().colwise() *= options.obs_weights.array();
                    G_work.array().rowwise() *= options.obs_weights.transpose().array();
                }
            }

            LmmEngine engine(lmm_opts);
            bool low_rank_accepted = false;
            LmmEngine::LowRankDiagnostics low_rank_diagnostics;
            std::string fallback = "none";
            bool null_model_converged = false;
            std::unique_ptr<KinshipOperator> kinship;
            bool use_matrix_free = options.fastlmm && options.fastlmm_matrix_free;
            if (use_matrix_free) {
                // LOCO: drop the current chromosome's variants from the operator.
                std::vector<int> loco_variants;
                if (loco_active && !chrom_filter.empty()) {
                    loco_variants.reserve(grm_variant_indices.size());
                    for (std::size_t t = 0; t < grm_variant_indices.size(); ++t) {
                        if (grm_variant_chrom[t] != chrom_filter) {
                            loco_variants.push_back(grm_variant_indices[t]);
                        }
                    }
                } else {
                    loco_variants = grm_variant_indices;
                }
                if (loco_variants.empty()) {
                    throw std::runtime_error(
                        "Matrix-free LOCO left no GRM variants for chromosome " + chrom_filter);
                }
                kinship = createKinshipOperator(&plink_reader, keep.keep_indices, loco_variants);
            }
            const int system_dim = matrix_free_scan
                ? static_cast<int>(y_sub.size())
                : static_cast<int>(G_work.rows());
            const int requested_rank = options.fastlmm ? options.fastlmm_max_rank
                                                      : options.lmm_low_rank;
            const int effective_rank = std::max(1, std::min(requested_rank, system_dim - 1));
            if (options.lmm_low_rank > 0 || options.fastlmm) {
                low_rank_diagnostics.requested_rank = requested_rank;

                bool use_matrix_free = options.fastlmm && options.fastlmm_matrix_free;

                if (effective_rank > 0 &&
                    (!options.fastlmm || effective_rank >= options.fastlmm_min_rank)) {

                    if (options.fastlmm && options.fastlmm_matrix_free) {
                        if (!kinship) {
                            throw std::runtime_error("Matrix-free path requires a kinship operator");
                        }
                        low_rank_diagnostics = engine.prepareLowRankMatrixFree(
                            y_sub, X_sub, *kinship, effective_rank,
                            options.lmm_low_rank_oversample, options.lmm_low_rank_seed,
                            options.fastlmm_trace_threshold,
                            options.fastlmm_min_rank,
                            options.n_threads);
                    } else {
                        low_rank_diagnostics = engine.prepareLowRank(
                            y_sub, X_sub, G_work, effective_rank,
                            options.lmm_low_rank_oversample, options.lmm_low_rank_seed,
                            options.fastlmm ? options.fastlmm_trace_threshold : 0.0,
                            options.fastlmm_min_rank);
                    }

                    low_rank_accepted = low_rank_diagnostics.relative_frobenius_error <=
                                        options.lmm_low_rank_error_tol;
                    if (options.fastlmm && !low_rank_diagnostics.trace_threshold_met) {
                        if (options.fastlmm_no_fallback) {
                            std::cerr << "WARNING: Trace threshold not met, continuing with fixed rank="
                                     << low_rank_diagnostics.effective_rank << std::endl;
                        } else {
                            low_rank_accepted = false;
                            fallback = "trace_threshold";
                        }
                    } else if (!low_rank_accepted) {
                        fallback = "error_tolerance";
                    }
                } else {
                    fallback = "invalid_effective_rank";
                }
            }
            if (!low_rank_accepted) {
                if (options.fastlmm && options.fastlmm_matrix_free) {
                    MatrixFreeVce::Options vce_opts;
                    vce_opts.use_covariate_projection = options.use_covariate_projection;
                    vce_opts.null_solver = static_cast<MatrixFreeVce::NullSolver>(options.lmm_null_solver);
                    vce_opts.max_iter = options.reml_maxit;
                    vce_opts.tol = options.reml_tol;
                    vce_opts.rank = effective_rank;
                    vce_opts.oversample = options.lmm_low_rank_oversample;
                    vce_opts.seed = options.lmm_low_rank_seed;
                    vce_opts.trace_threshold = 0.0;
                    vce_opts.min_rank = options.fastlmm_min_rank;
                    vce_opts.n_threads = options.n_threads;
                    MatrixFreeVce vce(vce_opts);
                    vce.prepare(y_sub, X_sub, *kinship);
                    vce.runNullModel();
                    null_model_converged = vce.converged();
                    engine.prepareFromMatrixFreeVce(vce);
                } else {
                    EigenVce::Options vce_opts;
                    vce_opts.use_covariate_projection = options.use_covariate_projection;
                    vce_opts.null_solver = static_cast<EigenVce::NullSolver>(options.lmm_null_solver);
                    vce_opts.max_iter = options.reml_maxit;
                    vce_opts.tol = options.reml_tol;
                    EigenVce vce(vce_opts);
                    vce.prepare(y_sub, X_sub, G_work);
                    vce.runNullModel();
                    null_model_converged = vce.converged();
                    engine.prepareFromVce(vce);
                }
            } else {
                // Low-rank preparation throws unless its bounded REML optimization succeeds.
                null_model_converged = true;
            }
            additive_null_converged.push_back(null_model_converged);
            std::cout << "  Null Model VCs: V(G)=" << engine.getSigma2g()
                      << ", V(e)=" << engine.getSigma2e() << std::endl;
            ValidationReport runtime_validation;
            validate_variance_components(runtime_validation,
                                         {engine.getSigma2g(), engine.getSigma2e()},
                                         "CosmicScanner null model",
                                         ValidationStage::Runtime);
            std::cout << runtime_validation.summary(false);
            runtime_validation.throw_if_fail();
            validation_report_.merge(runtime_validation);
            if (options.lmm_low_rank > 0 || options.fastlmm) {
                std::ofstream meta(options.out_prefix + ".low_rank.meta.tsv", std::ios::app);
                meta << (chrom_filter.empty() ? "GLOBAL" : chrom_filter)
                     << "\trandomized_subspace_iteration_tail_mean_woodbury"
                      << "\t" << low_rank_diagnostics.requested_rank
                     << "\t" << low_rank_diagnostics.effective_rank
                     << "\t" << low_rank_diagnostics.sketch_size
                     << "\t" << std::setprecision(12) << low_rank_diagnostics.residual_eigenvalue
                     << "\t" << low_rank_diagnostics.relative_frobenius_error
                     << "\t" << options.lmm_low_rank_error_tol
                     << "\t" << low_rank_diagnostics.decomposition_seconds
                      << "\t" << low_rank_diagnostics.estimated_peak_bytes
                      << "\t" << (24LL * system_dim * system_dim)
                      << "\t" << low_rank_diagnostics.effective_rank
                      << "\t" << low_rank_diagnostics.total_rank
                      << "\t" << low_rank_diagnostics.trace_fraction
                      << "\t" << low_rank_diagnostics.rank_selection
                      << "\t" << fallback << "\n";
            }
            log_stage_time(chrom_filter.empty() ? "null model" : ("null model chr " + chrom_filter), null_start);

            auto scan_start = std::chrono::steady_clock::now();
            std::vector<GwasResult> results;
            if (is_generic) {
                if (options.eb_lmm) throw std::runtime_error("EB-LMM is not supported for HapMap/text input");
                engine.scan(*generic_reader, results, keep.keep_indices, chrom_filter);
            } else if (is_vcf) {
                if (options.eb_lmm) engine.scanEB(vcf_reader, results, keep.keep_indices, chrom_filter);
                else engine.scan(vcf_reader, results, keep.keep_indices, chrom_filter);
            } else if (is_pgen) {
                if (options.eb_lmm) engine.scanEB(pgen_reader, results, keep.keep_indices, chrom_filter);
                else engine.scan(pgen_reader, results, keep.keep_indices, chrom_filter);
            } else if (is_bgen) {
                if (options.eb_lmm) engine.scanEB(bgen_reader, results, keep.keep_indices, chrom_filter);
                else engine.scan(bgen_reader, results, keep.keep_indices, chrom_filter);
            } else {
                if (options.eb_lmm) engine.scanEB(plink_reader, results, keep.keep_indices, chrom_filter);
                else engine.scan(plink_reader, results, keep.keep_indices, chrom_filter);
            }
            appendResults(out_file, results);

            // Write canonical sumstats
            std::string sumstats_file = options.out_prefix + ".cosmic.sumstats";
            appendCosmicSumstats(sumstats_file, results, "LMM");

            log_stage_time(chrom_filter.empty() ? "SNP scan" : ("SNP scan chr " + chrom_filter), scan_start);
        };

        if (!use_loco) {
            std::cout << (options.no_loco ? "LOCO disabled via flag. " : "Only 1 chromosome found. ");
            std::cout << "Using standard LMM (no LOCO)...\n";
            run_scan_with_grm(matrix_free_scan ? Eigen::MatrixXd()
                                               : DataFilter::subsetSquare(G_global, keep),
                              options.scan_chr);
        } else if (matrix_free_scan) {
            std::cout << "Running matrix-free LOCO Scan on " << chrom_labels.size()
                      << " chromosomes (no dense GRM)...\n";
            for (const std::string& chrom : chrom_labels) {
                if (!options.scan_chr.empty() && chrom != options.scan_chr) continue;
                std::cout << "Processing Chromosome " << chrom << "...\n";
                run_scan_with_grm(Eigen::MatrixXd(), chrom);
            }
        } else {
            std::cout << "Running LOCO Scan on " << split.chrom_labels.size()
                      << " chromosomes using cached split GRMs...\n";
            const int n_sub = keep.n_keep;
            Eigen::MatrixXd G_global_sub = DataFilter::subsetSquare(G_global, keep);
            for (size_t c = 0; c < split.chrom_labels.size(); ++c) {
                const std::string& chrom = split.chrom_labels[c];
                if (!options.scan_chr.empty() && chrom != options.scan_chr) continue;
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

        flushAssocResults(out_file);
        if (additive_null_converged.empty()) {
            throw std::runtime_error("No chromosome matched the requested additive scan filters");
        }
        const bool all_converged = std::all_of(additive_null_converged.begin(),
                                               additive_null_converged.end(),
                                               [](bool value) { return value; });
        gwas_output::appendCommonMetadata(options.out_prefix + ".scan.meta.txt", "additive",
                                          keep.n_keep, use_loco,
                                          options.use_covariate_projection, all_converged);
        std::cout << "Results saved to " << out_file << std::endl;
        if (options.eb_lmm) {
            std::string summary_file = options.out_prefix + ".eb.summary.txt";
            writeEbSummary(summary_file);
            std::cout << "EB summary saved to " << summary_file << std::endl;
        }
        // Write Manhattan/QQ plot TSVs
        if ((options.manhattan_tsv || options.qq_tsv) && !plot_results_.empty()) {
            writePlotTsv(options.out_prefix, plot_results_);
            if (options.manhattan_tsv) std::cout << "Manhattan TSV saved to " << options.out_prefix << ".manhattan.tsv\n";
            if (options.qq_tsv) std::cout << "QQ TSV saved to " << options.out_prefix << ".qq.tsv\n";
            plot_results_.clear();
        }
        log_stage_time("PLINK/VCF/text optimized GWAS total", total_start);
        return;
    }
    }

    // BGEN/PGEN path (original)
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
    // Header + data written in one pass by flushAssocResults

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
    lmm_opts.profile_scale = options.profile_scale;
    lmm_opts.gemma_refit = options.lmm_gemma_refit;
    lmm_opts.min_maf = options.min_maf;
    lmm_opts.min_hwe = options.min_hwe;
    lmm_opts.snp_model = options.snp_model;
    lmm_opts.use_covariate_projection = options.use_covariate_projection;
    lmm_opts.use_spa = options.use_spa;
    lmm_opts.extract_set = extract_ids_;
    lmm_opts.exclude_set = exclude_ids_;
    lmm_opts.scan_from_bp = options.scan_from_bp;
    lmm_opts.scan_to_bp = options.scan_to_bp;
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
        additive_null_converged.push_back(vce.converged());

        std::cout << "  Null Model VCs: V(G)=" << vce.getSigma2g() << ", V(e)=" << vce.getSigma2e() << std::endl;
        ValidationReport runtime_validation;
        validate_variance_components(runtime_validation,
                                     {vce.getSigma2g(), vce.getSigma2e()},
                                     "CosmicScanner null model",
                                     ValidationStage::Runtime);
        std::cout << runtime_validation.summary(false);
        runtime_validation.throw_if_fail();
        validation_report_.merge(runtime_validation);

        LmmEngine engine(lmm_opts);
        engine.prepareFromVce(vce);

        std::vector<GwasResult> results;
        if (is_bgen) {
            bgen_reader.reset();
            if (options.eb_lmm) engine.scanEB(bgen_reader, results, keep.keep_indices, options.scan_chr);
            else engine.scan(bgen_reader, results, keep.keep_indices, options.scan_chr);
        } else if (is_pgen) {
            pgen_reader.reset();
            if (options.eb_lmm) engine.scanEB(pgen_reader, results, keep.keep_indices, options.scan_chr);
            else engine.scan(pgen_reader, results, keep.keep_indices, options.scan_chr);
        } else {
            if (options.eb_lmm) engine.scanEB(plink_reader, results, keep.keep_indices, options.scan_chr);
            else engine.scan(plink_reader, results, keep.keep_indices, options.scan_chr);
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
            additive_null_converged.push_back(vce.converged());

            std::cout << "  Null Model VCs: V(G)=" << vce.getSigma2g() << ", V(e)=" << vce.getSigma2e() << std::endl;
            ValidationReport runtime_validation;
            validate_variance_components(runtime_validation,
                                         {vce.getSigma2g(), vce.getSigma2e()},
                                         "CosmicScanner null model",
                                         ValidationStage::Runtime);
            std::cout << runtime_validation.summary(false);
            runtime_validation.throw_if_fail();
            validation_report_.merge(runtime_validation);

            // Run LMM for this chromosome
            LmmEngine engine(lmm_opts);
            engine.prepareFromVce(vce);

            // Scan
            std::vector<GwasResult> chrom_results;
            if (is_bgen) { bgen_reader.reset(); if (options.eb_lmm) engine.scanEB(bgen_reader, chrom_results, keep.keep_indices, chrom); else engine.scan(bgen_reader, chrom_results, keep.keep_indices, chrom); }
            else if (is_pgen) { pgen_reader.reset(); if (options.eb_lmm) engine.scanEB(pgen_reader, chrom_results, keep.keep_indices, chrom); else engine.scan(pgen_reader, chrom_results, keep.keep_indices, chrom); }
            else { if (options.eb_lmm) engine.scanEB(plink_reader, chrom_results, keep.keep_indices, chrom); else engine.scan(plink_reader, chrom_results, keep.keep_indices, chrom); }

            appendResults(out_file, chrom_results);

            // Write canonical sumstats
            std::string sumstats_file = options.out_prefix + ".cosmic.sumstats";
            appendCosmicSumstats(sumstats_file, chrom_results, "LMM");
        }
    }
    flushAssocResults(out_file);
    if (additive_null_converged.empty()) {
        throw std::runtime_error("No chromosome matched the requested additive scan filters");
    }
    const bool all_converged = std::all_of(additive_null_converged.begin(),
                                           additive_null_converged.end(),
                                           [](bool value) { return value; });
    gwas_output::appendCommonMetadata(options.out_prefix + ".scan.meta.txt", "additive",
                                      keep.n_keep, use_loco,
                                      options.use_covariate_projection, all_converged);
    std::cout << "Results saved to " << out_file << std::endl;
    if (options.eb_lmm) {
        std::string summary_file = options.out_prefix + ".eb.summary.txt";
        writeEbSummary(summary_file);
        std::cout << "EB summary saved to " << summary_file << std::endl;
    }
}


} // namespace cosmic
