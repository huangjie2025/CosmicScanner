#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <memory>
#include <omp.h>

#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"
#include "hapmap_reader.h"
#include "text_geno_reader.h"
#include "npy_reader.h"
#include "ped_map_reader.h"
#include "data_loader.h"
#include "models/repeatability_model.h"
#include "models/multivariate_model.h"
#include "models/rrm_model.h"
#include "models/univariate_model.h"
#include "models/maternal_model.h"
#include "models/single_step_model.h"
#include "models/animal_model.h"
#include "models/common_environment_model.h"
#include "models/threshold_model.h"
#include "models/count_model.h"
#include "models/ordinal_model.h"
#include "models/survival_model.h"
#include "models/gxe_model.h"
#include "models/gemma_bslmm_model.h"
#include "models/gemma_vc_model.h"
#include "input/input_preparer.h"
#include "gemma_calccor.h"
#include "gemma_predictor.h"
#include "gemma_summary_tools.h"
#include "gwas_app.h"
#include "analysis_dataset.h"
#include "pheno_reader.h"
#include "run_manifest.h"
#include "solver_planner.h"
#include "genetics_workflow.h"
#include "association/cojo/cojo_engine.h"
#include "association/meta/meta_engine.h"
#include "association/fine_mapping/fine_mapping.h"

extern "C" {
#if defined(__APPLE__)
    void openblas_set_num_threads(int num_threads) __attribute__((weak_import));
#else
    void openblas_set_num_threads(int num_threads) __attribute__((weak));
#endif
}

#ifndef COSMICSCANNER_VERSION
#define COSMICSCANNER_VERSION "unknown"
#endif

static std::vector<std::string> split_csv_list(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ',')) {
        size_t b = 0;
        while (b < token.size() && std::isspace(static_cast<unsigned char>(token[b]))) ++b;
        size_t e = token.size();
        while (e > b && std::isspace(static_cast<unsigned char>(token[e - 1]))) --e;
        if (e > b) out.push_back(token.substr(b, e - b));
    }
    return out;
}

static std::string join_csv_list(const std::vector<std::string>& vals) {
    std::string out;
    for (const auto& v : vals) {
        if (!out.empty()) out += ",";
        out += v;
    }
    return out;
}

static void collect_prepared_rrm_covariates(const std::string& prepared_dir,
                                            std::string& tcovar_names,
                                            std::string& tclass_names) {
    if (!tcovar_names.empty() || !tclass_names.empty()) return;
    std::ifstream in(prepared_dir + "/covariate_audit.tsv");
    if (!in) return;
    std::string line;
    std::getline(in, line);
    std::vector<std::string> q, c;
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        std::string name, level, type;
        std::getline(ss, name, '\t');
        std::getline(ss, level, '\t');
        std::getline(ss, type, '\t');
        if (level != "record") continue;
        if (type == "categorical") c.push_back(name);
        else if (type == "quantitative") q.push_back(name);
    }
    tcovar_names = join_csv_list(q);
    tclass_names = join_csv_list(c);
}

static std::vector<int> subject_multivariate_columns(const std::string& pheno_file) {
    std::ifstream in(pheno_file);
    std::string header;
    std::getline(in, header);
    std::stringstream ss(header);
    std::vector<std::string> cols;
    std::string tok;
    while (ss >> tok) cols.push_back(tok);
    std::vector<int> out;
    for (int i = 3; i <= (int)cols.size(); ++i) out.push_back(i - 2);
    return out;
}

static void set_blas_threads(int num_threads) {
    if (openblas_set_num_threads) {
        openblas_set_num_threads(num_threads);
    }
}

static std::string planner_model_label(const cosmic::GwasConfig& cfg) {
    if (cfg.logistic_glmm) return "scanner-logistic-glmm";
    if (!cfg.logistic_test.empty()) return "scanner-logistic";
    if (cfg.saige_style) return "scanner-saige-style";
    if (cfg.farmcpu) return "scanner-farmcpu";
    if (cfg.blink) return "scanner-blink";
    if (cfg.firth_logistic) return "scanner-firth-logistic";
    if (cfg.conditional_joint) return "scanner-conditional-joint";
    if (!cfg.condition_snp_file.empty()) return "scanner-conditional-gwas";
    if (cfg.single_step_model) return "scanner-ssgwas";
    if (cfg.repeatability_model) return "scanner-repeatability-gwas";
    if (cfg.rrm_model) return "scanner-rrm-gwas";
    if (cfg.multivariate) return "scanner-multivariate-gwas";
    if (cfg.animal_model) return "scanner-animal-gwas";
    if (cfg.maternal_model) return "scanner-maternal-gwas";
    if (cfg.common_env_model) return "scanner-common-environment-gwas";
    if (cfg.linear_model) return "scanner-linear-gwas";
    return "scanner-lmm-gwas";
}

static bool is_standard_lmm_contract_path(const cosmic::GwasConfig& cfg) {
    return !cfg.rrm_model &&
           !cfg.repeatability_model &&
           !cfg.multivariate &&
           !cfg.maternal_model &&
           !cfg.single_step_model &&
           !cfg.animal_model &&
           !cfg.common_env_model &&
           !cfg.threshold_model &&
           !cfg.count_model &&
           !cfg.ordinal_model &&
           !cfg.survival_model &&
           !cfg.gxe_model &&
           !cfg.linear_model &&
           !cfg.eb_lmm &&
           !cfg.biobank_mode &&
           !cfg.run_skat &&
           !cfg.run_burden &&
             !cfg.farmcpu && !cfg.blink && !cfg.firth_logistic && cfg.logistic_test.empty() &&
             !cfg.logistic_glmm && cfg.condition_snp_file.empty();
}

static cosmic::SolverPlan make_solver_plan(const cosmic::GwasConfig& cfg,
                                           int n_samples,
                                           int n_markers) {
    cosmic::SolverRequest request;
    request.model_type = planner_model_label(cfg);
    request.n_samples = n_samples > 0 ? static_cast<std::size_t>(n_samples) : 0;
    request.n_markers = n_markers > 0 ? static_cast<std::size_t>(n_markers) : 0;
    request.system_dimension = request.n_samples;
    request.threads = cfg.num_threads > 0 ? cfg.num_threads : 1;
    request.tolerance = cfg.single_step_model ? cfg.ssgwas_pcg_tol : cfg.reml_tol;
    request.diagnostic_only = true;
    if (cfg.memory_budget_specified) {
        request.memory_budget = cosmic::make_memory_budget(cfg.memory_budget_bytes,
                                                           cosmic::MemoryBudgetSource::UserSpecified);
    }

    if (cfg.single_step_model) {
        request.requested_solver_label = cfg.ssgwas_solver;
        request.requested_solver = cfg.ssgwas_solver == "ldlt"
            ? cosmic::PlannerSolverKind::SparseDirect
            : cosmic::planner_solver_from_string(cfg.ssgwas_solver);
        request.current_solver_label = cfg.ssgwas_solver;
        request.matrix_free_operator_available = cfg.ssgwas_solver == "pcg";
        request.apy_available = cfg.ssgwas_apy_core_size > 0;
        request.apy_core_size = cfg.ssgwas_apy_core_size > 0 ? static_cast<std::size_t>(cfg.ssgwas_apy_core_size) : 0;
        request.matrix_sparsity = cosmic::PlannerSparsity::Sparse;
        request.sparse_matrix_available = cfg.ssgwas_solver == "ldlt";
    } else if (cfg.repeatability_model) {
        request.requested_solver_label = cfg.repeatability_mode;
        request.requested_solver = cosmic::planner_solver_from_string(cfg.repeatability_mode);
        request.current_solver_label = cfg.repeatability_mode;
        request.stcg_available = cfg.repeatability_mode == "stcg";
        request.matrix_free_operator_available = request.stcg_available;
    } else {
        request.requested_solver_label = "auto";
        request.requested_solver = cosmic::PlannerSolverKind::Auto;
        request.current_solver_label = cfg.lmm_null_solver;
    }

    if (cfg.single_step_model) {
        // ssGWAS builds the sparse MME inside the model; nnz is unknown at this diagnostic entry.
    } else if (cfg.make_grm_sparse) {
        request.matrix_sparsity = cosmic::PlannerSparsity::Sparse;
        request.sparse_matrix_available = true;
    } else {
        request.matrix_sparsity = cosmic::PlannerSparsity::Dense;
    }

    return cosmic::plan_solver(request);
}

static void log_solver_planner_diagnostics(const cosmic::SolverPlan& plan) {
    std::cout << cosmic::format_solver_plan(plan);
}

static cosmic::JsonValue argv_manifest(const std::vector<std::string>& argv) {
    cosmic::JsonValue out = cosmic::JsonValue::array();
    for (const auto& item : argv) out.push(cosmic::JsonValue::string(item));
    return out;
}

static cosmic::JsonValue scanner_inputs_manifest(const cosmic::GwasConfig& cfg) {
    cosmic::JsonValue inputs = cosmic::JsonValue::object();
    if (!cfg.bfile_prefix.empty()) {
        inputs.set("bfile_prefix", cosmic::JsonValue::string(cfg.bfile_prefix));
        inputs.set("bed", cosmic::input_file_manifest(cfg.bfile_prefix + ".bed"));
        inputs.set("bim", cosmic::input_file_manifest(cfg.bfile_prefix + ".bim"));
        inputs.set("fam", cosmic::input_file_manifest(cfg.bfile_prefix + ".fam"));
    }
    if (!cfg.pfile_prefix.empty()) {
        inputs.set("pfile_prefix", cosmic::JsonValue::string(cfg.pfile_prefix));
        inputs.set("pgen", cosmic::input_file_manifest(cfg.pfile_prefix + ".pgen"));
        inputs.set("pvar", cosmic::input_file_manifest(cfg.pfile_prefix + ".pvar"));
        inputs.set("psam", cosmic::input_file_manifest(cfg.pfile_prefix + ".psam"));
    }
    if (!cfg.pgen_file.empty()) inputs.set("pgen", cosmic::input_file_manifest(cfg.pgen_file));
    if (!cfg.pvar_file.empty()) inputs.set("pvar", cosmic::input_file_manifest(cfg.pvar_file));
    if (!cfg.psam_file.empty()) inputs.set("psam", cosmic::input_file_manifest(cfg.psam_file));
    if (!cfg.bgen_file.empty()) inputs.set("bgen", cosmic::input_file_manifest(cfg.bgen_file));
    if (!cfg.sample_file.empty()) inputs.set("sample", cosmic::input_file_manifest(cfg.sample_file));
    if (!cfg.pheno_file.empty()) inputs.set("phenotype", cosmic::input_file_manifest(cfg.pheno_file));
    if (!cfg.qcovar_file.empty()) inputs.set("qcovar", cosmic::input_file_manifest(cfg.qcovar_file));
    if (!cfg.covar_file.empty()) inputs.set("covar", cosmic::input_file_manifest(cfg.covar_file));
    if (!cfg.ped_file.empty()) inputs.set("pedigree", cosmic::input_file_manifest(cfg.ped_file));
    return inputs;
}

static void write_scanner_manifest(const cosmic::GwasConfig& cfg,
                                   const cosmic::AnalysisDataset& dataset,
                                   const cosmic::SolverPlan& solver_plan,
                                   const cosmic::ValidationReport& validation,
                                   std::chrono::system_clock::time_point start_wall,
                                   std::chrono::system_clock::time_point end_wall,
                                   double wall_seconds,
                                   std::uint64_t peak_rss_bytes) {
    const std::string manifest_path = cfg.out_prefix + ".run_manifest.json";
    cosmic::JsonValue root = cosmic::JsonValue::object();
    root.set("schema_version", cosmic::JsonValue::string("1"));
    cosmic::JsonValue software = cosmic::JsonValue::object();
    software.set("component", cosmic::JsonValue::string("CosmicScanner"));
    software.set("component_version", cosmic::JsonValue::string(COSMICSCANNER_VERSION));
    software.set("panomics_version", cosmic::JsonValue::string("unknown"));
    software.set("git_commit", cosmic::JsonValue::string("unknown"));
    software.set("build_type", cosmic::JsonValue::string("unknown"));
    root.set("software", std::move(software));

    cosmic::JsonValue run = cosmic::JsonValue::object();
    run.set("status", cosmic::JsonValue::string("completed"));
    run.set("timestamp_start", cosmic::JsonValue::string(cosmic::iso8601_utc(start_wall)));
    run.set("timestamp_end", cosmic::JsonValue::string(cosmic::iso8601_utc(end_wall)));
    run.set("argv", argv_manifest(cfg.command_argv));
    run.set("threads", cosmic::JsonValue::integer(cfg.num_threads > 0 ? cfg.num_threads : 1));
    if (cfg.bslmm_seed != 0) run.set("random_seed", cosmic::JsonValue::integer(cfg.bslmm_seed));
    else run.set("random_seed", cosmic::JsonValue::null());
    root.set("run", std::move(run));

    root.set("inputs", scanner_inputs_manifest(cfg));
    root.set("data", cosmic::analysis_dataset_manifest(dataset));

    cosmic::JsonValue model = cosmic::JsonValue::object();
    model.set("type", cosmic::JsonValue::string("standard_single_trait_lmm"));
    model.set("trait", cosmic::JsonValue::string(cfg.input_trait.empty() ? "trait" : cfg.input_trait));
    model.set("scan_tier", cosmic::JsonValue::string(cfg.scan_tier.empty() ? "legacy_default" : cfg.scan_tier));
    model.set("lmm_null_solver", cosmic::JsonValue::string(cfg.lmm_null_solver));
    model.set("grm_algorithm", cosmic::JsonValue::integer(cfg.grm_alg));
    root.set("model", std::move(model));
    root.set("solver", cosmic::solver_plan_manifest(solver_plan));
    root.set("validation", cosmic::validation_report_manifest(validation));

    cosmic::JsonValue perf = cosmic::JsonValue::object();
    perf.set("wall_time_seconds", cosmic::JsonValue::number(wall_seconds));
    perf.set("peak_rss_bytes", peak_rss_bytes ? cosmic::JsonValue::uinteger(peak_rss_bytes) : cosmic::JsonValue::null());
    root.set("performance", std::move(perf));

    cosmic::JsonValue outputs = cosmic::JsonValue::array();
    outputs.push(cosmic::output_record_manifest("association", cfg.out_prefix + ".assoc.txt"));
    outputs.push(cosmic::output_record_manifest("scan_metadata", cfg.out_prefix + ".scan.meta.txt"));
    cosmic::JsonValue manifest_out = cosmic::JsonValue::object();
    manifest_out.set("type", cosmic::JsonValue::string("manifest"));
    manifest_out.set("path", cosmic::JsonValue::string(manifest_path));
    manifest_out.set("exists", cosmic::JsonValue::boolean(true));
    outputs.push(std::move(manifest_out));
    root.set("outputs", std::move(outputs));

    cosmic::write_json_atomic(manifest_path, root);
}

static void log_analysis_dataset_shadow(const cosmic::GwasConfig& cfg,
                                        const std::vector<std::string>& fids,
                                        const std::vector<std::string>& iids,
                                        const cosmic::PlinkReader& plink_reader,
                                        const cosmic::BgenReader& bgen_reader,
                                        const cosmic::PgenReader& pgen_reader,
                                        const cosmic::VcfReader& vcf_reader,
                                        bool is_bgen,
                                        bool is_pgen,
                                        bool is_vcf) {
    try {
        cosmic::AnalysisDataset dataset;
        std::vector<cosmic::SampleKey> genotype_keys;
        genotype_keys.reserve(iids.size());
        for (std::size_t i = 0; i < iids.size(); ++i) {
            cosmic::SampleKey key;
            key.fid = i < fids.size() ? fids[i] : std::string();
            key.iid = iids[i];
            genotype_keys.push_back(std::move(key));
        }

        cosmic::GenotypeSourceMetadata genotype_metadata;
        genotype_metadata.sample_order = genotype_keys;
        genotype_metadata.n_samples = genotype_keys.size();
        if (is_vcf) {
            genotype_metadata.format = cosmic::GenotypeStorageFormat::Vcf;
            genotype_metadata.source_label = cfg.vcf_file;
            genotype_metadata.n_markers = vcf_reader.getNumVariants();
        } else if (is_bgen) {
            genotype_metadata.format = cosmic::GenotypeStorageFormat::Bgen;
            genotype_metadata.source_label = cfg.bgen_file;
            genotype_metadata.n_markers = bgen_reader.getNumVariants();
        } else if (is_pgen) {
            genotype_metadata.format = cosmic::GenotypeStorageFormat::PlinkPgen;
            genotype_metadata.source_label = cfg.pgen_file.empty() ? cfg.pfile_prefix : cfg.pgen_file;
            genotype_metadata.n_markers = pgen_reader.getNumVariants();
        } else {
            genotype_metadata.format = cosmic::GenotypeStorageFormat::PlinkBed;
            genotype_metadata.source_label = cfg.bfile_prefix;
            genotype_metadata.n_markers = plink_reader.getNumSnps();
        }

        cosmic::GenotypeSourceCapabilities capabilities;
        capabilities.supports_streaming = true;
        capabilities.supports_random_access = !is_bgen;
        capabilities.supports_prefetch = is_bgen;
        capabilities.supports_sample_subset = false;
        capabilities.materializes_full_matrix = false;
        genotype_metadata.capabilities = capabilities;
        dataset.setGenotypeSource(std::make_shared<cosmic::MetadataOnlyGenotypeSource>(
            std::move(genotype_metadata)));

        std::vector<cosmic::SampleKey> phenotype_keys;
        cosmic::PhenotypeTable phenotype;
        bool have_shadow_phenotype = false;
        if (!cfg.pheno_file.empty() && !cfg.repeatability_model && !cfg.rrm_model && !cfg.rrm_unbalance) {
            cosmic::PhenoReader reader;
            reader.load(cfg.pheno_file, cfg.mpheno);
            Eigen::VectorXd y;
            reader.align(fids, iids, y);
            phenotype.row_samples = genotype_keys;
            phenotype.trait_names = {cfg.input_trait.empty() ? std::string("trait") : cfg.input_trait};
            phenotype.values.reserve(static_cast<std::size_t>(y.size()));
            phenotype.missing_mask.reserve(static_cast<std::size_t>(y.size()));
            for (int i = 0; i < y.size(); ++i) {
                phenotype.values.push_back(y(i));
                const bool missing = std::isnan(y(i));
                phenotype.missing_mask.push_back(missing);
                if (!missing && static_cast<std::size_t>(i) < genotype_keys.size()) {
                    phenotype_keys.push_back(genotype_keys[static_cast<std::size_t>(i)]);
                }
            }
            have_shadow_phenotype = true;
        }

        dataset.setSampleRegistry(cosmic::SampleRegistry::fromGenotypeOrder(genotype_keys, phenotype_keys));
        if (have_shadow_phenotype) dataset.setPhenotypeTable(std::move(phenotype));
        std::cout << dataset.diagnosticSummary();
    } catch (const std::exception& e) {
        std::cout << "AnalysisDataset shadow mode\n"
                  << "Status               : metadata unavailable (" << e.what() << ")\n";
    }
}

using namespace std;
using namespace cosmic;

#include "gwas_app.h"

namespace cosmic {

GwasApp::GwasApp(const GwasConfig& config) : cfg(config) {}

int GwasApp::run() {
    try {
        if ((cfg.skat_o || cfg.acat_v) && (cfg.linear_model || cfg.biobank_mode || cfg.multivariate ||
            cfg.repeatability_model || cfg.rrm_model || cfg.single_step_model || cfg.animal_model ||
            cfg.maternal_model || cfg.common_env_model || cfg.threshold_model || cfg.count_model ||
            cfg.ordinal_model || cfg.survival_model || cfg.gxe_model))
            throw std::runtime_error("New set tests require the standard single-trait Gaussian null route");
        const bool multilocus = cfg.farmcpu || cfg.blink || cfg.firth_logistic ||
                                !cfg.logistic_test.empty() || cfg.logistic_glmm ||
                                !cfg.condition_snp_file.empty();
        if (multilocus && (cfg.linear_model || cfg.biobank_mode || cfg.multivariate || cfg.repeatability_model ||
            cfg.rrm_model || cfg.single_step_model || cfg.animal_model || cfg.maternal_model ||
            cfg.common_env_model || cfg.threshold_model || cfg.count_model || cfg.ordinal_model ||
            cfg.survival_model || cfg.gxe_model || cfg.run_skat || cfg.run_burden || cfg.use_spa ||
            cfg.eb_lmm || cfg.fast_lmm || cfg.lmm_gemma_mode != 0 || cfg.lmm_gemma_refit ||
            !cfg.scan_tier.empty() || cfg.out_ssf))
            throw std::runtime_error("FarmCPU, BLINK, Firth logistic, and conditional GWAS are experimental single-trait routes and cannot be combined with another model/test tier");
        if (cfg.saige_style && (cfg.linear_model || cfg.multivariate || cfg.repeatability_model ||
            cfg.rrm_model || cfg.single_step_model || cfg.animal_model || cfg.maternal_model ||
            cfg.common_env_model || cfg.threshold_model || cfg.count_model || cfg.ordinal_model ||
            cfg.survival_model || cfg.gxe_model || cfg.run_skat || cfg.run_burden || cfg.eb_lmm ||
            cfg.fast_lmm || cfg.lmm_gemma_mode != 0 || cfg.lmm_gemma_refit ||
            !cfg.scan_tier.empty() || cfg.out_ssf || multilocus))
            throw std::runtime_error("--saige-style cannot be combined with another model/test tier");
        if (prepareGeneticsWorkflow(cfg)) return 0;
        
        // Handle format conversion
        if (!cfg.convert_from.empty() && !cfg.convert_to.empty()) {
            std::cout << "--- Genotype Format Conversion ---\n";
            std::cout << "Converting from " << cfg.convert_from << " to " << cfg.convert_to << "\n";
            
            PlinkReader plink_reader;
            BgenReader bgen_reader;
            PgenReader pgen_reader;
            VcfReader vcf_reader;
            std::unique_ptr<SequentialGenotypeMatrix> generic_reader;
            
            // Load source genotype
            DataLoader::GenotypeMeta meta;
            if (!cfg.hmp_file.empty() || !cfg.geno_matrix_file.empty() || !cfg.npy_file.empty() || !cfg.ped_file.empty() || !cfg.bfile_prefix.empty() || !cfg.pfile_prefix.empty() || !cfg.bgen_file.empty() || !cfg.vcf_file.empty()) {
                if (!cfg.bfile_prefix.empty()) {
                    meta = DataLoader::loadGenotype(cfg.bfile_prefix, cfg.bgen_file, cfg.sample_file,
                                                    cfg.pfile_prefix, cfg.pgen_file, cfg.pvar_file, cfg.psam_file,
                                                    cfg.vcf_file, cfg.vcf_tbi_file, cfg.vcf_sample_file,
                                                    plink_reader, bgen_reader, pgen_reader, vcf_reader);
                } else if (!cfg.pfile_prefix.empty()) {
                    meta = DataLoader::loadGenotype(cfg.bfile_prefix, cfg.bgen_file, cfg.sample_file,
                                                    cfg.pfile_prefix, cfg.pgen_file, cfg.pvar_file, cfg.psam_file,
                                                    cfg.vcf_file, cfg.vcf_tbi_file, cfg.vcf_sample_file,
                                                    plink_reader, bgen_reader, pgen_reader, vcf_reader);
                } else if (!cfg.bgen_file.empty()) {
                    meta = DataLoader::loadGenotype(cfg.bfile_prefix, cfg.bgen_file, cfg.sample_file,
                                                    cfg.pfile_prefix, cfg.pgen_file, cfg.pvar_file, cfg.psam_file,
                                                    cfg.vcf_file, cfg.vcf_tbi_file, cfg.vcf_sample_file,
                                                    plink_reader, bgen_reader, pgen_reader, vcf_reader);
                } else if (!cfg.vcf_file.empty()) {
                    meta = DataLoader::loadGenotype(cfg.bfile_prefix, cfg.bgen_file, cfg.sample_file,
                                                    cfg.pfile_prefix, cfg.pgen_file, cfg.pvar_file, cfg.psam_file,
                                                    cfg.vcf_file, cfg.vcf_tbi_file, cfg.vcf_sample_file,
                                                    plink_reader, bgen_reader, pgen_reader, vcf_reader);
                } else if (!cfg.hmp_file.empty()) {
                    auto reader = std::make_unique<HapMapReader>();
                    reader->open(cfg.hmp_file);
                    generic_reader = std::move(reader);
                    meta.n_samples = static_cast<int>(generic_reader->getNumSamples());
                    meta.n_snps = static_cast<int>(generic_reader->getNumVariants());
                } else if (!cfg.geno_matrix_file.empty()) {
                    auto reader = std::make_unique<TextGenoReader>();
                    const auto orientation = cfg.geno_orientation == "samples-by-variants"
                        ? TextGenoOrientation::SamplesByVariants
                        : TextGenoOrientation::VariantsBySamples;
                    reader->open(cfg.geno_matrix_file, orientation, cfg.geno_var_file, cfg.geno_sample_file);
                    generic_reader = std::move(reader);
                    meta.n_samples = static_cast<int>(generic_reader->getNumSamples());
                    meta.n_snps = static_cast<int>(generic_reader->getNumVariants());
                } else if (!cfg.npy_file.empty()) {
                    auto reader = std::make_unique<NpyReader>();
                    const auto orientation = cfg.npy_orientation == "samples-by-variants"
                        ? NpyOrientation::SamplesByVariants
                        : NpyOrientation::VariantsBySamples;
                    reader->open(cfg.npy_file, orientation, cfg.npy_var_file, cfg.npy_sample_file);
                    generic_reader = std::move(reader);
                    meta.n_samples = static_cast<int>(generic_reader->getNumSamples());
                    meta.n_snps = static_cast<int>(generic_reader->getNumVariants());
                } else if (!cfg.ped_file.empty()) {
                    auto reader = std::make_unique<PedMapReader>();
                    reader->open(cfg.ped_file, cfg.map_file);
                    generic_reader = std::move(reader);
                    meta.n_samples = static_cast<int>(generic_reader->getNumSamples());
                    meta.n_snps = static_cast<int>(generic_reader->getNumVariants());
                }
            }
            
            // Perform conversion
            std::string out_prefix = cfg.out_prefix;
            if (out_prefix == "cosmic_gwas") out_prefix = "converted";
            
            if (cfg.convert_to == "text") {
                if (generic_reader) {
                    generic_reader->exportToText(out_prefix + ".geno", out_prefix + ".var.info", out_prefix + ".sample.info");
                } else if (meta.is_bgen) {
                    bgen_reader.exportToText(out_prefix + ".geno", out_prefix + ".var.info", out_prefix + ".sample.info");
                } else if (meta.is_pgen) {
                    pgen_reader.exportToText(out_prefix + ".geno", out_prefix + ".var.info", out_prefix + ".sample.info");
                } else {
                    plink_reader.exportToText(out_prefix + ".geno", out_prefix + ".var.info", out_prefix + ".sample.info");
                }
                std::cout << "Converted to text format: " << out_prefix << ".geno/.var.info/.sample.info\n";
            } else if (cfg.convert_to == "npy") {
                if (generic_reader) {
                    generic_reader->exportToNpy(out_prefix + ".npy", out_prefix + ".var.info", out_prefix + ".sample.info");
                } else if (meta.is_bgen) {
                    bgen_reader.exportToNpy(out_prefix + ".npy", out_prefix + ".var.info", out_prefix + ".sample.info");
                } else if (meta.is_pgen) {
                    pgen_reader.exportToNpy(out_prefix + ".npy", out_prefix + ".var.info", out_prefix + ".sample.info");
                } else {
                    plink_reader.exportToNpy(out_prefix + ".npy", out_prefix + ".var.info", out_prefix + ".sample.info");
                }
                std::cout << "Converted to NPY format: " << out_prefix << ".npy/.var.info/.sample.info\n";
            } else if (cfg.convert_to == "ped") {
                if (generic_reader) {
                    generic_reader->exportToPedMap(out_prefix + ".ped", out_prefix + ".map");
                } else if (meta.is_bgen) {
                    bgen_reader.exportToPedMap(out_prefix + ".ped", out_prefix + ".map");
                } else if (meta.is_pgen) {
                    pgen_reader.exportToPedMap(out_prefix + ".ped", out_prefix + ".map");
                } else {
                    plink_reader.exportToPedMap(out_prefix + ".ped", out_prefix + ".map");
                }
                std::cout << "Converted to PED/MAP format: " << out_prefix << ".ped/.map\n";
            } else if (cfg.convert_to == "vcf") {
                // VCF export requires phased data; basic implementation
                if (generic_reader) {
                    generic_reader->exportToVcf(out_prefix + ".vcf");
                } else if (meta.is_bgen) {
                    bgen_reader.exportToVcf(out_prefix + ".vcf");
                } else if (meta.is_pgen) {
                    pgen_reader.exportToVcf(out_prefix + ".vcf");
                } else {
                    plink_reader.exportToVcf(out_prefix + ".vcf");
                }
                std::cout << "Converted to VCF format: " << out_prefix << ".vcf\n";
            } else {
                throw std::runtime_error("Unsupported target format: " + cfg.convert_to + " (supported: text, npy, ped, vcf)");
            }
            
            std::cout << "Format conversion completed successfully.\n";
            return 0;
        }
        
        // Handle batch runner
        if (!cfg.batch_models.empty()) {
            std::cout << "--- Batch Model Runner ---\n";
            
            std::vector<std::string> models;
            std::stringstream ss(cfg.batch_models);
            std::string model;
            while (std::getline(ss, model, ',')) {
                // Trim whitespace
                model.erase(0, model.find_first_not_of(" \t"));
                model.erase(model.find_last_not_of(" \t") + 1);
                if (!model.empty()) models.push_back(model);
            }
            
            std::string out_prefix = cfg.batch_out_prefix;
            if (out_prefix.empty()) out_prefix = cfg.out_prefix;
            if (out_prefix == "cosmic_gwas") out_prefix = "batch";
            
            int failed_models = 0;
            for (const auto& m : models) {
                std::cout << "\n=== Running model: " << m << " ===\n";
                
                GwasConfig batch_cfg = cfg;
                batch_cfg.out_prefix = out_prefix + "." + m;
                batch_cfg.batch_models.clear();  // prevent recursive batch entry
                
                // Set model-specific flags
                if (m == "lm") {
                    batch_cfg.linear_model = true;
                    batch_cfg.linear_model_mode = 4;
                    batch_cfg.run_gwas = true;
                } else if (m == "lmm") {
                    batch_cfg.run_gwas = true;
                    batch_cfg.lmm_gemma_mode = 4;
                } else if (m == "fastlmm") {
                    batch_cfg.run_gwas = true;
                    batch_cfg.fastlmm = true;
                } else if (m == "fastlmm-matrix-free") {
                    batch_cfg.run_gwas = true;
                    batch_cfg.fastlmm = true;
                    batch_cfg.fastlmm_matrix_free = true;
                } else if (m == "saige") {
                    batch_cfg.run_gwas = true;
                    batch_cfg.saige_style = true;
                    batch_cfg.biobank_mode = true;
                } else if (m == "farmcpu") {
                    batch_cfg.run_gwas = true;
                    batch_cfg.farmcpu = true;
                } else if (m == "blink") {
                    batch_cfg.run_gwas = true;
                    batch_cfg.blink = true;
                } else {
                    std::cerr << "Error: Unknown model '" << m << "', skipping\n";
                    ++failed_models;
                    continue;
                }
                
                GwasApp batch_app(batch_cfg);
                const int status = batch_app.run();
                if (status != 0) {
                    ++failed_models;
                    std::cerr << "Batch model '" << m << "' failed with status " << status << "\n";
                }
            }
            
            if (models.empty()) {
                std::cerr << "Error: No batch models were specified.\n";
                return 1;
            }
            if (failed_models > 0) {
                std::cerr << "Batch run failed: " << failed_models << "/" << models.size()
                          << " models failed.\n";
                return 1;
            }
            std::cout << "\nBatch run completed: " << models.size() << " models succeeded.\n";
            return 0;
        }
        
        auto start_time = chrono::high_resolution_clock::now();
        auto start_wall = chrono::system_clock::now();
        cout << "CosmicScanner v" << COSMICSCANNER_VERSION << "\n";

        if (cfg.fine_map) {
            std::cout << "--- Native Conditional Fine-Mapping ---\n";
            const auto region = association::fine_mapping::parseRegion(cfg.finemap_region);
            association::cojo::CoJoOptions ld_options;
            ld_options.maf_filter = 0.0;
            ld_options.collinearity = 2.0;
            ld_options.region_chrom = region.chrom;
            ld_options.region_start = region.start;
            ld_options.region_end = region.end;
            ld_options.canonical_summary_only = true;
            association::cojo::CoJoEngine ld_engine(ld_options);

            PlinkReader plink_reader;
            BgenReader bgen_reader;
            PgenReader pgen_reader;
            VcfReader vcf_reader;
            auto meta = DataLoader::loadGenotype(cfg.bfile_prefix, cfg.bgen_file, cfg.sample_file,
                                                 cfg.pfile_prefix, cfg.pgen_file, cfg.pvar_file, cfg.psam_file,
                                                 cfg.vcf_file, cfg.vcf_tbi_file, cfg.vcf_sample_file,
                                                 plink_reader, bgen_reader, pgen_reader, vcf_reader);
            ld_engine.loadLDReferenceFromGenotype(cfg.finemap_summary,
                                                  &plink_reader, &bgen_reader, &pgen_reader, &vcf_reader,
                                                  meta.is_bgen, meta.is_pgen, meta.is_vcf);

            association::fine_mapping::FineMappingOptions options;
            options.prior_variance = cfg.finemap_prior_var;
            options.credible_level = cfg.credible_level;
            options.max_signals = cfg.finemap_max_signals;
            options.signal_p_threshold = cfg.finemap_signal_p;
            association::fine_mapping::FineMappingEngine engine(options);
            auto result = engine.run(ld_engine.getSummaryStats(), ld_engine.getLDReference());
            const std::string ld_source = !cfg.bfile_prefix.empty() ? "BED" : meta.is_pgen ? "PGEN" :
                                          meta.is_bgen ? "BGEN" : "VCF";
            association::fine_mapping::writeFineMappingResults(result, options, cfg.out_prefix, ld_source);
            std::cout << "Fine-mapping completed: variants=" << result.number_variants
                      << ", signals=" << result.number_signals
                      << ", LD=" << result.ld_dimension << "x" << result.ld_dimension
                      << ", runtime=" << result.runtime_seconds << "s\n";
            return 0;
        }

        // Handle COJO and meta-analysis before loading genotype data
        if (cfg.cojo || cfg.cojo_slct || cfg.cojo_joint) {
            std::cout << "--- COJO-style Conditional Joint Analysis ---\n";
            association::cojo::CoJoOptions cojo_opts;
            if (cfg.cojo) cojo_opts.mode = association::cojo::CoJoMode::CONDITIONAL;
            else if (cfg.cojo_slct) cojo_opts.mode = association::cojo::CoJoMode::STEPWISE;
            else cojo_opts.mode = association::cojo::CoJoMode::JOINT;
            cojo_opts.threshold = cfg.cojo_p;
            cojo_opts.collinearity = cfg.cojo_r2;
            cojo_opts.window_size = cfg.cojo_window;
            cojo_opts.max_iter = cfg.cojo_max_iter;
            association::cojo::CoJoEngine engine(cojo_opts);
            // COJO requires a genotype/LD reference; identity-LD fallback is not permitted
            if (!cfg.bfile_prefix.empty() || !cfg.bgen_file.empty() || !cfg.pfile_prefix.empty() ||
                !cfg.pgen_file.empty() || !cfg.vcf_file.empty()) {
                // We have genotype data - use it to build real LD matrix
                // Note: we need to load genotype data first
                PlinkReader plink_reader;
                BgenReader bgen_reader;
                PgenReader pgen_reader;
                VcfReader vcf_reader;
                auto meta = DataLoader::loadGenotype(cfg.bfile_prefix, cfg.bgen_file, cfg.sample_file,
                                                     cfg.pfile_prefix, cfg.pgen_file, cfg.pvar_file, cfg.psam_file,
                                                     cfg.vcf_file, cfg.vcf_tbi_file, cfg.vcf_sample_file,
                                                     plink_reader, bgen_reader, pgen_reader, vcf_reader);
                bool is_bgen = meta.is_bgen;
                bool is_pgen = meta.is_pgen;
                bool is_vcf = meta.is_vcf;
                engine.loadLDReferenceFromGenotype(cfg.cojo_summary,
                                                   &plink_reader, &bgen_reader, &pgen_reader, &vcf_reader,
                                                   is_bgen, is_pgen, is_vcf);
            } else {
                engine.loadSummaryStats(cfg.cojo_summary);
                engine.loadLDReference(cfg.cojo_ld_ref);
            }
            if (!cfg.cojo_cond_file.empty()) engine.loadConditionSnps(cfg.cojo_cond_file);
            auto results = engine.run();
            association::cojo::writeCoJoResults(results, cfg.out_prefix, true);
            return 0;
        }

        if (cfg.meta_analysis) {
            std::cout << "--- Summary-Statistics Meta-Analysis ---\n";
            association::meta::MetaAnalysisOptions meta_opts;
            if (cfg.meta_method == "fixed") meta_opts.method = association::meta::MetaAnalysisMethod::FIXED_EFFECT;
            else if (cfg.meta_method == "random-dl") meta_opts.method = association::meta::MetaAnalysisMethod::RANDOM_EFFECT_DL;
            else meta_opts.method = association::meta::MetaAnalysisMethod::RANDOM_EFFECT_REML;
            meta_opts.heterogeneity_threshold = cfg.meta_het_p;
            meta_opts.min_studies = cfg.meta_min_studies;
            meta_opts.use_gc_correction = cfg.meta_gc;
            meta_opts.gc_lambda = cfg.meta_gc_lambda;
            association::meta::MetaAnalysisEngine engine(meta_opts);
            if (!cfg.meta_studies_file.empty()) {
                std::ifstream mf(cfg.meta_studies_file);
                if (!mf) throw std::runtime_error("Cannot open meta-analysis study list: " + cfg.meta_studies_file);
                std::string study_file;
                int study_idx = 0;
                while (std::getline(mf, study_file)) {
                    const auto first = study_file.find_first_not_of(" \t\r");
                    if (first == std::string::npos || study_file[first] == '#') continue;
                    const auto last = study_file.find_last_not_of(" \t\r");
                    study_file = study_file.substr(first, last - first + 1);
                    engine.addStudyFromFile("study" + std::to_string(++study_idx), study_file);
                }
            }
            for (const auto& f : cfg.meta_study_files) {
                engine.addStudyFromFile("study" + std::to_string(engine.getStudies().size() + 1), f);
            }
            if (engine.getStudies().empty()) throw std::runtime_error("Meta-analysis study list is empty");
            auto results = engine.run();
            association::meta::MetaWriterOptions writer_opts;
            writer_opts.write_metal = (cfg.meta_out_format == "metal");
            writer_opts.write_simple = (cfg.meta_out_format == "simple");
            association::meta::writeMetaResults(results, cfg.out_prefix, writer_opts);
            return 0;
        }

        PlinkReader plink_reader;
        BgenReader bgen_reader;
        PgenReader pgen_reader;
        VcfReader vcf_reader;
        std::unique_ptr<SequentialGenotypeMatrix> generic_reader;
DataLoader::GenotypeMeta meta;
            const bool ped_is_pedigree = (cfg.single_step_model || cfg.animal_model || cfg.make_agh);
            if (!cfg.hmp_file.empty() || !cfg.geno_matrix_file.empty() || !cfg.npy_file.empty() || (!ped_is_pedigree && !cfg.ped_file.empty())) {
            if (!cfg.hmp_file.empty()) {
                auto reader = std::make_unique<HapMapReader>();
                reader->open(cfg.hmp_file);
                generic_reader = std::move(reader);
            } else if (!cfg.geno_matrix_file.empty()) {
                auto reader = std::make_unique<TextGenoReader>();
                const auto orientation = cfg.geno_orientation == "samples-by-variants"
                    ? TextGenoOrientation::SamplesByVariants
                    : TextGenoOrientation::VariantsBySamples;
                reader->open(cfg.geno_matrix_file, orientation, cfg.geno_var_file, cfg.geno_sample_file);
                generic_reader = std::move(reader);
            } else if (!cfg.npy_file.empty()) {
                auto reader = std::make_unique<NpyReader>();
                const auto orientation = cfg.npy_orientation == "samples-by-variants"
                    ? NpyOrientation::SamplesByVariants
                    : NpyOrientation::VariantsBySamples;
                reader->open(cfg.npy_file, orientation, cfg.npy_var_file, cfg.npy_sample_file);
                generic_reader = std::move(reader);
            } else {
                auto reader = std::make_unique<PedMapReader>();
                reader->open(cfg.ped_file, cfg.map_file);
                generic_reader = std::move(reader);
            }
            meta.n_samples = static_cast<int>(generic_reader->getNumSamples());
            meta.n_snps = static_cast<int>(generic_reader->getNumVariants());
            for (const auto& sample_id : generic_reader->getSampleIds()) {
                const size_t separator = sample_id.find(':');
                meta.fids.push_back(separator == std::string::npos ? sample_id : sample_id.substr(0, separator));
                meta.iids.push_back(separator == std::string::npos ? sample_id : sample_id.substr(separator + 1));
            }
            std::cout << "Loaded " << meta.n_samples << " individuals and " << meta.n_snps
                      << " variants from " << (!cfg.hmp_file.empty() ? "HapMap" : 
                                              (!cfg.geno_matrix_file.empty() ? "text genotype matrix" : 
                                               (!cfg.npy_file.empty() ? "NPY" : "PED/MAP"))) << ".\n";
        } else {
            meta = DataLoader::loadGenotype(cfg.bfile_prefix, cfg.bgen_file, cfg.sample_file,
                                            cfg.pfile_prefix, cfg.pgen_file, cfg.pvar_file, cfg.psam_file,
                                            cfg.vcf_file, cfg.vcf_tbi_file, cfg.vcf_sample_file,
                                            plink_reader, bgen_reader, pgen_reader, vcf_reader);
        }
        bool is_bgen = meta.is_bgen;
        bool is_pgen = meta.is_pgen;
        bool is_vcf = meta.is_vcf;
        int n_samples = meta.n_samples;
        vector<string> fids = std::move(meta.fids);
        vector<string> iids = std::move(meta.iids);

        if (generic_reader && (cfg.multivariate || cfg.repeatability_model || cfg.rrm_model ||
            cfg.maternal_model || cfg.common_env_model || cfg.single_step_model || cfg.animal_model ||
            cfg.threshold_model || cfg.count_model || cfg.survival_model || cfg.gxe_model ||
            cfg.bslmm_mode > 0 || cfg.vc_mode > 0 || cfg.gemma_eigen || cfg.gemma_calccor ||
            cfg.cojo || cfg.cojo_slct || cfg.cojo_joint)) {
            throw std::runtime_error("HapMap/text genotype input currently supports dense GRM and univariate LM/LMM/SAIGE only");
        }

        if (!cfg.input_file.empty()) {
            InputPrepareOptions prep_opts;
            prep_opts.input_file = cfg.input_file;
            prep_opts.input_format = cfg.input_format;
            prep_opts.input_model = cfg.input_model;
            prep_opts.output_dir = cfg.prepared_input_dir.empty() ? (cfg.out_prefix + ".prepared_input") : cfg.prepared_input_dir;
            prep_opts.id_col = cfg.id_col;
            prep_opts.iid_col = cfg.iid_col;
            prep_opts.time_col = cfg.time_col;
            prep_opts.trait_cols = split_csv_list(cfg.trait_cols);
            prep_opts.input_trait = cfg.input_trait;
            prep_opts.record_qcovars = split_csv_list(cfg.record_qcovars);
            prep_opts.record_ccovars = split_csv_list(cfg.record_ccovars);
            prep_opts.subject_qcovars = split_csv_list(cfg.subject_qcovars);
            prep_opts.subject_ccovars = split_csv_list(cfg.subject_ccovars);
            prep_opts.time_parser = cfg.time_parser;
            prep_opts.time_unique_policy = cfg.time_unique_policy;
            prep_opts.time_index_multiplier = cfg.time_index_multiplier;
            prep_opts.preview_schema = cfg.preview_schema;
            prep_opts.infer_columns = cfg.infer_columns;
            prep_opts.wide_to_long = cfg.wide_to_long;
            auto prepared = InputPreparer::prepare(prep_opts, fids, iids);
            cfg.prepared_input_dir = prepared.output_dir;
            cfg.tpheno = 1;
            cfg.mpheno = 2;
            cfg.tpheno_set = true;
            if (cfg.input_model == "rrm") {
                cfg.pheno_file = prepared.rrm_pheno_file;
                cfg.long_pheno_file.clear();
                cfg.qcovar_file.clear();
                cfg.covar_file.clear();
                cfg.rrm_unbalance = true;
                cfg.tcovar_names = join_csv_list(prep_opts.record_qcovars);
                cfg.tclass_names = join_csv_list(prep_opts.record_ccovars);
                cfg.rrm_model = true;
            } else if (cfg.input_model == "lmm") {
                cfg.pheno_file = prepared.subject_pheno_file;
                cfg.qcovar_file = prepared.subject_qcovar_file;
                cfg.covar_file = prepared.subject_covar_file;
                cfg.mpheno = 1;
            } else if (cfg.input_model == "multivariate") {
                cfg.pheno_file = prepared.multivariate_pheno_file;
                cfg.qcovar_file = prepared.subject_qcovar_file;
                cfg.covar_file = prepared.subject_covar_file;
                cfg.mpheno_list.clear();
                for (int i = 1; i <= (int)prep_opts.trait_cols.size(); ++i) cfg.mpheno_list.push_back(i);
                cfg.multivariate = true;
            } else {
                cfg.pheno_file = prepared.pheno_file;
                cfg.qcovar_file = prepared.qcovar_file;
                cfg.covar_file = prepared.covar_file;
                cfg.repeatability_model = true;
            }
            std::cout << "Prepared input written to " << prepared.output_dir
                      << " (matched records=" << prepared.matched_records
                      << ", unmatched records=" << prepared.unmatched_records
                      << ", selected trait=" << prepared.selected_trait << ")\n";
            if (cfg.preview_schema) {
                std::cout << "--preview-schema set; exiting after schema preview.\n";
                return 0;
            }
            if (cfg.prepare_input_only) {
                std::cout << "--prepare-input-only set; exiting after input preparation.\n";
                return 0;
            }
            } else if (!cfg.prepared_input_dir.empty() && cfg.pheno_file.empty() && cfg.long_pheno_file.empty()) {
            auto prepared = InputPreparer::loadPrepared(cfg.prepared_input_dir, cfg.input_trait, cfg.input_model);
            cfg.tpheno = 1;
            cfg.mpheno = 2;
            cfg.tpheno_set = true;
            if (cfg.input_model == "rrm") {
                cfg.pheno_file = prepared.rrm_pheno_file;
                cfg.long_pheno_file.clear();
                cfg.qcovar_file.clear();
                cfg.covar_file.clear();
                cfg.rrm_unbalance = true;
                collect_prepared_rrm_covariates(cfg.prepared_input_dir, cfg.tcovar_names, cfg.tclass_names);
                cfg.rrm_model = true;
            } else if (cfg.input_model == "lmm") {
                cfg.pheno_file = prepared.subject_pheno_file;
                cfg.qcovar_file = prepared.subject_qcovar_file;
                cfg.covar_file = prepared.subject_covar_file;
                cfg.mpheno = 1;
            } else if (cfg.input_model == "multivariate") {
                cfg.pheno_file = prepared.multivariate_pheno_file;
                cfg.qcovar_file = prepared.subject_qcovar_file;
                cfg.covar_file = prepared.subject_covar_file;
                cfg.mpheno_list = subject_multivariate_columns(prepared.multivariate_pheno_file);
                cfg.multivariate = true;
            } else {
                cfg.pheno_file = prepared.pheno_file;
                cfg.qcovar_file = prepared.qcovar_file;
                cfg.covar_file = prepared.covar_file;
                cfg.repeatability_model = true;
            }
            std::cout << "Using prepared input from " << prepared.output_dir
                      << " (selected trait=" << prepared.selected_trait << ")\n";
        }

        int standalone_count = 0;
        standalone_count += (cfg.make_grm || cfg.make_grm_sparse) ? 1 : 0;
        standalone_count += cfg.gemma_s_mode > 0 ? 1 : 0;
        standalone_count += cfg.gemma_q_mode > 0 ? 1 : 0;
        standalone_count += cfg.gemma_eigen ? 1 : 0;
        standalone_count += cfg.gemma_calccor ? 1 : 0;
        standalone_count += cfg.vc_mode > 0 ? 1 : 0;
        standalone_count += cfg.ci_mode > 0 ? 1 : 0;
        standalone_count += !cfg.predict_effect_file.empty() ? 1 : 0;
        standalone_count += cfg.bslmm_mode > 0 ? 1 : 0;
        if (standalone_count > 1) {
            throw std::runtime_error(
                "Multiple GEMMA-compatible standalone routes detected. "
                "Use only one of --make-grm/-gk, --gs/-gs, --gq/-gq, --gemma-eigen/-eigen, "
                "--vc/-vc, --ci/-ci, --calccor/-calccor, --predict-effect/-epm, or --bslmm/-bslmm per run.");
        }

        if (cfg.make_grm || cfg.make_grm_sparse) {
            DataLoader::runMakeGrm(cfg, plink_reader, bgen_reader, pgen_reader, vcf_reader,
                                   is_bgen, is_pgen, is_vcf, generic_reader.get());
            return 0;
        }

        if (cfg.gemma_s_mode > 0) {
            GemmaSummaryTools tools(cfg);
            tools.runGs(plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
        }

        if (cfg.gemma_q_mode > 0) {
            GemmaSummaryTools tools(cfg);
            tools.runGq();
        }

        if (cfg.gemma_eigen) {
            DataLoader::runGemmaEigen(cfg, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
        }

        if (cfg.gemma_calccor) {
            GemmaCalccor calccor(cfg);
            calccor.run(plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
        }

        if (cfg.vc_mode > 0) {
            if (!cfg.study_prefix.empty() || !cfg.ref_prefix.empty() ||
                !cfg.multi_study_file.empty() || !cfg.multi_ref_file.empty()) {
                GemmaSummaryTools tools(cfg);
                tools.runSummaryVc();
            } else if (!cfg.beta_file.empty() || !cfg.cor_file.empty()) {
                GemmaSummaryTools tools(cfg);
                tools.runBetaCorVc();
            } else {
                GemmaVcModel model(GemmaVcModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader,
                              is_bgen, is_pgen, is_vcf);
                model.run();
            }
        }

        if (cfg.ci_mode > 0) {
            GemmaSummaryTools tools(cfg);
            tools.runSummaryCi();
        }

        if (!cfg.predict_effect_file.empty()) {
            GemmaPredictOptions pred_opts;
            pred_opts.effect_file = cfg.predict_effect_file;
            pred_opts.output_file = cfg.predict_score_file.empty() ? (cfg.out_prefix + ".predict.txt") : cfg.predict_score_file;
            pred_opts.snp_col = cfg.predict_snp_col;
            pred_opts.effect_col = cfg.predict_effect_col;
            GemmaPredictor predictor(pred_opts);
            predictor.run(fids, iids, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
        }

        if (cfg.bslmm_mode > 0) {
            GemmaBslmmModel model(GemmaBslmmModel::Options::fromConfig(cfg));
            model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
            model.run();
        }

        if (cfg.run_gwas) {
            cout << (cfg.make_agh ? "Building A/G/H matrices...\n" : "Running GWAS...\n");

            // 检测互斥模型标志冲突
            int model_count = 0;
            model_count += cfg.rrm_model;
            model_count += cfg.repeatability_model;
            model_count += cfg.multivariate;
            model_count += cfg.maternal_model;
            model_count += cfg.single_step_model;
            model_count += cfg.animal_model;
            model_count += cfg.common_env_model;
            model_count += cfg.threshold_model;
            model_count += cfg.count_model;
            model_count += cfg.ordinal_model;
            model_count += cfg.survival_model;
            model_count += cfg.gxe_model;
            model_count += cfg.linear_model;
            model_count += cfg.farmcpu;
            model_count += cfg.blink;
            model_count += cfg.firth_logistic;
            model_count += cfg.saige_style;
            model_count += !cfg.logistic_test.empty();
            model_count += cfg.logistic_glmm;
            model_count += !cfg.condition_snp_file.empty();
            model_count += cfg.cojo + cfg.cojo_slct + cfg.cojo_joint;
            model_count += cfg.meta_analysis;
            if (model_count > 1) {
                throw std::runtime_error(
                    "Multiple model flags detected. Please use only one of:\n"
                    "  --long-pheno/--rrm-order, --repeatability-model, --multivariate,\n"
                    "  --maternal-model, --single-step, --animal-model, --common-env-model, --threshold-model,\n"
                    "  --count-model, --ordinal-model, --survival-model, --gxe-model, --lm,\n"
                    "  --farmcpu, --blink, --firth-logistic, --logistic, --logistic-glmm, --saige-style, --condition-snp");
            }
            if (!cfg.scan_tier.empty() &&
                (cfg.rrm_model || cfg.maternal_model || cfg.single_step_model ||
                 cfg.animal_model ||
                 cfg.threshold_model || cfg.count_model || cfg.ordinal_model ||
                 cfg.survival_model || cfg.gxe_model || cfg.linear_model)) {
                throw std::runtime_error(
                    "--scan-tier currently supports single-trait, multivariate, and repeatability Gaussian LMM GWAS only");
            }

            auto solver_plan = make_solver_plan(cfg, n_samples, meta.n_snps);
            log_solver_planner_diagnostics(solver_plan);
            if (!is_standard_lmm_contract_path(cfg)) {
                log_analysis_dataset_shadow(cfg, fids, iids, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
            }

            if (cfg.rrm_model) {
                RrmModel model(RrmModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.repeatability_model) {
                RepeatabilityModel model(RepeatabilityModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.common_env_model) {
                std::cout << "--- Common-Environment GWAS Model ---\n";
                CommonEnvironmentModel model(CommonEnvironmentModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.multivariate) {
                MultivariateModel model(MultivariateModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.maternal_model) {
                std::cout << "--- Maternal Effects GWAS Model ---\n";
                MaternalModel model(MaternalModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.single_step_model) {
                std::cout << "--- Single-Step GWAS (ssGWAS) Model ---\n";
                SingleStepModel model(SingleStepModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.animal_model) {
                std::cout << "--- Pedigree Animal Model (Simple Animal Model) ---\n";
                AnimalModel model(AnimalModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.threshold_model) {
                std::cout << "--- Threshold/Probit GWAS Model ---\n";
                ThresholdModel model(ThresholdModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.count_model) {
                std::cout << "--- Count/Binomial PSR-GWAS Model ---\n";
                CountModel model(CountModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.ordinal_model) {
                std::cout << "--- Ordinal PSR-GWAS Model ---\n";
                OrdinalModel model(OrdinalModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.survival_model) {
                std::cout << "--- Survival GWAS Model ---\n";
                SurvivalModel model(SurvivalModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.gxe_model) {
                std::cout << "--- GxE Reaction Norm GWAS Model ---\n";
                GxEModel model(GxEModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader, is_bgen, is_pgen, is_vcf);
                model.run();

            } else if (cfg.cojo || cfg.cojo_slct || cfg.cojo_joint) {
                std::cout << "--- COJO-style Conditional Joint Analysis ---\n";
                association::cojo::CoJoOptions cojo_opts;
                if (cfg.cojo) cojo_opts.mode = association::cojo::CoJoMode::CONDITIONAL;
                else if (cfg.cojo_slct) cojo_opts.mode = association::cojo::CoJoMode::STEPWISE;
                else cojo_opts.mode = association::cojo::CoJoMode::JOINT;
                cojo_opts.threshold = cfg.cojo_p;
                cojo_opts.collinearity = cfg.cojo_r2;
                cojo_opts.window_size = cfg.cojo_window;
                cojo_opts.max_iter = cfg.cojo_max_iter;
                association::cojo::CoJoEngine engine(cojo_opts);
                // Use already-loaded genotype data to build real LD matrix
                engine.loadLDReferenceFromGenotype(cfg.cojo_summary,
                                                   &plink_reader, &bgen_reader, &pgen_reader, &vcf_reader,
                                                   is_bgen, is_pgen, is_vcf);
                if (!cfg.cojo_cond_file.empty()) engine.loadConditionSnps(cfg.cojo_cond_file);
                auto results = engine.run();
                association::cojo::writeCoJoResults(results, cfg.out_prefix, true);

            } else if (cfg.meta_analysis) {
                std::cout << "--- Summary-Statistics Meta-Analysis ---\n";
                association::meta::MetaAnalysisOptions meta_opts;
                if (cfg.meta_method == "fixed") meta_opts.method = association::meta::MetaAnalysisMethod::FIXED_EFFECT;
                else if (cfg.meta_method == "random-dl") meta_opts.method = association::meta::MetaAnalysisMethod::RANDOM_EFFECT_DL;
                else meta_opts.method = association::meta::MetaAnalysisMethod::RANDOM_EFFECT_REML;
                meta_opts.heterogeneity_threshold = cfg.meta_het_p;
                meta_opts.min_studies = cfg.meta_min_studies;
                meta_opts.use_gc_correction = cfg.meta_gc;
                meta_opts.gc_lambda = cfg.meta_gc_lambda;
                association::meta::MetaAnalysisEngine engine(meta_opts);
                if (!cfg.meta_studies_file.empty()) {
                    std::ifstream mf(cfg.meta_studies_file);
                    std::string study_file;
                    int study_idx = 0;
                    while (mf >> study_file) {
                        engine.addStudyFromFile("study" + std::to_string(++study_idx), study_file);
                    }
                }
                for (const auto& f : cfg.meta_study_files) {
                    engine.addStudyFromFile("study" + std::to_string(engine.getStudies().size() + 1), f);
                }
                auto results = engine.run();
                association::meta::MetaWriterOptions writer_opts;
                writer_opts.write_metal = (cfg.meta_out_format == "metal");
                writer_opts.write_simple = (cfg.meta_out_format == "simple");
                association::meta::writeMetaResults(results, cfg.out_prefix, writer_opts);

            } else {
                UnivariateModel model(UnivariateModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, vcf_reader,
                              is_bgen, is_pgen, is_vcf, generic_reader.get());
                model.run();
                if (is_standard_lmm_contract_path(cfg) && model.hasAnalysisDataset()) {
                    auto manifest_end_wall = chrono::system_clock::now();
                    chrono::duration<double> manifest_elapsed = manifest_end_wall - start_wall;
                    write_scanner_manifest(cfg, model.analysisDataset(), solver_plan, model.validationReport(),
                                           start_wall, manifest_end_wall, manifest_elapsed.count(),
                                           cosmic::peak_rss_bytes());
                    cout << "Run manifest saved to " << cfg.out_prefix << ".run_manifest.json\n";
                }
            }
        }

        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end_time - start_time;
        cout << "\nTotal computation time: " << elapsed.count() << " seconds\n";

        const auto peak_bytes = cosmic::peak_rss_bytes();
        if (peak_bytes) {
            const double peak_mem_mb = static_cast<double>(peak_bytes) / (1024.0 * 1024.0);
            cout << "Peak memory usage: " << peak_mem_mb << " MB\n";
        }

        cout << "CosmicScanner finished successfully.\n";
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

} // end namespace cosmic
