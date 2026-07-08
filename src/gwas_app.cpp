#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <cctype>
#include <fstream>
#include <sstream>
#include <sys/resource.h>
#include <omp.h>

#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "data_loader.h"
#include "models/repeatability_model.h"
#include "models/multivariate_model.h"
#include "models/rrm_model.h"
#include "models/univariate_model.h"
#include "models/maternal_model.h"
#include "models/single_step_model.h"
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

extern "C" {
    void openblas_set_num_threads(int num_threads) __attribute__((weak));
}

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

using namespace std;
using namespace cosmic;

#include "gwas_app.h"

namespace cosmic {

GwasApp::GwasApp(const GwasConfig& config) : cfg(config) {}

int GwasApp::run() {
    try {
        auto start_time = chrono::high_resolution_clock::now();
        #ifndef COSMICSCANNER_VERSION
        #define COSMICSCANNER_VERSION "unknown"
        #endif
        cout << "CosmicScanner v" << COSMICSCANNER_VERSION << "\n";

        PlinkReader plink_reader;
        BgenReader bgen_reader;
        PgenReader pgen_reader;
        auto meta = DataLoader::loadGenotype(cfg.bfile_prefix, cfg.bgen_file, cfg.sample_file,
                                             cfg.pfile_prefix, cfg.pgen_file, cfg.pvar_file, cfg.psam_file,
                                             plink_reader, bgen_reader, pgen_reader);
        bool is_bgen = meta.is_bgen;
        bool is_pgen = meta.is_pgen;
        int n_samples = meta.n_samples;
        vector<string> fids = std::move(meta.fids);
        vector<string> iids = std::move(meta.iids);

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
            DataLoader::runMakeGrm(cfg, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
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
            calccor.run(plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
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
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
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
            predictor.run(fids, iids, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
        }

        if (cfg.bslmm_mode > 0) {
            GemmaBslmmModel model(GemmaBslmmModel::Options::fromConfig(cfg));
            model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
            model.run();
        }

        if (cfg.run_gwas) {
            cout << "Running GWAS...\n";

            // 检测互斥模型标志冲突
            int model_count = 0;
            model_count += cfg.rrm_model;
            model_count += cfg.repeatability_model;
            model_count += cfg.multivariate;
            model_count += cfg.maternal_model;
            model_count += cfg.single_step_model;
            model_count += cfg.threshold_model;
            model_count += cfg.count_model;
            model_count += cfg.ordinal_model;
            model_count += cfg.survival_model;
            model_count += cfg.gxe_model;
            model_count += cfg.linear_model;
            if (model_count > 1) {
                throw std::runtime_error(
                    "Multiple model flags detected. Please use only one of:\n"
                    "  --long-pheno/--rrm-order, --repeatability-model, --multivariate,\n"
                    "  --maternal-model, --single-step, --threshold-model,\n"
                    "  --count-model, --ordinal-model, --survival-model, --gxe-model, --lm");
            }

            if (cfg.rrm_model) {
                RrmModel model(RrmModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
                model.run();

            } else if (cfg.repeatability_model) {
                RepeatabilityModel model(RepeatabilityModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
                model.run();

            } else if (cfg.multivariate) {
                MultivariateModel model(MultivariateModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
                model.run();

            } else if (cfg.maternal_model) {
                std::cout << "--- Maternal Effects GWAS Model ---\n";
                MaternalModel model(MaternalModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
                model.run();

            } else if (cfg.single_step_model) {
                std::cout << "--- Single-Step GWAS (ssGWAS) Model ---\n";
                SingleStepModel model(SingleStepModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
                model.run();

            } else if (cfg.threshold_model) {
                std::cout << "--- Threshold/Probit GWAS Model ---\n";
                ThresholdModel model(ThresholdModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
                model.run();

            } else if (cfg.count_model) {
                std::cout << "--- Count/Binomial PSR-GWAS Model ---\n";
                CountModel model(CountModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
                model.run();

            } else if (cfg.ordinal_model) {
                std::cout << "--- Ordinal PSR-GWAS Model ---\n";
                OrdinalModel model(OrdinalModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
                model.run();

            } else if (cfg.survival_model) {
                std::cout << "--- Survival GWAS Model ---\n";
                SurvivalModel model(SurvivalModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
                model.run();

            } else if (cfg.gxe_model) {
                std::cout << "--- GxE Reaction Norm GWAS Model ---\n";
                GxEModel model(GxEModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
                model.run();

            } else {
                UnivariateModel model(UnivariateModel::Options::fromConfig(cfg));
                model.prepare(fids, iids, n_samples, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
                model.run();
            }
        }

        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end_time - start_time;
        cout << "\nTotal computation time: " << elapsed.count() << " seconds\n";

        struct rusage usage;
        if (getrusage(RUSAGE_SELF, &usage) == 0) {
            double peak_mem_mb = usage.ru_maxrss / 1024.0;
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
