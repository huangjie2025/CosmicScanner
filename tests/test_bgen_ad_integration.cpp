// Explicit internal harness, never a production CosmicScanner CLI route.
#include "gwas_scanner.h"
#include "bgen_ad_integration.h"
#include "gwas_output_schema.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <iostream>

static Eigen::MatrixXd matrix(const std::string& path) {
    std::ifstream input(path); int n = 0, columns = 0;
    input >> n >> columns;
    if (!input || n <= 0 || columns <= 0) throw std::runtime_error("harness matrix header");
    Eigen::MatrixXd result(n, columns); std::string cell;
    for (int i = 0; i < n; ++i) for (int j = 0; j < columns; ++j) {
        if (!(input >> cell)) throw std::runtime_error("harness matrix truncated");
        result(i, j) = std::stod(cell);
    }
    return result;
}

static void dump(const std::string& path, const Eigen::MatrixXd& value) {
    std::ofstream output(path, std::ios::binary);
    if (!output) throw std::runtime_error("harness dump output");
    output.write(reinterpret_cast<const char*>(value.data()), value.size() * sizeof(double));
}

int main(int argc, char** argv) {
    try {
        if (argc != 6) throw std::invalid_argument("internal harness: bgen phenotypes covariates output_directory options.tsv");
        const Eigen::MatrixXd phenotypes = matrix(argv[2]);
        Eigen::MatrixXf covariates;
        if (std::string(argv[3]) == "-") covariates = Eigen::MatrixXf::Ones(phenotypes.rows(), 1);
        else covariates = matrix(argv[3]).cast<float>();
        std::map<std::string, std::string> config;
        std::ifstream input(argv[5]); std::string key, value;
        while (input >> key >> value) config[key] = value;
        const auto get = [&](const std::string& name, const std::string& fallback) {
            const auto entry = config.find(name); return entry == config.end() ? fallback : entry->second;
        };
        std::filesystem::create_directories(argv[4]);
        cosmic::BgenReader reader; reader.open(argv[1]); reader.initialize();
        for (int rep = 0; rep < phenotypes.cols(); ++rep) {
            cosmic::GwasScanner::Options options;
            options.gwas_model_md = "additive_dominance";
            options.out_prefix = std::string(argv[4]) + "/rep" + std::to_string(rep);
            options.no_loco = get("loco", "true") == "false";
            options.reml_maxit = std::stoi(get("maxit", "200"));
            options.reml_tol = std::stod(get("tol", "1e-5"));
            options.min_maf = std::stod(get("maf", ".01"));
            options.min_hwe = std::stod(get("hwe", "0"));
            options.bgen_dom_max_missing = std::stod(get("max_missing", ".05"));
            options.bgen_dom_min_info = std::stod(get("min_info_d", ".01"));
            options.bgen_dom_entropy_warning = std::stod(get("entropy_warning", ".5"));
            options.scan_chr = get("chromosome", "");
            options.scan_from_bp = std::stoll(get("from_bp", "-1"));
            options.scan_to_bp = std::stoll(get("to_bp", "-1"));
            options.extract_file = get("extract", ""); options.exclude_file = get("exclude", "");
            cosmic::BgenAdResearchOptions research;
            research.population_file = get("population", ""); research.snp_weight_file = get("weights", "");
            const Eigen::VectorXd y = phenotypes.col(rep);
            auto keep = cosmic::DataFilter::filterSamples(y, covariates);
            if (config.count("keep")) {
                keep.keep_indices.clear(); std::stringstream stream(config["keep"]); std::string cell;
                while (std::getline(stream, cell, ',')) keep.keep_indices.push_back(std::stoi(cell));
                keep.n_keep = keep.keep_indices.size();
            }
            cosmic::BgenAdIntegrationDiagnostics diagnostics;
            cosmic::GwasScanner scanner(options);
            scanner.runBgenAdResearchScan(y, covariates, reader, keep, research, &diagnostics);
            if (!diagnostics.kernels) throw std::runtime_error("integration kernel diagnostics missing");
            const auto pair = diagnostics.kernels->global();
            dump(options.out_prefix + ".ka.bin", pair.additive);
            dump(options.out_prefix + ".kd.bin", pair.dominance);
            for (const auto& null : diagnostics.nulls) if (null.label != "global") {
                const auto outside = diagnostics.kernels->loco(null.label);
                dump(options.out_prefix + ".loco." + null.label + ".ka.bin", outside.additive);
                dump(options.out_prefix + ".loco." + null.label + ".kd.bin", outside.dominance);
            }
            // Controlled-null references distinguish encoding/scan parity from
            // the known legacy-vs-robust optimizer accuracy difference. Native
            // algorithms and existing unified row writer are reused unchanged.
            if (config.count("native_bed") && config.count("native_pgen")) {
                cosmic::PlinkReader native_bed; native_bed.load(config["native_bed"]);
                cosmic::PgenReader native_pgen; native_pgen.open(config["native_pgen"]);
                const auto ys = cosmic::DataFilter::subset(y, keep);
                const auto xs = cosmic::DataFilter::subsetRows(covariates, keep);
                for (const auto& kind : {std::string("bed"), std::string("pgen")}) {
                    std::ofstream reference(options.out_prefix + ".reference." + kind + ".ad.txt");
                    reference << std::setprecision(9) << cosmic::gwas_output::kUnifiedHeader;
                    for (const auto& null : diagnostics.nulls) {
                        const auto background = null.label == "global" ? diagnostics.kernels->global() : diagnostics.kernels->loco(null.label);
                        cosmic::EigenVceAD::Options fixed_options; fixed_options.verbose = false;
                        cosmic::EigenVceAD fixed(fixed_options);
                        fixed.prepare(ys, xs, background.additive, background.dominance);
                        fixed.fixedComponents(null.fit.sigma_a2, null.fit.sigma_d2, null.fit.sigma_e2);
                        cosmic::LmmEngineAD native; native.prepareFromVce(fixed);
                        std::vector<cosmic::AdResult> results;
                        if (kind == "bed") native.scan(native_bed, results, keep.keep_indices, null.label == "global" ? "" : null.label);
                        else native.scan(native_pgen, results, keep.keep_indices, null.label == "global" ? "" : null.label);
                        for (const auto& result : results) if (!result.rank_deficient)
                            cosmic::gwas_output::writeAdRow(reference, result);
                    }
                }
            }
            std::ofstream ids(options.out_prefix + ".analysis.ids.txt");
            for (const auto& id : diagnostics.kernels->sample_ids) ids << id << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Internal BGEN AD integration failure: " << error.what() << '\n';
        return 1;
    }
}
