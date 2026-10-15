#include "bgen_ad_scan_engine.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

static Eigen::MatrixXd readMatrixLib(const std::string& path) {
    std::ifstream input(path);
    int rows = 0, columns = 0;
    input >> rows >> columns;
    if (!input || rows < 1 || columns < 1) throw std::runtime_error("matrix header");
    Eigen::MatrixXd out(rows, columns);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < columns; ++j)
            if (!(input >> out(i, j))) throw std::runtime_error("matrix parse");
    return out;
}

static void dumpLib(const std::string& path, const Eigen::MatrixXd& value) {
    std::ofstream output(path, std::ios::binary);
    if (!output) throw std::runtime_error("matrix output");
    output.write(reinterpret_cast<const char*>(value.data()), value.size() * sizeof(double));
}

int main(int argc, char** argv) {
    try {
        if (argc != 11)
            throw std::invalid_argument(
                "test: bgen phenotypes covariates|- output scan_chr excluded_chr fixed|fitted a2 d2 e2");
        const std::string bgen = argv[1], phenotypePath = argv[2], covariatePath = argv[3],
                           output = argv[4], scanChr = argv[5], excludedChr = argv[6],
                           mode = argv[7];
        if (mode != "fixed" && mode != "fitted") throw std::invalid_argument("mode");
        const std::string excluded = excludedChr == "-" ? "" : excludedChr;
        const Eigen::MatrixXd phenotypes = readMatrixLib(phenotypePath);
        Eigen::MatrixXf covariates(phenotypes.rows(), 1);
        if (std::string(covariatePath) == "-") covariates.setOnes();
        else covariates = readMatrixLib(covariatePath).cast<float>();
        cosmic::BgenReader reader;
        reader.open(bgen);
        reader.initialize();
        cosmic::BgenAdKernelOptions kernelOptions;
        kernelOptions.max_missing = .05;
        kernelOptions.min_info_d = .01;
        const auto kernels = cosmic::BgenAdNullBackend::build(reader, kernelOptions);
        if (kernels.global().additive.rows() != phenotypes.rows())
            throw std::runtime_error("kernel/phenotype mismatch");
        const auto background = excluded.empty() ? kernels.global() : kernels.loco(excluded);
        dumpLib(output + ".ka.bin", background.additive);
        dumpLib(output + ".kd.bin", background.dominance);
        std::ofstream markers(output + ".markers.txt");
        for (const auto& id : kernels.accepted_variant_ids) markers << id << '\n';
        cosmic::BgenAdScanEngine::Options scanOptions;
        scanOptions.chromosome = scanChr;
        scanOptions.block_size = 256;
        scanOptions.max_missing = .05;
        scanOptions.min_info_d = .01;
        cosmic::BgenAdScanEngine engine(scanOptions);
        std::ofstream associations(output + ".internal.tsv");
        associations << std::setprecision(17)
                     << "replicate\tchr\tsnp\ta1\ta2\tfreq\tn\tbeta_a\tse_a\tp_add\tbeta_d\tse_d\tp_dom\tp_joint\tw_joint\n";
        std::ofstream nulls(output + ".null.tsv");
        nulls << std::setprecision(17)
              << "replicate\tsigma_a2\tsigma_d2\tsigma_e2\tconverged\tborder\tgradient\tkkt\tfits_in_scan\ttested\trank_deficient\tqc_filtered\ttarget_filtered\n";
        for (int rep = 0; rep < phenotypes.cols(); ++rep) {
            const Eigen::VectorXd y = phenotypes.col(rep);
            if (mode == "fixed") {
                const double a2 = std::stod(argv[8]), d2 = std::stod(argv[9]), e2 = std::stod(argv[10]);
                engine.prepareFixed(kernels, y, covariates, a2, d2, e2, excluded);
            } else {
                cosmic::EigenVceAD::Options options;
                options.verbose = false;
                options.max_em_iter = 200;
                options.tol = 1e-5;
                engine.fitNull(kernels, y, covariates, options, excluded);
            }
            const auto before = engine.nullFitCount();
            std::vector<cosmic::AdResult> results;
            engine.scan(reader, results);
            if (engine.nullFitCount() != before && mode == "fixed")
                throw std::runtime_error("fixed scan refit null");
            for (const auto& r : results) {
                if (r.rank_deficient) continue;
                associations << rep << '\t' << r.chrom << '\t' << r.snp_id << '\t' << r.a1 << '\t'
                              << r.a2 << '\t' << r.freq << '\t' << r.n_obs << '\t' << r.beta_add
                              << '\t' << r.se_add << '\t' << r.p_add << '\t' << r.beta_dom << '\t'
                              << r.se_dom << '\t' << r.p_dom << '\t' << r.p_joint << '\t'
                              << r.joint_chisq << '\n';
            }
            const auto& stats = engine.statistics();
            const bool wroteFit = engine.nullResult().has_value();
            if (wroteFit) {
                const auto& null = engine.nullResult().value();
                nulls << rep << '\t' << null.sigma_a2 << '\t' << null.sigma_d2 << '\t'
                      << null.sigma_e2 << '\t' << null.converged << '\t' << null.boundary_status
                      << '\t' << null.gradient_norm << '\t' << null.kkt_residual << '\t'
                      << engine.nullFitCount() - before << '\t' << stats.tested << '\t'
                      << stats.rank_deficient << '\t' << stats.qc_filtered << '\t'
                      << stats.target_filtered << '\n';
            } else {
                nulls << rep << "\t" << argv[8] << "\t" << argv[9] << "\t" << argv[10]
                      << "\ttrue\tFIXED\t0\t0\t" << engine.nullFitCount() - before << '\t'
                      << stats.tested << '\t' << stats.rank_deficient << '\t'
                      << stats.qc_filtered << '\t' << stats.target_filtered << '\n';
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "calibration driver failure: " << error.what() << '\n';
        return 1;
    }
}
