// Test-only E3.2a driver. No GwasScanner/CLI/output integration.
#include "bgen_ad_scan_engine.h"
#include "genetics/genotype_encoding.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

static Eigen::MatrixXd readMatrix(const std::string& path) {
    std::ifstream input(path);
    int n = 0, m = 0;
    input >> n >> m;
    if (!input || n < 1 || m < 1) throw std::runtime_error("test matrix header");
    Eigen::MatrixXd result(n, m);
    for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j) input >> result(i, j);
    if (!input) throw std::runtime_error("test matrix truncated");
    return result;
}

static void dumpMatrix(const std::string& path, const Eigen::MatrixXd& matrix) {
    std::ofstream output(path, std::ios::binary);
    output.write(reinterpret_cast<const char*>(matrix.data()), matrix.size() * sizeof(double));
}

static void writeRows(std::ostream& output, int rep, const std::vector<cosmic::AdResult>& results) {
    for (const auto& r : results) {
        if (r.rank_deficient) continue;
        output << rep << '\t' << r.chrom << '\t' << r.snp_id << '\t' << r.a1 << '\t' << r.a2
               << '\t' << r.freq << '\t' << r.n_obs << '\t' << r.beta_add << '\t' << r.se_add
               << '\t' << r.p_add << '\t' << r.beta_dom << '\t' << r.se_dom << '\t' << r.p_dom
               << '\t' << r.p_joint << '\t' << r.joint_chisq << '\n';
    }
}

int main(int argc, char** argv) {
    try {
        if (argc != 9) throw std::invalid_argument("test: bgen bed|- pgen|- phenotypes output fixed|fitted|parity indices|- covariates|-");
        const std::string bgen_path = argv[1], bed_path = argv[2], pgen_path = argv[3];
        const std::string output = argv[5], mode = argv[6];
        if (mode != "fixed" && mode != "fitted" && mode != "parity") throw std::invalid_argument("test scan mode");
        const Eigen::MatrixXd phenotypes = readMatrix(argv[4]);
        Eigen::MatrixXf covariates = std::string(argv[8]) == "-"
            ? Eigen::MatrixXf::Ones(phenotypes.rows(), 1)
            : readMatrix(argv[8]).cast<float>().eval();
        std::vector<int> indices;
        if (std::string(argv[7]) != "-") {
            std::stringstream stream(argv[7]); std::string entry;
            while (std::getline(stream, entry, ',')) indices.push_back(std::stoi(entry));
        }
        cosmic::BgenReader reader;
        reader.open(bgen_path); reader.initialize();
        cosmic::BgenAdKernelOptions kernel_options;
        // Keep the causal marker out of the calibration background, not out
        // of the scan. Shared KA/KD membership still holds.
        kernel_options.exclude.insert("rs000016");
        const auto kernels = cosmic::BgenAdNullBackend::build(reader, kernel_options, indices);
        const auto pair = kernels.global();
        if (pair.additive.rows() != phenotypes.rows()) throw std::runtime_error("phenotype/kernel sample size mismatch");
        dumpMatrix(output + ".ka.bin", pair.additive);
        dumpMatrix(output + ".kd.bin", pair.dominance);
        cosmic::BgenAdScanEngine::Options engine_options;
        engine_options.block_size = 37;
        cosmic::BgenAdScanEngine engine(engine_options);
        std::vector<cosmic::AdResult> results;
        bool unprepared_rejected = false;
        try { engine.scan(reader, results); } catch (const std::logic_error&) { unprepared_rejected = true; }
        if (!unprepared_rejected) throw std::runtime_error("scan accepted an unprepared null");

        cosmic::PlinkReader bed;
        cosmic::PgenReader pgen;
        if (mode == "parity") {
            bed.load(bed_path); pgen.open(pgen_path);
            if (!pgen.getRecordCapabilities().hardcall_only) throw std::runtime_error("parity PGEN not hardcall-only");
            // Read probabilities again and compare each fixed column to native
            // BED/PGEN prepared hardcall columns on exactly the selected samples.
            reader.reset(); pgen.reset();
            double max_encoding = 0.0;
            for (uint32_t j = 0; j < reader.getNumVariants(); ++j) {
                std::string c, id, a1, a2, pc, pid, ref, alt; uint32_t position, pp;
                cosmic::BgenGenotypeProbs probabilities;
                Eigen::VectorXf bd, pg;
                if (!reader.readVariantProbs(c, id, position, a1, a2, probabilities) ||
                    !pgen.readVariant(pc, pid, pp, ref, alt, pg)) throw std::runtime_error("parity variant ordering");
                bed.getSnpGenotypes(j, bd);
                const auto bgen = cosmic::BgenAdScanEngine::encode(probabilities, kernels.sample_indices);
                for (const auto* values : {&bd, &pg}) {
                    Eigen::VectorXd raw(phenotypes.rows());
                    for (int i = 0; i < raw.size(); ++i) raw[i] = (*values)[kernels.sample_indices[i]];
                    if (bgen.observed == 0) continue;
                    const auto a = cosmic::genetics::prepareAssociationMarker(raw, "additive", .01, 0);
                    const auto d = cosmic::genetics::prepareAssociationMarker(raw, "heterozygote", .01, 0);
                    max_encoding = std::max(max_encoding, (a.values - bgen.additive).cwiseAbs().maxCoeff());
                    max_encoding = std::max(max_encoding, (d.values - bgen.dominance).cwiseAbs().maxCoeff());
                }
            }
            if (max_encoding > 1e-12) throw std::runtime_error("one-hot fixed-marker column parity failed");
            std::ofstream encoding(output + ".encoding.tsv");
            encoding << std::setprecision(17) << "max_column_error\n" << max_encoding << '\n';
        }
        std::ofstream association(output + ".internal.tsv"), native_bed(output + ".bed.tsv"), native_pgen(output + ".pgen.tsv");
        const std::string header = "replicate\tchr\tsnp\ta1\ta2\tfreq\tn\tbeta_a\tse_a\tp_add\tbeta_d\tse_d\tp_dom\tp_joint\tw_joint\n";
        association << std::setprecision(17) << header;
        native_bed << std::setprecision(17) << header;
        native_pgen << std::setprecision(17) << header;
        std::ofstream nulls(output + ".null.tsv");
        nulls << std::setprecision(17) << "replicate\tsigma_a2\tsigma_d2\tsigma_e2\tconverged\tfits_in_scan\ttested\trank_deficient\tqc_filtered\tprojection_error\n";
        for (int rep = 0; rep < phenotypes.cols(); ++rep) {
            const Eigen::VectorXd y = phenotypes.col(rep);
            double a = .25, d = .20, e = .55;
            if (mode == "fitted") {
                cosmic::EigenVceAD::Options options; options.verbose = false; options.max_em_iter = 200; options.tol = 1e-5;
                engine.fitNull(kernels, y, covariates, options);
                a = engine.nullResult()->sigma_a2; d = engine.nullResult()->sigma_d2; e = engine.nullResult()->sigma_e2;
                bool rejected = false;
                try { engine.updateFixedPhenotype(y); } catch (const std::logic_error&) { rejected = true; }
                if (!rejected) throw std::runtime_error("a fitted null was reused for a different phenotype");
            } else if (rep == 0) engine.prepareFixed(kernels, y, covariates, a, d, e);
            else engine.updateFixedPhenotype(y);
            const auto before_scan = engine.nullFitCount();
            engine.scan(reader, results);
            if (engine.nullFitCount() != before_scan) throw std::runtime_error("per-SNP null refitting");
            const auto& statistics = engine.statistics();
            if (statistics.tested + statistics.rank_deficient + statistics.qc_filtered + statistics.target_filtered != reader.getNumVariants()) {
                throw std::runtime_error("variant accounting");
            }
            writeRows(association, rep, results);
            // Check the rebuilt workspace, not just variance component values.
            const Eigen::MatrixXd inverse = (a * pair.additive + d * pair.dominance + e * Eigen::MatrixXd::Identity(y.size(), y.size())).inverse();
            const Eigen::MatrixXd c = covariates.cast<double>();
            const Eigen::MatrixXd projection = inverse - inverse * c * (c.transpose() * inverse * c).inverse() * c.transpose() * inverse;
            const Eigen::MatrixXd rotated_projection = engine.rotation() * engine.varianceEigenvalues().cwiseInverse().asDiagonal() * engine.rotation().transpose();
            const double error = (projection - rotated_projection).cwiseAbs().maxCoeff();
            if (error > 1e-10) throw std::runtime_error("robust/fixed V workspace does not match GLS projection");
            nulls << rep << '\t' << a << '\t' << d << '\t' << e << "\ttrue\t"
                  << engine.nullFitCount() - before_scan << '\t' << statistics.tested << '\t' << statistics.rank_deficient
                  << '\t' << statistics.qc_filtered << '\t' << error << '\n';
            if (mode == "parity") {
                cosmic::EigenVceAD::Options fixed_options; fixed_options.verbose = false;
                cosmic::EigenVceAD fixed(fixed_options);
                fixed.prepare(y, covariates, pair.additive, pair.dominance); fixed.fixedComponents(a, d, e);
                cosmic::LmmEngineAD native; native.prepareFromVce(fixed);
                std::vector<cosmic::AdResult> native_results;
                native.scan(bed, native_results, kernels.sample_indices); writeRows(native_bed, rep, native_results);
                native.scan(pgen, native_results, kernels.sample_indices); writeRows(native_pgen, rep, native_results);
                if (rep == 0) {
                    auto check_options = engine_options;
                    check_options.block_size = 1;
                    cosmic::BgenAdScanEngine scalar_blocks(check_options);
                    scalar_blocks.prepareFixed(kernels, y, covariates, a, d, e);
                    std::vector<cosmic::AdResult> scalar_results;
                    scalar_blocks.scan(reader, scalar_results);
                    for (size_t j = 0; j < results.size(); ++j) {
                        if (scalar_results[j].rank_deficient != results[j].rank_deficient) throw std::runtime_error("block membership parity");
                        if (!results[j].rank_deficient &&
                            (std::abs(scalar_results[j].beta_add - results[j].beta_add) > 1e-10 ||
                             std::abs(scalar_results[j].beta_dom - results[j].beta_dom) > 1e-10 ||
                             std::abs(scalar_results[j].p_joint - results[j].p_joint) > 1e-10)) {
                            throw std::runtime_error("block-size invariance");
                        }
                    }
                    check_options.chromosome = "1";
                    check_options.from_bp = 500; check_options.to_bp = 2000;
                    check_options.extract = {"rs000005", "rs000016", "rs000080"};
                    check_options.exclude = {"rs000005"};
                    cosmic::BgenAdScanEngine filtered(check_options);
                    filtered.prepareFixed(kernels, y, covariates, a, d, e);
                    std::vector<cosmic::AdResult> filtered_results;
                    filtered.scan(reader, filtered_results);
                    if (filtered.statistics().tested != 1 || filtered_results[15].rank_deficient ||
                        std::abs(filtered_results[15].beta_dom - results[15].beta_dom) > 1e-10) {
                        throw std::runtime_error("target filters changed the fixed null or target result");
                    }
                    auto wrong_samples = kernels;
                    std::swap(wrong_samples.sample_ids[0], wrong_samples.sample_ids[1]);
                    cosmic::BgenAdScanEngine wrong;
                    wrong.prepareFixed(wrong_samples, y, covariates, a, d, e);
                    bool sample_rejected = false;
                    try { wrong.scan(reader, filtered_results); } catch (const std::invalid_argument&) { sample_rejected = true; }
                    if (!sample_rejected) throw std::runtime_error("null/genotype sample mismatch accepted");
                    cosmic::BgenAdScanEngine uncertified;
                    cosmic::EigenVceAD::Options zero_budget; zero_budget.verbose = false; zero_budget.max_em_iter = 0;
                    bool null_rejected = false;
                    try { uncertified.fitNull(kernels, y, covariates, zero_budget); }
                    catch (const std::runtime_error&) { null_rejected = true; }
                    if (!null_rejected) throw std::runtime_error("uncertified null accepted by scan engine");
                }
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "BGEN AD scan test failure: " << error.what() << '\n';
        return 1;
    }
}
