// E3.1 test-only driver. Not a CosmicScanner GWAS/CLI entrypoint.
#include "bgen_ad_null.h"
#include "grm_builder.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void writeMatrix(const std::string& path, const Eigen::MatrixXd& matrix) {
    std::ofstream output(path, std::ios::binary);
    if (!output) throw std::runtime_error("Cannot write test matrix " + path);
    output.write(reinterpret_cast<const char*>(matrix.data()), matrix.size() * sizeof(double));
}

int main(int argc, char** argv) {
    try {
        if (argc != 9) throw std::invalid_argument(
            "test driver: format input output indices|- chromosome|- min_info_d phenotype|- maxit");
        const std::string format = argv[1], input = argv[2], output = argv[3];
        std::vector<int> indices;
        if (std::string(argv[4]) != "-") {
            std::stringstream stream(argv[4]);
            std::string cell;
            while (std::getline(stream, cell, ',')) indices.push_back(std::stoi(cell));
        }
        const std::string chromosome = argv[5];
        const double min_info_d = std::stod(argv[6]);
        cosmic::BgenAdKernelPair kernels;
        double loco_error = 0.0;
        long long rejected_info = 0;
        long long rejected_missing = 0;
        std::vector<std::string> accepted_ids;
        if (format == "bgen" || format == "bgen-legacy") {
            cosmic::BgenReader reader;
            reader.open(input);
            reader.initialize();
            cosmic::BgenAdKernelOptions options;
            options.min_info_d = min_info_d;
            options.block_size = 37; // exercise partial blocks and chromosome grouping
            if (chromosome != "-") options.include_chromosomes.insert(chromosome);
            const auto result = cosmic::BgenAdNullBackend::build(reader, options, indices);
            kernels = result.global();
            accepted_ids = result.accepted_variant_ids;
            rejected_info = result.rejected_info;
            rejected_missing = result.rejected_missing;
            for (const auto& entry : result.chromosome_raw) {
                if (entry.second.markers == result.global_raw.markers) {
                    bool rejected = false;
                    try { result.loco(entry.first); } catch (const std::runtime_error&) { rejected = true; }
                    if (!rejected) throw std::runtime_error("empty LOCO outside not rejected");
                    continue;
                }
                const auto subtraction = result.loco(entry.first);
                auto direct_options = options;
                direct_options.exclude_chromosomes.insert(entry.first);
                const auto direct = cosmic::BgenAdNullBackend::build(reader, direct_options, indices).global();
                loco_error = std::max(loco_error,
                    (subtraction.additive - direct.additive).cwiseAbs().maxCoeff());
                loco_error = std::max(loco_error,
                    (subtraction.dominance - direct.dominance).cwiseAbs().maxCoeff());
            }
            if (loco_error >= 1e-10) throw std::runtime_error("dual LOCO subtraction/rebuild failure");
            bool duplicate_rejected = false;
            try { cosmic::BgenAdNullBackend::build(reader, options, {0, 0}); }
            catch (const std::invalid_argument&) { duplicate_rejected = true; }
            if (!duplicate_rejected) throw std::runtime_error("duplicate sample indices accepted");
            bool bounds_rejected = false;
            try { cosmic::BgenAdNullBackend::build(reader, options, {0, static_cast<int>(reader.getNumSamples())}); }
            catch (const std::invalid_argument&) { bounds_rejected = true; }
            if (!bounds_rejected) throw std::runtime_error("out-of-range sample indices accepted");
            bool unknown_rejected = false;
            try { result.loco("absent_chromosome"); }
            catch (const std::invalid_argument&) { unknown_rejected = true; }
            if (!unknown_rejected) throw std::runtime_error("unknown LOCO chromosome accepted");
        } else {
            cosmic::GrmBuilder::Options options;
            if (chromosome != "-") options.chrom_filter = chromosome;
            cosmic::PlinkReader bed;
            cosmic::PgenReader pgen;
            if (format == "bed") {
                if (!indices.empty()) throw std::invalid_argument("BED reference must be materialized subset");
                bed.load(input);
            } else if (format == "pgen") {
                pgen.open(input);
                if (!pgen.getRecordCapabilities().hardcall_only) throw std::runtime_error("non-hardcall PGEN");
            } else throw std::invalid_argument("unknown test format");
            options.algorithm = 1;
            cosmic::GrmBuilder additive(options);
            if (format == "bed") additive.compute(bed); else additive.compute(pgen, indices);
            options.algorithm = 4;
            cosmic::GrmBuilder dominance(options);
            if (format == "bed") dominance.compute(bed); else dominance.compute(pgen, indices);
            kernels.additive = additive.getGrm();
            kernels.dominance = dominance.getGrm();
            kernels.denominator_a = additive.getSum2pq();
            kernels.denominator_d = dominance.getSum2pq();
            kernels.markers = additive.getTotalSnps();
            if (kernels.markers != dominance.getTotalSnps()) throw std::runtime_error("hardcall marker sets differ");
        }

        writeMatrix(output + ".ka.bin", kernels.additive);
        writeMatrix(output + ".kd.bin", kernels.dominance);
        std::ofstream metadata(output + ".meta.tsv");
        metadata << std::setprecision(17) << "Item\tValue\nN\t" << kernels.additive.rows()
                 << "\nMARKERS\t" << kernels.markers
                 << "\nDENOMINATOR_A\t" << kernels.denominator_a
                 << "\nDENOMINATOR_D\t" << kernels.denominator_d
                 << "\nLOCO_MAX_ERROR\t" << loco_error
                 << "\nREJECTED_INFO\t" << rejected_info
                 << "\nREJECTED_MISSING\t" << rejected_missing << '\n';
        std::ofstream marker_file(output + ".markers.txt");
        for (const auto& id : accepted_ids) marker_file << id << '\n';

        if (std::string(argv[7]) != "-") {
            std::ifstream phenotype(argv[7]);
            int n = 0, replicates = 0;
            phenotype >> n >> replicates;
            if (!phenotype || n != kernels.additive.rows() || replicates < 1) {
                throw std::invalid_argument("phenotype test matrix dimensions");
            }
            Eigen::MatrixXd values(n, replicates);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < replicates; ++j) phenotype >> values(i, j);
            if (!phenotype) throw std::runtime_error("truncated phenotype test matrix");
            Eigen::MatrixXf covariates = Eigen::MatrixXf::Ones(n, 1);
            cosmic::EigenVceAD::Options options;
            options.verbose = false;
            options.max_em_iter = std::stoi(argv[8]);
            options.tol = 1e-5;
            bool aliased_rejected = false;
            try {
                cosmic::BgenAdNullBackend::fit(kernels, values.col(0),
                    Eigen::MatrixXf::Ones(n, 2), options);
            } catch (const std::invalid_argument&) { aliased_rejected = true; }
            if (!aliased_rejected) throw std::runtime_error("aliased covariate columns accepted");
            std::ofstream estimates(output + ".vce.tsv");
            estimates << std::setprecision(17)
                      << "replicate\tsigma_a2\tsigma_d2\tsigma_e2\tlogL\tconverged\titerations"
                      << "\tlikelihood_delta\tgradient_norm\tkkt_residual\traw_gradient_norm\tboundary_status"
                      << "\tselected_start\tgradient_tolerance\tlikelihood_tolerance\tinformation_rank\tsigma_d2_local_se\n";
            std::ofstream candidates(output + ".candidates.tsv");
            candidates << std::setprecision(17)
                       << "replicate\tmodel\tstart\tsigma_a2\tsigma_d2\tsigma_e2\tlogL\titerations"
                       << "\tlikelihood_delta\tgradient_norm\tkkt_residual\tface_converged\tkkt_valid\tboundary_status\n";
            for (int rep = 0; rep < replicates; ++rep) {
                cosmic::BgenAdNullResult fit;
                if (format == "bgen") {
                    fit = cosmic::BgenAdNullBackend::fit(kernels, values.col(rep), covariates, options);
                } else {
                    // Native BED/PGEN and explicit before-control use precisely
                    // the unchanged legacy EigenVceAD solver.
                    cosmic::EigenVceAD legacy(options);
                    legacy.prepare(values.col(rep), covariates, kernels.additive, kernels.dominance);
                    if (!legacy.runNullModel()) throw std::runtime_error("legacy test null failed");
                    fit.sigma_a2 = legacy.getSigma2a();
                    fit.sigma_d2 = legacy.getSigma2d();
                    fit.sigma_e2 = legacy.getSigma2e();
                    fit.log_likelihood = legacy.getLogLikelihood();
                    fit.converged = legacy.converged();
                    fit.iterations = legacy.lastIterations();
                    fit.boundary_status = "LEGACY_NOT_RECORDED";
                    fit.selected_start = "legacy";
                }
                estimates << rep << '\t' << fit.sigma_a2 << '\t' << fit.sigma_d2 << '\t'
                          << fit.sigma_e2 << '\t' << fit.log_likelihood << '\t'
                          << (fit.converged ? "true" : "false") << '\t' << fit.iterations
                          << '\t' << fit.likelihood_delta << '\t' << fit.gradient_norm << '\t' << fit.kkt_residual
                          << '\t' << fit.raw_gradient_norm << '\t' << fit.boundary_status << '\t' << fit.selected_start
                          << '\t' << fit.gradient_tolerance << '\t' << fit.likelihood_tolerance
                          << '\t' << fit.information_rank << '\t' << fit.sigma_d2_local_se << '\n';
                for (const auto& candidate : fit.candidates) {
                    candidates << rep << '\t' << candidate.model << '\t' << candidate.start
                               << '\t' << candidate.sigma_a2 << '\t' << candidate.sigma_d2 << '\t' << candidate.sigma_e2
                               << '\t' << candidate.log_likelihood << '\t' << candidate.iterations
                               << '\t' << candidate.likelihood_delta << '\t' << candidate.gradient_norm << '\t' << candidate.kkt_residual
                               << '\t' << (candidate.face_converged ? "true" : "false")
                               << '\t' << (candidate.kkt_valid ? "true" : "false") << '\t' << candidate.boundary_status << '\n';
                }
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "BGEN AD null test failure: " << error.what() << '\n';
        return 1;
    }
}
