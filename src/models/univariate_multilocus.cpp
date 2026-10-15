#include "models/univariate_model.h"

#include "association/multilocus/multilocus_engine.h"
#include "association/generalized/firth_logistic.h"
#include "grm_builder.h"
#include "lmm_engine.h"
#include "qc_utils.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace cosmic {
namespace {

std::unordered_set<std::string> readIdSet(const std::string& path) {
    std::unordered_set<std::string> ids;
    if (path.empty()) return ids;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open SNP list: " + path);
    std::string id, extra;
    while (std::getline(input, id)) {
        const auto first = id.find_first_not_of(" \t\r");
        if (first == std::string::npos) continue;
        const auto last = id.find_last_not_of(" \t\r");
        id = id.substr(first, last - first + 1);
        if (id.find_first_of(" \t") != std::string::npos)
            throw std::invalid_argument("SNP list must contain one ID per line: " + path);
        if (!ids.insert(id).second) throw std::invalid_argument("Duplicate SNP ID in list: " + id);
    }
    if (ids.empty()) throw std::invalid_argument("SNP list is empty: " + path);
    return ids;
}

std::vector<std::string> readIdList(const std::string& path) {
    const auto ids = readIdSet(path);
    std::ifstream input(path);
    std::vector<std::string> ordered;
    std::string id;
    while (std::getline(input, id)) {
        const auto first = id.find_first_not_of(" \t\r");
        if (first == std::string::npos) continue;
        const auto last = id.find_last_not_of(" \t\r");
        ordered.push_back(id.substr(first, last - first + 1));
    }
    if (ordered.size() != ids.size()) throw std::logic_error("SNP list ordering failed");
    return ordered;
}

bool inRegion(const PlinkBimInfo& variant, const UnivariateModel::Options& options) {
    return (options.scan_chr.empty() || variant.chrom == options.scan_chr) &&
           (options.scan_from_bp < 0 || variant.bp_pos >= options.scan_from_bp) &&
           (options.scan_to_bp < 0 || variant.bp_pos <= options.scan_to_bp);
}

} // namespace

void UnivariateModel::runMultiLocusModel() {
    using association::MultiLocusEngine;
    using association::MultiLocusMethod;
    const bool conditional = !options_.condition_snp_file.empty();
    const auto extract = readIdSet(options_.extract_file);
    const auto exclude = readIdSet(options_.exclude_file);
    const auto condition = conditional ? readIdList(options_.condition_snp_file)
                                       : std::vector<std::string>{};
    const size_t total_markers = is_pgen_ ? pgen_reader_->getNumVariants()
                               : is_bgen_ ? bgen_reader_->getNumVariants()
                                          : plink_reader_->getNumSnps();
    constexpr long double memory_limit = 1024.0L * 1024.0L * 1024.0L;
    const size_t estimated_markers = extract.empty() ? total_markers : std::min(total_markers, extract.size());
    const long double estimate = static_cast<long double>(keep_.n_keep) * estimated_markers * sizeof(double);
    if (estimate > memory_limit)
        throw std::runtime_error("Experimental multi-locus dense genotype matrix exceeds 1 GiB; reduce markers with --extract or use a smaller dataset");

    GenotypeBlock block;
    std::vector<PlinkBimInfo> metadata;
    std::vector<Eigen::VectorXd> columns;
    columns.reserve(total_markers);
    metadata.reserve(total_markers);
    auto accept = [&](const PlinkBimInfo& variant, const Eigen::VectorXf& dosage) {
        if (!inRegion(variant, options_)) return;
        if (!extract.empty() && !extract.count(variant.id)) return;
        if (exclude.count(variant.id)) return;
        Eigen::VectorXd genotype(keep_.n_keep);
        double sum = 0.0;
        int called = 0;
        long long hom0 = 0, het = 0, hom2 = 0;
        bool complete_hardcall = true;
        for (int i = 0; i < keep_.n_keep; ++i) {
            const double value = dosage[keep_.keep_indices[i]];
            genotype[i] = value;
            if (!std::isfinite(value)) {
                complete_hardcall = false;
                continue;
            }
            if (value < 0.0 || value > 2.0) throw std::runtime_error("Invalid genotype dosage: " + variant.id);
            sum += value;
            ++called;
            const double rounded = std::round(value);
            if (std::abs(value - rounded) > 1e-6) complete_hardcall = false;
            else if (rounded < 0.5) ++hom0;
            else if (rounded < 1.5) ++het;
            else ++hom2;
        }
        if (called == 0) return;
        const double mean = sum / called;
        const double maf = std::min(mean / 2.0, 1.0 - mean / 2.0);
        if (maf < options_.min_maf) return;
        if (options_.min_hwe > 0.0 && complete_hardcall && calculateHWE(hom0, het, hom2) < options_.min_hwe) return;
        for (int i = 0; i < genotype.size(); ++i) if (!std::isfinite(genotype[i])) genotype[i] = mean;
        if ((genotype.array() - genotype.mean()).matrix().squaredNorm() <= 1e-12) return;
        columns.push_back(std::move(genotype));
        metadata.push_back(variant);
        block.variant_ids.push_back(variant.id);
    };

    if (is_pgen_) {
        pgen_reader_->reset();
        Eigen::VectorXf dosage;
        std::string chrom, id, ref, alt;
        uint32_t position = 0;
        while (pgen_reader_->readVariant(chrom, id, position, ref, alt, dosage)) {
            accept(PlinkBimInfo{chrom, id, 0.0, static_cast<long long>(position), alt, ref}, dosage);
        }
    } else if (is_bgen_) {
        bgen_reader_->reset();
        Eigen::VectorXf dosage;
        std::string chrom, id, a1, a2;
        uint32_t position = 0;
        while (bgen_reader_->readVariant(chrom, id, position, a1, a2, dosage)) {
            accept(PlinkBimInfo{chrom, id, 0.0, static_cast<long long>(position), a1, a2}, dosage);
        }
    } else {
        const auto& bim = plink_reader_->getBimInfo();
        Eigen::VectorXf dosage;
        for (size_t j = 0; j < total_markers; ++j) {
            plink_reader_->getSnpGenotypes(j, dosage);
            accept(bim[j], dosage);
        }
    }
    if (columns.empty()) throw std::runtime_error("No testable variants for multi-locus analysis");
    for (const auto& id : condition) {
        if (std::find(block.variant_ids.begin(), block.variant_ids.end(), id) == block.variant_ids.end())
            throw std::runtime_error("Condition SNP is absent or filtered: " + id);
    }
    block.dosages.resize(keep_.n_keep, columns.size());
    for (size_t j = 0; j < columns.size(); ++j) block.dosages.col(j) = columns[j];

    Eigen::MatrixXd phenotype(keep_.n_keep, 1);
    Eigen::MatrixXd covariates(keep_.n_keep, X_.cols());
    for (int i = 0; i < keep_.n_keep; ++i) {
        phenotype(i, 0) = y_[keep_.keep_indices[i]];
        covariates.row(i) = X_.row(keep_.keep_indices[i]).cast<double>();
    }
    std::vector<Eigen::Index> identity(keep_.n_keep);
    for (int i = 0; i < keep_.n_keep; ++i) identity[i] = i;
    ScanContext context{phenotype, covariates, identity};

    if (options_.firth_logistic) {
        association::FirthLogisticTest test;
        test.prepare(context);
        test.scanBlock(block);
        const std::string output = options_.out_prefix + ".cosmic.assoc";
        std::ofstream out(output);
        if (!out) throw std::runtime_error("Cannot write Firth logistic association output: " + output);
        out << "MODEL\tCHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tITERATION\tPSEUDO_QTN\n";
        out << std::setprecision(17);
        const auto& results = test.results();
        int converged = 0;
        for (size_t j = 0; j < results.size(); ++j) {
            const auto& result = results[j];
            if (!result.converged) continue;
            ++converged;
            out << "FIRTH_LOGISTIC\t" << metadata[j].chrom << '\t' << metadata[j].id << '\t'
                << metadata[j].bp_pos << '\t' << metadata[j].alt << '\t' << metadata[j].ref << '\t'
                << columns[j].mean() / 2.0 << '\t' << keep_.n_keep << '\t' << result.beta << '\t'
                << result.se << '\t' << result.p_value << '\t' << result.iterations << "\t0\n";
        }
        if (converged == 0) throw std::runtime_error("Firth logistic produced no converged marker tests");
        std::ofstream meta(options_.out_prefix + ".generalized.meta.tsv");
        if (!meta) throw std::runtime_error("Cannot write Firth logistic metadata");
        meta << "key\tvalue\nmodel\tFIRTH_LOGISTIC\nstatus\texperimental\n"
             << "samples\t" << keep_.n_keep << "\nvariants_loaded\t" << results.size()
             << "\nvariants_converged\t" << converged << "\n";
        std::cout << "Experimental Firth logistic results written to " << output << '\n';
        return;
    }

    association::MultiLocusOptions engine_options;
    if (options_.farmcpu) {
        engine_options = association::farmCpuOptions(options_.farmcpu_max_iterations,
                                                      options_.farmcpu_pseudo_qtn);
        GrmBuilder::Options grm_options;
        grm_options.min_maf = static_cast<float>(options_.min_maf);
        grm_options.min_hwe = static_cast<float>(options_.min_hwe);
        grm_options.algorithm = options_.grm_alg;
        grm_options.thread_num = options_.n_threads;
        grm_options.extract_set = extract;
        grm_options.exclude_set = exclude;
        GrmBuilder builder(grm_options);
        if (is_pgen_) { pgen_reader_->reset(); builder.compute(*pgen_reader_); }
        else if (is_bgen_) { bgen_reader_->reset(); builder.compute(*bgen_reader_); }
        else builder.compute(*plink_reader_);
        Eigen::MatrixXd grm = DataFilter::subsetSquare(builder.getGrm(), keep_);
        EigenVce::Options vce_options;
        vce_options.use_covariate_projection = options_.use_covariate_projection;
        EigenVce vce(vce_options);
        vce.prepare(phenotype.col(0), covariates.cast<float>(), grm);
        vce.runNullModel();
        LmmEngine::Options lmm_options;
        lmm_options.thread_num = options_.n_threads;
        lmm_options.min_maf = options_.min_maf;
        lmm_options.min_hwe = options_.min_hwe;
        lmm_options.fast_mode = true;
        lmm_options.scan_from_bp = options_.scan_from_bp;
        lmm_options.scan_to_bp = options_.scan_to_bp;
        lmm_options.extract_set = extract;
        lmm_options.exclude_set = exclude;
        LmmEngine lmm(lmm_options);
        lmm.prepareFromVce(vce);
        std::vector<GwasResult> initial;
        if (is_pgen_) { pgen_reader_->reset(); lmm.scan(*pgen_reader_, initial, keep_.keep_indices, options_.scan_chr); }
        else if (is_bgen_) { bgen_reader_->reset(); lmm.scan(*bgen_reader_, initial, keep_.keep_indices, options_.scan_chr); }
        else lmm.scan(*plink_reader_, initial, keep_.keep_indices, options_.scan_chr);
        std::unordered_map<std::string, double> p_by_id;
        for (const auto& result : initial) if (std::isfinite(result.p_value)) p_by_id[result.snp_id] = result.p_value;
        engine_options.initial_p_values.reserve(block.variant_ids.size());
        for (const auto& id : block.variant_ids) {
            const auto found = p_by_id.find(id);
            engine_options.initial_p_values.push_back(found == p_by_id.end() ? 1.0 : found->second);
        }
    } else if (options_.blink) {
        engine_options = association::blinkOptions(options_.farmcpu_pseudo_qtn);
    } else {
        engine_options = options_.conditional_joint ? association::conditionalJointOptions(condition)
                                                    : association::conditionalOptions(condition);
    }

    MultiLocusEngine engine(engine_options);
    engine.prepare(context);
    engine.scanBlock(block);
    const std::string model = options_.farmcpu ? "farmcpu" : options_.blink ? "blink" :
                              options_.conditional_joint ? "conditional_joint" : "conditional";
    const std::string output = options_.out_prefix +
        (options_.conditional_joint ? ".conditional_joint.assoc" :
         conditional ? ".conditional.assoc" : ".cosmic.assoc");
    std::ofstream out(output);
    if (!out) throw std::runtime_error("Cannot write multi-locus association output: " + output);
    out << "MODEL\tCHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tITERATION\tPSEUDO_QTN\n";
    out << std::setprecision(17);
    const auto& results = engine.results();
    for (size_t j = 0; j < results.size(); ++j) {
        const auto& result = results[j];
        if (!std::isfinite(result.p_value) || result.se <= 0.0) continue;
        out << model << '\t' << metadata[j].chrom << '\t' << metadata[j].id << '\t'
            << metadata[j].bp_pos << '\t' << metadata[j].alt << '\t' << metadata[j].ref << '\t'
            << columns[j].mean() / 2.0 << '\t' << result.n_obs << '\t' << result.beta << '\t'
            << result.se << '\t' << result.p_value << '\t' << result.iteration << '\t'
            << (result.selected_pseudo_qtn ? 1 : 0) << '\n';
    }
    std::ofstream meta(options_.out_prefix + ".multilocus.meta.tsv");
    if (!meta) throw std::runtime_error("Cannot write multi-locus metadata");
    meta << "key\tvalue\nmodel\t" << model << "\nstatus\texperimental\n"
         << "samples\t" << keep_.n_keep << "\nvariants\t" << results.size() << "\n"
         << "selected_pseudo_qtn\t" << engine.selectedIndices().size() << "\n"
         << "dense_genotype_bytes\t" << static_cast<unsigned long long>(estimate) << "\n";
    if (conditional) {
        meta << "condition_file\t" << options_.condition_snp_file << "\ncondition_count\t"
             << condition.size() << "\n";
    } else {
        meta << "max_iterations\t" << engine_options.max_iterations << "\nmax_pseudo_qtn\t"
             << engine_options.max_pseudo_qtn << "\nld_r2\t" << engine_options.ld_r2 << "\n";
    }
    std::cout << "Experimental " << model << " association results written to " << output << '\n';
}

} // namespace cosmic
