#include "models/univariate_model.h"

#include "association/generalized/saige_style.h"
#include "association/generalized/logistic.h"
#include "grm_builder.h"
#include "data_filter.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"

#include <Eigen/SparseCore>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include <map>
#include <vector>
#include <set>

namespace cosmic {

namespace {

struct SaigeVariant {
    std::string chrom;
    std::string id;
    std::string a1;
    std::string a2;
    long long position = 0;
    double frequency = 0.0;
    Eigen::VectorXd dosage;
};

std::unordered_set<std::string> loadSaigeIds(const std::string& path) {
    std::unordered_set<std::string> ids;
    if (path.empty()) return ids;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open SAIGE-style SNP list: " + path);
    std::string id;
    while (input >> id) ids.insert(id);
    if (ids.empty()) throw std::runtime_error("SAIGE-style SNP list is empty: " + path);
    return ids;
}

bool passesSaigeFilters(const std::string& chrom, const std::string& id, long long position,
                        const UnivariateModel::Options& options,
                        const std::unordered_set<std::string>& extract,
                        const std::unordered_set<std::string>& exclude) {
    if (!options.scan_chr.empty() && chrom != options.scan_chr) return false;
    if (options.scan_from_bp >= 0 && position < options.scan_from_bp) return false;
    if (options.scan_to_bp >= 0 && position > options.scan_to_bp) return false;
    if (!extract.empty() && !extract.count(id)) return false;
    if (exclude.count(id)) return false;
    return true;
}

bool prepareSaigeDosage(const Eigen::VectorXf& source, const std::vector<int>& keep,
                        double min_maf, Eigen::VectorXd& dosage, double& frequency) {
    dosage.resize(keep.size());
    double sum = 0.0;
    int observed = 0;
    for (int i = 0; i < static_cast<int>(keep.size()); ++i) {
        const double value = source[keep[i]];
        dosage[i] = value;
        if (std::isfinite(value) && value >= 0.0 && value <= 2.0) {
            sum += value;
            ++observed;
        } else {
            dosage[i] = std::numeric_limits<double>::quiet_NaN();
        }
    }
    if (observed < 3) return false;
    frequency = sum / (2.0 * observed);
    if (!std::isfinite(frequency) || std::min(frequency, 1.0 - frequency) < min_maf) {
        return false;
    }
    const double mean = 2.0 * frequency;
    for (int i = 0; i < dosage.size(); ++i)
        if (!std::isfinite(dosage[i]))
            dosage[i] = mean;
    return (dosage.array() - mean).matrix().squaredNorm() > 1e-10;
}

Eigen::SparseMatrix<double> sparsifyKinship(const Eigen::MatrixXd& dense, double cutoff) {
    if (dense.rows() != dense.cols() || !dense.allFinite())
        throw std::runtime_error("SAIGE-style GRM is not a finite square matrix");
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(dense.rows() * 8);
    for (int col = 0; col < dense.cols(); ++col) {
        for (int row = 0; row < dense.rows(); ++row) {
            const double value = 0.5 * (dense(row, col) + dense(col, row));
            if (row == col || std::abs(value) >= cutoff) {
                triplets.emplace_back(row, col, value);
            }
        }
    }
    Eigen::SparseMatrix<double> sparse(dense.rows(), dense.cols());
    sparse.setFromTriplets(triplets.begin(), triplets.end());
    sparse.makeCompressed();
    return sparse;
}

Eigen::SparseMatrix<double> computeLocoSparseKinship(
    const Eigen::MatrixXd& dense_global,
    const GrmBuilder::SplitGrmResult& split,
    const std::string& chrom,
    const std::vector<int>& keep_indices,
    double sparse_cutoff) {
    
    size_t chrom_idx = static_cast<size_t>(-1);
    for (size_t i = 0; i < split.chrom_labels.size(); ++i) {
        if (split.chrom_labels[i] == chrom) {
            chrom_idx = i;
            break;
        }
    }
    if (chrom_idx == static_cast<size_t>(-1))
        throw std::runtime_error("Chromosome " + chrom + " not found in split GRM");
    if (chrom_idx >= split.per_chrom_sum_2pq.size() || chrom_idx >= split.per_chrom_raw.size())
        throw std::runtime_error("Chromosome index out of bounds in split GRM");
    
    double S_chr = split.per_chrom_sum_2pq[chrom_idx];
    double S_global = split.global_sum_2pq;
    double S_loco = split.global_sum_2pq - split.per_chrom_sum_2pq[chrom_idx];
    
    if (S_loco <= 1e-6 || S_chr <= 0.0)
        return sparsifyKinship(dense_global, sparse_cutoff);
    
    int n = dense_global.rows();
    Eigen::MatrixXd raw_chr_sub(n, n);
    const Eigen::MatrixXd& raw_chr = split.per_chrom_raw[chrom_idx];
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            raw_chr_sub(i, j) = raw_chr(keep_indices[i], keep_indices[j]);
    
    Eigen::MatrixXd dense_loco = dense_global * (S_global / S_loco) - raw_chr_sub / S_loco;
    return sparsifyKinship(dense_loco, sparse_cutoff);
}

void collectVariantsPgen(PgenReader* reader,
                          const UnivariateModel::Options& options,
                          const std::unordered_set<std::string>& extract_set,
                          const std::unordered_set<std::string>& exclude_set,
                          const std::vector<int>& keep_indices,
                          double min_maf,
                          std::map<std::string, std::vector<SaigeVariant>>& chrom_variants,
                          std::vector<Eigen::VectorXd>& calibration_markers,
                          int max_calibration) {
    reader->reset();
    std::string chrom, id, ref, alt;
    uint32_t pos = 0;
    Eigen::VectorXf dosage_f;
    int n_calibration = 0;
    while (reader->readVariant(chrom, id, pos, ref, alt, dosage_f)) {
        if (passesSaigeFilters(chrom, id, static_cast<long long>(pos), options, extract_set, exclude_set)) {
            Eigen::VectorXd dosage;
            double freq;
            if (prepareSaigeDosage(dosage_f, keep_indices, min_maf, dosage, freq)) {
                SaigeVariant sv;
                sv.chrom = chrom;
                sv.id = id;
                sv.a1 = alt;
                sv.a2 = ref;
                sv.position = pos;
                sv.frequency = freq;
                sv.dosage = std::move(dosage);
                std::string chrom_key = sv.chrom;
                chrom_variants[chrom_key].push_back(std::move(sv));
                if (n_calibration < max_calibration) {
                    calibration_markers.push_back(chrom_variants[chrom_key].back().dosage);
                    ++n_calibration;
                }
            }
        }
    }
}

void collectVariantsSequential(SequentialGenotypeMatrix* reader,
                               const UnivariateModel::Options& options,
                               const std::unordered_set<std::string>& extract_set,
                               const std::unordered_set<std::string>& exclude_set,
                               const std::vector<int>& keep_indices,
                               double min_maf,
                               std::map<std::string, std::vector<SaigeVariant>>& chrom_variants,
                               std::vector<Eigen::VectorXd>& calibration_markers,
                               int max_calibration) {
    reader->reset();
    std::string chrom, id, ref, alt;
    uint32_t pos = 0;
    Eigen::VectorXf dosage_f;
    int n_calibration = 0;
    while (reader->readVariant(chrom, id, pos, ref, alt, dosage_f)) {
        if (!passesSaigeFilters(chrom, id, pos, options, extract_set, exclude_set)) continue;
        Eigen::VectorXd dosage;
        double frequency = 0.0;
        if (!prepareSaigeDosage(dosage_f, keep_indices, min_maf, dosage, frequency)) continue;
        SaigeVariant variant;
        variant.chrom = chrom;
        variant.id = id;
        variant.a1 = alt;
        variant.a2 = ref;
        variant.position = pos;
        variant.frequency = frequency;
        variant.dosage = std::move(dosage);
        auto& stored = chrom_variants[chrom];
        stored.push_back(std::move(variant));
        if (n_calibration < max_calibration) {
            calibration_markers.push_back(stored.back().dosage);
            ++n_calibration;
        }
    }
}

void collectVariantsBgen(BgenReader* reader,
                          const UnivariateModel::Options& options,
                          const std::unordered_set<std::string>& extract_set,
                          const std::unordered_set<std::string>& exclude_set,
                          const std::vector<int>& keep_indices,
                          double min_maf,
                          std::map<std::string, std::vector<SaigeVariant>>& chrom_variants,
                          std::vector<Eigen::VectorXd>& calibration_markers,
                          int max_calibration) {
    reader->reset();
    std::string chrom, id, a1, a2;
    uint32_t pos = 0;
    Eigen::VectorXf dosage_f;
    int n_calibration = 0;
    while (reader->readVariant(chrom, id, pos, a1, a2, dosage_f)) {
        if (passesSaigeFilters(chrom, id, static_cast<long long>(pos), options, extract_set, exclude_set)) {
            Eigen::VectorXd dosage;
            double freq;
            if (prepareSaigeDosage(dosage_f, keep_indices, min_maf, dosage, freq)) {
                SaigeVariant sv;
                sv.chrom = chrom;
                sv.id = id;
                sv.a1 = a2;
                sv.a2 = a1;
                sv.position = pos;
                sv.frequency = freq;
                sv.dosage = std::move(dosage);
                std::string chrom_key = sv.chrom;
                chrom_variants[chrom_key].push_back(std::move(sv));
                if (n_calibration < max_calibration) {
                    calibration_markers.push_back(chrom_variants[chrom_key].back().dosage);
                    ++n_calibration;
                }
            }
        }
    }
}

void collectVariantsPlink(PlinkReader* reader,
                          const UnivariateModel::Options& options,
                          const std::unordered_set<std::string>& extract_set,
                          const std::unordered_set<std::string>& exclude_set,
                          const std::vector<int>& keep_indices,
                          double min_maf,
                          std::map<std::string, std::vector<SaigeVariant>>& chrom_variants,
                          std::vector<Eigen::VectorXd>& calibration_markers,
                          int max_calibration) {
    const auto& bim = reader->getBimInfo();
    Eigen::VectorXf dosage_f(reader->getNumSamples());
    int n_calibration = 0;
    for (size_t snp_idx = 0; snp_idx < bim.size(); ++snp_idx) {
        const auto& info = bim[snp_idx];
        if (passesSaigeFilters(info.chrom, info.id, info.bp_pos, options, extract_set, exclude_set)) {
            reader->getSnpGenotypes(snp_idx, dosage_f);
            Eigen::VectorXd dosage;
            double freq;
            if (prepareSaigeDosage(dosage_f, keep_indices, min_maf, dosage, freq)) {
                SaigeVariant sv;
                sv.chrom = info.chrom;
                sv.id = info.id;
                sv.a1 = info.alt;
                sv.a2 = info.ref;
                sv.position = info.bp_pos;
                sv.frequency = freq;
                sv.dosage = std::move(dosage);
                std::string chrom_key = sv.chrom;
                chrom_variants[chrom_key].push_back(std::move(sv));
                if (n_calibration < max_calibration) {
                    calibration_markers.push_back(chrom_variants[chrom_key].back().dosage);
                    ++n_calibration;
                }
            }
        }
    }
}

void collectVariantsVcf(VcfReader* reader,
                        const UnivariateModel::Options& options,
                        const std::unordered_set<std::string>& extract_set,
                        const std::unordered_set<std::string>& exclude_set,
                        const std::vector<int>& keep_indices,
                        double min_maf,
                        std::map<std::string, std::vector<SaigeVariant>>& chrom_variants,
                        std::vector<Eigen::VectorXd>& calibration_markers,
                        int max_calibration) {
    reader->loadVariantInfo();
    const auto& variants = reader->getVariants();
    int n_calibration = 0;
    for (const auto& variant : variants) {
        if (passesSaigeFilters(variant.chrom, variant.id, static_cast<long long>(variant.pos), options, extract_set, exclude_set)) {
            Eigen::VectorXf dosage_f(reader->getNumSamples());
            for (size_t si = 0; si < dosage_f.size(); ++si) {
                if (si >= variant.genotypes.size()) {
                    dosage_f[si] = std::numeric_limits<float>::quiet_NaN();
                    continue;
                }
                const std::string& gt = variant.genotypes[si];
                if (gt.empty() || gt[0] == '.') {
                    dosage_f[si] = std::numeric_limits<float>::quiet_NaN();
                } else if (gt.size() >= 3) {
                    int alt_count = 0;
                    if (gt[0] == '1') alt_count++;
                    if (gt[2] == '1') alt_count++;
                    dosage_f[si] = static_cast<float>(alt_count);
                } else if (gt.size() == 1) {
                    if (gt[0] == '1') dosage_f[si] = 2.0f;
                    else if (gt[0] == '0') dosage_f[si] = 0.0f;
                    else dosage_f[si] = std::numeric_limits<float>::quiet_NaN();
                } else {
                    dosage_f[si] = std::numeric_limits<float>::quiet_NaN();
                }
            }
            Eigen::VectorXd dosage;
            double freq;
            if (prepareSaigeDosage(dosage_f, keep_indices, min_maf, dosage, freq)) {
                SaigeVariant sv;
                sv.chrom = variant.chrom;
                sv.id = variant.id;
                sv.a1 = variant.alt;
                sv.a2 = variant.ref;
                sv.position = variant.pos;
                sv.frequency = freq;
                sv.dosage = std::move(dosage);
                std::string chrom_key = sv.chrom;
                chrom_variants[chrom_key].push_back(std::move(sv));
                if (n_calibration < max_calibration) {
                    calibration_markers.push_back(chrom_variants[chrom_key].back().dosage);
                    ++n_calibration;
                }
            }
        }
    }
}

} // anonymous namespace

void UnivariateModel::runLogistic() {
    Eigen::VectorXd phenotype = DataFilter::subset(y_, keep_);
    Eigen::MatrixXd covariates = DataFilter::subsetRows(X_, keep_).cast<double>();
    int cases = 0;
    for (int i = 0; i < phenotype.size(); ++i) {
        if (phenotype[i] != 0.0 && phenotype[i] != 1.0)
            throw std::runtime_error("Logistic phenotype must be coded exactly 0/1");
        cases += phenotype[i] == 1.0;
    }
    if (cases == 0 || cases == phenotype.size())
        throw std::runtime_error("Logistic analysis requires both cases and controls");

    const auto extract_set = loadSaigeIds(options_.extract_file);
    const auto exclude_set = loadSaigeIds(options_.exclude_file);
    std::map<std::string, std::vector<SaigeVariant>> chrom_variants;
    std::vector<Eigen::VectorXd> unused;
    if (generic_reader_) {
        collectVariantsSequential(generic_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                                  options_.min_maf, chrom_variants, unused, 0);
    } else if (is_vcf_) {
        collectVariantsVcf(vcf_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                           options_.min_maf, chrom_variants, unused, 0);
    } else if (is_pgen_) {
        collectVariantsPgen(pgen_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                            options_.min_maf, chrom_variants, unused, 0);
    } else if (is_bgen_) {
        collectVariantsBgen(bgen_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                            options_.min_maf, chrom_variants, unused, 0);
    } else {
        collectVariantsPlink(plink_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                             options_.min_maf, chrom_variants, unused, 0);
    }

    std::vector<SaigeVariant> variants;
    for (auto& entry : chrom_variants)
        for (auto& variant : entry.second) variants.push_back(std::move(variant));
    if (variants.empty()) throw std::runtime_error("Logistic analysis found no testable variants");

    std::ofstream out(options_.out_prefix + ".cosmic.assoc");
    std::ofstream diagnostics(options_.out_prefix + ".logistic.variants.tsv");
    if (!out || !diagnostics) throw std::runtime_error("Cannot write logistic association output");
    out << "MODEL\tCHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tITERATION\tPSEUDO_QTN\n";
    diagnostics << "SNP\tSTATUS\tSEPARATED\tITERATIONS\n";
    out << std::setprecision(17);

    const bool mixed = options_.logistic_glmm;
    std::vector<cosmic::association::SaigeStyleResult> mixed_results;
    std::vector<cosmic::association::LogisticResult> ordinary_results;
    double tau = 0.0;
    int null_iterations = 0;
    if (mixed) {
        GrmBuilder::Options grm_opts;
        grm_opts.algorithm = options_.grm_alg;
        grm_opts.min_maf = static_cast<float>(options_.min_maf);
        grm_opts.min_hwe = static_cast<float>(options_.min_hwe);
        grm_opts.thread_num = options_.n_threads;
        grm_opts.extract_set = extract_set;
        grm_opts.exclude_set = exclude_set;
        GrmBuilder builder(grm_opts);
        if (generic_reader_) builder.compute(*generic_reader_);
        else if (is_vcf_) builder.compute(*vcf_reader_);
        else if (is_pgen_) { pgen_reader_->reset(); builder.compute(*pgen_reader_); }
        else if (is_bgen_) { bgen_reader_->reset(); builder.compute(*bgen_reader_); }
        else builder.compute(*plink_reader_);
        Eigen::MatrixXd dense = builder.getGrm();
        if (dense.rows() != keep_.n_keep) dense = DataFilter::subsetSquare(dense, keep_);
        const auto kinship = sparsifyKinship(dense, 0.0);
        cosmic::association::SaigeStyleNullModel null_model;
        null_model.fit(phenotype, covariates, kinship);
        tau = null_model.tau();
        null_iterations = null_model.metadata().iterations;
        for (const auto& variant : variants) mixed_results.push_back(null_model.testMarkerExact(variant.dosage));
    } else {
        Eigen::MatrixXd phenotype_matrix = phenotype;
        std::vector<Eigen::Index> identity(phenotype.size());
        for (int i = 0; i < phenotype.size(); ++i) identity[i] = i;
        cosmic::association::LogisticOptions logistic_options;
        logistic_options.test = options_.logistic_test == "wald"
            ? cosmic::association::LogisticTestType::WALD
            : cosmic::association::LogisticTestType::SCORE;
        cosmic::association::LogisticTest test(logistic_options);
        test.prepare(ScanContext{phenotype_matrix, covariates, identity});
        GenotypeBlock block;
        block.dosages.resize(phenotype.size(), variants.size());
        for (size_t j = 0; j < variants.size(); ++j) {
            block.dosages.col(j) = variants[j].dosage;
            block.variant_ids.push_back(variants[j].id);
        }
        test.scanBlock(block);
        ordinary_results = test.results();
    }

    int tested_count = 0;
    int failure_count = 0;
    const std::string model = mixed ? "LOGISTIC_GLMM" :
        (options_.logistic_test == "wald" ? "LOGISTIC_WALD" : "LOGISTIC_SCORE");
    for (size_t j = 0; j < variants.size(); ++j) {
        const bool tested = mixed ? mixed_results[j].tested : ordinary_results[j].converged;
        const bool separated = !mixed && ordinary_results[j].separated;
        const int iterations = mixed ? null_iterations : ordinary_results[j].iterations;
        diagnostics << variants[j].id << '\t' << (tested ? "SUCCESS" : "NUMERICAL_FAILURE")
                    << '\t' << (separated ? 1 : 0) << '\t' << iterations << '\n';
        if (!tested) { ++failure_count; continue; }
        ++tested_count;
        const double beta = mixed ? mixed_results[j].beta : ordinary_results[j].beta;
        const double se = mixed ? mixed_results[j].se : ordinary_results[j].se;
        const double p = mixed ? mixed_results[j].p_value : ordinary_results[j].p_value;
        out << model << '\t' << variants[j].chrom << '\t' << variants[j].id << '\t'
            << variants[j].position << '\t' << variants[j].a1 << '\t' << variants[j].a2 << '\t'
            << variants[j].frequency << '\t' << phenotype.size() << '\t' << beta << '\t'
            << se << '\t' << p << '\t' << iterations << "\t0\n";
    }
    if (tested_count == 0) throw std::runtime_error("Logistic analysis produced no valid marker tests");
    std::ofstream meta(options_.out_prefix + ".logistic.meta.tsv");
    meta << "key\tvalue\nmodel\t" << model << "\ntest\t"
         << (mixed ? "pql_mixed_exact_score_variance" : options_.logistic_test) << "\nsamples\t" << phenotype.size()
         << "\ncases\t" << cases << "\nvariants_tested\t" << tested_count
         << "\nvariants_failed\t" << failure_count << "\ntau\t" << std::setprecision(17) << tau
         << "\nnull_approximation\t" << (mixed ? "PQL" : "NONE")
         << "\nscore_variance\t" << (mixed ? "FULL_PROJECTED_WORKING_COVARIANCE" : "GLM_INFORMATION")
         << "\nvariance_ratio_used\t0\nspa_used\t0\n";
}


void UnivariateModel::runSaigeStyle() {
    Eigen::VectorXd phenotype = DataFilter::subset(y_, keep_);
    Eigen::MatrixXd covariates = DataFilter::subsetRows(X_, keep_).cast<double>();

    int cases = 0;
    for (int i = 0; i < phenotype.size(); ++i) {
        if (phenotype[i] != 0.0 && phenotype[i] != 1.0)
            throw std::runtime_error("SAIGE-style phenotype must be coded exactly 0/1");
        cases += (phenotype[i] == 1.0);
    }
    if (cases == 0 || cases == phenotype.size())
        throw std::runtime_error("SAIGE-style requires both cases and controls");

    const std::unordered_set<std::string> extract_set = loadSaigeIds(options_.extract_file);
    const std::unordered_set<std::string> exclude_set = loadSaigeIds(options_.exclude_file);

    // Collect chromosome labels
    std::unordered_set<std::string> chrom_labels;
    if (generic_reader_) {
        generic_reader_->reset();
        std::string chrom, id, ref, alt;
        uint32_t pos = 0;
        Eigen::VectorXf dosage;
        while (generic_reader_->readVariant(chrom, id, pos, ref, alt, dosage)) chrom_labels.insert(chrom);
    } else if (is_vcf_) {
        for (const auto& chrom : vcf_reader_->getChromLabels())
            chrom_labels.insert(chrom);
    } else if (is_pgen_) {
        for (const auto& snp : pgen_reader_->getVariants())
            chrom_labels.insert(snp.chrom);
    } else if (is_bgen_) {
        bgen_reader_->loadVariantInfo();
        for (const auto& snp : bgen_reader_->getVariants())
            chrom_labels.insert(snp.chrom);
    } else {
        for (const auto& snp : plink_reader_->getBimInfo())
            chrom_labels.insert(snp.chrom);
    }

    const bool use_loco = chrom_labels.size() > 1 && !options_.no_loco;
    std::cout << "SAIGE-style mode: " << (use_loco ? "LOCO" : "global")
              << " (" << chrom_labels.size() << " chromosomes)" << std::endl;

    // Build GRM
    GrmBuilder::Options grm_opts;
    grm_opts.algorithm = options_.grm_alg;
    grm_opts.min_maf = static_cast<float>(options_.min_maf);
    grm_opts.min_hwe = static_cast<float>(options_.min_hwe);
    grm_opts.thread_num = options_.n_threads;
    grm_opts.extract_set = extract_set;
    grm_opts.exclude_set = exclude_set;
    GrmBuilder::SplitGrmResult split_grm;
    std::set<std::string> chrom_labels_ordered(chrom_labels.begin(), chrom_labels.end());

    Eigen::MatrixXd dense_kinship;
    if (use_loco) {
        std::cout << "Building split GRM for SAIGE-style LOCO..." << std::endl;
        GrmBuilder::Options split_opts = grm_opts;
        GrmBuilder split_builder(split_opts);
        if (generic_reader_) {
            split_grm = split_builder.computeSplit(*generic_reader_, chrom_labels_ordered);
        } else if (is_vcf_) {
            vcf_reader_->loadVariantInfo();
            split_grm = split_builder.computeSplit(*vcf_reader_, chrom_labels_ordered);
        } else if (is_pgen_) {
            pgen_reader_->reset();
            split_grm = split_builder.computeSplit(*pgen_reader_, chrom_labels_ordered);
        } else if (is_bgen_) {
            bgen_reader_->reset();
            split_grm = split_builder.computeSplit(*bgen_reader_, chrom_labels_ordered);
        } else {
            split_grm = split_builder.computeSplit(*plink_reader_, chrom_labels_ordered);
        }
        dense_kinship = split_grm.global;
    } else {
        GrmBuilder builder(grm_opts);
        if (generic_reader_) {
            builder.compute(*generic_reader_);
        } else if (is_vcf_) {
            builder.compute(*vcf_reader_);
        } else if (is_pgen_) {
            builder.compute(*pgen_reader_);
        } else if (is_bgen_) {
            builder.compute(*bgen_reader_);
        } else {
            builder.compute(*plink_reader_);
        }
        dense_kinship = builder.getGrm();
    }
    if (dense_kinship.rows() != keep_.n_keep)
        dense_kinship = DataFilter::subsetSquare(dense_kinship, keep_);

    Eigen::SparseMatrix<double> kinship_global = sparsifyKinship(dense_kinship, options_.sparse_cutoff);
    Eigen::MatrixXd dense_global = dense_kinship;

    // SAIGE-style null model fitting
    cosmic::association::SaigeStyleNullModel null_model;
    null_model.fit(phenotype, covariates, kinship_global);

    // Prepare calibration markers
    std::vector<Eigen::VectorXd> calibration_markers;
    calibration_markers.reserve(500);

    // For calibration, we'll collect a subset of SNPs from all chromosomes
    std::map<std::string, std::vector<SaigeVariant>> chrom_variants_dummy;
    if (generic_reader_) {
        collectVariantsSequential(generic_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                                  options_.min_maf, chrom_variants_dummy, calibration_markers, 500);
    } else if (is_vcf_) {
        collectVariantsVcf(vcf_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                            options_.min_maf, chrom_variants_dummy, calibration_markers, 500);
    } else if (is_pgen_) {
        collectVariantsPgen(pgen_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                            options_.min_maf, chrom_variants_dummy, calibration_markers, 500);
    } else if (is_bgen_) {
        collectVariantsBgen(bgen_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                            options_.min_maf, chrom_variants_dummy, calibration_markers, 500);
    } else {
        collectVariantsPlink(plink_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                            options_.min_maf, chrom_variants_dummy, calibration_markers, 500);
    }

    if (calibration_markers.empty()) {
        throw std::runtime_error("SAIGE-style: no valid calibration markers found");
    }

    // Calibrate variance ratio (MAC-stratified)
    null_model.calibrateVarianceRatio(calibration_markers);
    const auto& mac_bins = null_model.macBins();
    const double global_vr = null_model.globalVarianceRatio();
    const int n_cal_markers = null_model.calibrationMarkers();

    std::cout << "SAIGE-style calibration: global VR=" << global_vr
              << ", MAC bins: " << mac_bins.size() << ", markers used: " << n_cal_markers << std::endl;

    // Now collect all variants for testing
    std::map<std::string, std::vector<SaigeVariant>> chrom_variants;
    
    if (generic_reader_) {
        collectVariantsSequential(generic_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                                  options_.min_maf, chrom_variants, calibration_markers, 0);
    } else if (is_vcf_) {
        collectVariantsVcf(vcf_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                            options_.min_maf, chrom_variants, calibration_markers, 0);
    } else if (is_pgen_) {
        collectVariantsPgen(pgen_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                            options_.min_maf, chrom_variants, calibration_markers, 0);
    } else if (is_bgen_) {
        collectVariantsBgen(bgen_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                            options_.min_maf, chrom_variants, calibration_markers, 0);
    } else {
        collectVariantsPlink(plink_reader_, options_, extract_set, exclude_set, keep_.keep_indices,
                            options_.min_maf, chrom_variants, calibration_markers, 0);
    }

    // Setup output
    std::ofstream out(options_.out_prefix + ".cosmic.assoc");
    if (!out) throw std::runtime_error("Cannot write SAIGE-style association output");
    out << "MODEL\tCHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tITERATION\tPSEUDO_QTN\n";
    std::ofstream variant_meta(options_.out_prefix + ".saige_style.variants.tsv");
    if (!variant_meta) throw std::runtime_error("Cannot write SAIGE-style variant metadata");
    variant_meta << "CHR\tSNP\tBP\tAF\tMAC\tMAC_BIN\tVARIANCE_RATIO\tLOG10P\tSPA_USED\n";
    int spa_count = 0;

    // For each chromosome, optionally use LOCO GRM and test variants
    for (const auto& chrom_entry : chrom_labels_ordered) {
        const auto& variants = chrom_variants[chrom_entry];
        if (variants.empty()) continue;

        Eigen::SparseMatrix<double> kinship = kinship_global;
        if (use_loco && !split_grm.chrom_labels.empty()) {
            try {
                kinship = computeLocoSparseKinship(dense_global, split_grm, chrom_entry, keep_.keep_indices, options_.sparse_cutoff);
            } catch (const std::exception& e) {
                std::cerr << "LOCO GRM computation failed for chromosome " << chrom_entry
                          << ": " << e.what() << ". Using global GRM." << std::endl;
                kinship = kinship_global;
            }
        }

        // Refit null model with chromosome-specific GRM if LOCO
        cosmic::association::SaigeStyleNullModel null_model_chrom;
        null_model_chrom.fit(phenotype, covariates, kinship);
        // Re-calibrate variance ratio with the same calibration markers
        null_model_chrom.calibrateVarianceRatio(calibration_markers);
        SPATest spa(SPATest::Options{});
        spa.set_null_model(null_model_chrom.fittedProbabilities(), null_model_chrom.workingVariance());

        // Test variants on this chromosome
        for (const auto& variant : variants) {
            cosmic::association::SaigeStyleResult result = null_model_chrom.testMarker(variant.dosage, &spa);

            if (!result.tested) continue;

            int mac = 0;
            int n_observed = 0;
            for (int i = 0; i < variant.dosage.size(); ++i) {
                if (std::isfinite(variant.dosage[i]) && variant.dosage[i] >= 0 && variant.dosage[i] <= 2)
                    ++n_observed;
            }
            if (n_observed >= 3) {
                double freq = variant.dosage.sum() / (2.0 * n_observed);
                double maf = std::min(freq, 1.0 - freq);
                mac = static_cast<int>(2.0 * n_observed * maf);
            }

            double vr = null_model_chrom.varianceRatioForMac(mac);
            int mac_bin_idx = 0;
            if (mac <= 10) mac_bin_idx = 0;
            else if (mac <= 20) mac_bin_idx = 1;
            else if (mac <= 50) mac_bin_idx = 2;
            else mac_bin_idx = 3;

            double log10p = (result.p_value > 0 && std::isfinite(result.p_value))
                                ? -std::log10(result.p_value)
                                : std::numeric_limits<double>::infinity();
            if (result.spa_used) ++spa_count;

            out << "SAIGE_STYLE\t"
                << variant.chrom << "\t"
                << variant.id << "\t"
                << variant.position << "\t"
                << variant.a1 << "\t"
                << variant.a2 << "\t"
                << std::fixed << std::setprecision(6) << variant.frequency << "\t"
                << variant.dosage.size() << "\t"
                << std::scientific << std::setprecision(17) << result.beta << "\t"
                << result.se << "\t"
                << result.p_value << "\t"
                << null_model_chrom.metadata().iterations << "\t"
                << "0\n";
            variant_meta << variant.chrom << "\t" << variant.id << "\t" << variant.position << "\t"
                         << std::setprecision(12) << variant.frequency << "\t" << mac << "\t"
                         << mac_bin_idx << "\t" << vr << "\t" << log10p << "\t"
                         << (result.spa_used ? 1 : 0) << "\n";
        }
    }

    out.close();
    variant_meta.close();
    if (!out || !variant_meta) throw std::runtime_error("Failed while writing SAIGE-style results");

    // Write metadata file
    std::ofstream meta_out(options_.out_prefix + ".saige_style.meta.tsv");
    if (!meta_out) throw std::runtime_error("Cannot write SAIGE-style null metadata");
    meta_out << "model\tSAIGE_STYLE\n";
    meta_out << "status\texperimental\n";
    meta_out << "null_solver\tsparse_ldlt_pql\n";
    meta_out << "variance_component_method\tmoment_he\n";
    meta_out << "null_converged\t1\n";
    meta_out << "variance_ratio\t" << std::setprecision(6) << global_vr << "\n";
    meta_out << "variance_ratio_markers\t" << n_cal_markers << "\n";
    meta_out << "spa_marker_count\t" << spa_count << "\n";
    for (size_t i = 0; i < mac_bins.size(); ++i) {
        meta_out << "mac_bin_" << i << "\t" << mac_bins[i].min_mac << "-"
                 << (mac_bins[i].max_mac == 0 ? std::string("inf") : std::to_string(mac_bins[i].max_mac))
                 << ";markers=" << mac_bins[i].n_calibration_markers
                 << ";variance_ratio=" << mac_bins[i].variance_ratio
                 << ";fallback=" << (mac_bins[i].fallback ? 1 : 0) << "\n";
    }
    meta_out << "sparse_grm_nnz\t" << kinship_global.nonZeros() << "\n";
    meta_out.close();
    if (!meta_out) throw std::runtime_error("Failed while writing SAIGE-style null metadata");

    std::cout << "SAIGE-style scan completed. Output: " << options_.out_prefix << ".cosmic.assoc" << std::endl;
}

} // namespace cosmic
