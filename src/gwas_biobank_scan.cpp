#include "gwas_scanner.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace cosmic {

Eigen::MatrixXd orthonormal_covariates(const Eigen::MatrixXd& C);
std::unordered_set<std::string> load_snp_id_set(const std::string& path);
void log_stage_time(const std::string& label, const std::chrono::steady_clock::time_point& start);
std::vector<std::string> sort_chrom_labels(const std::set<std::string>& chroms_set);

void write_biobank_predictions(const std::string& path,
                               const PlinkReader& reader,
                               const std::vector<int>& keep_indices,
                               const std::vector<Eigen::VectorXd>& preds);
void write_biobank_predictions(const std::string& path,
                               const std::vector<std::string>& sample_ids,
                               const std::vector<int>& keep_indices,
                               const std::vector<Eigen::VectorXd>& preds);
std::vector<Eigen::VectorXd> read_biobank_predictions(const std::string& path,
                                                      const PlinkReader& reader,
                                                      const std::vector<int>& keep_indices);
std::vector<Eigen::VectorXd> read_biobank_predictions_by_sample_id(const std::string& path,
                                                                   const std::vector<std::string>& sample_ids,
                                                                   const std::vector<int>& keep_indices);
void scan_pgen_biobank_linear_parallel(const PgenReader& reader,
                                       uint32_t start_idx,
                                       uint32_t count,
                                       const Eigen::VectorXd& y,
                                       const Eigen::MatrixXf& X,
                                       const Eigen::VectorXd& loco_pred,
                                       bool include_loco,
                                       const std::vector<int>& keep_indices,
                                       const std::string& chrom_filter,
                                       long long from_bp,
                                       long long to_bp,
                                       double min_maf,
                                       const std::unordered_set<std::string>& extract_ids,
                                       const std::unordered_set<std::string>& exclude_ids,
                                       int n_threads,
                                       std::vector<GwasResult>& results);
void scan_bgen_biobank_linear_parallel(const BgenReader& reader,
                                       uint32_t start_idx,
                                       uint32_t count,
                                       const Eigen::VectorXd& y,
                                       const Eigen::MatrixXf& X,
                                       const Eigen::VectorXd& loco_pred,
                                       bool include_loco,
                                       const std::vector<int>& keep_indices,
                                       const std::string& chrom_filter,
                                       long long from_bp,
                                       long long to_bp,
                                       double min_maf,
                                       const std::unordered_set<std::string>& extract_ids,
                                       const std::unordered_set<std::string>& exclude_ids,
                                       int n_threads,
                                       std::vector<GwasResult>& results);
void scan_plink_biobank_linear(const PlinkReader& reader,
                               const Eigen::VectorXd& y,
                               const Eigen::MatrixXf& X,
                               const Eigen::VectorXd& loco_pred,
                               bool include_loco,
                               const std::vector<int>& keep_indices,
                               const std::string& chrom_filter,
                               long long from_bp,
                               long long to_bp,
                               double min_maf,
                               const std::unordered_set<std::string>& extract_ids,
                               const std::unordered_set<std::string>& exclude_ids,
                               int n_threads,
                               std::vector<GwasResult>& results);

void GwasScanner::runBiobankScan(
    const Eigen::VectorXd& y,
    const Eigen::MatrixXf& X,
    PlinkReader& plink_reader,
    BgenReader& bgen_reader,
    PgenReader& pgen_reader,
    bool is_bgen,
    bool is_pgen,
    const DataFilter::KeepInfo& keep
) {
    std::cout << "\n=========================================\n";
    std::cout << "  Biobank Mode (REGENIE-style Ridge)     \n";
    std::cout << "=========================================\n";

    if (options.regenie_step < 0 || options.regenie_step > 2) {
        throw std::runtime_error("--regenie-step must be 0, 1, or 2.");
    }
    if (options.regenie_step != 2 && (is_bgen || is_pgen) && options.regenie_model_snps.empty()) {
        throw std::runtime_error("Biobank Step 1 with PGEN/BGEN requires --regenie-model-snps <file>.");
    }
    if (options.regenie_step == 2 && options.regenie_use_pred.empty()) {
        throw std::runtime_error("--regenie-step 2 requires --use-pred <file>.");
    }
    if (options.regenie_step == 1 && options.regenie_write_pred.empty()) {
        std::cout << "  Warning: --regenie-step 1 requested without --write-pred; predictions will not be saved.\n";
    }

    Eigen::VectorXd y_sub = DataFilter::subset(y, keep);
    Eigen::MatrixXf X_sub = DataFilter::subsetRows(X, keep);
    int n = keep.n_keep;
    int m = is_bgen ? bgen_reader.getNumVariants() : (is_pgen ? pgen_reader.getNumVariants() : plink_reader.getNumSnps());
    Eigen::VectorXd y_step1 = y_sub;
    Eigen::MatrixXd Q_cov = orthonormal_covariates(X_sub.cast<double>());
    if (Q_cov.cols() > 0) y_step1 -= Q_cov * (Q_cov.transpose() * y_step1);
    std::vector<Eigen::VectorXd> loco_predictions(23, Eigen::VectorXd::Zero(n));

    std::set<std::string> observed_chroms;
    if (is_bgen) {
        bgen_reader.loadVariantInfo();
        for (const auto& v : bgen_reader.getVariants()) observed_chroms.insert(v.chrom);
    } else if (is_pgen) {
        for (const auto& v : pgen_reader.getVariants()) observed_chroms.insert(v.chrom);
    } else {
        for (const auto& snp : plink_reader.getBimInfo()) observed_chroms.insert(snp.chrom);
    }

    if (!options.regenie_use_pred.empty()) {
        std::cout << "  Loading biobank prediction cache: " << options.regenie_use_pred << std::endl;
        if (is_bgen) {
            loco_predictions = read_biobank_predictions_by_sample_id(options.regenie_use_pred, bgen_reader.getSampleIds(), keep.keep_indices);
        } else if (is_pgen) {
            loco_predictions = read_biobank_predictions_by_sample_id(options.regenie_use_pred, pgen_reader.getSampleIds(), keep.keep_indices);
        } else {
            loco_predictions = read_biobank_predictions(options.regenie_use_pred, plink_reader, keep.keep_indices);
        }
    } else {
        RegenieStep1::Options step1_opts;
        step1_opts.block_size = options.regenie_block_size;
        step1_opts.use_loco = !options.no_loco;
        step1_opts.n_threads = options.n_threads;

        auto model_snps = load_snp_id_set(options.regenie_model_snps);
        if (!model_snps.empty()) {
            std::cout << "  Step 1 model SNP list: " << model_snps.size() << " SNP IDs\n";
        }

        int candidate_model_snps = 0;
        const auto& bim = plink_reader.getBimInfo();
        if (is_bgen) {
            for (const auto& v : bgen_reader.getVariants())
                if (model_snps.empty() || model_snps.find(v.id) != model_snps.end()) candidate_model_snps++;
        } else if (is_pgen) {
            for (const auto& v : pgen_reader.getVariants())
                if (model_snps.empty() || model_snps.find(v.id) != model_snps.end()) candidate_model_snps++;
        } else {
            for (const auto& snp : plink_reader.getBimInfo())
                if (model_snps.empty() || model_snps.find(snp.id) != model_snps.end()) candidate_model_snps++;
        }
        step1_opts.n_model_snps = candidate_model_snps;
        RegenieStep1 regenie_step1(step1_opts);

        std::cout << "  Loading Genotypes in blocks for Step 1..." << std::endl;

        int model_snps_used = 0;
        int current_block_snps = 0;
        Eigen::MatrixXd X_block(n, step1_opts.block_size);
        int current_block_chrom = 0;

        if (is_bgen) {
            std::string chrom, rsid, ref, alt;
            uint32_t pos;
            Eigen::VectorXf dosage_f;

            bgen_reader.reset();
            while (bgen_reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f)) {
                if (!model_snps.empty() && model_snps.find(rsid) == model_snps.end()) continue;
                if (!extract_ids_.empty() && extract_ids_.find(rsid) == extract_ids_.end()) continue;
                if (!exclude_ids_.empty() && exclude_ids_.find(rsid) != exclude_ids_.end()) continue;

                int chr_int = 0;
                try { chr_int = std::stoi(chrom); } catch (...) { chr_int = 0; }

                if (current_block_snps == 0) current_block_chrom = chr_int;
                if (current_block_snps > 0 && chr_int != current_block_chrom) {
                    regenie_step1.add_block(X_block.leftCols(current_block_snps), y_step1, current_block_chrom);
                    current_block_snps = 0;
                    current_block_chrom = chr_int;
                }

                float sum = 0.0f; int count = 0;
                for (int i = 0; i < n; ++i) {
                    float val = dosage_f[keep.keep_indices[i]];
                    if (!std::isnan(val)) { sum += val; count++; }
                }
                float mean = count > 0 ? sum / count : 0.0f;
                for (int i = 0; i < n; ++i) {
                    float val = dosage_f[keep.keep_indices[i]];
                    X_block(i, current_block_snps) = std::isnan(val) ? mean : val;
                }

                current_block_snps++;
                model_snps_used++;

                if (current_block_snps == step1_opts.block_size) {
                    regenie_step1.add_block(X_block.leftCols(current_block_snps), y_step1, current_block_chrom);
                    current_block_snps = 0;
                }
                if (!model_snps.empty() && model_snps_used >= candidate_model_snps) break;
            }
        } else if (is_pgen) {
            std::string chrom, rsid, ref, alt;
            uint32_t pos;
            Eigen::VectorXf dosage_f;

            pgen_reader.reset();
            while (pgen_reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f)) {
                if (!model_snps.empty() && model_snps.find(rsid) == model_snps.end()) continue;
                if (!extract_ids_.empty() && extract_ids_.find(rsid) == extract_ids_.end()) continue;
                if (!exclude_ids_.empty() && exclude_ids_.find(rsid) != exclude_ids_.end()) continue;

                int chr_int = 0;
                try { chr_int = std::stoi(chrom); } catch (...) { chr_int = 0; }

                if (current_block_snps == 0) current_block_chrom = chr_int;
                if (current_block_snps > 0 && chr_int != current_block_chrom) {
                    regenie_step1.add_block(X_block.leftCols(current_block_snps), y_step1, current_block_chrom);
                    current_block_snps = 0;
                    current_block_chrom = chr_int;
                }

                float sum = 0.0f; int count = 0;
                for (int i = 0; i < n; ++i) {
                    float val = dosage_f[keep.keep_indices[i]];
                    if (!std::isnan(val)) { sum += val; count++; }
                }
                float mean = count > 0 ? sum / count : 0.0f;
                for (int i = 0; i < n; ++i) {
                    float val = dosage_f[keep.keep_indices[i]];
                    X_block(i, current_block_snps) = std::isnan(val) ? mean : val;
                }

                current_block_snps++;
                model_snps_used++;

                if (current_block_snps == step1_opts.block_size) {
                    regenie_step1.add_block(X_block.leftCols(current_block_snps), y_step1, current_block_chrom);
                    current_block_snps = 0;
                }
                if (!model_snps.empty() && model_snps_used >= candidate_model_snps) break;
            }
        } else {
            int scanned_snps = 0;
            std::string chrom;
            Eigen::VectorXf dosage_f;

            while (scanned_snps < m) {
                chrom = bim[scanned_snps].chrom;
                if (!model_snps.empty() && model_snps.find(bim[scanned_snps].id) == model_snps.end()) {
                    scanned_snps++;
                    continue;
                }
                if (!extract_ids_.empty() && extract_ids_.find(bim[scanned_snps].id) == extract_ids_.end()) {
                    scanned_snps++;
                    continue;
                }
                if (!exclude_ids_.empty() && exclude_ids_.find(bim[scanned_snps].id) != exclude_ids_.end()) {
                    scanned_snps++;
                    continue;
                }

                Eigen::VectorXd d_d(plink_reader.getNumSamples());
                plink_reader.readSnp(scanned_snps, d_d);
                dosage_f = d_d.cast<float>();

                int chr_int = 0;
                try { chr_int = std::stoi(chrom); } catch (...) { chr_int = 0; }

                if (current_block_snps == 0) current_block_chrom = chr_int;
                if (current_block_snps > 0 && chr_int != current_block_chrom) {
                    regenie_step1.add_block(X_block.leftCols(current_block_snps), y_step1, current_block_chrom);
                    current_block_snps = 0;
                    current_block_chrom = chr_int;
                }

                float sum = 0.0f; int count = 0;
                for (int i = 0; i < n; ++i) {
                    float val = dosage_f[keep.keep_indices[i]];
                    if (!std::isnan(val)) { sum += val; count++; }
                }
                float mean = count > 0 ? sum / count : 0.0f;
                for (int i = 0; i < n; ++i) {
                    float val = dosage_f[keep.keep_indices[i]];
                    X_block(i, current_block_snps) = std::isnan(val) ? mean : val;
                }

                current_block_snps++;
                model_snps_used++;
                scanned_snps++;

                if (current_block_snps == step1_opts.block_size) {
                    regenie_step1.add_block(X_block.leftCols(current_block_snps), y_step1, current_block_chrom);
                    current_block_snps = 0;
                }
                if (!model_snps.empty() && model_snps_used >= candidate_model_snps) break;
            }
        }

        if (current_block_snps > 0) {
            regenie_step1.add_block(X_block.leftCols(current_block_snps), y_step1, current_block_chrom);
        }
        if (model_snps_used == 0) throw std::runtime_error("No SNPs available for biobank Step 1 model fitting.");
        std::cout << "  Step 1 model SNPs used: " << model_snps_used << std::endl;

        regenie_step1.fit_level1_stacking(y_step1);
        loco_predictions[0] = regenie_step1.get_loco_prediction(0);
        for (int chr = 1; chr <= 22; ++chr) {
            Eigen::VectorXd pred = regenie_step1.get_loco_prediction(chr);
            if (pred.size() == n) loco_predictions[chr] = pred;
        }
        if (is_bgen) {
            write_biobank_predictions(options.regenie_write_pred, bgen_reader.getSampleIds(), keep.keep_indices, loco_predictions);
        } else if (is_pgen) {
            write_biobank_predictions(options.regenie_write_pred, pgen_reader.getSampleIds(), keep.keep_indices, loco_predictions);
        } else {
            write_biobank_predictions(options.regenie_write_pred, plink_reader, keep.keep_indices, loco_predictions);
        }
        if (!options.regenie_write_pred.empty()) {
            std::cout << "  Wrote biobank prediction cache: " << options.regenie_write_pred << std::endl;
        }
    }

    if (options.regenie_step == 1) {
        std::cout << "  Biobank Step 1 complete. Skipping association scan due to --regenie-step 1.\n";
        return;
    }

    std::cout << "\n[Step 2] Running Single-SNP Association Tests...\n";
    if (!options.scan_chr.empty() || options.scan_from_bp >= 0 || options.scan_to_bp >= 0) {
        std::cout << "  Step 2 scan window: chr=" << (options.scan_chr.empty() ? "ALL" : options.scan_chr)
                  << " from=" << (options.scan_from_bp >= 0 ? std::to_string(options.scan_from_bp) : "START")
                  << " to=" << (options.scan_to_bp >= 0 ? std::to_string(options.scan_to_bp) : "END") << "\n";
    }
    std::string out_file = options.out_prefix + ".biobank.assoc.txt";
    // Header + data written in one pass by flushAssocResults

    if (options.no_loco) {
        Eigen::VectorXd loco_pred = loco_predictions[0];
        auto scan_start = std::chrono::steady_clock::now();
        std::vector<GwasResult> results;
        if (is_bgen) {
            bgen_reader.loadVariantInfo();
            scan_bgen_biobank_linear_parallel(bgen_reader, 0, bgen_reader.getNumVariants(),
                                             y_sub, X_sub, loco_pred, true,
                                             keep.keep_indices, options.scan_chr, options.scan_from_bp,
                                             options.scan_to_bp, options.min_maf,
                                             extract_ids_, exclude_ids_, options.n_threads, results);
        } else if (is_pgen) {
            scan_pgen_biobank_linear_parallel(pgen_reader, 0, pgen_reader.getNumVariants(),
                                              y_sub, X_sub, loco_pred, true,
                                              keep.keep_indices, options.scan_chr, options.scan_from_bp,
                                              options.scan_to_bp, options.min_maf,
                                              extract_ids_, exclude_ids_, options.n_threads, results);
        } else {
            scan_plink_biobank_linear(plink_reader, y_sub, X_sub, loco_pred, true,
                                      keep.keep_indices, options.scan_chr, options.scan_from_bp,
                                      options.scan_to_bp, options.min_maf,
                                      extract_ids_, exclude_ids_, options.n_threads, results);
        }
        appendResults(out_file, results);
        log_stage_time("biobank Step2 scan", scan_start);
    } else {
        std::vector<std::string> sorted_chroms = sort_chrom_labels(observed_chroms);
        for (const std::string& chrom_label : sorted_chroms) {
            if (!options.scan_chr.empty() && chrom_label != options.scan_chr) continue;
            int chr = 0;
            try { chr = std::stoi(chrom_label); } catch (...) { chr = 0; }
            if (chr <= 0 || chr > 22) continue;
            Eigen::VectorXd loco_pred = loco_predictions[chr];
            if (loco_pred.size() != n) loco_pred = Eigen::VectorXd::Zero(n);

            std::cout << "  Scanning Chromosome " << chrom_label << " using LOCO predictions...\n";
            auto scan_start = std::chrono::steady_clock::now();
            std::vector<GwasResult> results;
            if (is_bgen) {
                bgen_reader.loadVariantInfo();
                auto [start_idx, count] = bgen_reader.getChromRange(chrom_label);
                if (count > 0) {
                    scan_bgen_biobank_linear_parallel(bgen_reader, start_idx, count,
                                                     y_sub, X_sub, loco_pred, true,
                                                     keep.keep_indices, chrom_label, options.scan_from_bp,
                                                     options.scan_to_bp, options.min_maf,
                                                     extract_ids_, exclude_ids_, options.n_threads, results);
                }
            } else if (is_pgen) {
                auto [start_idx, count] = pgen_reader.getChromRange(chrom_label);
                if (count > 0) {
                    scan_pgen_biobank_linear_parallel(pgen_reader, start_idx, count,
                                                      y_sub, X_sub, loco_pred, true,
                                                      keep.keep_indices, chrom_label, options.scan_from_bp,
                                                      options.scan_to_bp, options.min_maf,
                                                      extract_ids_, exclude_ids_, options.n_threads, results);
                }
            } else {
                scan_plink_biobank_linear(plink_reader, y_sub, X_sub, loco_pred, true,
                                          keep.keep_indices, chrom_label, options.scan_from_bp,
                                          options.scan_to_bp, options.min_maf,
                                          extract_ids_, exclude_ids_, options.n_threads, results);
            }
            appendResults(out_file, results);
            log_stage_time("biobank Step2 scan chr " + chrom_label, scan_start);
        }
    }

    flushAssocResults(out_file);
    std::cout << "  Biobank GWAS Scan complete. Results saved to: " << out_file << std::endl;
}

} // namespace cosmic
