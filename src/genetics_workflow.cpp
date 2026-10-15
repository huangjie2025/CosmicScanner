#include "genetics_workflow.h"
#include "config.h"
#include "genetics/ld/ld_calculator.h"
#include "genetics/population/pca.h"
#include "genetics/qc/qc_filters.h"
#include "grm_builder.h"
#include "plink_reader.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cosmic {
namespace {
using Output = std::unique_ptr<FILE, decltype(&std::fclose)>;

Output openOutput(const std::string& path) {
    Output out(std::fopen(path.c_str(), "wx"), &std::fclose);
    if (!out) throw std::runtime_error("Cannot exclusively create genetics output: " + path);
    return out;
}

void finishOutput(Output& out) {
    const bool failed = std::ferror(out.get()) != 0;
    const int status = std::fclose(out.release());
    if (failed || status != 0) throw std::runtime_error("Failed writing genetics output");
}
} // namespace

bool prepareGeneticsWorkflow(GwasConfig& cfg) {
    using namespace genetics;
    const bool clump = !cfg.ld_clump_file.empty();
    const bool ld = cfg.ld_prune || clump;
    const bool standalone = cfg.pca || ld || cfg.cojo || cfg.cojo_slct || cfg.cojo_joint || cfg.meta_analysis || cfg.fine_map;
    const bool qc = cfg.qc_requested || cfg.geno != 1.0 || cfg.mind != 1.0 ||
                    !cfg.keep_file.empty() || !cfg.remove_file.empty();
    // Legacy maf/hwe/extract/exclude alone must not opt into physical subsetting.
    if (!standalone && !qc) return false;

    // COJO and meta-analysis are standalone summary-stats operations that don't require genotype data
    if (cfg.cojo || cfg.cojo_slct || cfg.cojo_joint || cfg.meta_analysis || cfg.fine_map) {
        return false; // Handled separately in main.cpp
    }

    if (cfg.bfile_prefix.empty() || !cfg.bgen_file.empty() || !cfg.sample_file.empty() ||
        !cfg.pfile_prefix.empty() || !cfg.pgen_file.empty() || !cfg.pvar_file.empty() ||
        !cfg.psam_file.empty())
        throw std::invalid_argument("Genetics QC/PCA/LD requires BED/BIM/FAM via --bfile only; BGEN/PGEN are unsupported");
    const bool old_standalone = cfg.make_grm || cfg.make_grm_sparse || cfg.gemma_s_mode != 0 ||
        cfg.gemma_q_mode != 0 || cfg.gemma_eigen || cfg.gemma_calccor || cfg.vc_mode != 0 ||
        cfg.ci_mode != 0 || !cfg.predict_effect_file.empty() || cfg.bslmm_mode != 0 ||
        cfg.make_agh || cfg.write_agh || cfg.write_ped_kernel || cfg.prepare_input_only ||
        cfg.preview_schema || !cfg.write_template_config.empty() ||
        cfg.cojo || cfg.cojo_slct || cfg.cojo_joint || cfg.meta_analysis;
    const int models = int(cfg.rrm_model) + int(cfg.repeatability_model) + int(cfg.multivariate) +
        int(cfg.maternal_model) + int(cfg.single_step_model) + int(cfg.animal_model) +
        int(cfg.common_env_model) + int(cfg.threshold_model) + int(cfg.count_model) +
        int(cfg.ordinal_model) + int(cfg.survival_model) + int(cfg.gxe_model) + int(cfg.linear_model) +
        int(cfg.farmcpu) + int(cfg.blink) + int(cfg.firth_logistic) + int(cfg.saige_style) +
        int(!cfg.logistic_test.empty()) + int(cfg.logistic_glmm) +
        int(!cfg.condition_snp_file.empty());
    if (int(cfg.pca) + int(cfg.ld_prune) + int(clump) > 1 || old_standalone || models > 1 ||
        (standalone && (cfg.run_gwas || models != 0 || cfg.fast_lmm || cfg.lmm_gemma_mode != 0 ||
                        cfg.eb_lmm || cfg.biobank_mode || cfg.run_skat || cfg.run_burden ||
                        cfg.use_spa || !cfg.scan_tier.empty())))
        throw std::invalid_argument("Conflicting genetics tasks: PCA, LD prune, LD clump, and existing standalone/GWAS routes cannot be combined; QC may precede a Gaussian GWAS");
    if (cfg.single_step_model || !cfg.genotyped_list.empty() || !cfg.ssgwas_apy_core_file.empty())
        throw std::invalid_argument("Genotype QC is unsupported with single-step/genotyped-list sample semantics");
    if (cfg.threshold_model || cfg.count_model || cfg.ordinal_model || cfg.survival_model ||
        cfg.use_spa || cfg.run_skat || cfg.run_burden || cfg.firth_logistic || cfg.saige_style ||
        !cfg.logistic_test.empty() || cfg.logistic_glmm)
        throw std::invalid_argument("Genotype QC downstream support is limited to Gaussian models");
    if (!cfg.snp_weight_file.empty() || !cfg.pop_class_file.empty() ||
        !cfg.repeatability_grm_cache_prefix.empty() || !cfg.repeatability_weight_file.empty() ||
        !cfg.regenie_model_snps.empty() || !cfg.regenie_use_pred.empty() ||
        !cfg.regenie_write_pred.empty() || cfg.regenie_step != 0 ||
        !cfg.multivariate_emmvgwas_G_file.empty() || !cfg.multivariate_emmvgwas_R_file.empty() ||
        !cfg.rrm_compatible_varcom.empty() || !cfg.category_file.empty() ||
        !cfg.multi_category_file.empty() || !cfg.category_weight_file.empty() ||
        !cfg.snp_weight_summary_file.empty() || !cfg.beta_file.empty() || !cfg.cor_file.empty() ||
        !cfg.study_prefix.empty() || !cfg.ref_prefix.empty() || !cfg.multi_study_file.empty() ||
        !cfg.multi_ref_file.empty())
        throw std::invalid_argument("Genotype QC cannot safely reuse caches, external kernels/predictions, or sample/variant-indexed auxiliary files; regenerate them on the QC subset separately");
    if (standalone && (cfg.snp_model != "additive" || !cfg.scan_chr.empty() ||
                       cfg.scan_from_bp >= 0 || cfg.scan_to_bp >= 0))
        throw std::invalid_argument("PCA/LD uses additive BED dosages and QC extract/exclude lists, not SNP model or GWAS region options");
    const auto probability = [](double x, double maximum, const char* name) {
        if (!std::isfinite(x) || x < 0.0 || x > maximum)
            throw std::invalid_argument(std::string("Invalid genetics option: ") + name);
    };
    probability(cfg.geno, 1.0, "geno");
    probability(cfg.mind, 1.0, "mind");
    probability(cfg.min_maf, 0.5, "maf");
    probability(cfg.min_hwe, 1.0, "hwe");
    if (ld) {
        probability(cfg.ld_r2, 1.0, "ld-r2");
        if (cfg.ld_window < 0 || cfg.ld_step < 1)
            throw std::invalid_argument("LD requires window >= 0 bp and step >= 1");
        if (clump) probability(cfg.clump_p, 1.0, "clump-p");
    }
    if (cfg.pca && (cfg.pca_k < 1 || cfg.grm_alg < 0 || cfg.grm_alg > 5))
        throw std::invalid_argument("PCA requires k >= 1 and a supported GRM algorithm (0..5)");
    if (cfg.out_prefix.empty() || std::filesystem::path(cfg.out_prefix).filename().empty())
        throw std::invalid_argument("Genetics output prefix needs a filename");
    std::vector<std::string> output_paths{cfg.out_prefix + ".genetics.meta.tsv"};
    if (cfg.pca) {
        output_paths.push_back(cfg.out_prefix + ".eigenvec");
        output_paths.push_back(cfg.out_prefix + ".eigenval");
    } else if (cfg.ld_prune) {
        output_paths.push_back(cfg.out_prefix + ".prune.in");
        output_paths.push_back(cfg.out_prefix + ".prune.out");
    } else if (clump) {
        output_paths.push_back(cfg.out_prefix + ".clump.tsv");
        output_paths.push_back(cfg.out_prefix + ".lead.tsv");
    }
    for (const auto& path : output_paths) {
        if (std::filesystem::exists(std::filesystem::symlink_status(path)))
            throw std::runtime_error("Genetics output already exists: " + path);
    }

    struct SummaryRow { double p; std::string chrom; long long pos; };
    std::unordered_map<std::string, SummaryRow> summary;
    if (clump) {
        std::ifstream in(cfg.ld_clump_file);
        if (!in) throw std::runtime_error("Cannot open LD clump summary: " + cfg.ld_clump_file);
        std::string line, token;
        if (!std::getline(in, line)) throw std::invalid_argument("Empty LD clump summary");
        std::istringstream header(line);
        std::vector<std::string> columns;
        int snp_col = -1, p_col = -1, chr_col = -1, bp_col = -1;
        std::unordered_set<std::string> names;
        while (header >> token) {
            const int index = static_cast<int>(columns.size());
            columns.push_back(token);
            // Accepted aliases: SNP/P/CHR/BP or snp/pvalue/chr/pos. Extra columns
            // are allowed, but every row must match the complete header width.
            std::string canonical = token;
            if (token == "snp") canonical = "SNP";
            if (token == "pvalue") canonical = "P";
            if (token == "chr") canonical = "CHR";
            if (token == "pos") canonical = "BP";
            if (!names.insert(canonical).second)
                throw std::invalid_argument("Duplicate LD summary column: " + token);
            if (canonical == "SNP") snp_col = index;
            if (canonical == "P") p_col = index;
            if (canonical == "CHR") chr_col = index;
            if (canonical == "BP") bp_col = index;
        }
        if (snp_col < 0 || p_col < 0 || chr_col < 0 || bp_col < 0)
            throw std::invalid_argument("LD summary requires SNP P CHR BP (aliases: snp pvalue chr pos)");
        size_t line_number = 1;
        while (std::getline(in, line)) {
            ++line_number;
            std::istringstream row(line);
            std::vector<std::string> fields;
            while (row >> token) fields.push_back(token);
            if (fields.empty()) continue;
            if (fields.size() != columns.size())
                throw std::invalid_argument("LD summary column count mismatch at line " + std::to_string(line_number));
            double p;
            long long pos;
            try {
                size_t used_p = 0, used_bp = 0;
                p = std::stod(fields[p_col], &used_p);
                pos = std::stoll(fields[bp_col], &used_bp);
                if (used_p != fields[p_col].size() || used_bp != fields[bp_col].size() ||
                    !std::isfinite(p) || p < 0.0 || p > 1.0 || pos < 0)
                    throw std::invalid_argument("range");
            } catch (const std::exception&) {
                throw std::invalid_argument("Invalid LD summary P/BP at line " + std::to_string(line_number));
            }
            if (!summary.emplace(fields[snp_col], SummaryRow{p, fields[chr_col], pos}).second)
                throw std::invalid_argument("Duplicate LD summary SNP ID: " + fields[snp_col]);
        }
        if (in.bad() || summary.empty()) throw std::invalid_argument("Empty or unreadable LD clump summary");
    }
    if (ld) {
        // Validate even QC-excluded matching IDs, before materializing any output.
        std::ifstream in(cfg.bfile_prefix + ".bim");
        if (!in) throw std::runtime_error("Cannot open LD BIM metadata");
        std::unordered_set<std::string> ids;
        std::string line;
        size_t overlap = 0;
        while (std::getline(in, line)) {
            std::istringstream row(line);
            std::string chr, id, a1, a2, extra;
            double cm;
            long long bp;
            if (!(row >> chr)) continue;
            if (!(row >> id >> cm >> bp >> a1 >> a2) || (row >> extra) || !std::isfinite(cm) || bp < 0)
                throw std::invalid_argument("Invalid six-column LD BIM metadata");
            if (!ids.insert(id).second) throw std::invalid_argument("Duplicate LD BIM SNP ID: " + id);
            const auto it = summary.find(id);
            if (it != summary.end()) {
                if (it->second.chrom != chr || it->second.pos != bp)
                    throw std::invalid_argument("LD summary CHR/BP mismatch for SNP: " + id);
                ++overlap;
            }
        }
        if (in.bad() || ids.empty()) throw std::invalid_argument("Empty or unreadable LD BIM metadata");
        if (clump && overlap == 0) throw std::invalid_argument("LD summary has no BIM overlap");
    }

    const std::string original_prefix = cfg.bfile_prefix;
    const auto filtered = QCFilter::preprocess(original_prefix, cfg.out_prefix + ".qc",
        QCOptions{cfg.geno, cfg.mind, cfg.min_maf, cfg.min_hwe, cfg.keep_file,
                  cfg.remove_file, cfg.extract_file, cfg.exclude_file});
    const size_t n = filtered.samples.size(), m = filtered.variants.size();
    size_t analyzed = m, missing_summary = 0, leads_count = 0, members_count = 0;
    long double ld_workspace_bytes = 0;
    long long grm_snps = 0;
    if (cfg.pca) {
        if (static_cast<size_t>(cfg.pca_k) > n)
            throw std::invalid_argument("PCA k exceeds the number of QC-surviving samples");
        std::cout << "PCA: dense GRM/eigendecomposition uses O(n^2) memory; the LD 512 MiB guard does not apply.\n";
        PlinkReader reader;
        reader.load(filtered.prefix);
        GrmBuilder::Options options;
        // QC already applied the thresholds on this physical sample/variant set.
        options.min_maf = 0.0f;
        options.min_hwe = 0.0f;
        options.algorithm = cfg.grm_alg;
        options.force_slow = cfg.force_slow;
        options.thread_num = cfg.num_threads > 0 ? cfg.num_threads : 1;
        GrmBuilder builder(options);
        builder.compute(reader);
        grm_snps = builder.getTotalSnps();
        if (grm_snps <= 0) throw std::runtime_error("PCA has no usable GRM variants after QC");
        const auto result = PCA::computeFromGrm(builder.getGrm(), cfg.pca_k);
        auto vectors = openOutput(cfg.out_prefix + ".eigenvec");
        auto values = openOutput(cfg.out_prefix + ".eigenval");
        std::fprintf(vectors.get(), "FID\tIID");
        for (int j = 0; j < cfg.pca_k; ++j) std::fprintf(vectors.get(), "\tPC%d", j + 1);
        std::fprintf(vectors.get(), "\n");
        for (size_t i = 0; i < n; ++i) {
            const auto& sample = reader.getFamInfo()[i];
            std::fprintf(vectors.get(), "%s\t%s", sample.fid.c_str(), sample.iid.c_str());
            for (int j = 0; j < cfg.pca_k; ++j)
                std::fprintf(vectors.get(), "\t%.17g", result.scores(static_cast<Eigen::Index>(i), j));
            std::fprintf(vectors.get(), "\n");
        }
        for (int j = 0; j < cfg.pca_k; ++j) std::fprintf(values.get(), "%.17g\n", result.eigenvalues[j]);
        finishOutput(vectors);
        finishOutput(values);
    } else if (ld) {
        if (n < 2 || n > static_cast<size_t>(std::numeric_limits<int>::max()) ||
            m > static_cast<size_t>(std::numeric_limits<int>::max()))
            throw std::invalid_argument("LD requires >= 2 samples and int-sized dimensions");
        PlinkReader reader;
        reader.load(filtered.prefix);
        const auto& bim = reader.getBimInfo();
        std::vector<int> order;
        for (size_t j = 0; j < m; ++j) {
            if (!clump || summary.count(bim[j].id)) order.push_back(static_cast<int>(j));
        }
        analyzed = order.size();
        missing_summary = m - analyzed;
        if (order.empty()) throw std::invalid_argument("LD summary has no usable overlap after QC");
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            return bim[a].chrom < bim[b].chrom ||
                   (bim[a].chrom == bim[b].chrom && bim[a].bp_pos < bim[b].bp_pos);
        });
        // Two dense buffers plus vector/index slack. The moved API argument
        // normally uses one buffer; this is deliberately conservative, not RSS.
        constexpr size_t limit = 512ULL * 1024 * 1024;
        ld_workspace_bytes = 2.0L * sizeof(double) * n * analyzed +
                             16.0L * n + 256.0L * m;
        if (ld_workspace_bytes > limit)
            throw std::runtime_error("LD dense workspace exceeds the conservative 512 MiB guard; reduce samples/variants with QC keep/extract. This guard excludes BED mmap, metadata strings, summary tables, and process overhead");
        std::cout << "LD workspace estimate: " << ld_workspace_bytes
                  << " bytes (512 MiB limit; excludes BED mmap, metadata/summary strings and process overhead).\n";
        Eigen::MatrixXd dosages(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(analyzed));
        Eigen::VectorXd p(static_cast<Eigen::Index>(analyzed));
        Eigen::VectorXf column;
        std::vector<std::string> chrom;
        std::vector<long long> pos;
        for (size_t j = 0; j < analyzed; ++j) {
            const int original = order[j];
            reader.getSnpGenotypes(original, column);
            float minimum = 2.0f, maximum = 0.0f;
            size_t called = 0;
            for (Eigen::Index i = 0; i < column.size(); ++i) {
                if (std::isnan(column[i])) continue;
                if (!std::isfinite(column[i]) || column[i] < 0.0f || column[i] > 2.0f)
                    throw std::runtime_error("Invalid LD dosage for SNP: " + bim[original].id);
                minimum = std::min(minimum, column[i]);
                maximum = std::max(maximum, column[i]);
                ++called;
            }
            if (called == 0 || minimum == maximum)
                throw std::invalid_argument("LD rejects all-missing/constant dosage SNP: " + bim[original].id + "; exclude it and rerun with a fresh output prefix");
            dosages.col(static_cast<Eigen::Index>(j)) = column.cast<double>();
            chrom.push_back(bim[original].chrom);
            pos.push_back(bim[original].bp_pos);
            if (clump) p[static_cast<Eigen::Index>(j)] = summary.at(bim[original].id).p;
        }
        if (cfg.ld_prune) {
            const auto kept = LDPruner::prune(std::move(dosages), chrom, pos,
                                              cfg.ld_window, cfg.ld_step, cfg.ld_r2);
            std::vector<bool> retained(m, false);
            for (int index : kept) retained[order[index]] = true;
            leads_count = kept.size();
            auto in = openOutput(cfg.out_prefix + ".prune.in");
            auto out = openOutput(cfg.out_prefix + ".prune.out");
            for (size_t j = 0; j < m; ++j)
                std::fprintf(retained[j] ? in.get() : out.get(), "%s\n", bim[j].id.c_str());
            finishOutput(in);
            finishOutput(out);
        } else {
            const auto lead = LDClumper::clump(std::move(dosages), chrom, pos,
                cfg.ld_window, cfg.ld_step, cfg.ld_r2, p, cfg.clump_p);
            std::vector<int> original_to_sorted(m, -1);
            for (size_t j = 0; j < analyzed; ++j) original_to_sorted[order[j]] = static_cast<int>(j);
            auto members = openOutput(cfg.out_prefix + ".clump.tsv");
            auto leads = openOutput(cfg.out_prefix + ".lead.tsv");
            std::fprintf(members.get(), "SNP\tP\tCHR\tBP\tLEAD_ID\n");
            std::fprintf(leads.get(), "SNP\tP\tCHR\tBP\tLEAD_ID\n");
            for (size_t original = 0; original < m; ++original) {
                const int j = original_to_sorted[original];
                if (j < 0 || lead[j] < 0) continue;
                const auto& variant = bim[original];
                const auto& lead_id = bim[order[lead[j]]].id;
                std::fprintf(members.get(), "%s\t%.17g\t%s\t%lld\t%s\n", variant.id.c_str(),
                    p[j], variant.chrom.c_str(), variant.bp_pos, lead_id.c_str());
                ++members_count;
                if (lead[j] == j) {
                    std::fprintf(leads.get(), "%s\t%.17g\t%s\t%lld\t%s\n", variant.id.c_str(),
                        p[j], variant.chrom.c_str(), variant.bp_pos, lead_id.c_str());
                    ++leads_count;
                }
            }
            finishOutput(members);
            finishOutput(leads);
        }
    }
    auto metadata = openOutput(output_paths.front());
    std::fprintf(metadata.get(), "key\tvalue\nworkflow\t%s\ninput_prefix\t%s\nqc_prefix\t%s\n"
        "samples_qc\t%zu\nvariants_qc\t%zu\nvariants_analyzed\t%zu\n"
        "geno\t%.17g\nmind\t%.17g\nmaf\t%.17g\nhwe\t%.17g\n",
        cfg.pca ? "pca" : cfg.ld_prune ? "ld_prune" : clump ? "ld_clump" : cfg.run_gwas ? "qc_then_gwas" : "qc_only",
        original_prefix.c_str(), filtered.prefix.c_str(), n, m, analyzed,
        cfg.geno, cfg.mind, cfg.min_maf, cfg.min_hwe);
    if (cfg.pca)
        std::fprintf(metadata.get(), "pca_k\t%d\ngrm_algorithm\t%d\ngrm_variants_used\t%lld\n"
            "pca_scores\tunit_eigenvectors_descending_eigenvalues\npca_memory\tdense_O(n^2);no_LD_memory_guard\n",
            cfg.pca_k, cfg.grm_alg, grm_snps);
    if (ld)
        std::fprintf(metadata.get(), "ld_window_bp\t%d\nld_step\t%d\nld_r2\t%.17g\n"
            "ld_workspace_estimated_bytes\t%.0Lf\nld_workspace_limit_bytes\t536870912\n"
            "ld_memory_scope\ttwo_dense_buffers_plus_vector_index_slack;excludes_BED_mmap_metadata_strings_summary_tables_process_overhead\n"
            "ld_sort\tlexicographic_CHR_then_BP_stable_ties\noutput_order\toriginal_BIM\n"
            "missing_dosages\tcolumn_mean_imputed\nconstant_dosages\trejected\nretained_or_leads\t%zu\n",
            cfg.ld_window, cfg.ld_step, cfg.ld_r2, ld_workspace_bytes, leads_count);
    if (clump)
        std::fprintf(metadata.get(), "clump_summary\t%s\nclump_p\t%.17g\nvariants_missing_summary\t%zu\n"
            "clump_members\t%zu\nunassigned_overlap\t%zu\n"
            "clump_omissions\tmissing_summary_and_unassigned_variants;leads_included_as_members\n",
            cfg.ld_clump_file.c_str(), cfg.clump_p, missing_summary, members_count, analyzed - members_count);
    finishOutput(metadata);
    std::cout << "Genetics workflow complete; QC BED prefix: " << filtered.prefix
              << "; metadata: " << output_paths.front() << '\n';
    if (!standalone && cfg.run_gwas) {
        cfg.bfile_prefix = filtered.prefix;
        return false;
    }
    return true;
}
} // namespace cosmic
