#pragma once
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unordered_set>
#include <omp.h>
#include <Eigen/Dense>
#include "covar_reader.h"
#include "pheno_reader.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"
#include "grm_builder.h"
#include "config.h"

namespace cosmic {
namespace DataLoader {

struct GenotypeMeta {
    bool is_bgen = false;
    bool is_pgen = false;
    bool is_vcf = false;
    int n_samples = 0;
    int n_snps = 0;
    std::vector<std::string> fids;
    std::vector<std::string> iids;
};

inline GenotypeMeta loadGenotype(
    const std::string& bfile_prefix,
    const std::string& bgen_file,
    const std::string& sample_file,
    const std::string& pfile_prefix,
    const std::string& pgen_file,
    const std::string& pvar_file,
    const std::string& psam_file,
    const std::string& vcf_file,
    const std::string& vcf_tbi_file,
    const std::string& vcf_sample_file,
    PlinkReader& plink_reader,
    BgenReader& bgen_reader,
    PgenReader& pgen_reader,
    VcfReader& vcf_reader)
{
    GenotypeMeta meta;
    meta.is_bgen = !bgen_file.empty();
    meta.is_pgen = (!pfile_prefix.empty() || !pgen_file.empty());
    meta.is_vcf = !vcf_file.empty();

    if (meta.is_vcf) {
        std::cout << "Loading VCF data from: " << vcf_file << "\n";
        vcf_reader.open(vcf_file, vcf_tbi_file);
        vcf_reader.reset();
        meta.n_samples = vcf_reader.getNumSamples();
        meta.n_snps = vcf_reader.getNumVariants();

        if (!vcf_sample_file.empty()) {
            std::cout << "Loading VCF sample file: " << vcf_sample_file << "\n";
            std::ifstream sfs(vcf_sample_file);
            std::string line;
            int cnt = 0;
            while (std::getline(sfs, line)) {
                if (line.empty() || line[0] == '#') continue;
                std::stringstream ss(line);
                std::string fid, iid;
                ss >> fid >> iid;
                meta.fids.push_back(fid);
                meta.iids.push_back(iid);
                cnt++;
            }
            if (cnt != meta.n_samples) {
                std::cerr << "Warning: VCF sample file count (" << cnt << ") != VCF sample count (" << meta.n_samples << ").\n";
            }
        } else {
            const auto& ids = vcf_reader.getSampleIds();
            for (const auto& id : ids) {
                // Without a sample mapping file, treat the full VCF sample ID as IID with FID="0"
                // Only split on '_' or ':' if user provides --vcf-sample mapping file
                meta.fids.push_back("0");
                meta.iids.push_back(id);
            }
        }
    } else if (meta.is_bgen) {
        std::cout << "Loading BGEN data from: " << bgen_file << "\n";
        bgen_reader.open(bgen_file);
        meta.n_samples = bgen_reader.initialize();
        meta.n_snps = bgen_reader.getNumVariants();

        if (!sample_file.empty()) {
            std::cout << "Loading sample file: " << sample_file << "\n";
            std::ifstream sfs(sample_file);
            std::string line;
            std::getline(sfs, line);
            std::getline(sfs, line);
            int cnt = 0;
            while (std::getline(sfs, line)) {
                if (line.empty()) continue;
                std::stringstream ss(line);
                std::string fid, iid;
                ss >> fid >> iid;
                meta.fids.push_back(fid);
                meta.iids.push_back(iid);
                cnt++;
            }
            if (cnt != meta.n_samples) {
                std::cerr << "Warning: Sample file count (" << cnt << ") != BGEN sample count (" << meta.n_samples << ").\n";
            }
        } else {
            const auto& ids = bgen_reader.getSampleIds();
            for (const auto& id : ids) {
                size_t sep = id.find('_');
                if (sep == std::string::npos) sep = id.find(':');
                if (sep != std::string::npos) {
                    meta.fids.push_back(id.substr(0, sep));
                    meta.iids.push_back(id.substr(sep + 1));
                } else {
                    meta.fids.push_back(id);
                    meta.iids.push_back(id);
                }
            }
        }
    } else if (meta.is_pgen) {
        std::cout << "Loading PGEN data...\n";
        std::string pgen = pgen_file;
        std::string pvar = pvar_file;
        std::string psam = psam_file;
        if (!pfile_prefix.empty()) {
            if (pgen.empty()) pgen = pfile_prefix + ".pgen";
            if (pvar.empty()) pvar = pfile_prefix + ".pvar";
            if (psam.empty()) psam = pfile_prefix + ".psam";
        }
        pgen_reader.open(pgen, pvar, psam);
        meta.n_samples = pgen_reader.getNumSamples();
        meta.n_snps = pgen_reader.getNumVariants();

        bool loaded_psam_ids = false;
        if (!psam.empty()) {
            std::ifstream psam_fs(psam);
            std::string header;
            if (std::getline(psam_fs, header)) {
                std::vector<std::string> columns;
                std::stringstream hs(header);
                std::string col;
                while (hs >> col) {
                    if (!col.empty() && col[0] == '#') col.erase(0, 1);
                    columns.push_back(col);
                }
                int fid_col = -1;
                int iid_col = -1;
                for (int i = 0; i < static_cast<int>(columns.size()); ++i) {
                    if (columns[i] == "FID") fid_col = i;
                    if (columns[i] == "IID") iid_col = i;
                }

                if (iid_col >= 0) {
                    std::vector<std::string> psam_fids;
                    std::vector<std::string> psam_iids;
                    std::string line;
                    while (std::getline(psam_fs, line)) {
                        if (line.empty()) continue;
                        std::stringstream ss(line);
                        std::vector<std::string> fields;
                        std::string field;
                        while (ss >> field) fields.push_back(field);
                        if (iid_col >= static_cast<int>(fields.size())) continue;
                        std::string iid = fields[iid_col];
                        std::string fid = (fid_col >= 0 && fid_col < static_cast<int>(fields.size())) ? fields[fid_col] : iid;
                        psam_fids.push_back(fid);
                        psam_iids.push_back(iid);
                    }
                    if (static_cast<int>(psam_iids.size()) == meta.n_samples) {
                        meta.fids = std::move(psam_fids);
                        meta.iids = std::move(psam_iids);
                        loaded_psam_ids = true;
                    } else if (!psam_iids.empty()) {
                        std::cerr << "Warning: PSAM sample count (" << psam_iids.size()
                                  << ") != PGEN sample count (" << meta.n_samples
                                  << "). Falling back to PGEN sample IDs.\n";
                    }
                }
            }
        }

        if (!loaded_psam_ids) {
            const auto& ids = pgen_reader.getSampleIds();
            for (const auto& id : ids) {
                size_t sep = id.find(':');
                if (sep != std::string::npos) {
                    meta.fids.push_back(id.substr(0, sep));
                    meta.iids.push_back(id.substr(sep + 1));
                } else {
                    meta.fids.push_back(id);
                    meta.iids.push_back(id);
                }
            }
        }
    } else {
        std::cout << "Loading PLINK data from: " << bfile_prefix << ".*\n";
        plink_reader.load(bfile_prefix);
        meta.n_samples = plink_reader.getNumSamples();
        meta.n_snps = plink_reader.getNumSnps();

        const auto& fam = plink_reader.getFamInfo();
        for (const auto& f : fam) {
            meta.fids.push_back(f.fid);
            meta.iids.push_back(f.iid);
        }
    }

    std::cout << "Loaded " << meta.n_samples << " individuals and " << meta.n_snps << " SNPs.\n";
    return meta;
}

inline Eigen::MatrixXf loadCovariates(
    const std::string& qcovar_file,
    const std::string& covar_file,
    const std::vector<std::string>& fids,
    const std::vector<std::string>& iids,
    int n_samples)
{
    Eigen::MatrixXf covs_q, covs_c;
    bool has_q = false, has_c = false;

    if (!qcovar_file.empty()) {
        CovarReader cov_reader;
        cov_reader.load(qcovar_file);
        cov_reader.align(fids, iids, covs_q);
        has_q = true;
    }
    if (!covar_file.empty()) {
        CovarReader cov_reader;
        cov_reader.loadCategorical(covar_file);
        cov_reader.align(fids, iids, covs_c);
        has_c = true;
    }

    int total_cols = 1;
    if (has_q) total_cols += covs_q.cols();
    if (has_c) total_cols += covs_c.cols();

    Eigen::MatrixXf X(n_samples, total_cols);
    X.col(0).setOnes();
    int current_col = 1;
    if (has_q) {
        X.middleCols(current_col, covs_q.cols()) = covs_q;
        current_col += covs_q.cols();
    }
    if (has_c) {
        X.middleCols(current_col, covs_c.cols()) = covs_c;
    }

    if (has_q || has_c) {
        std::cout << "Using " << total_cols << " fixed-effect columns (including intercept).\n";
    }
    return X;
}

inline Eigen::VectorXd loadPhenotype(
    const std::string& pheno_file,
    int mpheno,
    const std::vector<std::string>& fids,
    const std::vector<std::string>& iids,
    int n_samples,
    PlinkReader& plink_reader,
    bool is_bgen,
    bool is_pgen)
{
    Eigen::VectorXd y(n_samples);
    if (pheno_file.empty()) {
        if (is_bgen || is_pgen) {
            throw std::runtime_error("BGEN/PGEN input requires --pheno file.");
        }
        const auto& fam = plink_reader.getFamInfo();
        for (int i = 0; i < n_samples; ++i) y(i) = fam[i].phenotype;
    } else {
        PhenoReader pheno_reader;
        pheno_reader.load(pheno_file, mpheno);
        pheno_reader.align(fids, iids, y);
    }
    return y;
}

inline void runMakeGrm(const GwasConfig& cfg, PlinkReader& plink_reader, BgenReader& bgen_reader, PgenReader& pgen_reader, VcfReader& vcf_reader, bool is_bgen, bool is_pgen, bool is_vcf, SequentialGenotypeMatrix* generic_reader = nullptr) {
    std::cout << (cfg.make_grm_sparse ? "Calculating Sparse GRM..." : "Calculating GRM...") << "\n";
    GrmBuilder::Options opts;
    opts.min_maf = (float)cfg.min_maf;
    opts.min_hwe = (float)cfg.min_hwe;
    opts.bgen_dom_max_missing = cfg.bgen_dom_max_missing;
    opts.bgen_dom_min_info = cfg.bgen_dom_min_info;
    opts.bgen_dom_entropy_warning = cfg.bgen_dom_entropy_warning;
    opts.algorithm = cfg.grm_alg;
    opts.force_slow = cfg.force_slow;
    opts.snp_weight_file = cfg.snp_weight_file;
    opts.pop_class_file = cfg.pop_class_file;
    auto load_ids = [](const std::string& path) {
        std::unordered_set<std::string> ids;
        if (path.empty()) return ids;
        std::ifstream input(path);
        if (!input) throw std::runtime_error("Cannot open SNP list: " + path);
        std::string id;
        while (input >> id) ids.insert(id);
        if (ids.empty()) throw std::runtime_error("SNP list is empty: " + path);
        return ids;
    };
    opts.extract_set = load_ids(cfg.extract_file);
    opts.exclude_set = load_ids(cfg.exclude_file);
    if (cfg.make_grm_sparse) {
        opts.sparse = true;
        opts.sparse_threshold = (float)cfg.sparse_cutoff;
        opts.sparse_out_prefix = cfg.out_prefix;
        std::cout << "  Sparse Threshold: " << cfg.sparse_cutoff << "\n";
    }
    #ifdef _OPENMP
    opts.thread_num = omp_get_max_threads();
    #endif

    GrmBuilder builder(opts);
    if (generic_reader) {
        if (cfg.make_grm_sparse) throw std::runtime_error("Sparse GRM is not supported for HapMap/text input");
        builder.compute(*generic_reader);
    } else if (is_vcf) {
        if (cfg.make_grm_sparse) builder.computeSparse(vcf_reader);
        else builder.compute(vcf_reader);
    } else if (is_bgen) {
        if (cfg.make_grm_sparse) builder.computeSparse(bgen_reader);
        else builder.compute(bgen_reader);
    } else if (is_pgen) {
        if (cfg.make_grm_sparse) builder.computeSparse(pgen_reader);
        else builder.compute(pgen_reader);
    } else {
        builder.compute(plink_reader);
    }

    if (!cfg.make_grm_sparse) {
        builder.save(cfg.out_prefix);
    }
}

inline void runGemmaEigen(const GwasConfig& cfg, PlinkReader& plink_reader, BgenReader& bgen_reader, PgenReader& pgen_reader, bool is_bgen, bool is_pgen) {
    std::cout << "Calculating GEMMA-style eigen decomposition...\n";
    GrmBuilder::Options opts;
    opts.min_maf = (float)cfg.min_maf;
    opts.min_hwe = (float)cfg.min_hwe;
    opts.algorithm = cfg.grm_alg;
    opts.force_slow = cfg.force_slow;
    opts.snp_weight_file = cfg.snp_weight_file;
    opts.pop_class_file = cfg.pop_class_file;
    #ifdef _OPENMP
    opts.thread_num = omp_get_max_threads();
    #endif

    GrmBuilder builder(opts);
    if (is_bgen) builder.compute(bgen_reader);
    else if (is_pgen) builder.compute(pgen_reader);
    else builder.compute(plink_reader);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(builder.getGrm());
    if (es.info() != Eigen::Success) throw std::runtime_error("GRM eigen decomposition failed");

    std::ofstream fd(cfg.out_prefix + ".eigenD.txt");
    std::ofstream fu(cfg.out_prefix + ".eigenU.txt");
    if (!fd || !fu) throw std::runtime_error("Cannot write GEMMA-style eigen output");
    fd << std::scientific << std::setprecision(12);
    fu << std::scientific << std::setprecision(12);
    const auto& eval = es.eigenvalues();
    const auto& evec = es.eigenvectors();
    for (int i = 0; i < eval.size(); ++i) fd << eval[i] << "\n";
    for (int i = 0; i < evec.rows(); ++i) {
        for (int j = 0; j < evec.cols(); ++j) {
            if (j) fu << "\t";
            fu << evec(i, j);
        }
        fu << "\n";
    }
    std::cout << "GEMMA-style eigen outputs written to " << cfg.out_prefix
              << ".eigenD.txt and " << cfg.out_prefix << ".eigenU.txt\n";
}

} // namespace DataLoader

using DataLoader::loadGenotype;
using DataLoader::loadCovariates;
using DataLoader::loadPhenotype;
using DataLoader::runMakeGrm;
using DataLoader::runGemmaEigen;

} // namespace cosmic
