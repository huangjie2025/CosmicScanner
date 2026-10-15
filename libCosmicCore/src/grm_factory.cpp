#include "grm_factory.h"
#include "logger.h"
#include <omp.h>

namespace cosmic {

std::vector<std::string> GrmFactory::extractIds(PlinkReader& reader) {
    std::vector<std::string> ids;
    const auto& fam = reader.getFamInfo();
    ids.reserve(fam.size());
    for (const auto& f : fam) ids.push_back(f.iid.empty() ? f.fid : f.iid);
    return ids;
}

std::vector<std::string> GrmFactory::extractIds(PgenReader& reader) {
    return reader.getSampleIds();
}

std::vector<std::string> GrmFactory::extractIds(BgenReader& reader) {
    const auto& sids = reader.getSampleIds();
    std::vector<std::string> ids;
    ids.reserve(sids.size());
    for (const auto& id : sids) ids.push_back(id);
    return ids;
}

GrmBuildResult GrmFactory::buildGrm(const GrmBuildConfig& cfg,
                                    const std::string& bfile,
                                    const std::string& pfile,
                                    PlinkReader* plink,
                                    PgenReader* pgen,
                                    BgenReader* bgen) {
    GrmBuilder::Options opts;
    opts.min_maf = cfg.min_maf;
    opts.min_hwe = cfg.min_hwe;
    opts.algorithm = cfg.algorithm;
    opts.compute_homo_hete = cfg.compute_homo_hete;
    opts.snp_weight_file = cfg.snp_weight_file;
    opts.pop_class_file = cfg.pop_class_file;
    opts.thread_num = cfg.threads;
    if (cfg.block_size > 0) opts.block_size = cfg.block_size;

    // Use pre-loaded reader if available, otherwise load from file
    if (bgen) {
        bgen->reset();
        return doCompute(opts, *bgen);
    } else if (pgen) {
        pgen->reset();
        return doCompute(opts, *pgen);
    } else if (plink) {
        return doCompute(opts, *plink);
    } else if (!bfile.empty()) {
        PlinkReader reader;
        reader.load(bfile);
        return doCompute(opts, reader);
    } else if (!pfile.empty()) {
        PgenReader reader;
        reader.open(pfile);
        return doCompute(opts, reader);
    } else {
        throw std::runtime_error("GrmFactory: no genotype source provided.");
    }
}

GrmBuildResult GrmFactory::buildChrGrm(const GrmBuildConfig& base_cfg,
                                       const std::string& chrom,
                                       PlinkReader* plink,
                                       PgenReader* pgen,
                                       BgenReader* bgen) {
    GrmBuildConfig chr_cfg = base_cfg;
    chr_cfg.chrom_filter = chrom;
    return buildGrm(chr_cfg, "", "", plink, pgen, bgen);
}

} // namespace cosmic
