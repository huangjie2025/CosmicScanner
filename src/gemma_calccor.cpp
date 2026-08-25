#include "gemma_calccor.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_set>

#include "qc_utils.h"

namespace cosmic {
namespace {

std::unordered_set<std::string> load_id_set(const std::string& path) {
    std::unordered_set<std::string> ids;
    if (path.empty()) return ids;
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open SNP list: " + path);
    std::string id;
    while (in >> id) ids.insert(id);
    return ids;
}

bool is_missing(double x) {
    return !std::isfinite(x);
}

double dot_product(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    return a.dot(b);
}

} // namespace

GemmaCalccor::GemmaCalccor(const GwasConfig& cfg) : cfg_(cfg) {}

bool GemmaCalccor::passesIdFilters(const std::string& rs) const {
    static thread_local std::string cached_extract_path;
    static thread_local std::string cached_exclude_path;
    static thread_local std::unordered_set<std::string> extract;
    static thread_local std::unordered_set<std::string> exclude;

    if (cached_extract_path != cfg_.extract_file) {
        extract = load_id_set(cfg_.extract_file);
        cached_extract_path = cfg_.extract_file;
    }
    if (cached_exclude_path != cfg_.exclude_file) {
        exclude = load_id_set(cfg_.exclude_file);
        cached_exclude_path = cfg_.exclude_file;
    }

    if (!extract.empty() && !extract.count(rs)) return false;
    if (!exclude.empty() && exclude.count(rs)) return false;
    return true;
}

bool GemmaCalccor::passesRegionFilters(const std::string& chr, long long bp) const {
    if (!cfg_.scan_chr.empty() && chr != cfg_.scan_chr) return false;
    if (cfg_.scan_from_bp >= 0 && bp >= 0 && bp < cfg_.scan_from_bp) return false;
    if (cfg_.scan_to_bp >= 0 && bp >= 0 && bp > cfg_.scan_to_bp) return false;
    return true;
}

bool GemmaCalccor::normalizeVariant(const std::string& chr,
                                    const std::string& rs,
                                    long long bp,
                                    double cm,
                                    const std::string& allele1,
                                    const std::string& allele0,
                                    const Eigen::VectorXd& dosage,
                                    VariantRecord& rec) const {
    if (!passesIdFilters(rs) || !passesRegionFilters(chr, bp)) return false;

    const int n = static_cast<int>(dosage.size());
    int n_obs = 0;
    double sum = 0.0;
    long long hom0 = 0;
    long long het = 0;
    long long hom2 = 0;
    bool hard_calls = true;
    for (int i = 0; i < n; ++i) {
        const double g = dosage(i);
        if (is_missing(g)) continue;
        sum += g;
        ++n_obs;
        const double rounded = std::round(g);
        if (std::abs(g - rounded) <= 1e-6 && rounded >= 0.0 && rounded <= 2.0) {
            if (rounded < 0.5) ++hom0;
            else if (rounded < 1.5) ++het;
            else ++hom2;
        } else {
            hard_calls = false;
        }
    }
    if (n_obs == 0) return false;

    const double mean = sum / static_cast<double>(n_obs);
    const double af = mean / 2.0;
    const double maf = std::min(af, 1.0 - af);
    if (maf + 1e-12 < cfg_.min_maf) return false;
    if (cfg_.min_hwe > 0.0) {
        if (hard_calls && n_obs == n) {
            if (calculateHWE(hom0, het, hom2) < cfg_.min_hwe) return false;
        } else {
            static bool warned = false;
            if (!warned) {
                std::cerr << "Warning: --calccor applies --hwe only to fully observed hard-call genotypes; "
                          << "dosage/missing variants are not HWE-filtered.\n";
                warned = true;
            }
        }
    }

    Eigen::VectorXd centered(n);
    centered.setZero();
    double ss = 0.0;
    for (int i = 0; i < n; ++i) {
        if (is_missing(dosage(i))) {
            centered(i) = 0.0;
        } else {
            centered(i) = dosage(i) - mean;
            ss += centered(i) * centered(i);
        }
    }

    rec.chr = chr;
    rec.rs = rs;
    rec.bp = bp;
    rec.cm = cm;
    rec.allele1 = allele1;
    rec.allele0 = allele0;
    rec.n_miss = n - n_obs;
    rec.n_obs = n_obs;
    rec.af = maf;
    rec.var = ss / static_cast<double>(n);
    rec.centered = std::move(centered);
    return true;
}

std::vector<GemmaCalccor::VariantRecord> GemmaCalccor::loadVariants(
    PlinkReader& plink_reader,
    BgenReader& bgen_reader,
    PgenReader& pgen_reader,
    bool is_bgen,
    bool is_pgen) const {
    std::vector<VariantRecord> variants;

    if (is_pgen) {
        pgen_reader.reset();
        variants.reserve(pgen_reader.getNumVariants());
        std::string chr, rs, ref, alt;
        uint32_t pos = 0;
        Eigen::VectorXf dosage_f;
        while (pgen_reader.readVariant(chr, rs, pos, ref, alt, dosage_f)) {
            Eigen::VectorXd dosage = dosage_f.cast<double>();
            VariantRecord rec;
            if (normalizeVariant(chr, rs, pos, -9.0, alt, ref, dosage, rec)) {
                variants.push_back(std::move(rec));
            }
        }
    } else if (is_bgen) {
        bgen_reader.reset();
        variants.reserve(bgen_reader.getNumVariants());
        std::string chr, rs, a1, a2;
        uint32_t pos = 0;
        Eigen::VectorXf dosage_f;
        while (bgen_reader.readVariant(chr, rs, pos, a1, a2, dosage_f)) {
            Eigen::VectorXd dosage = dosage_f.cast<double>();
            VariantRecord rec;
            if (normalizeVariant(chr, rs, pos, -9.0, a2, a1, dosage, rec)) {
                variants.push_back(std::move(rec));
            }
        }
    } else {
        const auto& bim = plink_reader.getBimInfo();
        variants.reserve(bim.size());
        for (size_t i = 0; i < bim.size(); ++i) {
            Eigen::VectorXd dosage;
            plink_reader.readSnp(static_cast<int>(i), dosage);
            VariantRecord rec;
            if (normalizeVariant(bim[i].chrom, bim[i].id, bim[i].bp_pos, bim[i].cm_pos,
                                 bim[i].alt, bim[i].ref, dosage, rec)) {
                variants.push_back(std::move(rec));
            }
        }
    }

    return variants;
}

void GemmaCalccor::writeCor(const std::vector<VariantRecord>& variants) const {
    std::ofstream out(cfg_.out_prefix + ".cor.txt");
    if (!out) throw std::runtime_error("Cannot write GEMMA-style correlation output");

    out << "chr\trs\tps\tn_mis\tn_obs\tallele1\tallele0\taf\twindow_size\tvar\tcor\n";

    for (size_t i = 0; i < variants.size(); ++i) {
        const auto& vi = variants[i];
        std::vector<double> cors;
        const double vi_ss = dot_product(vi.centered, vi.centered);

        for (size_t j = i + 1; j < variants.size(); ++j) {
            const auto& vj = variants[j];
            if (vj.chr != vi.chr) break;
            if (cfg_.calccor_window_cm != 0.0 &&
                (vi.cm < 0.0 || vj.cm < 0.0 || (vj.cm - vi.cm) >= cfg_.calccor_window_cm)) {
                break;
            }
            if (cfg_.calccor_window_bp != 0 &&
                (vi.bp < 0 || vj.bp < 0 || (vj.bp - vi.bp) >= cfg_.calccor_window_bp)) {
                break;
            }
            if (cfg_.calccor_window_ns != 0 &&
                static_cast<int>(cors.size()) >= cfg_.calccor_window_ns) {
                break;
            }

            const double vj_ss = dot_product(vj.centered, vj.centered);
            double r = std::numeric_limits<double>::quiet_NaN();
            if (vi_ss > 0.0 && vj_ss > 0.0) {
                r = dot_product(vi.centered, vj.centered) / std::sqrt(vi_ss * vj_ss);
            }
            cors.push_back(r);
        }

        out << vi.chr << '\t' << vi.rs << '\t' << vi.bp << '\t'
            << vi.n_miss << '\t' << vi.n_obs << '\t'
            << vi.allele1 << '\t' << vi.allele0 << '\t'
            << std::fixed << std::setprecision(3) << vi.af << '\t'
            << cors.size() << '\t'
            << std::scientific << std::setprecision(6) << vi.var << '\t';

        if (cors.empty()) {
            out << "NA";
        } else {
            for (size_t k = 0; k < cors.size(); ++k) {
                if (k) out << ',';
                if (std::isfinite(cors[k])) out << cors[k];
                else out << "nan";
            }
        }
        out << '\n';
    }
}

void GemmaCalccor::run(PlinkReader& plink_reader,
                       BgenReader& bgen_reader,
                       PgenReader& pgen_reader,
                       bool is_bgen,
                       bool is_pgen) {
    if (cfg_.calccor_window_bp < 0) {
        throw std::runtime_error("--windowbp must be >= 0");
    }
    if (cfg_.calccor_window_cm < 0.0) {
        throw std::runtime_error("--windowcm must be >= 0");
    }
    if (cfg_.calccor_window_ns < 0) {
        throw std::runtime_error("--windowns must be >= 0");
    }
    if (cfg_.calccor_window_cm != 0.0 && (is_bgen || is_pgen)) {
        std::cerr << "Warning: --windowcm is only meaningful for PLINK BIM cM positions; "
                  << "BGEN/PGEN metadata do not currently provide cM, so no cM neighbours will be emitted.\n";
    }

    std::cout << "Calculating GEMMA-style SNP correlation windows...\n";
    auto variants = loadVariants(plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen);
    if (variants.empty()) throw std::runtime_error("--calccor found no usable SNPs after filters");
    writeCor(variants);
    std::cout << "GEMMA-style correlation output written to " << cfg_.out_prefix
              << ".cor.txt (" << variants.size() << " SNPs)\n";
}

} // namespace cosmic
