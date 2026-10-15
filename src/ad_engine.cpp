#include "ad_engine.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include "genetics/genotype_encoding.h"

namespace cosmic {

namespace {

// Chi-square survival function (copied from lmm_engine.cpp for self-containment).
static double gamma_q(double a, double x) {
    if (x < 0.0 || a <= 0.0) return 1.0;
    if (x < a + 1.0) {
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int i = 0; i < 100; ++i) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::abs(del) < std::abs(sum) * 1e-7) {
                return 1.0 - sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
            }
        }
    } else {
        double b = x + 1.0 - a;
        double c = 1.0 / 1e-30;
        double d = 1.0 / b;
        double h = d;
        for (int i = 1; i <= 100; ++i) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (std::abs(d) < 1e-30) d = 1e-30;
            c = b + an / c;
            if (std::abs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (std::abs(del - 1.0) < 1e-7) {
                return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
            }
        }
    }
    return 0.0;
}

static double chisq_prob(double stat, int df) {
    if (stat <= 0) return 1.0;
    return gamma_q(df / 2.0, stat / 2.0);
}

constexpr size_t kBlockSize = 512;

} // namespace

void LmmEngineAD::prepareFromVce(const EigenVceAD& vce) {
    sigma2_a = vce.getSigma2a();
    sigma2_d = vce.getSigma2d();
    sigma2_e = vce.getSigma2e();
    D_V = vce.getEigenvalues();
    W_proj = vce.getProjectedW();
    Uty = vce.getUty();
    UtX = vce.getUtX();
}

LmmEngineAD::ScanWorkspace LmmEngineAD::prepareWorkspace() const {
    const int n_eig = static_cast<int>(D_V.size());
    ScanWorkspace workspace;
    workspace.v_inv_diag.resize(n_eig);
    for (int i = 0; i < n_eig; ++i) {
        const double d = std::max(D_V[i], 1e-20);
        workspace.v_inv_diag[i] = 1.0 / d;
    }

    // Null fixed-effect quantities (UtX is 0-column for the projected AD null;
    // kept general so the engine also works for a non-projected null).
    Eigen::MatrixXd Xt_Vinv_X_inv(UtX.cols(), UtX.cols());
    if (UtX.cols() > 0) {
        Eigen::MatrixXd weighted_UtX = UtX.array().colwise() * workspace.v_inv_diag.array();
        Eigen::LLT<Eigen::MatrixXd> lltX(UtX.transpose() * weighted_UtX);
        if (lltX.info() == Eigen::Success) {
            Xt_Vinv_X_inv = lltX.solve(
                Eigen::MatrixXd::Identity(UtX.cols(), UtX.cols()));
        } else {
            Xt_Vinv_X_inv = (UtX.transpose() * weighted_UtX).inverse();
        }
    }

    Eigen::VectorXd Xt_Vinv_y = Eigen::VectorXd::Zero(UtX.cols());
    if (UtX.cols() > 0) {
        for (int i = 0; i < n_eig; ++i) {
            Xt_Vinv_y += workspace.v_inv_diag[i] * (UtX.row(i).transpose() * Uty[i]);
        }
    }
    Eigen::VectorXd beta(UtX.cols());
    if (UtX.cols() > 0) beta = Xt_Vinv_X_inv * Xt_Vinv_y;

    Eigen::VectorXd resid = Uty;
    if (UtX.cols() > 0) resid -= UtX * beta;
    workspace.p_y = resid.cwiseProduct(workspace.v_inv_diag);
    return workspace;
}

void LmmEngineAD::evaluateBlock(const ScanWorkspace& workspace,
                                const Eigen::MatrixXd& ga_block,
                                const Eigen::MatrixXd& gd_block,
                                const std::vector<int>& snp_indices,
                                const std::vector<double>& freqs,
                                const std::vector<VariantMeta>& meta,
                                std::vector<AdResult>& results,
                                long long& tested,
                                long long& rank_deficient) const {
    const int b = static_cast<int>(snp_indices.size());
    if (b == 0) return;
    const Eigen::MatrixXd Ga_use = ga_block.leftCols(b);
    const Eigen::MatrixXd Gd_use = gd_block.leftCols(b);

    const Eigen::MatrixXd UtGa = W_proj.transpose() * Ga_use;
    const Eigen::MatrixXd UtGd = W_proj.transpose() * Gd_use;

    const Eigen::VectorXd score_a = UtGa.transpose() * workspace.p_y;
    const Eigen::VectorXd score_d = UtGd.transpose() * workspace.p_y;

    const Eigen::MatrixXd qA = UtGa.array().colwise() * workspace.v_inv_diag.array();
    const Eigen::MatrixXd qD = UtGd.array().colwise() * workspace.v_inv_diag.array();
    const Eigen::VectorXd A11 = (UtGa.array() * qA.array()).colwise().sum().transpose();
    const Eigen::VectorXd A22 = (UtGd.array() * qD.array()).colwise().sum().transpose();
    const Eigen::MatrixXd A12_full = UtGd.transpose() * qA;

    for (int col = 0; col < b; ++col) {
        const int j = snp_indices[col];
        const double i11 = A11(col);
        const double i22 = A22(col);
        const double i12 = A12_full(col, col);

        AdResult& res = results[j];
        res.chrom = meta[static_cast<size_t>(j)].chrom;
        res.snp_id = meta[static_cast<size_t>(j)].snp_id;
        res.pos = meta[static_cast<size_t>(j)].pos;
        res.a1 = meta[static_cast<size_t>(j)].a1;
        res.a2 = meta[static_cast<size_t>(j)].a2;
        res.freq = freqs[col];
        res.n_obs = static_cast<int>(ga_block.rows());

        const double det = i11 * i22 - i12 * i12;
        const double scale = std::max(i11 * i22, std::numeric_limits<double>::min());
        if (!(i11 > 1e-12) || !(i22 > 1e-12) || !(det > 1e-12 * scale) ||
            !std::isfinite(det)) {
            ++rank_deficient;
            res.rank_deficient = true;
            continue;
        }

        const double s1 = score_a(col);
        const double s2 = score_d(col);
        if (!std::isfinite(s1) || !std::isfinite(s2)) {
            ++rank_deficient;
            res.rank_deficient = true;
            continue;
        }

        const double beta_a = (i22 * s1 - i12 * s2) / det;
        const double beta_d = (-i12 * s1 + i11 * s2) / det;
        const double var_a = i22 / det;
        const double var_d = i11 / det;
        if (!(var_a > 0.0) || !(var_d > 0.0) || !std::isfinite(beta_a) ||
            !std::isfinite(beta_d)) {
            ++rank_deficient;
            res.rank_deficient = true;
            continue;
        }

        const double chisq_add = beta_a * beta_a / var_a;
        const double chisq_dom = beta_d * beta_d / var_d;
        const double chisq_joint = (i22 * s1 * s1 - 2.0 * i12 * s1 * s2 + i11 * s2 * s2) / det;

        res.beta_add = beta_a;
        res.se_add = std::sqrt(var_a);
        res.p_add = chisq_prob(chisq_add, 1);
        res.beta_dom = beta_d;
        res.se_dom = std::sqrt(var_d);
        res.p_dom = chisq_prob(chisq_dom, 1);
        res.log10_p_dom = (res.p_dom > 0.0) ? -std::log10(res.p_dom) : 300.0;
        res.joint_chisq = chisq_joint;
        res.p_joint = chisq_prob(chisq_joint, 2);
        res.log10_p_joint = (res.p_joint > 0.0) ? -std::log10(res.p_joint) : 300.0;
        res.rank_deficient = false;
        ++tested;
    }
}

void LmmEngineAD::scan(const PlinkReader& reader, std::vector<AdResult>& results,
                       const std::vector<int>& keep_indices,
                       const std::string& chrom_filter) {
    const int m = static_cast<int>(reader.getNumSnps());
    const int n = static_cast<int>(keep_indices.size());
    results.resize(m);

    const ScanWorkspace workspace = prepareWorkspace();

    const auto& bim = reader.getBimInfo();
    const uint8_t* raw_data = reader.getData();
    const size_t bytes_per_snp = reader.getBytesPerSnp();

    std::vector<VariantMeta> meta(m);
    for (int j = 0; j < m; ++j) {
        meta[static_cast<size_t>(j)].chrom = bim[j].chrom;
        meta[static_cast<size_t>(j)].snp_id = bim[j].id;
        meta[static_cast<size_t>(j)].a1 = bim[j].alt;
        meta[static_cast<size_t>(j)].a2 = bim[j].ref;
        meta[static_cast<size_t>(j)].pos = bim[j].bp_pos;
    }

    Eigen::MatrixXd Ga_block(n, kBlockSize);
    Eigen::MatrixXd Gd_block(n, kBlockSize);
    std::vector<int> snp_indices;
    std::vector<double> freqs;
    snp_indices.reserve(kBlockSize);
    freqs.reserve(kBlockSize);

    long long rank_deficient = 0;
    long long tested = 0;

    const auto flush_block = [&]() {
        evaluateBlock(workspace, Ga_block, Gd_block, snp_indices, freqs, meta,
                      results, tested, rank_deficient);
        snp_indices.clear();
        freqs.clear();
    };

    for (int j = 0; j < m; ++j) {
        if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;
        if (options.scan_from_bp >= 0 && bim[j].bp_pos < options.scan_from_bp) continue;
        if (options.scan_to_bp >= 0 && bim[j].bp_pos > options.scan_to_bp) continue;
        if (!options.extract_set.empty() && options.extract_set.find(bim[j].id) == options.extract_set.end()) continue;
        if (!options.exclude_set.empty() && options.exclude_set.find(bim[j].id) != options.exclude_set.end()) continue;
        if (snp_indices.size() == kBlockSize) flush_block();

        Eigen::VectorXd marker_raw(n);
        const uint8_t* snp_data = raw_data + 3 + static_cast<size_t>(j) * bytes_per_snp;
        for (int i = 0; i < n; ++i) {
            const int sample_idx = keep_indices[i];
            const uint8_t code = (snp_data[sample_idx >> 2] >> ((sample_idx & 3) * 2)) & 3;
            if (code == 1) {
                marker_raw[i] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }
            marker_raw[i] = (code == 0) ? 2.0 : (code == 2 ? 1.0 : 0.0);
        }

        const auto prepared_add = genetics::prepareAssociationMarker(
            marker_raw, "additive", options.min_maf, options.min_hwe);
        if (!prepared_add.passes_maf || !prepared_add.passes_hwe) continue;
        // Dominance design column requires hardcalls; throws otherwise.
        const auto prepared_dom = genetics::prepareAssociationMarker(
            marker_raw, "heterozygote", options.min_maf, options.min_hwe);

        const int col = static_cast<int>(snp_indices.size());
        Ga_block.col(col) = prepared_add.values;
        Gd_block.col(col) = prepared_dom.values;
        snp_indices.push_back(j);
        freqs.push_back(prepared_add.allele_frequency);
    }
    flush_block();

    std::cout << "  [AD scan] tested=" << tested << ", rank-deficient/excluded="
              << rank_deficient << std::endl;
}

void LmmEngineAD::scan(PgenReader& reader, std::vector<AdResult>& results,
                       const std::vector<int>& keep_indices,
                       const std::string& chrom_filter) {
    const int m = static_cast<int>(reader.getNumVariants());
    const int n = static_cast<int>(keep_indices.size());
    results.resize(m);

    const ScanWorkspace workspace = prepareWorkspace();

    std::vector<VariantMeta> meta(m);
    Eigen::MatrixXd Ga_block(n, kBlockSize);
    Eigen::MatrixXd Gd_block(n, kBlockSize);
    std::vector<int> snp_indices;
    std::vector<double> freqs;
    snp_indices.reserve(kBlockSize);
    freqs.reserve(kBlockSize);

    long long rank_deficient = 0;
    long long tested = 0;

    const auto flush_block = [&]() {
        evaluateBlock(workspace, Ga_block, Gd_block, snp_indices, freqs, meta,
                      results, tested, rank_deficient);
        snp_indices.clear();
        freqs.clear();
    };

    // PGEN decodes exact hardcall allele counts per variant: 0 = hom A2,
    // 1 = het, 2 = hom A1 (counting the ALT/A1 allele), NaN = missing.
    // These are hardcalls, not expected dosages, so the existing PGEN additive
    // definition and the BED dominance encoding both apply unchanged.
    reader.reset();
    std::string chrom, rsid, ref, alt;
    uint32_t pos = 0;
    Eigen::VectorXf dosage;

    for (int j = 0; j < m; ++j) {
        if (!reader.readVariant(chrom, rsid, pos, ref, alt, dosage)) break;
        meta[static_cast<size_t>(j)].chrom = chrom;
        meta[static_cast<size_t>(j)].snp_id = rsid;
        meta[static_cast<size_t>(j)].a1 = alt;
        meta[static_cast<size_t>(j)].a2 = ref;
        meta[static_cast<size_t>(j)].pos = pos;

        if (!chrom_filter.empty() && chrom != chrom_filter) continue;
        if (options.scan_from_bp >= 0 &&
            static_cast<long long>(pos) < options.scan_from_bp) continue;
        if (options.scan_to_bp >= 0 &&
            static_cast<long long>(pos) > options.scan_to_bp) continue;
        if (!options.extract_set.empty() && options.extract_set.find(rsid) == options.extract_set.end()) continue;
        if (!options.exclude_set.empty() && options.exclude_set.find(rsid) != options.exclude_set.end()) continue;
        if (snp_indices.size() == kBlockSize) flush_block();

        Eigen::VectorXd marker_raw(n);
        for (int i = 0; i < n; ++i) {
            const float value = dosage[keep_indices[i]];
            marker_raw[i] = std::isnan(value)
                ? std::numeric_limits<double>::quiet_NaN()
                : static_cast<double>(value);
        }

        const auto prepared_add = genetics::prepareAssociationMarker(
            marker_raw, "additive", options.min_maf, options.min_hwe);
        if (!prepared_add.passes_maf || !prepared_add.passes_hwe) continue;
        const auto prepared_dom = genetics::prepareAssociationMarker(
            marker_raw, "heterozygote", options.min_maf, options.min_hwe);

        const int col = static_cast<int>(snp_indices.size());
        Ga_block.col(col) = prepared_add.values;
        Gd_block.col(col) = prepared_dom.values;
        snp_indices.push_back(j);
        freqs.push_back(prepared_add.allele_frequency);
    }
    flush_block();

    std::cout << "  [AD scan] tested=" << tested << ", rank-deficient/excluded="
              << rank_deficient << std::endl;
}

} // namespace cosmic
