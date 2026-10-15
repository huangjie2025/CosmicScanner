#include "models/multivariate_model.h"
#include "models/multivariate_helpers.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace cosmic {

std::ofstream mv_open_output_file(const std::string& path, std::ios_base::openmode mode) {
    std::ofstream ofs(path, mode);
    if (!ofs) throw std::runtime_error("Failed to open output file: " + path);
    return ofs;
}

void mv_validate_output_parent(const std::string& out_prefix) {
    std::filesystem::path parent = std::filesystem::path(out_prefix).parent_path();
    if (!parent.empty() && !std::filesystem::exists(parent)) {
        throw std::runtime_error("Output directory does not exist: " + parent.string());
    }
}

void mv_write_vbeta_header(std::ostream& os, int n_traits) {
    for (int i = 0; i < n_traits; ++i) {
        for (int j = i; j < n_traits; ++j) os << "\tVBETA_" << (i + 1) << "_" << (j + 1);
    }
}

void mv_write_vbeta_values(std::ostream& os, const Eigen::MatrixXd& vbeta, int n_traits) {
    for (int i = 0; i < n_traits; ++i) {
        for (int j = i; j < n_traits; ++j) {
            double value = (i < vbeta.rows() && j < vbeta.cols()) ? vbeta(i, j) : std::numeric_limits<double>::quiet_NaN();
            os << "\t" << value;
        }
    }
}

void MultivariateModel::writeResults(const std::string& out_file, const std::vector<MvGwasResult>& results) {
    int n_traits = Y_.cols();
    auto ofs = mv_open_output_file(out_file, std::ios::out | std::ios::app);
    for (const auto& res : results) {
        if (std::isnan(res.p_value)) continue;
        ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
            << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
            << res.wald_stat << "\t" << res.p_value;
        for (int t = 0; t < n_traits; ++t) {
            double beta = res.beta(t);
            const double variance = res.se(t, t);
            if (!std::isfinite(variance) || variance < 0.0) {
                throw std::runtime_error("mvLMM covariance TRUE_SINGULARITY for SNP " + res.snp_id +
                                         ": invalid trait variance");
            }
            double se = std::sqrt(variance);
            ofs << "\t" << beta << "\t" << se;
        }
        if (options_.output_vbeta) mv_write_vbeta_values(ofs, res.se, n_traits);
        if (options_.use_eb) {
            double eb_confidence = std::isfinite(res.eb_confidence) ? res.eb_confidence : 0.0;
            eb_confidence = std::max(0.0, std::min((double)n_traits, eb_confidence));
            eb_effective_tests_ += eb_confidence;
            eb_tested_markers_++;
            ofs << "\t" << res.eb_marker_lambda << "\t" << eb_confidence;
        }
        ofs << "\n";
    }
}

void MultivariateModel::writeEmmvgwasOutputs(const MvEigenVce& vce,
                                             const std::vector<MvGwasResult>& results) const {
    if (!options_.use_emmvgwas_lambda) throw std::runtime_error("--mv-emmvgwas-output requires --mv-compatible emmvgwas");

    const int n_traits = Y_.cols();
    auto write_matrix_csv = [](const std::string& path, const Eigen::MatrixXd& M) {
        auto ofs = mv_open_output_file(path);
        ofs << std::setprecision(12);
        for (int j = 0; j < M.cols(); ++j) ofs << (j ? "," : "") << "V" << (j + 1);
        ofs << "\n";
        for (int i = 0; i < M.rows(); ++i) {
            for (int j = 0; j < M.cols(); ++j) ofs << (j ? "," : "") << M(i, j);
            ofs << "\n";
        }
    };

    write_matrix_csv(options_.out_prefix + "G.csv", vce.getEmmvgwasG());
    write_matrix_csv(options_.out_prefix + "R.csv", vce.getEmmvgwasR());
    {
        auto ofs = mv_open_output_file(options_.out_prefix + "beta.csv");
        ofs << std::setprecision(12) << "V1\n";
        Eigen::VectorXd beta = vce.getEmmvgwasBeta();
        for (int i = 0; i < beta.size(); ++i) ofs << beta(i) << "\n";
    }
    {
        auto ofs = mv_open_output_file(options_.out_prefix + "_association.csv");
        ofs << std::setprecision(12);
        ofs << "snp";
        for (int t = 0; t < n_traits; ++t) ofs << ",trait" << (t + 1);
        ofs << ",wald,p,log10p\n";
        int snp_index = 1;
        for (const auto& res : results) {
            if (!std::isfinite(res.p_value)) { ++snp_index; continue; }
            ofs << snp_index;
            for (int t = 0; t < n_traits; ++t) ofs << "," << res.beta(t);
            double log10p = res.p_value > 0.0 ? -std::log10(res.p_value) : std::numeric_limits<double>::infinity();
            ofs << "," << res.wald_stat << "," << res.p_value << "," << log10p << "\n";
            ++snp_index;
        }
    }
    {
        auto ofs = mv_open_output_file(options_.out_prefix + "_association_with_ids.csv");
        ofs << std::setprecision(12);
        ofs << "snp,SNP_ID,CHR,BP,A1,A2";
        for (int t = 0; t < n_traits; ++t) ofs << ",trait" << (t + 1);
        ofs << ",wald,p,log10p\n";
        int snp_index = 1;
        for (const auto& res : results) {
            if (!std::isfinite(res.p_value)) { ++snp_index; continue; }
            double log10p = res.p_value > 0.0 ? -std::log10(res.p_value) : std::numeric_limits<double>::infinity();
            ofs << snp_index << "," << res.snp_id << "," << res.chrom << "," << res.pos
                << "," << res.a1 << "," << res.a2;
            for (int t = 0; t < n_traits; ++t) ofs << "," << res.beta(t);
            ofs << "," << res.wald_stat << "," << res.p_value << "," << log10p << "\n";
            ++snp_index;
        }
    }
}

void MultivariateModel::writeEbSummary(const std::string& out_file) const {
    std::ofstream ofs(out_file);
    if (!ofs) throw std::runtime_error("Cannot write multivariate EB summary: " + out_file);
    ofs << std::setprecision(12);
    double m_eff = std::max(eb_effective_tests_, 1.0);
    ofs << "Metric\tValue\n";
    ofs << "method\tEB-mvLMM random-SNP scan (Phi=lambda*Ve)\n";
    ofs << "tested_markers\t" << eb_tested_markers_ << "\n";
    ofs << "traits\t" << Y_.cols() << "\n";
    ofs << "effective_tests\t" << eb_effective_tests_ << "\n";
    ofs << "effective_tests_per_marker\t" << (eb_tested_markers_ > 0 ? eb_effective_tests_ / (double)eb_tested_markers_ : 0.0) << "\n";
    ofs << "bonferroni_0.05\t" << (eb_tested_markers_ > 0 ? 0.05 / (double)eb_tested_markers_ : 0.0) << "\n";
    ofs << "eb_mbonf_0.05\t" << (0.05 / m_eff) << "\n";
    ofs << "note\tValidated summary-level EB shrinkage auxiliary scan; report standard mvLMM results and ordinary Bonferroni/FDR alongside EB-mBonf.\n";
}

} // namespace cosmic
