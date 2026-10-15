#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "dominance_engine.h"
#include "eigen_vce_d.h"
#include "grm_builder.h"
#include "plink_reader.h"

using namespace cosmic;

namespace {

int failures = 0;

#define CHECK(condition, message)                                                  \
    do {                                                                           \
        if (!(condition)) {                                                         \
            ++failures;                                                             \
            std::cerr << "FAIL: " << message << "\n";                            \
        }                                                                           \
    } while (0)

void writeBed(const std::string& path, const std::vector<int>& genotypes,
              int n, int m) {
    std::ofstream output(path, std::ios::binary);
    const unsigned char magic[] = {0x6c, 0x1b, 0x01};
    output.write(reinterpret_cast<const char*>(magic), 3);
    for (int j = 0; j < m; ++j) {
        std::vector<unsigned char> bytes((n + 3) / 4, 0);
        for (int i = 0; i < n; ++i) {
            const int genotype = genotypes[static_cast<size_t>(j) * n + i];
            const unsigned char code = genotype < 0 ? 1 :
                (genotype == 0 ? 3 : (genotype == 1 ? 2 : 0));
            bytes[i / 4] |= code << (2 * (i % 4));
        }
        output.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    }
}

std::string writeFixture(const std::filesystem::path& directory,
                         std::vector<int>& genotypes, int n, int m) {
    std::filesystem::create_directories(directory);
    std::mt19937 rng(24680);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    genotypes.assign(static_cast<size_t>(n) * m, 0);
    for (int j = 0; j < m; ++j) {
        const double p = 0.15 + 0.3 * (j % 5) / 4.0;
        for (int i = 0; i < n; ++i) {
            const double u = uniform(rng);
            genotypes[static_cast<size_t>(j) * n + i] =
                u < p * p ? 2 : (u < p * p + 2.0 * p * (1.0 - p) ? 1 : 0);
        }
    }
    for (int i = 0; i < n; ++i) {
        genotypes[i] = 0;                         // SNP 0: monomorphic.
        genotypes[static_cast<size_t>(n) + i] = 1; // SNP 1: zero centered dominance column.
    }
    for (int i = 0; i < 8; ++i) {
        genotypes[static_cast<size_t>(2) * n + i] = -1; // SNP 2: missing calls.
    }
    for (int i = 0; i < n; ++i) {
        genotypes[static_cast<size_t>(3) * n + i] = i % 3; // Stable causal marker.
    }

    const std::string prefix = (directory / "dominance_scan").string();
    writeBed(prefix + ".bed", genotypes, n, m);
    {
        std::ofstream bim(prefix + ".bim");
        for (int j = 0; j < m; ++j) {
            bim << (j < m / 2 ? 1 : 2) << "\tsnp" << j << "\t0\t"
                << (1000 + j * 10) << "\tA\tG\n";
        }
    }
    {
        std::ofstream fam(prefix + ".fam");
        for (int i = 0; i < n; ++i) {
            fam << "F" << i << "\tI" << i << "\t0\t0\t0\t-9\n";
        }
    }
    return prefix;
}

Eigen::MatrixXd manualVitezica(const std::vector<int>& genotypes, int n, int m,
                               double min_maf) {
    Eigen::MatrixXd raw = Eigen::MatrixXd::Zero(n, n);
    double denominator = 0.0;
    for (int j = 0; j < m; ++j) {
        int observed = 0;
        double dosage_sum = 0.0;
        for (int i = 0; i < n; ++i) {
            const int genotype = genotypes[static_cast<size_t>(j) * n + i];
            if (genotype >= 0) {
                dosage_sum += genotype;
                ++observed;
            }
        }
        if (observed == 0) continue;
        const double p = dosage_sum / (2.0 * observed);
        const double q = 1.0 - p;
        if (std::min(p, q) < min_maf || std::min(p, q) < 1e-8) continue;

        Eigen::VectorXd z(n);
        for (int i = 0; i < n; ++i) {
            const int genotype = genotypes[static_cast<size_t>(j) * n + i];
            z[i] = genotype < 0 ? 0.0 :
                (genotype == 0 ? -2.0 * p * p :
                 (genotype == 1 ? 2.0 * p * q : -2.0 * q * q));
        }
        raw.noalias() += z * z.transpose();
        denominator += 4.0 * p * p * q * q;
    }
    return raw / denominator;
}

Eigen::VectorXd preparedDominance(const std::vector<int>& genotypes, int marker,
                                  int n) {
    Eigen::VectorXd result(n);
    double observed_sum = 0.0;
    int observed = 0;
    for (int i = 0; i < n; ++i) {
        const int genotype = genotypes[static_cast<size_t>(marker) * n + i];
        if (genotype >= 0) {
            observed_sum += genotype == 1 ? 1.0 : 0.0;
            ++observed;
        }
    }
    const double mean = observed_sum / observed;
    for (int i = 0; i < n; ++i) {
        const int genotype = genotypes[static_cast<size_t>(marker) * n + i];
        result[i] = genotype < 0 ? 0.0 : (genotype == 1 ? 1.0 - mean : -mean);
    }
    return result;
}

void directGls(const Eigen::VectorXd& y, const Eigen::MatrixXd& X,
               const Eigen::MatrixXd& K, const Eigen::VectorXd& marker,
               double sigma_d2, double sigma_e2,
               double& beta, double& se) {
    Eigen::MatrixXd V = sigma_d2 * K;
    V.diagonal().array() += sigma_e2;
    Eigen::LDLT<Eigen::MatrixXd> ldlt(V);
    const Eigen::MatrixXd V_inv_X = ldlt.solve(X);
    const Eigen::VectorXd V_inv_y = ldlt.solve(y);
    const Eigen::VectorXd V_inv_g = ldlt.solve(marker);
    const Eigen::MatrixXd middle = (X.transpose() * V_inv_X).inverse();
    const Eigen::VectorXd py = V_inv_y - V_inv_X * middle * (X.transpose() * V_inv_y);
    const Eigen::VectorXd pg = V_inv_g - V_inv_X * middle * (X.transpose() * V_inv_g);
    const double information = marker.dot(pg);
    beta = marker.dot(py) / information;
    se = std::sqrt(1.0 / information);
}

} // namespace

int main() {
    constexpr int n = 96;
    constexpr int m = 24;
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() / "cosmic_dominance_scan_unit";
    std::filesystem::remove_all(directory);

    std::vector<int> genotypes;
    const std::string prefix = writeFixture(directory, genotypes, n, m);
    PlinkReader reader;
    reader.load(prefix);

    GrmBuilder::Options grm_options;
    grm_options.algorithm = 4;
    grm_options.min_maf = 0.01f;
    grm_options.write_gcta_format = false;
    GrmBuilder builder(grm_options);
    builder.compute(reader);
    const Eigen::MatrixXd K_D = builder.getGrm();
    const Eigen::MatrixXd expected_K_D = manualVitezica(genotypes, n, m, 0.01);
    const double kd_max_diff = (K_D - expected_K_D).cwiseAbs().maxCoeff();
    std::cout << "algorithm4_max_abs_diff=" << kd_max_diff << "\n";
    CHECK(kd_max_diff < 1e-6,
          "production algorithm 4 must equal Vitezica coding normalized by sum 4p^2q^2");

    Eigen::MatrixXf X(n, 2);
    X.col(0).setOnes();
    for (int i = 0; i < n; ++i) X(i, 1) = static_cast<float>((i % 7) - 3);
    const Eigen::VectorXd causal = preparedDominance(genotypes, 3, n);
    std::mt19937 rng(13579);
    std::normal_distribution<double> normal(0.0, 0.35);
    Eigen::VectorXd y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = 1.25 * causal[i] + 0.15 * X(i, 1) + normal(rng);
    }

    const std::vector<int> keep = [&]() {
        std::vector<int> indices(n);
        for (int i = 0; i < n; ++i) indices[i] = i;
        return indices;
    }();

    const auto run_scan = [&](double sigma_d2, double sigma_e2) {
        EigenVceD::Options options;
        options.verbose = false;
        EigenVceD vce(options);
        vce.prepare(y, X, K_D);
        vce.fixedComponents(sigma_d2, sigma_e2);
        LmmEngineD engine;
        engine.prepareFromVce(vce);
        std::vector<DomResult> results;
        engine.scan(reader, results, keep);
        return results;
    };

    for (const auto& components : std::vector<std::pair<double, double>>{{0.3, 0.7}, {0.0, 1.0}}) {
        const auto results = run_scan(components.first, components.second);
        CHECK(results.size() == m, "scan result vector must preserve BIM indexing");
        CHECK(std::isnan(results[0].p_dom), "monomorphic marker must be filtered");
        CHECK(std::isnan(results[1].p_dom), "constant heterozygote indicator must be rank deficient");
        CHECK(std::isfinite(results[2].p_dom), "marker with missing calls must be mean-imputed and tested");
        CHECK(std::isfinite(results[3].p_dom), "causal marker must be tested");

        double expected_beta = 0.0;
        double expected_se = 0.0;
        directGls(y, X.cast<double>(), K_D, causal, components.first,
                  components.second, expected_beta, expected_se);
        CHECK(std::abs(results[3].beta_dom - expected_beta) < 1e-10,
              "scan beta must match direct GLS");
        CHECK(std::abs(results[3].se_dom - expected_se) < 1e-10,
              "scan SE must match direct GLS");
        std::cout << "sigma_d2=" << components.first << " causal_beta="
                  << results[3].beta_dom << " causal_se=" << results[3].se_dom
                  << " causal_p=" << results[3].p_dom << "\n";
        CHECK(results[3].p_dom < 0.01, "pure dominance QTL must be detected");
    }

    std::filesystem::remove_all(directory);
    if (failures != 0) {
        std::cerr << failures << " dominance scan assertion(s) failed\n";
        return 1;
    }
    std::cout << "Dominance scan production-reference tests passed\n";
    return 0;
}
