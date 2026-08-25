#include "rare_variant_test.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <Eigen/Eigenvalues>

namespace cosmic {

namespace {

double gamma_q_rv(double a, double x) {
    if (x < 0.0 || a <= 0.0) return 1.0;
    if (x == 0.0) return 1.0;
    constexpr int kMaxIter = 200;
    constexpr double kEps = 1e-12;
    constexpr double kTiny = 1e-300;

    if (x < a + 1.0) {
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int i = 0; i < kMaxIter; ++i) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::abs(del) < std::abs(sum) * kEps) {
                double p = sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
                return std::min(1.0, std::max(0.0, 1.0 - p));
            }
        }
    } else {
        double b = x + 1.0 - a;
        double c = 1.0 / kTiny;
        double d = 1.0 / std::max(b, kTiny);
        double h = d;
        for (int i = 1; i <= kMaxIter; ++i) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (std::abs(d) < kTiny) d = kTiny;
            c = b + an / c;
            if (std::abs(c) < kTiny) c = kTiny;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (std::abs(del - 1.0) < kEps) {
                double q = std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
                return std::min(1.0, std::max(0.0, q));
            }
        }
    }
    return 1.0;
}

} // namespace

// --- SetFileReader Implementation ---
std::vector<VariantSet> SetFileReader::read(const std::string& filename) {
    std::vector<VariantSet> sets;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open SetFile: " + filename);
    }

    std::unordered_map<std::string, std::vector<std::string>> set_map;
    std::string line, set_id, snp_id;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (iss >> set_id >> snp_id) {
            set_map[set_id].push_back(snp_id);
        }
    }

    for (const auto& kv : set_map) {
        VariantSet vs;
        vs.set_id = kv.first;
        vs.snp_ids = kv.second;
        sets.push_back(vs);
    }

    std::cout << "[CosmicScanner] Loaded " << sets.size() << " gene/region sets from " << filename << std::endl;
    return sets;
}

// --- RareVariantTest Implementation ---
RareVariantTest::RareVariantTest(Method method) : method_(method), beta_a_(1.0), beta_b_(25.0) {}

void RareVariantTest::setBetaWeights(double a, double b) {
    beta_a_ = a;
    beta_b_ = b;
}

Eigen::VectorXd RareVariantTest::calculateWeights(const Eigen::VectorXd& MAF) const {
    // Weight = BetaPDF(MAF, a, b)
    int k = MAF.size();
    Eigen::VectorXd W(k);
    for (int i = 0; i < k; ++i) {
        double f = MAF(i);
        if (f < 1e-6) f = 1e-6; // avoid log(0)
        if (f > 0.5) f = 0.5;   // force minor
        
        // Beta unnormalized weight: f^(a-1) * (1-f)^(b-1)
        double w = std::pow(f, beta_a_ - 1.0) * std::pow(1.0 - f, beta_b_ - 1.0);
        W(i) = w;
    }
    return W;
}

double RareVariantTest::mixtureChiSqPval(double Q, const Eigen::VectorXd& lambda) const {
    // Satterthwaite Approximation: c * ChiSq(df)
    // E[Q] = sum(lambda)
    // Var[Q] = 2 * sum(lambda^2)
    
    double mean_Q = lambda.sum();
    double var_Q = 2.0 * lambda.squaredNorm();
    
    if (var_Q <= 0.0) return 1.0;
    
    double df = (2.0 * mean_Q * mean_Q) / var_Q;
    double c = var_Q / (2.0 * mean_Q);
    
    double scaled_Q = Q / c;
    
    if (scaled_Q <= 0) return 1.0;
    return gamma_q_rv(0.5 * df, 0.5 * scaled_Q);
}

RareVariantTest::Result RareVariantTest::runTest(const std::string& set_id, 
                                                 const Eigen::MatrixXd& X_set, 
                                                 const Eigen::VectorXd& score_vec, 
                                                 const Eigen::MatrixXd& P_matrix, 
                                                 const Eigen::VectorXd& MAF) {
    // Calculate Weights W = diag(w_1 ... w_k)
    Eigen::VectorXd W = calculateWeights(MAF);
    
    if (method_ == Method::BURDEN) {
        return runBurden(set_id, X_set, score_vec, P_matrix, W);
    } else if (method_ == Method::SKAT) {
        return runSKAT(set_id, X_set, score_vec, P_matrix, W);
    }
    
    // SKAT-O skeleton fallback
    return runSKAT(set_id, X_set, score_vec, P_matrix, W);
}

RareVariantTest::Result RareVariantTest::runBurden(const std::string& set_id, 
                                                   const Eigen::MatrixXd& X_set, 
                                                   const Eigen::VectorXd& score_vec, 
                                                   const Eigen::MatrixXd& P_matrix, 
                                                   const Eigen::VectorXd& W) const {
    Result res;
    res.set_id = set_id;
    res.n_snps = X_set.cols();
    
    // Burden score: S = sum(w_i * U_i) = W' * U
    double burden_score = W.dot(score_vec);
    
    // Burden Variance: V = W' * X' * P * X * W
    Eigen::VectorXd XW = X_set * W;
    double var_burden = XW.transpose() * P_matrix * XW;
    
    if (var_burden <= 0.0) return res;
    
    // Z-test
    double z = burden_score / std::sqrt(var_burden);
    res.stat = z * z;
    res.p_value = std::erfc(std::abs(z) / std::sqrt(2.0));
    
    return res;
}

RareVariantTest::Result RareVariantTest::runSKAT(const std::string& set_id, 
                                                 const Eigen::MatrixXd& X_set, 
                                                 const Eigen::VectorXd& score_vec, 
                                                 const Eigen::MatrixXd& P_matrix, 
                                                 const Eigen::VectorXd& W) const {
    Result res;
    res.set_id = set_id;
    res.n_snps = X_set.cols();
    
    // SKAT Statistic: Q = sum(w_i^2 * U_i^2) = U' W W U
    Eigen::VectorXd WU = W.cwiseProduct(score_vec);
    double Q = WU.squaredNorm();
    res.stat = Q;
    
    // Asymptotic Distribution: Q ~ sum(lambda_i * ChiSq_1)
    // where lambda_i are eigenvalues of V = W * X' * P * X * W
    // V is K x K (K = number of SNPs)
    
    // Construct V
    Eigen::MatrixXd X_W = X_set * W.asDiagonal(); // N x K
    Eigen::MatrixXd V = X_W.transpose() * P_matrix * X_W; // K x K
    
    // Compute eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(V);
    Eigen::VectorXd lambda = es.eigenvalues();
    
    // Filter out negative or near-zero eigenvalues due to numerical instability
    Eigen::VectorXd positive_lambda(lambda.size());
    int valid_count = 0;
    for (int i = 0; i < lambda.size(); ++i) {
        if (lambda(i) > 1e-8) {
            positive_lambda(valid_count++) = lambda(i);
        }
    }
    
    if (valid_count == 0) return res;
    
    positive_lambda.conservativeResize(valid_count);
    
    // Get P-value from mixture of Chi-squares
    res.p_value = mixtureChiSqPval(Q, positive_lambda);
    
    return res;
}

} // namespace cosmic
