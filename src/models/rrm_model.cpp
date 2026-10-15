#include "models/rrm_model.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <unordered_map>
#include <omp.h>
#include "pheno_reader.h"
#include "data_loader.h"
#include "gwas_scanner.h"

namespace cosmic {

RrmModel::Options RrmModel::Options::fromConfig(const GwasConfig& cfg) {
    Options opts;
    opts.out_prefix = cfg.out_prefix;
    opts.min_maf = cfg.min_maf;
    opts.min_hwe = cfg.min_hwe;
    opts.grm_alg = cfg.grm_alg;
    opts.n_threads = 1;
    #ifdef _OPENMP
    opts.n_threads = omp_get_max_threads();
    #endif
    opts.no_loco = cfg.no_loco;
    opts.fast_lmm = cfg.fast_lmm;
    opts.use_spa = cfg.use_spa;
    opts.rrm_order = cfg.rrm_order;
    opts.rrm_curve_grid = cfg.rrm_curve_grid;
    opts.rrm_top = cfg.rrm_top;
    opts.rrm_ai_iter = cfg.rrm_ai_iter;
    opts.rrm_sig_threshold = cfg.rrm_sig_threshold;
    opts.rrm_loco_global_null = cfg.rrm_loco_global_null;
    opts.rrm_basis = cfg.rrm_basis;
    opts.rrm_gwas_method = cfg.rrm_gwas_method;
    opts.tpheno = cfg.tpheno;
    opts.mpheno = cfg.mpheno;
    opts.long_pheno_file = cfg.long_pheno_file;
    opts.pheno_file = cfg.pheno_file;
    opts.bfile_prefix = cfg.bfile_prefix;
    opts.rrm_unbalance = cfg.rrm_unbalance;
    opts.rrm_compatible = cfg.rrm_compatible;
    opts.rrm_compatible_varcom = cfg.rrm_compatible_varcom;
    opts.mpheno_list = cfg.mpheno_list;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.tpoint_col_name = cfg.tpoint_col_name;
    opts.tclass_names = cfg.tclass_names;
    opts.tcovar_names = cfg.tcovar_names;
    return opts;
}

RrmModel::RrmModel(const Options& opts) : options_(opts) {}

void RrmModel::prepare(const std::vector<std::string>& fids,
                              const std::vector<std::string>& iids,
                              int n_samples,
                              PlinkReader& plink_reader,
                              BgenReader& bgen_reader,
                              PgenReader& pgen_reader,
                              VcfReader& vcf_reader,
                              bool is_bgen, bool is_pgen, bool is_vcf) {
    fids_ = fids;
    iids_ = iids;
    n_samples_ = n_samples;
    plink_reader_ = &plink_reader;
    bgen_reader_ = &bgen_reader;
    pgen_reader_ = &pgen_reader;
    vcf_reader_ = &vcf_reader;
    is_bgen_ = is_bgen;
    is_pgen_ = is_pgen;
    is_vcf_ = is_vcf;

    if (options_.rrm_compatible == "gma-fixreg") {
        if (options_.pheno_file.empty()) throw std::runtime_error("--rrm-compatible gma-fixreg requires --pheno.");
        if (options_.mpheno_list.empty()) throw std::runtime_error("--rrm-compatible gma-fixreg requires --mpheno-list.");
        if (options_.rrm_compatible_varcom.empty()) throw std::runtime_error("--rrm-compatible gma-fixreg requires --rrm-compatible-varcom.");
        if (is_bgen_ || is_pgen_) throw std::runtime_error("--rrm-compatible gma-fixreg currently supports PLINK --bfile only.");
        std::cout << "Loading balanced phenotype matrix for GMA-compatible fixreg RRM...\n";
        PhenoReader pheno_reader;
        pheno_reader.loadMatrix(options_.pheno_file, options_.mpheno_list);
        size_t matched = pheno_reader.alignMatrix(fids_, iids_, Y_wide_);
        std::cout << "Matched " << matched << " complete balanced phenotype records.\n";
        keep_.keep_indices.clear();
        for (int i = 0; i < n_samples_; ++i) {
            if (i < Y_wide_.rows() && Y_wide_.row(i).allFinite()) keep_.keep_indices.push_back(i);
        }
        keep_.n_keep = (int)keep_.keep_indices.size();
        if (keep_.n_keep == 0) throw std::runtime_error("No complete phenotype records for GMA-compatible fixreg mode.");
        return;
    }

    if (options_.rrm_compatible == "gmat-unbalance") {
        options_.rrm_unbalance = true;
        options_.rrm_gwas_method = "coefficient";
    } else if (!options_.rrm_compatible.empty()) {
        throw std::runtime_error("Unsupported --rrm-compatible mode: " + options_.rrm_compatible);
    }

    if (options_.rrm_unbalance) {
        if (options_.pheno_file.empty()) throw std::runtime_error("--rrm-unbalance requires --pheno file.");
        if (options_.tpoint_col_name.empty() && options_.tpheno < 1) throw std::runtime_error("--rrm-unbalance requires --tpoint <col_name> or --tpheno <col>.");
        std::cout << "Loading unbalanced longitudinal data for Random Regression Model...\n";
        int tcol = options_.tpheno;
        if (!options_.tpoint_col_name.empty() && tcol < 2) {
            std::ifstream ifs(options_.pheno_file);
            std::string header;
            std::getline(ifs, header);
            std::stringstream hss(header);
            std::vector<std::string> cols;
            std::string tok;
            while (hss >> tok) cols.push_back(tok);
            for (size_t i = 2; i < cols.size(); ++i) {
                if (cols[i] == options_.tpoint_col_name) { tcol = (int)i + 1; break; }
            }
            if (tcol < 2) throw std::runtime_error("--tpoint column '" + options_.tpoint_col_name + "' not found in phenotype file.");
        }
        if (options_.mpheno < 1) {
            throw std::runtime_error("--mpheno must specify phenotype value column (1-indexed, after FID IID).");
        }
        PhenoReader pheno_reader;
        records_ = pheno_reader.loadLongitudinalRecords(options_.pheno_file, tcol, options_.mpheno);
        X_base_ = DataLoader::loadCovariates(options_.qcovar_file, options_.covar_file, fids_, iids_, n_samples_);
        keep_.keep_indices.resize(n_samples_);
        for (int i = 0; i < n_samples_; ++i) keep_.keep_indices[i] = i;
        keep_.n_keep = n_samples_;
        return;
    }

    if (options_.long_pheno_file.empty()) throw std::runtime_error("RRM model requires --long-pheno file.");

    std::cout << "Loading longitudinal phenotype data for Random Regression Model...\n";
    PhenoReader pheno_reader;
    records_ = pheno_reader.loadLongitudinalRecords(options_.long_pheno_file, options_.tpheno, options_.mpheno);

    X_base_ = DataLoader::loadCovariates(options_.qcovar_file, options_.covar_file, fids_, iids_, n_samples_);
    std::cout << "Note: only time-invariant covariates are supported in the current RRM GWAS stage.\n";

    keep_.keep_indices.resize(n_samples_);
    for (int i = 0; i < n_samples_; ++i) keep_.keep_indices[i] = i;
    keep_.n_keep = n_samples_;
}

void RrmModel::run() {
    if (options_.rrm_compatible == "gma-fixreg") {
        runGmaFixregCompatible();
        return;
    }

    if (options_.rrm_unbalance) {
        runUnbalanced();
        return;
    }

    GwasScanner::Options opts;
    opts.out_prefix = options_.out_prefix;
    opts.min_maf = options_.min_maf;
    opts.min_hwe = options_.min_hwe;
    opts.fast_lmm = options_.fast_lmm;
    opts.use_spa = options_.use_spa;
    opts.no_loco = options_.no_loco;
    opts.grm_alg = options_.grm_alg;
    opts.rrm_curve_grid = options_.rrm_curve_grid;
    opts.rrm_top = options_.rrm_top;
    opts.rrm_ai_iter = options_.rrm_ai_iter;
    opts.rrm_sig_threshold = options_.rrm_sig_threshold;
    opts.rrm_loco_global_null = options_.rrm_loco_global_null;
    opts.rrm_basis = options_.rrm_basis;
    #ifdef _OPENMP
    opts.n_threads = options_.n_threads;
    #endif

    GwasScanner scanner(opts);
scanner.runRrmScan(records_, options_.rrm_order, X_base_, fids_, iids_,
                        *plink_reader_, *bgen_reader_, *pgen_reader_, *vcf_reader_,
                        is_bgen_, is_pgen_, is_vcf_, keep_);
}

static double rrm_gamma_q(double a, double x) {
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

static double rrm_chisq_prob(double stat, int df) {
    if (stat <= 0.0) return 1.0;
    return rrm_gamma_q(df / 2.0, stat / 2.0);
}

static double factorial_d(int n) {
    double out = 1.0;
    for (int i = 2; i <= n; ++i) out *= i;
    return out;
}

static Eigen::MatrixXd gma_legendre(int n_time, int order) {
    Eigen::MatrixXd L(n_time, order + 1);
    const double tmin = 1.0;
    const double tmax = (double)n_time;
    for (int i = 0; i < n_time; ++i) {
        double t = (n_time == 1) ? 0.0 : 2.0 * ((double)i + 1.0 - tmin) / (tmax - tmin) - 1.0;
        for (int k = 0; k <= order; ++k) {
            int c = k / 2;
            double p = 0.0;
            for (int r = 0; r <= c; ++r) {
                p += std::sqrt((2.0 * k + 1.0) / 2.0) * std::pow(0.5, k) *
                     (std::pow(-1.0, r) * factorial_d(2 * k - 2 * r) /
                      (factorial_d(r) * factorial_d(k - r) * factorial_d(k - 2 * r))) *
                     std::pow(t, k - 2 * r);
            }
            L(i, k) = p;
        }
    }
    return L;
}

static std::vector<double> read_gma_varcom(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open GMA varcom file: " + path);
    std::string line;
    std::getline(ifs, line);
    std::vector<double> vals;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string a, b, c;
        double v;
        if (ss >> a >> b >> c >> v) vals.push_back(v);
    }
    if (vals.empty()) throw std::runtime_error("No variance components read from: " + path);
    return vals;
}

static Eigen::MatrixXd lower_tri_cov(const std::vector<double>& vals, int offset, int dim) {
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(dim, dim);
    int idx = offset;
    for (int j = 0; j < dim; ++j) {
        for (int k = 0; k <= j; ++k) {
            C(j, k) = vals[idx++];
            C(k, j) = C(j, k);
        }
    }
    return C;
}

static bool gma_covariances(const Eigen::VectorXd& varcom, int d,
                            Eigen::MatrixXd& cov_add,
                            Eigen::MatrixXd& cov_per,
                            double& ve) {
    int tri = d * (d + 1) / 2;
    cov_add = Eigen::MatrixXd::Zero(d, d);
    cov_per = Eigen::MatrixXd::Zero(d, d);
    int idx = 0;
    for (int j = 0; j < d; ++j) {
        for (int k = 0; k <= j; ++k) {
            cov_add(j, k) = varcom(idx++);
            cov_add(k, j) = cov_add(j, k);
        }
    }
    for (int j = 0; j < d; ++j) {
        for (int k = 0; k <= j; ++k) {
            cov_per(j, k) = varcom(idx++);
            cov_per(k, j) = cov_per(j, k);
        }
    }
    ve = varcom(2 * tri);
    if (ve < 0.0) return false;
    return Eigen::LLT<Eigen::MatrixXd>(cov_add).info() == Eigen::Success &&
           Eigen::LLT<Eigen::MatrixXd>(cov_per).info() == Eigen::Success;
}

static Eigen::MatrixXd gma_tit_matrix(const Eigen::MatrixXd& L, int d, int var_idx) {
    const int tri = d * (d + 1) / 2;
    if (var_idx == 2 * tri) return Eigen::MatrixXd::Identity(L.rows(), L.rows());
    int local = var_idx % tri;
    int idx = 0;
    int row = 0, col = 0;
    for (int j = 0; j < d; ++j) {
        for (int k = 0; k <= j; ++k) {
            if (idx == local) { row = j; col = k; }
            ++idx;
        }
    }
    Eigen::MatrixXd one = Eigen::MatrixXd::Zero(d, d);
    one(row, col) = 1.0;
    one(col, row) = 1.0;
    return L * one * L.transpose();
}

static Eigen::MatrixXd gma_em_matrix_eigen(const Eigen::VectorXd& varcom, int d, int ran_df, int num_record) {
    int tri = d * (d + 1) / 2;
    Eigen::MatrixXd cov_add, cov_per;
    double ve = 0.0;
    gma_covariances(varcom, d, cov_add, cov_per, ve);
    std::vector<std::pair<int, int>> ind;
    for (int j = 0; j < d; ++j) {
        for (int k = 0; k <= j; ++k) ind.emplace_back(j, k);
    }
    Eigen::MatrixXd em = Eigen::MatrixXd::Zero(2 * tri + 1, 2 * tri + 1);
    auto fill_block = [&](const Eigen::MatrixXd& C, int offset) {
        Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(tri, tri);
        for (int j = 0; j < tri; ++j) {
            for (int k = 0; k <= j; ++k) {
                int a0 = ind[j].first, a1 = ind[j].second;
                int b0 = ind[k].first, b1 = ind[k].second;
                tmp(k, j) = tmp(j, k) = (C(a0, b0) * C(a1, b1) + C(a0, b1) * C(a1, b0)) / (2.0 * ran_df);
            }
        }
        em.block(offset, offset, tri, tri) = tmp.inverse();
    };
    fill_block(cov_add, 0);
    fill_block(cov_per, tri);
    em(2 * tri, 2 * tri) = num_record / (ve * ve);
    return 0.5 * em;
}

static bool gma_fd_ai_update(const std::vector<Eigen::MatrixXd>& X,
                             const Eigen::MatrixXd& Yrot,
                             const Eigen::MatrixXd& L,
                             const Eigen::VectorXd& eval,
                             const Eigen::VectorXd& varcom,
                             Eigen::VectorXd& delta) {
    const int n = (int)X.size();
    const int t = (int)L.rows();
    const int d = (int)L.cols();
    const int p = (int)X[0].cols();
    const int nvar = (int)varcom.size();
    const int tri = d * (d + 1) / 2;
    Eigen::MatrixXd cov_add, cov_per;
    double ve = 0.0;
    if (!gma_covariances(varcom, d, cov_add, cov_per, ve)) return false;
    Eigen::MatrixXd LCaddLt = L * cov_add * L.transpose();
    Eigen::MatrixXd LCperLt = L * cov_per * L.transpose();
    std::vector<Eigen::MatrixXd> Vinv(n), xv(n);
    std::vector<Eigen::RowVectorXd> yv(n);
    Eigen::MatrixXd xvx_sum = Eigen::MatrixXd::Zero(p, p);
    Eigen::RowVectorXd yvx_sum = Eigen::RowVectorXd::Zero(p);
    for (int i = 0; i < n; ++i) {
        Eigen::MatrixXd V = eval(i) * LCaddLt + LCperLt + ve * Eigen::MatrixXd::Identity(t, t);
        Vinv[i] = V.inverse();
        xv[i] = X[i].transpose() * Vinv[i];
        yv[i] = Yrot.row(i) * Vinv[i];
        xvx_sum.noalias() += xv[i] * X[i];
        yvx_sum.noalias() += yv[i] * X[i];
    }
    Eigen::MatrixXd xvx_inv = xvx_sum.inverse();
    std::vector<Eigen::MatrixXd> tit(nvar), xvkvx(nvar);
    std::vector<Eigen::RowVectorXd> yvkvx(nvar);
    std::vector<double> yvkvy(nvar, 0.0);
    Eigen::VectorXd fd = Eigen::VectorXd::Zero(nvar);
    Eigen::MatrixXd ai = Eigen::MatrixXd::Zero(nvar, nvar);
    for (int m = 0; m < nvar; ++m) {
        tit[m] = gma_tit_matrix(L, d, m);
        double trvk = 0.0;
        xvkvx[m] = Eigen::MatrixXd::Zero(p, p);
        yvkvx[m] = Eigen::RowVectorXd::Zero(p);
        for (int i = 0; i < n; ++i) {
            double km = (m < tri) ? eval(i) : 1.0;
            Eigen::MatrixXd vt = Vinv[i] * tit[m];
            trvk += km * vt.trace();
            xvkvx[m].noalias() += km * xv[i] * tit[m] * xv[i].transpose();
            double yky = (yv[i] * tit[m] * yv[i].transpose())(0, 0);
            yvkvy[m] += km * yky;
            yvkvx[m].noalias() += km * yv[i] * tit[m] * xv[i].transpose();
        }
        double trpk = (xvx_inv.array() * xvkvx[m].array()).sum();
        double ypkpy2 = (yvx_sum * xvx_inv * yvkvx[m].transpose())(0, 0);
        double ypkpy3 = (yvx_sum * xvx_inv * xvkvx[m] * xvx_inv * yvx_sum.transpose())(0, 0);
        double ypkpy = yvkvy[m] - 2.0 * ypkpy2 + ypkpy3;
        fd(m) = -0.5 * (trvk - trpk) + 0.5 * ypkpy;
        for (int q = 0; q <= m; ++q) {
            double val1 = 0.0;
            Eigen::RowVectorXd yvkvkvx = Eigen::RowVectorXd::Zero(p);
            Eigen::RowVectorXd yvkvkvx2 = Eigen::RowVectorXd::Zero(p);
            Eigen::MatrixXd xvkvkvx2 = Eigen::MatrixXd::Zero(p, p);
            for (int i = 0; i < n; ++i) {
                double km = (m < tri) ? eval(i) : 1.0;
                double kq = (q < tri) ? eval(i) : 1.0;
                double kk = km * kq;
                val1 += kk * (yv[i] * tit[m] * Vinv[i] * tit[q] * yv[i].transpose())(0, 0);
                yvkvkvx.noalias() += kk * yv[i] * tit[m] * Vinv[i] * tit[q] * xv[i].transpose();
                yvkvkvx2.noalias() += kk * yv[i] * tit[q] * Vinv[i] * tit[m] * xv[i].transpose();
                xvkvkvx2.noalias() += kk * xv[i] * tit[m] * Vinv[i] * tit[q] * xv[i].transpose();
            }
            double val2 = (yvkvkvx * xvx_inv * yvx_sum.transpose())(0, 0);
            double val3 = (yvkvx[m] * xvx_inv * yvkvx[q].transpose())(0, 0);
            double val4 = (yvkvx[m] * xvx_inv * xvkvx[q] * xvx_inv * yvx_sum.transpose())(0, 0);
            double val5 = (yvx_sum * xvx_inv * yvkvkvx2.transpose())(0, 0);
            double val6 = (yvx_sum * xvx_inv * xvkvkvx2 * xvx_inv * yvx_sum.transpose())(0, 0);
            double val7 = (yvx_sum * xvx_inv * xvkvx[m] * xvx_inv * yvkvx[q].transpose())(0, 0);
            double val8 = (yvx_sum * xvx_inv * xvkvx[m] * xvx_inv * xvkvx[q] * xvx_inv * yvx_sum.transpose())(0, 0);
            ai(m, q) = 0.5 * (val1 - val2 - val3 + val4 - val5 + val6 + val7 - val8);
            ai(q, m) = ai(m, q);
        }
    }
    delta = ai.fullPivLu().solve(fd);
    Eigen::MatrixXd ca, cp;
    double cv = 0.0;
    if (delta.allFinite() && gma_covariances(varcom + delta, d, ca, cp, cv)) return true;

    Eigen::MatrixXd em = gma_em_matrix_eigen(varcom, d, n, n * t);
    for (int step = 1; step <= 100; ++step) {
        double gamma = step / 100.0;
        Eigen::MatrixXd blended = (1.0 - gamma) * ai + gamma * em;
        delta = blended.fullPivLu().solve(fd);
        if (delta.allFinite() && gma_covariances(varcom + delta, d, ca, cp, cv)) return true;
    }
    return false;
}

static Eigen::MatrixXd read_gma_grm(const std::string& path, const std::vector<std::string>& ids) {
    std::unordered_map<std::string, int> id_to_idx;
    for (int i = 0; i < (int)ids.size(); ++i) id_to_idx[ids[i]] = i;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(ids.size(), ids.size());
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open GMA GRM file: " + path);
    std::string a, b;
    double v;
    while (ifs >> a >> b >> v) {
        auto ia = id_to_idx.find(a);
        auto ib = id_to_idx.find(b);
        if (ia == id_to_idx.end() || ib == id_to_idx.end()) continue;
        K(ia->second, ib->second) = v;
    }
    K = K + K.transpose();
    K.diagonal() *= 0.5;
    return K;
}

void RrmModel::runGmaFixregCompatible() {
    const int n = keep_.n_keep;
    const int n_time = (int)options_.mpheno_list.size();
    const int order = options_.rrm_order;
    const int d = order + 1;
    if (n_time <= 0) throw std::runtime_error("GMA-compatible fixreg requires at least one time point.");
    std::vector<double> varcom = read_gma_varcom(options_.rrm_compatible_varcom);
    int tri = d * (d + 1) / 2;
    if ((int)varcom.size() != 2 * tri + 1) {
        throw std::runtime_error("GMA varcom dimension does not match --rrm-order.");
    }
    Eigen::MatrixXd cov_add = lower_tri_cov(varcom, 0, d);
    Eigen::MatrixXd cov_per = lower_tri_cov(varcom, tri, d);
    double ve = varcom.back();

    std::vector<std::string> ids;
    Eigen::MatrixXd Y(n, n_time);
    for (int r = 0; r < n; ++r) {
        int idx = keep_.keep_indices[r];
        ids.push_back(iids_[idx]);
        Y.row(r) = Y_wide_.row(idx);
    }

    std::string grm_path = options_.bfile_prefix + ".grm";
    Eigen::MatrixXd K = read_gma_grm(grm_path, ids);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K);
    if (es.info() != Eigen::Success) throw std::runtime_error("Eigen decomposition failed for GMA GRM.");
    Eigen::VectorXd eval = es.eigenvalues();
    Eigen::MatrixXd U = es.eigenvectors();
    Eigen::MatrixXd Yrot = U.transpose() * Y;
    Eigen::MatrixXd L = gma_legendre(n_time, order);
    Eigen::VectorXd ones_rot = U.transpose() * Eigen::VectorXd::Ones(n);
    const auto& bim = plink_reader_->getBimInfo();
    std::string out_file = options_.out_prefix + ".rrm.gma_fixreg.assoc.txt";
    std::ofstream ofs(out_file);
    if (!ofs) throw std::runtime_error("Cannot write output: " + out_file);
    ofs << "CHR\tSNP\tBP\tA1\tA2\tCHISQ\tP\tchi_val\tp_val";
    for (int k = 0; k < d; ++k) ofs << "\tBETA_L" << k;
    ofs << "\n";
    ofs << std::setprecision(12);

    Eigen::VectorXd dosage_all(plink_reader_->getNumSamples());
    for (size_t s = 0; s < plink_reader_->getNumSnps(); ++s) {
        plink_reader_->readSnp((int)s, dosage_all);
        Eigen::VectorXd g(n);
        double sum = 0.0;
        int count = 0;
        for (int r = 0; r < n; ++r) {
            double val = dosage_all[keep_.keep_indices[r]];
            if (std::isfinite(val)) { sum += val; count++; }
        }
        double mean = count > 0 ? sum / count : 0.0;
        for (int r = 0; r < n; ++r) {
            double val = dosage_all[keep_.keep_indices[r]];
            g(r) = std::isfinite(val) ? val : mean;
        }
        Eigen::VectorXd grot = U.transpose() * g;

        std::vector<Eigen::MatrixXd> X(n);
        for (int i = 0; i < n; ++i) {
            X[i].resize(n_time, 2 * d);
            X[i].leftCols(d) = ones_rot(i) * L;
            X[i].rightCols(d) = grot(i) * L;
        }

        Eigen::VectorXd var_update(varcom.size());
        for (int j = 0; j < (int)varcom.size(); ++j) var_update(j) = varcom[j];
        for (int iter = 0; iter < 10; ++iter) {
            Eigen::VectorXd delta;
            if (!gma_fd_ai_update(X, Yrot, L, eval, var_update, delta)) break;
            Eigen::VectorXd candidate = var_update + delta;
            Eigen::MatrixXd ca, cp;
            double cv = 0.0;
            if (!gma_covariances(candidate, d, ca, cp, cv)) break;
            var_update = candidate;
        }
        Eigen::MatrixXd cov_add_snp, cov_per_snp;
        double ve_snp = 0.0;
        if (!gma_covariances(var_update, d, cov_add_snp, cov_per_snp, ve_snp)) {
            cov_add_snp = cov_add;
            cov_per_snp = cov_per;
            ve_snp = ve;
        }
        Eigen::MatrixXd LCaddLt = L * cov_add_snp * L.transpose();
        Eigen::MatrixXd LCperLt = L * cov_per_snp * L.transpose();
        Eigen::MatrixXd xvx = Eigen::MatrixXd::Zero(2 * d, 2 * d);
        Eigen::VectorXd xvy = Eigen::VectorXd::Zero(2 * d);
        for (int i = 0; i < n; ++i) {
            Eigen::MatrixXd V = eval(i) * LCaddLt + LCperLt + ve_snp * Eigen::MatrixXd::Identity(n_time, n_time);
            Eigen::MatrixXd Vinv = V.inverse();
            xvx.noalias() += X[i].transpose() * Vinv * X[i];
            xvy.noalias() += X[i].transpose() * Vinv * Yrot.row(i).transpose();
        }
        Eigen::MatrixXd cov_beta = xvx.inverse();
        Eigen::VectorXd beta = cov_beta * xvy;
        Eigen::VectorXd snp_beta = beta.tail(d);
        Eigen::MatrixXd snp_cov = cov_beta.bottomRightCorner(d, d);
        double chisq = (snp_beta.transpose() * snp_cov.inverse() * snp_beta)(0, 0);
        double p = rrm_chisq_prob(chisq, d);
        const auto& snp = bim[s];
        ofs << snp.chrom << '\t' << snp.id << '\t' << snp.bp_pos << '\t' << snp.alt << '\t' << snp.ref
            << '\t' << chisq << '\t' << p << '\t' << chisq << '\t' << p;
        for (int k = 0; k < d; ++k) ofs << '\t' << snp_beta(k);
        ofs << '\n';
    }
    std::cout << "Wrote GMA-compatible fixreg association results: " << out_file << "\n";
}

void RrmModel::runUnbalanced() {
    if (records_.empty()) throw std::runtime_error("No unbalanced longitudinal records loaded.");
    int n_ind = keep_.n_keep;
    int n_obs = (int)records_.size();
    int n_basis = options_.rrm_order + 1;

    Eigen::MatrixXf X_tdep;
    std::vector<std::string> tdep_names;

    if (!options_.tclass_names.empty() || !options_.tcovar_names.empty()) {
        std::cout << "Parsing time-dependent covariates from phenotype file...\n";
        std::ifstream ifs(options_.pheno_file);
        if (!ifs) throw std::runtime_error("Cannot open phenotype file: " + options_.pheno_file);
        std::string header;
        std::getline(ifs, header);
        std::stringstream hss(header);
        std::vector<std::string> col_names;
        std::string token;
        while (hss >> token) col_names.push_back(token);
        if (col_names.size() < 4) throw std::runtime_error("Phenotype file must have at least 4 columns.");

        auto parse_col_list = [&](const std::string& spec, std::vector<int>& out_idx) {
            if (spec.empty()) return;
            std::stringstream ss(spec);
            std::string name;
            while (std::getline(ss, name, ',')) {
                auto it = std::find(col_names.begin(), col_names.end(), name);
                if (it == col_names.end())
                    throw std::runtime_error("Column '" + name + "' not found in phenotype file.");
                out_idx.push_back((int)(it - col_names.begin()));
            }
        };

        std::vector<int> tclass_cols, tcovar_cols;
        parse_col_list(options_.tclass_names, tclass_cols);
        parse_col_list(options_.tcovar_names, tcovar_cols);

        for (int c : tclass_cols) tdep_names.push_back(col_names[c] + "_class");
        for (int c : tcovar_cols) tdep_names.push_back(col_names[c]);

        int n_tdep = (int)tclass_cols.size() + (int)tcovar_cols.size();
        if (n_tdep > 0) {
            X_tdep.resize(n_obs, n_tdep);
            std::map<std::string, std::map<std::string, int>> class_levels;

            for (int ri = 0; ri < n_obs; ++ri) {
                std::string line;
                if (!std::getline(ifs, line)) break;
                std::stringstream lss(line);
                std::vector<std::string> tokens;
                std::string val;
                while (lss >> val) tokens.push_back(val);
                if (tokens.size() < col_names.size()) continue;
                for (int c : tclass_cols) {
                    if (c < (int)tokens.size()) {
                        auto& levels = class_levels[col_names[c]];
                        if (levels.find(tokens[c]) == levels.end()) {
                            levels[tokens[c]] = (int)levels.size() + 1;
                        }
                    }
                }
            }

            ifs.clear();
            ifs.seekg(0);
            std::getline(ifs, header);

            for (int ri = 0; ri < n_obs; ++ri) {
                std::string line;
                if (!std::getline(ifs, line)) break;
                std::stringstream lss(line);
                std::vector<std::string> tokens;
                std::string val;
                while (lss >> val) tokens.push_back(val);
                int col = 0;
                for (int c : tclass_cols) {
                    X_tdep(ri, col) = (double)(class_levels[col_names[c]][tokens[c]] - 1);
                    col++;
                }
                for (int c : tcovar_cols) {
                    double v = (c < (int)tokens.size()) ? std::stod(tokens[c]) : 0.0;
                    X_tdep(ri, col) = (float)v;
                    col++;
                }
            }
        }
    }

    GwasScanner::Options opts;
    opts.out_prefix = options_.out_prefix;
    opts.min_maf = options_.min_maf;
    opts.min_hwe = options_.min_hwe;
    opts.fast_lmm = options_.fast_lmm;
    opts.use_spa = options_.use_spa;
    opts.no_loco = options_.no_loco;
    opts.grm_alg = options_.grm_alg;
    opts.rrm_curve_grid = options_.rrm_curve_grid;
    opts.rrm_top = options_.rrm_top;
    opts.rrm_ai_iter = options_.rrm_ai_iter;
    opts.rrm_sig_threshold = options_.rrm_sig_threshold;
    opts.rrm_loco_global_null = options_.rrm_loco_global_null;
    opts.rrm_basis = options_.rrm_basis;
    opts.rrm_gwas_method = options_.rrm_gwas_method;
#ifdef _OPENMP
    opts.n_threads = options_.n_threads;
#endif

    GwasScanner scanner(opts);
    scanner.runRrmUnbalancedScan(records_, options_.rrm_order, X_base_, X_tdep, tdep_names,
                                 fids_, iids_, *plink_reader_, *bgen_reader_, *pgen_reader_,
                                 is_bgen_, is_pgen_, keep_);
}

} // namespace cosmic
