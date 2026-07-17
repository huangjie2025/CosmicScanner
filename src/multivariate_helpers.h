#pragma once

#include <fstream>
#include <ios>
#include <string>
#include <Eigen/Dense>
#include "multivariate_model.h"

namespace cosmic {

struct MvCrtFactors {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    bool ok = false;
};

std::ofstream mv_open_output_file(const std::string& path, std::ios_base::openmode mode = std::ios::out);
void mv_validate_output_parent(const std::string& out_prefix);
void mv_write_vbeta_header(std::ostream& os, int n_traits);
void mv_write_vbeta_values(std::ostream& os, const Eigen::MatrixXd& vbeta, int n_traits);
double mv_chisq_prob(double stat, int df);
bool mv_fit_gls(const Eigen::MatrixXd& UtY,
                const Eigen::MatrixXd& UtX,
                const Eigen::VectorXd& D,
                const Eigen::MatrixXd& Vg,
                const Eigen::MatrixXd& Ve,
                Eigen::VectorXd& beta_all,
                Eigen::MatrixXd& beta_cov,
                double* loglik = nullptr);
bool mv_em_refit_rotated(const Eigen::MatrixXd& UtY,
                         const Eigen::MatrixXd& UtX,
                         const Eigen::VectorXd& D,
                         int max_iter,
                         double tol,
                         Eigen::MatrixXd& Vg,
                         Eigen::MatrixXd& Ve,
                         const Eigen::MatrixXd* init_vg = nullptr,
                         const Eigen::MatrixXd* init_ve = nullptr,
                         int* iterations_out = nullptr,
                         double* update_norm_out = nullptr,
                         std::string* stop_reason_out = nullptr);
double mv_chisq_stat_from_p(double p, int df);
double mv_gc_correct_p(double p, double lambda, int df);
MvCrtFactors mv_compute_crt_factors(const Eigen::MatrixXd& UtX,
                                    const Eigen::VectorXd& D,
                                    const Eigen::MatrixXd& Vg,
                                    const Eigen::MatrixXd& Ve);
double mv_pcrt(int mode, int df, double p, const MvCrtFactors& crt);
bool mv_scan_filter_passes(const PlinkBimInfo& var, const MultivariateModel::Options& options);

} // namespace cosmic
