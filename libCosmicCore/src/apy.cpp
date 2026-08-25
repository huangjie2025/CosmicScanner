#include "apy.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>

namespace cosmic {
namespace apy {

static void require(bool ok, const std::string& msg) {
  if (!ok) throw std::runtime_error(msg);
}

static double dot_rows_scaled(const Matrix& Z, std::size_t i, std::size_t j, double scale) {
  require(Z.cols != 0, "dot_rows_scaled: empty Z");
  const std::size_t base_i = i * Z.cols;
  const std::size_t base_j = j * Z.cols;
  double s = 0.0;
  for (std::size_t k = 0; k < Z.cols; ++k) s += Z.a[base_i + k] * Z.a[base_j + k];
  return s * scale;
}

static double diag_row_scaled(const Matrix& Z, std::size_t i, double scale) {
  const std::size_t base = i * Z.cols;
  double s = 0.0;
  for (std::size_t k = 0; k < Z.cols; ++k) {
    const double x = Z.a[base + k];
    s += x * x;
  }
  return s * scale;
}

void CholeskyLower::factorize_in_place() {
  require(L.rows == L.cols, "cholesky: matrix must be square");
  const std::size_t n = L.rows;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      double s = L(i, j);
      for (std::size_t k = 0; k < j; ++k) s -= L(i, k) * L(j, k);
      if (i == j) {
        require(s > 0.0, "cholesky: matrix not SPD");
        L(i, j) = std::sqrt(s);
      } else {
        L(i, j) = s / L(j, j);
      }
    }
    for (std::size_t j = i + 1; j < n; ++j) L(i, j) = 0.0;
  }
}

std::vector<double> CholeskyLower::solve(const std::vector<double>& b) const {
  require(L.rows == L.cols, "cholesky solve: matrix must be square");
  require(b.size() == L.rows, "cholesky solve: dimension mismatch");
  const std::size_t n = L.rows;
  std::vector<double> y(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    double s = b[i];
    for (std::size_t k = 0; k < i; ++k) s -= L(i, k) * y[k];
    y[i] = s / L(i, i);
  }
  std::vector<double> x(n, 0.0);
  for (std::size_t ii = 0; ii < n; ++ii) {
    const std::size_t i = n - 1 - ii;
    double s = y[i];
    for (std::size_t k = i + 1; k < n; ++k) s -= L(k, i) * x[k];
    x[i] = s / L(i, i);
  }
  return x;
}

std::vector<int> make_core_status_random(std::size_t n, std::size_t ncore, std::uint64_t seed) {
  require(ncore > 0 && ncore <= n, "random core: invalid ncore");
  std::vector<std::size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::mt19937_64 rng(seed);
  std::shuffle(idx.begin(), idx.end(), rng);
  std::vector<int> s(n, 0);
  for (std::size_t i = 0; i < ncore; ++i) s[idx[i]] = 1;
  return s;
}

std::vector<int> make_core_status_pca_like(const Matrix& Z, std::size_t ncore, std::uint64_t seed, std::size_t proj_dim) {
  require(ncore > 0 && ncore <= Z.rows, "pca core: invalid ncore");
  require(proj_dim > 0, "pca core: invalid proj_dim");

  std::mt19937_64 rng(seed);
  std::normal_distribution<double> nd(0.0, 1.0);
  Matrix Omega(Z.cols, proj_dim, 0.0);
  for (double& x : Omega.a) x = nd(rng);

  std::vector<double> score(Z.rows, 0.0);
  for (std::size_t i = 0; i < Z.rows; ++i) {
    const std::size_t base = i * Z.cols;
    double row_score = 0.0;
    for (std::size_t k = 0; k < proj_dim; ++k) {
      double s = 0.0;
      for (std::size_t j = 0; j < Z.cols; ++j) s += Z.a[base + j] * Omega(j, k);
      row_score += s * s;
    }
    score[i] = row_score;
  }

  std::vector<std::size_t> idx(Z.rows);
  std::iota(idx.begin(), idx.end(), 0);
  std::partial_sort(idx.begin(), idx.begin() + static_cast<std::ptrdiff_t>(ncore), idx.end(),
                    [&](std::size_t a, std::size_t b) { return score[a] > score[b]; });
  std::vector<int> s(Z.rows, 0);
  for (std::size_t i = 0; i < ncore; ++i) s[idx[i]] = 1;
  return s;
}

APYModel APYModel::build_from_Z(const Matrix& Z, const std::vector<int>& core_status, double scale, double ridge) {
  require(Z.rows > 0 && Z.cols > 0, "build_from_Z: empty Z");
  require(core_status.size() == Z.rows, "build_from_Z: core_status length mismatch");
  require(std::isfinite(scale) && scale > 0.0, "build_from_Z: scale must be > 0");
  require(std::isfinite(ridge) && ridge >= 0.0, "build_from_Z: ridge must be >= 0");

  APYModel m;
  m.n = Z.rows;
  m.scale = scale;
  m.ridge = ridge;

  m.core_orig.reserve(m.n);
  m.non_orig.reserve(m.n);
  for (std::size_t i = 0; i < m.n; ++i) {
    if (core_status[i] == 1) m.core_orig.push_back(i);
    else m.non_orig.push_back(i);
  }
  require(!m.core_orig.empty(), "build_from_Z: core set is empty");
  m.ncore = m.core_orig.size();
  m.nnon = m.non_orig.size();

  m.apy_to_orig.reserve(m.n);
  m.apy_to_orig.insert(m.apy_to_orig.end(), m.core_orig.begin(), m.core_orig.end());
  m.apy_to_orig.insert(m.apy_to_orig.end(), m.non_orig.begin(), m.non_orig.end());
  m.orig_to_apy.assign(m.n, 0);
  for (std::size_t i = 0; i < m.n; ++i) m.orig_to_apy[m.apy_to_orig[i]] = i;

  m.Gcc = Matrix(m.ncore, m.ncore, 0.0);
  for (std::size_t i = 0; i < m.ncore; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      const double v = dot_rows_scaled(Z, m.core_orig[i], m.core_orig[j], m.scale);
      m.Gcc(i, j) = v;
      m.Gcc(j, i) = v;
    }
  }
  for (std::size_t i = 0; i < m.ncore; ++i) m.Gcc(i, i) += m.ridge;

  m.chol = CholeskyLower(m.Gcc);
  m.chol.factorize_in_place();

  m.X = Matrix(m.ncore, m.nnon, 0.0);
  m.d_inv.assign(m.nnon, 0.0);

  std::vector<double> gcn(m.ncore, 0.0);
  for (std::size_t j = 0; j < m.nnon; ++j) {
    const std::size_t non_i = m.non_orig[j];
    for (std::size_t i = 0; i < m.ncore; ++i) gcn[i] = dot_rows_scaled(Z, m.core_orig[i], non_i, m.scale);

    auto x = m.chol.solve(gcn);
    for (std::size_t i = 0; i < m.ncore; ++i) m.X(i, j) = x[i];

    double s = 0.0;
    for (std::size_t i = 0; i < m.ncore; ++i) s += gcn[i] * x[i];
    double gnn = diag_row_scaled(Z, non_i, m.scale) + m.ridge;
    double d = gnn - s;
    require(std::isfinite(d) && std::abs(d) > 0.0, "build_from_Z: APY residual variance is zero");
    m.d_inv[j] = 1.0 / d;
  }

  return m;
}

std::vector<double> APYModel::multiply_Ginv_APY(const std::vector<double>& v_orig_order) const {
  require(v_orig_order.size() == n, "multiply: v length mismatch");
  std::vector<double> v_apy(n, 0.0);
  for (std::size_t apy_i = 0; apy_i < n; ++apy_i) v_apy[apy_i] = v_orig_order[apy_to_orig[apy_i]];

  std::vector<double> v_c(v_apy.begin(), v_apy.begin() + static_cast<std::ptrdiff_t>(ncore));
  std::vector<double> v_n(v_apy.begin() + static_cast<std::ptrdiff_t>(ncore), v_apy.end());

  std::vector<double> Xt_vc(nnon, 0.0);
  for (std::size_t j = 0; j < nnon; ++j) {
    double s = 0.0;
    for (std::size_t i = 0; i < ncore; ++i) s += X(i, j) * v_c[i];
    Xt_vc[j] = s;
  }

  std::vector<double> b(nnon, 0.0);
  std::vector<double> c(nnon, 0.0);
  for (std::size_t j = 0; j < nnon; ++j) {
    b[j] = d_inv[j] * Xt_vc[j];
    c[j] = d_inv[j] * v_n[j];
  }

  auto y_c = chol.solve(v_c);
  for (std::size_t i = 0; i < ncore; ++i) {
    double s = 0.0;
    for (std::size_t j = 0; j < nnon; ++j) s += X(i, j) * (b[j] - c[j]);
    y_c[i] += s;
  }

  std::vector<double> y_n(nnon, 0.0);
  for (std::size_t j = 0; j < nnon; ++j) y_n[j] = c[j] - b[j];

  std::vector<double> y_apy(n, 0.0);
  for (std::size_t i = 0; i < ncore; ++i) y_apy[i] = y_c[i];
  for (std::size_t j = 0; j < nnon; ++j) y_apy[ncore + j] = y_n[j];

  std::vector<double> y_orig(n, 0.0);
  for (std::size_t apy_i = 0; apy_i < n; ++apy_i) y_orig[apy_to_orig[apy_i]] = y_apy[apy_i];
  return y_orig;
}

} // namespace apy
} // namespace cosmic
