#pragma once
#include <vector>
#include <cstddef>
#include <cstdint>

namespace cosmic {
namespace apy {

struct Matrix {
  std::size_t rows = 0;
  std::size_t cols = 0;
  std::vector<double> a;

  Matrix() = default;
  Matrix(std::size_t r, std::size_t c, double v = 0.0) : rows(r), cols(c), a(r * c, v) {}

  double& operator()(std::size_t r, std::size_t c) { return a[r * cols + c]; }
  double operator()(std::size_t r, std::size_t c) const { return a[r * cols + c]; }
};

struct CholeskyLower {
  Matrix L;
  CholeskyLower() = default;
  explicit CholeskyLower(Matrix A) : L(std::move(A)) {}
  void factorize_in_place();
  std::vector<double> solve(const std::vector<double>& b) const;
};

struct APYModel {
  std::size_t n = 0;
  std::size_t ncore = 0;
  std::size_t nnon = 0;
  double scale = 1.0;
  double ridge = 0.0;
  std::vector<std::size_t> apy_to_orig;
  std::vector<std::size_t> orig_to_apy;
  std::vector<std::size_t> core_orig;
  std::vector<std::size_t> non_orig;
  Matrix Gcc;
  CholeskyLower chol;
  Matrix X;
  std::vector<double> d_inv;

  static APYModel build_from_Z(const Matrix& Z, const std::vector<int>& core_status, double scale, double ridge);
  std::vector<double> multiply_Ginv_APY(const std::vector<double>& v_orig_order) const;
};

std::vector<int> make_core_status_random(std::size_t n, std::size_t ncore, std::uint64_t seed);
std::vector<int> make_core_status_pca_like(const Matrix& Z, std::size_t ncore, std::uint64_t seed, std::size_t proj_dim);

} // namespace apy
} // namespace cosmic
