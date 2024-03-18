#ifndef STMESH_UTILITY_HPP
#define STMESH_UTILITY_HPP

#include <array>
#include <concepts> // IWYU pragma: keep
#include <limits>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <internal_use_only/config.hpp>

namespace boost::geometry {}

namespace stmesh {
using FLOAT_T = stmesh::cmake::FLOATING_POINT_TYPE;
template <unsigned D> using VectorF = Eigen::Vector<FLOAT_T, static_cast<int>(D)>;
using Vector4F = VectorF<static_cast<unsigned>(4)>;
using Vector3F = VectorF<static_cast<unsigned>(3)>;

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace bg = boost::geometry;

template <typename Scalar, int Rows, int Cols> struct MatrixHash {
  size_t operator()(const Eigen::Matrix<Scalar, Rows, Cols> &matrix) const;
};

using Vector4FHash = MatrixHash<FLOAT_T, 4, 1>;
using Vector3FHash = MatrixHash<FLOAT_T, 3, 1>;

namespace detail {
template <std::floating_point T> T consteval sqrtNewtonRaphson(T x, T curr, T prev) {
  return curr == prev ? curr : sqrtNewtonRaphson(x, (curr + x / curr) / 2, curr);
}
} // namespace detail

template <std::floating_point T> T consteval sqrt(T x) {
  return x >= 0 && x < std::numeric_limits<T>::infinity() ? detail::sqrtNewtonRaphson(x, x, T{0})
                                                          : std::numeric_limits<T>::quiet_NaN();
}

[[nodiscard]] constexpr unsigned nChoosek(unsigned n, unsigned k) noexcept {
  if (k > n)
    return 0;
  if (k * 2 > n)
    k = n - k;
  if (k == 0)
    return 1;

  unsigned result = n;
  for (unsigned i = 2; i <= k; ++i) {
    result *= (n - i + 1);
    result /= i;
  }
  return result;
}

constexpr unsigned factorial(unsigned n) noexcept {
  unsigned result = 1;
  for (unsigned i = 2; i <= n; ++i)
    result *= i;
  return result;
}

template <int D>
[[nodiscard]] std::array<VectorF<static_cast<unsigned>(D)>, 1U << static_cast<unsigned>(D)>
allCorners(const Eigen::AlignedBox<FLOAT_T, D> &box);

template <int D, int N>
[[nodiscard]] Eigen::Matrix<FLOAT_T, D, D - N> kernel(const Eigen::Matrix<FLOAT_T, D, N> &matrix) noexcept;
} // namespace stmesh
#endif