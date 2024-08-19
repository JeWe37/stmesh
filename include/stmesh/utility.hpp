#ifndef STMESH_UTILITY_HPP
#define STMESH_UTILITY_HPP

#include <array>
#include <concepts> // IWYU pragma: keep
#include <cstddef>
#include <limits>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <fmt/core.h>
#include <fmt/ostream.h>

#include <internal_use_only/config.hpp>

namespace boost::geometry {}

namespace stmesh {
using FLOAT_T = stmesh::cmake::FLOATING_POINT_TYPE; ///< The floating point type used in the library
template <unsigned D>
using VectorF = Eigen::Vector<FLOAT_T, static_cast<int>(D)>; ///< A vector of dimension D of FLOAT_T
using Vector4F = VectorF<static_cast<unsigned>(4)>;          ///< A vector of dimension 4 of FLOAT_T
using Vector3F = VectorF<static_cast<unsigned>(3)>;          ///< A vector of dimension 3 of FLOAT_T

template <unsigned D>
using ArrayF = Eigen::Array<FLOAT_T, static_cast<int>(D), 1>; ///< An array of dimension D of FLOAT_T
using Array4F = ArrayF<static_cast<unsigned>(4)>;             ///< An array of dimension 4 of FLOAT_T
using Array3F = ArrayF<static_cast<unsigned>(3)>;             ///< An array of dimension 3 of FLOAT_T

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace bg = boost::geometry; ///< Alias for boost::geometry

/// A hash function for Eigen::Matrix
/**
 * A hash function for Eigen::Matrix. This is used to allow Eigen::Matrix to be used as a key in a std::unordered_map.
 *
 * @tparam Scalar The scalar type of the matrix
 * @tparam Rows The number of rows in the matrix
 * @tparam Cols The number of columns in the matrix
 */
template <typename Scalar, int Rows, int Cols> struct MatrixHash {
  /// The hash function
  /**
   * The hash function for Eigen::Matrix. Calculated by hashing each element of the matrix and combining the hashes.
   *
   * @param matrix The matrix to hash
   * @return The hash of the matrix
   */
  size_t operator()(const Eigen::Matrix<Scalar, Rows, Cols> &matrix) const;
};

using Vector4FHash = MatrixHash<FLOAT_T, 4, 1>; ///< A hash function for Vector4F
using Vector3FHash = MatrixHash<FLOAT_T, 3, 1>; ///< A hash function for Vector3F

namespace detail {
template <std::floating_point T> T consteval sqrtNewtonRaphson(T x, T curr, T prev) {
  return curr == prev ? curr : sqrtNewtonRaphson(x, (curr + x / curr) / 2, curr);
}
} // namespace detail

/// A sign function
/**
 * A sign function. This function returns the sign of a number.
 *
 * @tparam T The type of the number
 * @param val The number to get the sign of
 * @return The sign of val
 */
template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

/// A consteval square root function
/**
 * A consteval square root function. This function uses the Newton-Raphson method to calculate the square root of a
 * number. If the number is negative or infinity, the function returns NaN.
 * Should the std::sqrt function become constexpr, this function can be removed.
 *
 * @param x The number to calculate the square root of
 * @return The square root of x
 */
template <std::floating_point T> T consteval sqrt(T x) {
  return x >= 0 && x < std::numeric_limits<T>::infinity() ? detail::sqrtNewtonRaphson(x, x, T{0})
                                                          : std::numeric_limits<T>::quiet_NaN();
}

/// A constexpr power function
/**
 * A constexpr power function. This function calculates the power of a number.
 *
 * @param x The exponent
 * @param N The base
 * @return $N^x$
 */
[[nodiscard]] constexpr static size_t expN(size_t x, size_t N) noexcept {
  size_t result = 1;
  for (size_t i = 0; i < x; ++i)
    result *= N;
  return result;
}

/// A constexpr Choose function
/**
 * A constexpr Choose function. This function calculates the number of ways to choose k elements from a set of n
 * elements.
 *
 * @param n The number of elements in the set
 * @param k The number of elements to choose
 * @return The number of ways to choose k elements from a set of n elements
 */
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

/// A constexpr factorial function
/**
 * A constexpr factorial function. This function calculates the factorial of a number.
 *
 * @param n The number to calculate the factorial of
 * @return The factorial of n
 */
constexpr unsigned factorial(unsigned n) noexcept {
  unsigned result = 1;
  for (unsigned i = 2; i <= n; ++i)
    result *= i;
  return result;
}

/// A function to list all corners of a box
/**
 * A function to list all corners of a box. This function returns an array of vectors, each of which is a corner of the
 * box.
 *
 * @param box The box to list the corners of
 * @return An array of vectors, each of which is a corner of the box
 */
template <int D>
[[nodiscard]] std::array<VectorF<static_cast<unsigned>(D)>, 1U << static_cast<unsigned>(D)>
allCorners(const Eigen::AlignedBox<FLOAT_T, D> &box);

/// A function to calculate the kernel of a matrix
/**
 * A function to calculate the kernel of a matrix. This function returns a matrix whose columns are an orthonormal basis
 * for the kernel of the input matrix. Internally this function uses Eigen's SVD decomposition to calculate the kernel.
 *
 * @param matrix The matrix to calculate the kernel of
 * @return A matrix whose columns are an orthonormal basis for the kernel of the input matrix
 */
template <int D, int N>
[[nodiscard]] Eigen::Matrix<FLOAT_T, D, D - N> kernel(const Eigen::Matrix<FLOAT_T, D, N> &matrix) noexcept;
} // namespace stmesh

/// A formatter for Eigen::Matrix
/**
 * A formatter for Eigen::Matrix. This formatter allows Eigen::Matrix to be used with fmt::format.
 * The formatter uses the ostream formatter to format the Eigen::Matrix.
 *
 * @tparam T The type of the Eigen::Matrix
 */
template <typename T>
requires std::is_base_of_v<Eigen::DenseBase<T>, T>
struct fmt::formatter<T> : ostream_formatter {};
#endif