// NOLINTBEGIN(misc-include-cleaner)
#include <stmesh/utility.hpp>

#include <functional>

namespace stmesh {
template <typename Scalar, int Rows, int Cols>
size_t MatrixHash<Scalar, Rows, Cols>::operator()(const Eigen::Matrix<Scalar, Rows, Cols> &matrix) const {
  size_t seed = 0;
  for (Eigen::Index i = 0; i < matrix.size(); ++i) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    Scalar elem = matrix.data()[i];
    // NOLINTNEXTLINE(*-magic-numbers)
    seed ^= std::hash<Scalar>()(elem) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
  }
  return seed;
}

template struct MatrixHash<FLOAT_T, 4, 1>;
template struct MatrixHash<FLOAT_T, 3, 1>;

template <int D>
[[nodiscard]] std::array<VectorF<static_cast<unsigned>(D)>, 1U << static_cast<unsigned>(D)>
allCorners(const Eigen::AlignedBox<FLOAT_T, D> &box) {
  constexpr unsigned N = 1U << static_cast<unsigned>(D);
  std::array<VectorF<static_cast<unsigned>(D)>, N> result;
  for (size_t i = 0; i < N; ++i) {
    for (unsigned j = 0; j < D; ++j)
      result.at(i)[j] = (i & (1U << j)) ? box.max()[j] : box.min()[j];
  }
  return result;
}

template std::array<VectorF<4>, 16> allCorners(const Eigen::AlignedBox<FLOAT_T, 4> &box);

template <int D, int N>
[[nodiscard]] Eigen::Matrix<FLOAT_T, D, D - N> kernel(const Eigen::Matrix<FLOAT_T, D, N> &matrix) noexcept {
  const Eigen::JacobiSVD<Eigen::Matrix<FLOAT_T, D, N>> svd(matrix, Eigen::ComputeFullU);
  return svd.matrixU().rightCols(D - N);
}

template Eigen::Matrix<FLOAT_T, 4, 1> kernel(const Eigen::Matrix<FLOAT_T, 4, 3> &matrix);
} // namespace stmesh

// NOLINTEND(misc-include-cleaner)