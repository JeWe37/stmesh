#ifndef STMESH_GEOMETRIC_SIMPLEX_HPP
#define STMESH_GEOMETRIC_SIMPLEX_HPP

#include <array>
#include <concepts> // IWYU pragma: keep

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "sdf.hpp"
#include "utility.hpp"

namespace stmesh {
template <unsigned D, unsigned N = D + 1> class GeometricSimplex {
  Eigen::Matrix<FLOAT_T, static_cast<int>(D), static_cast<int>(N)> vertices_;

public:
  explicit GeometricSimplex(const std::ranges::sized_range auto &vertices) : vertices_(D, vertices.size()) {
    Eigen::Index i = 0;
    for (const VectorF<D> &vert : vertices)
      vertices_.col(i++) = vert;
  }

  explicit GeometricSimplex(const Eigen::Matrix<FLOAT_T, static_cast<int>(D), static_cast<int>(N)> &vertices);

  GeometricSimplex();

  [[nodiscard]] const Eigen::Matrix<FLOAT_T, static_cast<int>(D), static_cast<int>(N)> &vertices() const noexcept;

  [[nodiscard]] Eigen::Matrix<FLOAT_T, static_cast<int>(N + 1), static_cast<int>(N + 1)>
  caleyMengerMatrix() const noexcept;

  [[nodiscard]] FLOAT_T caleyMengerDeterminant() const noexcept;

  [[nodiscard]] FLOAT_T shortestEdgeLength() const noexcept;

  [[nodiscard]] FLOAT_T content() const noexcept;

  [[nodiscard]] HyperSphere<D> circumsphere() const noexcept;

  [[nodiscard]] FLOAT_T radiusEdgeRatio() const noexcept;

  [[nodiscard]] FLOAT_T quality() const noexcept;

  // [[nodiscard]] AffineSubspace<D, N - 1> affineSubspace() const noexcept { return {vertices_}; }

  [[nodiscard]] Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(D)> normalRay() const noexcept
  requires(D == N);

  template <unsigned M> [[nodiscard]] std::array<GeometricSimplex<D, M>, nChoosek(N, M)> subSimplices() const noexcept {
    std::array<GeometricSimplex<D, M>, nChoosek(N, M)> result;
    std::array<bool, N> bitmask{};
    std::fill_n(bitmask.begin(), M, true);

    size_t i = 0;
    do {
      Eigen::Matrix<FLOAT_T, static_cast<int>(D), static_cast<int>(M)> vertices_matrix(D, M);
      Eigen::Index k = 0;
      for (size_t j = 0; j < N; ++j) {
        if (bitmask[j])
          vertices_matrix.col(k++) = vertices_.col(static_cast<Eigen::Index>(j));
      }
      result[i++] = GeometricSimplex<D, M>(vertices_matrix);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    return result;
  }

  [[nodiscard]] Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(D)> toParameterizedLine() const noexcept
  requires(N == 2);

  [[nodiscard]] CGAL::Polyhedron_3<CGAL::Cartesian<FLOAT_T>>
  planeCut(const Eigen::Hyperplane<FLOAT_T, static_cast<int>(D)> &plane) const
      // NOLINTNEXTLINE(*-magic-numbers)
  requires(D == 4 && N == 5);

  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> boundingBox() const noexcept;

  template <unsigned M = N, unsigned L = 1> [[nodiscard]] auto allSubSimplicies() const noexcept {
    if constexpr (M == L)
      return std::make_tuple(subSimplices<M>());
    else
      // smallest to largest dimension
      return std::tuple_cat(allSubSimplicies<M - 1, L>(), std::make_tuple(subSimplices<M>()));
  }

  [[nodiscard]] bool wellShaped(FLOAT_T rho_bar, FLOAT_T tau_bar) const noexcept;

  [[nodiscard]] bool sliver(FLOAT_T rho_bar, FLOAT_T tau_bar) const noexcept;

  [[nodiscard]] unsigned sliverSimplex(FLOAT_T rho_bar, FLOAT_T tau_bar) const noexcept;

  [[nodiscard]] bool smallSliverSimplex(FLOAT_T rho_bar, FLOAT_T tau_bar, FLOAT_T max_radius) const noexcept;
};
} // namespace stmesh

#endif