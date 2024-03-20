#ifndef STMESH_GEOMETRIC_SIMPLEX_HPP
#define STMESH_GEOMETRIC_SIMPLEX_HPP

#include <array>
#include <concepts> // IWYU pragma: keep
#include <cstddef>
#include <utility>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "sdf.hpp"
#include "utility.hpp"

namespace stmesh {
/// A geometric simplex
/**
 * A geometric simplex. A geometric simplex is a generalization of a triangle to higher dimensions. It is defined by a
 * set of vertices, and can be used to calculate various properties of the simplex, such as the content, the
 * circumsphere, and the quality. The simplex may be embedded in any dimension, and may have any number of vertices.
 *
 * @tparam D The dimension of the space
 * @tparam N The number of vertices in the simplex
 */
template <unsigned D, unsigned N = D + 1> class GeometricSimplex {
  Eigen::Matrix<FLOAT_T, static_cast<int>(D), static_cast<int>(N)> vertices_;

public:
  /// Construct a geometric simplex
  /**
   * Construct a geometric simplex. The vertices of the simplex are given as a range of vertices. The vertices are
   * copied into the simplex.
   *
   * @param vertices The vertices of the simplex
   */
  explicit GeometricSimplex(const std::ranges::sized_range auto &vertices) : vertices_(D, vertices.size()) {
    Eigen::Index i = 0;
    for (const VectorF<D> &vert : vertices)
      vertices_.col(i++) = vert;
  }

  /// Construct a geometric simplex
  /**
   * Construct a geometric simplex. The vertices of the simplex are given as a matrix of vertices. The vertices are
   * copied into the simplex.
   *
   * @param vertices The vertices of the simplex
   */
  explicit GeometricSimplex(const Eigen::Matrix<FLOAT_T, static_cast<int>(D), static_cast<int>(N)> &vertices);

  /// Construct a geometric simplex
  /**
   * Construct a geometric simplex. Default constructor, the vertices are not initialized.
   */
  GeometricSimplex();

  /// Get the vertices of the simplex
  /**
   * Get the vertices of the simplex. The vertices are returned as a matrix of vertices.
   *
   * @return The vertices of the simplex
   */
  [[nodiscard]] const Eigen::Matrix<FLOAT_T, static_cast<int>(D), static_cast<int>(N)> &vertices() const noexcept;

  /// Compute the caley menger matrix
  /**
   * Compute the caley menger matrix. The caley menger matrix is a matrix that contains the squared distances between
   * the vertices of the simplex. The caley menger matrix is used to calculate the content of the simplex.
   *
   * @return The caley menger matrix
   */
  [[nodiscard]] Eigen::Matrix<FLOAT_T, static_cast<int>(N + 1), static_cast<int>(N + 1)>
  caleyMengerMatrix() const noexcept;

  /// Compute the caley menger determinant
  /**
   * Compute the caley menger determinant. The caley menger determinant is a value that is used to calculate the content
   * of the simplex. The caley menger determinant is the determinant of the caley menger matrix.
   *
   * @return The caley menger determinant
   */
  [[nodiscard]] FLOAT_T caleyMengerDeterminant() const noexcept;

  /// Compute the length of the shortest edge
  /**
   * Compute the length of the shortest edge. This is used to calculate the radius-edge ratio of the simplex and
   * quality.
   *
   * @return The length of the shortest edge
   */
  [[nodiscard]] FLOAT_T shortestEdgeLength() const noexcept;

  /// Compute the content of the simplex
  /**
   * Compute the content of the simplex. The content is for example the volume of a 3-simplex, or the area for a
   * 2-simplex. The dimension is determined by N, not D. The content is calculated using the caley menger determinant.
   *
   * @return The content of the simplex
   */
  [[nodiscard]] FLOAT_T content() const noexcept;

  /// Compute the circumsphere of the simplex
  /**
   * Compute the circumsphere of the simplex. The circumsphere is the smallest sphere that contains the simplex. The
   * circumsphere is calculated using CGAL's construction and is returned as a hyper sphere.
   *
   * @return The circumsphere of the simplex
   */
  [[nodiscard]] HyperSphere<D> circumsphere() const noexcept;

  /// Compute the radius-edge ratio of the simplex
  /**
   * Compute the radius-edge ratio of the simplex. The radius-edge ratio is the ratio of the radius of the circumsphere
   * to the length of the shortest edge.
   *
   * @return The radius-edge ratio of the simplex
   */
  [[nodiscard]] FLOAT_T radiusEdgeRatio() const noexcept;

  /// Compute the quality of the simplex
  /**
   * Compute the quality of the simplex. The quality is a measure of how well shaped the simplex is. The quality is
   * calculated as the content over the D-th power of the radius-edge ratio.
   *
   * @return The quality of the simplex
   */
  [[nodiscard]] FLOAT_T quality() const noexcept;

  /// Compute the normal ray of the simplex
  /**
   * Compute the normal ray of the simplex. The normal ray is a line that is normal to the simplex and passes through
   * the center of the circumsphere. The normal ray is only defined for D = N, i.e. a hypersurface.
   *
   * @return The normal ray of the simplex
   */
  [[nodiscard]] Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(D)> normalRay() const noexcept
  requires(D == N);

  /// Compute the sub simplices of the simplex
  /**
   * Compute the sub simplices of the simplex. The sub simplices are the simplices that are formed by removing one or
   * more vertices from the simplex. The sub simplices are returned as an array of simplices. Essentailly these are the
   * faces of the simplex of dimension M.
   *
   * @return The sub simplices of the simplex
   */
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

  /// Convert the simplex to an Eigen paramterized line
  /**
   * Convert the simplex to an Eigen paramterized line. This is only defined for N = 2, i.e. a line.
   *
   * @return The paramterized line and the length of the line
   */
  [[nodiscard]] std::pair<Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(D)>, FLOAT_T>
  toParameterizedLine() const noexcept
  requires(N == 2);

  /// Cut the simplex by a hyperplane
  /**
   * Cut the simplex by a hyperplane. The hyperplane is defined by a normal and a distance from the origin. The
   * cut is then returned as a polyhedron with up to 6 vertices. This is only implemented for D==4 and N==5, i.e.
   * cutting a 4-simplex in 4D to yield a polyhedron.
   * The hyperplane must have a normal pointing in the 4th dimension.
   *
   * @param plane The hyperplane to cut the simplex by
   * @return The polyhedron that is the result of the cut
   */
  [[nodiscard]] CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel>
  planeCut(const Eigen::Hyperplane<FLOAT_T, static_cast<int>(D)> &plane) const
  requires(D == 4 && N == 5);

  /// Get the bounding box of the simplex
  /**
   * Get the bounding box of the simplex. The bounding box is a minimal box that contains the simplex.
   *
   * @return The bounding box of the simplex
   */
  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> boundingBox() const noexcept;

  /// Get all sub simplicies of the simplex
  /**
   * Get all sub simplicies of the simplex. These are the sub simplices of any dimension from L to M.
   *
   * @return All sub simplicies of the simplex
   * @sa subSimplices
   */
  template <unsigned M = N, unsigned L = 1> [[nodiscard]] auto allSubSimplicies() const noexcept {
    if constexpr (M == L)
      return std::make_tuple(subSimplices<M>());
    else
      // smallest to largest dimension
      return std::tuple_cat(allSubSimplicies<M - 1, L>(), std::make_tuple(subSimplices<M>()));
  }

  /// Check if the simplex is well shaped
  /**
   * Check if the simplex is well shaped. A simplex is well shaped if the radius-edge ratio is less than rho_bar
   * and the quality is better than tau_bar.
   *
   * @param rho_bar The maximum radius-edge ratio
   * @param tau_bar The minimum quality
   * @return Whether the simplex is well shaped
   */
  [[nodiscard]] bool wellShaped(FLOAT_T rho_bar, FLOAT_T tau_bar) const noexcept;

  /// Check if the simplex is a sliver
  /**
   * Check if the simplex is a sliver. A simplex is a sliver if the radius-edge ratio is less than rho_bar, but the
   * quality is less than tau_bar. Additionally, all sub simplices are checked that they are well shaped.
   *
   * @param rho_bar The maximum radius-edge ratio
   * @param tau_bar The minimum quality
   * @return Whether the simplex is a sliver
   */
  [[nodiscard]] bool sliver(FLOAT_T rho_bar, FLOAT_T tau_bar) const noexcept;

  /// Check if the simplex is a sliver simplex
  /**
   * Check if the simplex is a sliver simplex. A simplex is a sliver simplex any of its sub simplices are slivers.
   *
   * @param rho_bar The maximum radius-edge ratio
   * @param tau_bar The minimum quality
   * @return Whether the simplex is a sliver simplex
   */
  [[nodiscard]] unsigned sliverSimplex(FLOAT_T rho_bar, FLOAT_T tau_bar) const noexcept;

  /// Check if the simplex is a small sliver simplex
  /**
   * Check if the simplex is a small sliver simplex. A simplex is a small sliver simplex if it is a sliver simplex and
   * the circumradius is less than max_radius.
   *
   * @param rho_bar The maximum radius-edge ratio
   * @param tau_bar The minimum quality
   * @param max_radius The maximum circumradius
   * @return Whether the simplex is a small sliver simplex
   */
  [[nodiscard]] bool smallSliverSimplex(FLOAT_T rho_bar, FLOAT_T tau_bar, FLOAT_T max_radius) const noexcept;
};
} // namespace stmesh

#endif