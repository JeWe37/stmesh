#ifndef STMESH_MESH_PROJECT_HPP
#define STMESH_MESH_PROJECT_HPP

#include <cstddef>
#include <utility>

#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/index/parameters.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "geometry_helpers.hpp"
#include "problem_types.hpp"
#include "triangulation.hpp"
#include "utility.hpp"

namespace stmesh {
/// A class to get a value at a point in the mesh
/**
 * A class to get a value at a point in the mesh. This class uses an rtree to efficiently query the vertices within a
 * radius of a point. The value is then interpolated based on the barycentric coordinates of the point in the simplex.
 */
class MeshProjector {
  using Tree = bg::index::rtree<std::pair<BGBox, size_t>,
                                // NOLINTNEXTLINE(*-magic-numbers)
                                bg::index::rstar<16>>;

  Tree tree_;
  const TriangulationFromMixdWithData *triangulation_;

public:
  /// A constructor from a triangulation
  /**
   * A constructor from a triangulation. This constructor creates a mesh projector from a triangulation. The simplices
   * of the triangulation are inserted into the rtree. Its data contains the values that will be projected.
   *
   * @param triangulation The triangulation to create the mesh projector from
   */
  explicit MeshProjector(const TriangulationFromMixdWithData *triangulation);

  /// Projects a point to the mesh
  /**
   * Projects a point to the mesh. This function projects a point to the mesh. The value at the point is interpolated
   * based on the barycentric coordinates of the point in the simplex.
   *
   * @param point The point to project
   * @return The value at the point
   */
  [[nodiscard]] Eigen::Vector<FLOAT_T, Eigen::Dynamic> project(const Vector4F &point) const;

  /// Gets the problem type of the mesh projector
  /**
   * Gets the problem type of the mesh projector. It originates from the triangulation.
   *
   * @return The problem type
   */
  [[nodiscard]] const ProblemType &problemType() const noexcept;
};
} // namespace stmesh

#endif // STMESH_MESH_PROJECT_HPP