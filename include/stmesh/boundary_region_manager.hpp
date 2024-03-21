#ifndef STMESH_BOUNDARY_REGION_MANAGER_HPP
#define STMESH_BOUNDARY_REGION_MANAGER_HPP

#include <cstddef>
#include <vector>

#include "sdf.hpp"
#include "stmesh/utility.hpp"

namespace stmesh {
/// Concept for a boundary region manager
/**
 * Concept for a boundary region manager. A boundary region manager is a class that manages boundary regions. Inside a
 * boundary region, boundary faces are assigned a specific boundary id, e.g. for use in a MIXD file.
 *
 * @tparam T The type of the boundary region manager
 */
template <typename T>
concept BoundaryRegionManager = requires(const T t, Vector4F vec) {
  { t.findBoundaryRegion(vec) } -> std::convertible_to<size_t>;
};

/// A boundary region manager for no boundary regions
/**
 * A boundary region manager for no boundary regions. This boundary region manager always returns 0 for any point.
 */
class NoopBoundaryManager {
public:
  /// Find the boundary region of a point
  /**
   * This function always returns 0. The parameter is unused.
   *
   * @return The index of the boundary region that the point is inside
   */
  [[nodiscard]] static size_t findBoundaryRegion(const Vector4F &) noexcept;
};

/// A boundary region manager for a hypercube boundary regions
/**
 * A boundary region manager for a hypercube boundary regions. This boundary region manager assigns boundary regions to
 * points based on the hypercube that the point is inside.
 */
class HypercubeBoundaryManager {
  std::vector<HyperCube4> boundary_regions_;

public:
  /// Add a boundary region to the boundary region manager
  /**
   * Add a boundary region to the boundary region manager. The boundary region is a hypercube that defines the boundary
   * region. Boundary regions added later have a higher index than boundary regions added earlier.
   *
   * @param boundary_region The boundary region to add
   * @return The index of the boundary region
   */
  size_t addBoundaryRegion(const HyperCube4 &boundary_region);

  /// Find the boundary region of a point
  /**
   * Find the boundary region of a point. The boundary region is the index of the boundary region that the point is
   * inside.
   *
   * @param point The point to find the boundary region of
   * @return The index of the boundary region that the point is inside
   */
  [[nodiscard]] size_t findBoundaryRegion(const Vector4F &point) const noexcept;
};
} // namespace stmesh

#endif