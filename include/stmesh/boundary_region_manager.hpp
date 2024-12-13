#ifndef STMESH_BOUNDARY_REGION_MANAGER_HPP
#define STMESH_BOUNDARY_REGION_MANAGER_HPP

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <vector>

#include "sdf.hpp"
#include "utility.hpp"
#include "writable_triangulation.hpp"

namespace stmesh {
/// Concept for a boundary region manager
/**
 * Concept for a boundary region manager. A boundary region manager is a class that manages boundary regions. Inside a
 * boundary region, boundary faces are assigned a specific boundary id, e.g. for use in a MIXD file.
 *
 * @tparam T The type of the boundary region manager
 */
template <typename T>
concept BoundaryRegionManager = requires(const T t, detail::DummyCell cell, size_t j) {
  { t.findBoundaryRegion(cell, j) } -> std::convertible_to<size_t>;
};

/// A boundary region manager for no boundary regions
/**
 * A boundary region manager for no boundary regions. This boundary region manager always returns 0 for any cell.
 */
class NoopBoundaryManager {
public:
  /// Find the boundary region of a point
  /**
   * This function always returns 0. The parameter is unused.
   *
   * @return Zero
   */
  [[nodiscard]] static size_t findBoundaryRegion(const WritableCell auto & /*unused*/, size_t /*unused*/) noexcept {
    return 0;
  }
};

static_assert(BoundaryRegionManager<NoopBoundaryManager>);

/// A boundary region manager for a hypercube boundary regions
/**
 * A boundary region manager for a hypercube boundary regions. This boundary region manager assigns boundary regions to
 * cell faces based on the hypercube that their verticies are inside of.
 */
class HypercubeBoundaryManager {
  std::vector<HyperCube4> boundary_regions_;

  [[nodiscard]] size_t pointBoundaryRegion(const Vector4F &point) const noexcept;

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

  /// Find the boundary region of a cell
  /**
   * Find the boundary region of a cell. The boundary region is the index of the highest-index boundary region that a
   * vertex of the cell face is inside of.
   *
   * @param cell The cell to find the boundary region of
   * @return The index of the boundary region that the cell is inside
   */
  [[nodiscard]] size_t findBoundaryRegion(const WritableCell auto &cell, size_t j) const noexcept {
    size_t max_index = 0;
    for (size_t i = 0; i < 5; ++i) {
      if (i != j)
        max_index = std::max(max_index,
                             pointBoundaryRegion(cell.geometricSimplex().vertices().col(static_cast<Eigen::Index>(i))));
    }
    return max_index;
  }
};

static_assert(BoundaryRegionManager<HypercubeBoundaryManager>);
} // namespace stmesh

#endif