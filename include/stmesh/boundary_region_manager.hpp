#ifndef STMESH_BOUNDARY_REGION_MANAGER_HPP
#define STMESH_BOUNDARY_REGION_MANAGER_HPP

#include <cstddef>
#include <vector>

#include "sdf.hpp"
#include "stmesh/utility.hpp"

namespace stmesh {
template <typename T>
concept BoundaryRegionManager = requires(const T t, Vector4F vec) {
  { t.findBoundaryRegion(vec) } -> std::convertible_to<size_t>;
};

class HypercubeBoundaryManager {
  std::vector<HyperCube4> boundary_regions_;

public:
  size_t addBoundaryRegion(const HyperCube4 &boundary_region);

  [[nodiscard]] size_t findBoundaryRegion(const Vector4F &point) const noexcept;
};
} // namespace stmesh

#endif