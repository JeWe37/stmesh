#ifndef STMESH_BOUNDARY_REGION_MANAGER_HPP
#define STMESH_BOUNDARY_REGION_MANAGER_HPP

#include <vector>

#include "sdf.hpp"
#include "stmesh/utility.hpp"

namespace stmesh {
class BoundaryRegionManager {
  std::vector<HyperCube4> boundary_regions_;

public:
  size_t addBoundaryRegion(const HyperCube4 &boundary_region);

  [[nodiscard]] size_t findBoundaryRegion(const Vector4F &point) const;
};
} // namespace stmesh

#endif