#include "stmesh/boundary_region_manager.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>

#include "stmesh/sdf.hpp"
#include "stmesh/utility.hpp"

namespace stmesh {
size_t HypercubeBoundaryManager::addBoundaryRegion(const HyperCube4 &boundary_region) {
  boundary_regions_.emplace_back(boundary_region);
  return boundary_regions_.size() + 1;
}

[[nodiscard]] size_t HypercubeBoundaryManager::findBoundaryRegion(const Vector4F &point) const noexcept {
  if (auto it = std::ranges::find_if(boundary_regions_,
                                     [&point](const auto &region) { return region.signedDistance(point) < FLOAT_T{}; });
      it != boundary_regions_.end())
    return static_cast<size_t>(std::distance(boundary_regions_.begin(), it) + 1);
  return 0;
}
} // namespace stmesh
