#include "stmesh/geometry_helpers.hpp"

#include "stmesh/utility.hpp"

namespace stmesh {
[[nodiscard]] auto pointFromVector(const Vector4F &vector) noexcept -> BGPoint {
  BGPoint point;
  point.set<0>(vector[0]);
  point.set<1>(vector[1]);
  point.set<2>(vector[2]);
  point.set<3>(vector[3]);
  return point;
}

[[nodiscard]] auto vectorFromPoint(const BGPoint &point) noexcept -> Vector4F {
  return {static_cast<FLOAT_T>(point.get<0>()), static_cast<FLOAT_T>(point.get<1>()),
          static_cast<FLOAT_T>(point.get<2>()), static_cast<FLOAT_T>(point.get<3>())};
}

[[nodiscard]] auto boxFromAABB(const Eigen::AlignedBox<FLOAT_T, 4> &aabb) noexcept -> BGBox {
  return {pointFromVector(aabb.min()), pointFromVector(aabb.max())};
}
} // namespace stmesh