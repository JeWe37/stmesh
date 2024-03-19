#ifndef STMESH_MARCHING_HYPERCUBES_HPP
#define STMESH_MARCHING_HYPERCUBES_HPP

#include <array>
#include <optional>

#include "utility.hpp"

namespace stmesh {
std::optional<Vector4F> surfaceRayIntersection(const Vector4F &start, const Vector4F &end,
                                               const std::array<bool, (1U << 4U)> &corner_values, bool *inside);
}; // namespace stmesh
#endif