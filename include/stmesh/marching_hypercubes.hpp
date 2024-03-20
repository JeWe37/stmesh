#ifndef STMESH_MARCHING_HYPERCUBES_HPP
#define STMESH_MARCHING_HYPERCUBES_HPP

#include <array>
#include <optional>

#include "utility.hpp"

namespace stmesh {
/// Computes the intersection of a ray with a surface defined by marching hypercubes
/**
 * Computes the intersection of a ray with a surface defined by marching hypercubes. The surface is defined by the
 * values of the corners of the hypercubes. The ray is defined by a start and end point. The intersection is computed
 * using the marching hypercubes algorithm.
 * Should the ray not intersect the surface, std::nullopt is returned. Otherwise, the intersection point is returned.
 * To check if the ray is inside the surface, a pointer to a boolean can be passed. If the ray is inside the surface,
 * the boolean will be set to true. If the pointer is nullptr, the inside check is not performed.
 *
 * @param start The start point of the ray
 * @param end The end point of the ray
 * @param corner_values The values of the corners of the hypercubes
 * @param inside A pointer to a boolean that will be set to true if the ray is inside the surface
 * @return The intersection point of the ray with the surface, or std::nullopt if there is no intersection
 */
[[nodiscard]] std::optional<Vector4F> surfaceRayIntersection(const Vector4F &start, const Vector4F &end,
                                                             const std::array<bool, (1U << 4U)> &corner_values,
                                                             bool *inside);
}; // namespace stmesh
#endif