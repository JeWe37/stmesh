#ifndef STMESH_GEOMETRY_HELPERS_HPP
#define STMESH_GEOMETRY_HELPERS_HPP

#include <boost/geometry.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/parameters.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "utility.hpp"

namespace stmesh {
namespace bg = boost::geometry; ///< Alias for boost::geometry

using BGPoint = bg::model::point<FLOAT_T, 4, bg::cs::cartesian>; ///< A boost::geometry point in 4D
using BGBox = bg::model::box<BGPoint>;                           ///< A boost::geometry box in 4D

/// Converts an Eigen vector to a boost::geometry point
/**
 * Converts an Eigen vector to a boost::geometry point. This function is used to convert Eigen vectors to
 * boost::geometry points, which can be used to query the rtree.
 *
 * @param vector The vector to convert
 * @return The point
 */
[[nodiscard]] BGPoint pointFromVector(const Vector4F &vector) noexcept;

/// Converts a boost::geometry point to an Eigen vector
/**
 * Converts a boost::geometry point to an Eigen vector. Useful for converting points from the rtree to Eigen
 * vectors.
 *
 * @param point The point to convert
 * @return The vector
 */
[[nodiscard]] Vector4F vectorFromPoint(const BGPoint &point) noexcept;

/// Converts an Eigen aligned box to a boost::geometry box
/**
 * Converts an Eigen aligned box to a boost::geometry box. Useful for converting bounding boxes to boost::geometry
 * boxes, which can be inserted into the rtree.
 *
 * @param aabb The aligned box to convert
 * @return The box
 */
[[nodiscard]] BGBox boxFromAABB(const Eigen::AlignedBox<FLOAT_T, 4> &aabb) noexcept;
} // namespace stmesh
#endif