#ifndef STMESH_SDF_HPP
#define STMESH_SDF_HPP

#include <concepts> // IWYU pragma: keep

#include <Eigen/Geometry>
#include <fmt/core.h>

#include "sdf_mixins.hpp"
#include "utility.hpp"

namespace stmesh {
/// A concept for a signed distance field
/**
 * A concept for a signed distance field. A signed distance field is a function that returns the distance to the
 * surface of a shape, the signed distance to the surface of a shape, the normal of the surface of a shape, and the
 * bounding box of the shape.
 *
 * @tparam T The type to check
 * @tparam D The dimension of the space
 */
template <typename T, unsigned D>
concept SignedDistanceField = requires(T t, VectorF<D> vec) {
  { t.distance(vec) } -> std::convertible_to<FLOAT_T>;
  { t.signedDistance(vec) } -> std::convertible_to<FLOAT_T>;
  { t.normal(vec) } -> std::convertible_to<VectorF<D>>;
  { t.boundingBox() } -> std::convertible_to<Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)>>;
};

/// A tag for exact signed distance fields
/**
 * A tag for exact signed distance fields. An exact signed distance field is a signed distance field that returns
 * exact values for the signed distance, rather than simply an upper bound. This is useful for algorithms that require
 * exact signed distance fields, such as finding the closest point on the surface of a shape.
 *
 * @tparam T The signed distnace field to declare as exact
 */
template <typename T>
requires SignedDistanceField<T, T::kDimension>
struct ExactSDFTag {
  static constexpr bool kValue = false; ///< Whether the type is an exact signed distance field
};

/// A concept for an exact signed distance field
/**
 * A concept for an exact signed distance field. An exact signed distance field is a signed distance field that returns
 * exact values for the signed distance, rather than simply an upper bound. This is useful for algorithms that require
 * exact signed distance fields, such as finding the closest point on the surface of a shape.
 * This concept checks if the ExactSDFTag is defined for the type.
 *
 * @tparam T The type to check
 */
template <typename T>
concept ExactSDF = requires { ExactSDFTag<T>::kValue; };

/// An SDF for a hyper sphere
/**
 * An SDF for a hyper sphere. Provides some additional functionality for a hyper sphere, such as scaling and
 * equality comparison, as well as getters.
 *
 * @tparam D The dimension of the space
 */
template <unsigned D>
class HyperSphere : public stmesh::DistanceMixin<HyperSphere<D>, D>, public SamplingMixin<HyperSphere<D>, D> {
  FLOAT_T radius_;
  VectorF<D> center_;

public:
  static inline constexpr unsigned kDimension = D; ///< The dimension of the space

  /// Construct a hyper sphere
  /**
   * Construct a hyper sphere. The hyper sphere is defined by a radius and a center.
   *
   * @param radius The radius of the hyper sphere
   * @param center The center of the hyper sphere
   */
  HyperSphere(FLOAT_T radius, const VectorF<D> &center) noexcept;

  /// Get the radius of the hyper sphere
  /**
   * Get the radius of the hyper sphere.
   *
   * @return The radius of the hyper sphere
   */
  [[nodiscard]] FLOAT_T radius() const noexcept;

  /// Scale the hyper sphere
  /**
   * Scale the hyper sphere. The radius of the hyper sphere is scaled by the factor.
   * The center of the hyper sphere is not changed.
   *
   * @param factor The factor to scale the hyper sphere by
   */
  void scale(FLOAT_T factor) noexcept;

  /// Get the center of the hyper sphere
  /**
   * Get the center of the hyper sphere.
   *
   * @return The center of the hyper sphere
   */
  [[nodiscard]] const VectorF<D> &center() const noexcept;

  /// Calculate the signed distance of a point from the hyper sphere
  /**
   * Calculate the signed distance of a point from the hyper sphere. The signed distance is the distance of the point
   * from the surface of the hyper sphere, with the sign indicating whether the point is inside(negative) or
   * outside(positive) the hyper sphere.
   *
   * @param point The point at which to calculate the signed distance
   * @return The signed distance of the point from the hyper sphere
   */
  [[nodiscard]] FLOAT_T signedDistance(const VectorF<D> &point) const noexcept;

  /// Calculate the normal of the hyper sphere at a point
  /**
   * Calculate the normal of the hyper sphere at a point. The normal is the direction from the point to the center of
   * the hyper sphere, or away from the center if the point is inside the hyper sphere.
   *
   * @param point The point at which to calculate the normal
   * @return The normal of the hyper sphere at the point
   */
  [[nodiscard]] VectorF<D> normal(const VectorF<D> &point) const noexcept;

  /// Get the bounding box of the hyper sphere
  /**
   * Get the bounding box of the hyper sphere. The bounding box is a minimal box that contains the hyper sphere.
   *
   * @return The bounding box of the hyper sphere
   */
  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> boundingBox() const noexcept;

  /// Check if the hyper sphere is equal to another hyper sphere
  /**
   * Check if the hyper sphere is equal to another hyper sphere. Two hyper spheres are equal if they have the same
   * radius and center.
   *
   * @param other The other hyper sphere to compare to
   * @return Whether the hyper sphere is equal to the other hyper sphere
   */
  bool operator==(const HyperSphere<D> &other) const noexcept;
};

template <unsigned D> struct ExactSDFTag<HyperSphere<D>> {
  static constexpr bool kValue = true;
};

using HyperSphere3 = HyperSphere<3>; ///< A hyper sphere in 3D
using HyperSphere4 = HyperSphere<4>; ///< A hyper sphere in 4D

/// An SDF for a hyper cube
/**
 * An SDF for a hyper cube. Only offers the basic functionality of an SDF, such as signed distance, normal, and
 * bounding box.
 *
 * @tparam D The dimension of the space
 */
template <unsigned D>
class HyperCube : public stmesh::DistanceMixin<HyperCube<D>, D>, public CentralDifferenceNormalMixin<HyperCube<D>, D> {
  VectorF<D> min_;
  VectorF<D> max_;

public:
  static inline constexpr unsigned kDimension = D; ///< The dimension of the space

  /// Construct a hyper cube
  /**
   * Construct a hyper cube. The hyper cube is defined by a minimum and maximum point.
   *
   * @param min The minimum point of the hyper cube
   * @param max The maximum point of the hyper cube
   */
  HyperCube(const VectorF<D> &min, const VectorF<D> &max) noexcept;

  /// Get the signed distance of a point from the hyper cube
  /**
   * Get the signed distance of a point from the hyper cube. The signed distance is the distance of the point from the
   * surface of the hyper cube, with the sign indicating whether the point is inside(negative) or outside(positive) the
   * hyper cube.
   *
   * @param point The point at which to calculate the signed distance
   * @return The signed distance of the point from the hyper cube
   */
  [[nodiscard]] FLOAT_T signedDistance(VectorF<D> point) const noexcept;

  /// Get the bounding box of the hyper cube
  /**
   * Get the bounding box of the hyper cube. The bounding box is simply the hyper cube itself.
   *
   * @return The bounding box of the hyper cube
   */
  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> boundingBox() const noexcept;

  friend class ::fmt::formatter<HyperCube<D>>;
};

template <unsigned D> struct ExactSDFTag<HyperCube<D>> {
  static constexpr bool kValue = true;
};

using HyperCube4 = HyperCube<4>; ///< A hyper cube in 4D
} // namespace stmesh

/// A formatter for stmesh::HyperCube<D>
/**
 * A formatter for stmesh::HyperCube<D>. This formatter allows stmesh::HyperCube<D> to be used with fmt::format.
 *
 * @tparam D The dimension of the space
 */
template <unsigned D> class fmt::formatter<stmesh::HyperCube<D>> {
public:
  /// Parse the format string
  constexpr static auto parse(format_parse_context &ctx) { return ctx.begin(); }
  /// Format the hyper cube
  /**
   * Format the hyper cube. The hyper cube is formatted as a string containing the minimum and maximum points of the
   * hyper cube.
   *
   * @param hypercube The hyper cube to format
   * @param ctx The format context
   * @return The iterator to the end of the formatted string
   */
  template <typename Context> static constexpr auto format(const stmesh::HyperCube<D> &hypercube, Context &ctx) {
    return fmt::format_to(ctx.out(), "(({}),({}))", hypercube.min_.transpose(), hypercube.max_.transpose());
  }
};

#endif