#ifndef STMESH_SDF_HPP
#define STMESH_SDF_HPP

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdlib>
#include <utility>

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
concept SignedDistanceField = SignedDistance<T, D> && requires(T t, VectorF<D> vec) {
  { t.distance(vec) } -> std::convertible_to<FLOAT_T>;
  { t.normal(vec) } -> std::convertible_to<VectorF<D>>;
  { t.boundingBox() } -> std::convertible_to<Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)>>;
};

/// A concept for an SDF of any dimension
/**
 * A concept for an SDF of any dimension. This concept checks if the type is a signed distance field of any dimension.
 *
 * @tparam T The type to check
 */
template <typename T>
concept SDF = SignedDistanceField<T, T::kDimension>;

/// A tag for exact signed distance fields
/**
 * A tag for exact signed distance fields. An exact signed distance field is a signed distance field that returns
 * exact values for the signed distance, rather than simply an upper bound. This is useful for algorithms that require
 * exact signed distance fields, such as finding the closest point on the surface of a shape.
 *
 * @tparam T The signed distance field to declare as exact
 */
template <SDF T> struct ExactSDFTag {
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
  static constexpr unsigned kDimension = D; ///< The dimension of the space

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
  static constexpr unsigned kDimension = D; ///< The dimension of the space

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

/// An SDF for an extruded lower-dimensional SDF
/**
 * An SDF for an extruded lower-dimensional SDF. Only offers the basic functionality of an SDF, such as signed distance,
 * normal, and bounding box.
 *
 * @tparam T The base SDF to extrude
 */
template <ExactSDF T> class ExtrudedSDF : public stmesh::DistanceMixin<ExtrudedSDF<T>, T::kDimension + 1> {
public:
  static constexpr unsigned kDimension = T::kDimension + 1; ///< The dimension of the space
private:
  T base_;
  Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(kDimension)> axis_;
  Eigen::Matrix<FLOAT_T, static_cast<int>(T::kDimension), static_cast<int>(kDimension)> projection_;
  Eigen::Matrix<FLOAT_T, static_cast<int>(kDimension), static_cast<int>(kDimension)> inverse_;

  std::pair<VectorF<kDimension>, FLOAT_T> project(const VectorF<kDimension> &point) const noexcept {
    FLOAT_T lambda = axis_.direction().dot(point - axis_.origin()) / axis_.direction().squaredNorm();
    VectorF<T::kDimension> projected = projection_ * (point - lambda * axis_.direction());
    FLOAT_T projected_distance = base_.signedDistance(projected);
    bool lambda_inside = lambda >= 0.0 && lambda <= 1.0;
    if (projected_distance <= 0.0) {
      FLOAT_T start_distance = lambda * axis_.direction().norm();
      FLOAT_T end_distance = (1.0 - lambda) * axis_.direction().norm();
      if (start_distance <= -projected_distance || end_distance <= -projected_distance) {
        if (start_distance <= end_distance)
          return {-axis_.direction().normalized(), -start_distance};
        return {axis_.direction().normalized(), -end_distance};
      }
    }
    // side is closest
    VectorF<T::kDimension> base_normal = base_.normal(projected);
    VectorF<T::kDimension> closest_point = projected - projected_distance * base_normal;
    VectorF<kDimension> unprojected =
        inverse_ * (VectorF<kDimension>() << closest_point, std::max(std::min(lambda, 1.0), 0.0)).finished();
    FLOAT_T sign = lambda_inside && projected_distance <= 0.0 ? -1.0 : 1.0;
    return {lambda_inside ? (projection_.transpose() * base_normal).eval()
                          : ((point - unprojected).normalized() * sign).eval(),
            (point - unprojected).norm() * sign};
  }

  static Eigen::Matrix<FLOAT_T, static_cast<int>(kDimension), static_cast<int>(kDimension)>
  getInverse(const VectorF<kDimension> &direction,
             const Eigen::Matrix<FLOAT_T, static_cast<int>(T::kDimension), static_cast<int>(kDimension)> &projection) {
    Eigen::Matrix<FLOAT_T, static_cast<int>(kDimension), static_cast<int>(kDimension)> mat;
    mat << projection, direction.transpose() / direction.squaredNorm();
    return mat.inverse();
  }

public:
  /// Construct an extruded SDF
  /**
   * Construct an extruded SDF. The extruded SDF is defined by a base SDF, and an axis along which to extrude the base
   * SDF.
   *
   * @param base The base SDF
   * @param axis The axis of the extrusion
   */
  ExtrudedSDF(const T &base, const Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(kDimension)> &axis) noexcept
      : base_(base), axis_(axis),
        projection_(axis.direction().transpose().fullPivLu().kernel().transpose().rowwise().normalized()),
        inverse_(getInverse(axis.direction(), projection_)) {}

  /// Get the signed distance of a point from the extruded SDF
  /**
   * Get the signed distance of a point from the extruded SDF. The signed distance is the distance of the point from the
   * surface of the extruded SDF, with the sign indicating whether the point is inside(negative) or outside(positive)
   * the extruded SDF.
   *
   * @param point The point at which to calculate the signed distance
   * @return The signed distance of the point from the extruded SDF
   */
  [[nodiscard]] FLOAT_T signedDistance(const VectorF<kDimension> &point) const noexcept {
    return project(point).second;
  }

  /// Get the normal of the extruded SDF at a point
  /**
   * Get the normal of the extruded SDF at a point. The normal is the direction from the point away from the surface of
   * the extruded SDF.
   *
   * @param point The point at which to calculate the normal
   * @return The normal of the extruded SDF at the point
   */
  [[nodiscard]] VectorF<kDimension> normal(const VectorF<kDimension> &point) const noexcept {
    return project(point).first;
  }

  /// Get the bounding box of the extruded SDF
  /**
   * Get the bounding box of the extruded SDF. Note the bounding box is only an approximation, as the extruded SDF is
   * not a simple shape. The true minimal bounding box is contained within it.
   *
   * @return The bounding box of the extruded SDF
   */
  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(kDimension)> boundingBox() const noexcept {
    Eigen::Matrix<FLOAT_T, static_cast<int>(kDimension), static_cast<int>(1U << kDimension)> points;
    Eigen::Index i = 0;
    VectorF<kDimension> appended;
    for (FLOAT_T lambda : {0.0, 1.0}) {
      for (const VectorF<T::kDimension> &corner : allCorners(base_.boundingBox())) {
        appended << corner, lambda;
        points.col(i++) = inverse_ * appended;
      }
    }
    return {points.rowwise().minCoeff(), points.rowwise().maxCoeff()};
  }
};

template <ExactSDF T> struct ExactSDFTag<ExtrudedSDF<T>> {
  static constexpr bool kValue = true;
};

/// An SDF for a cylinder
/**
 * An SDF for a cylinder. Only offers the basic functionality of an SDF, such as signed distance, normal, and
 * bounding box. Only works in 3D.
 */
class CylinderSDF : public stmesh::DistanceMixin<CylinderSDF, 3>, public stmesh::SamplingMixin<CylinderSDF, 3> {
  FLOAT_T radius_;
  Eigen::ParametrizedLine<FLOAT_T, 3> axis_;

public:
  static constexpr unsigned kDimension = 3; ///< The dimension of the space

  /// Construct a cylinder SDF
  /**
   * Construct a cylinder SDF. The cylinder SDF is defined by a radius and an axis.
   *
   * @param radius The radius of the cylinder
   * @param axis The axis of the cylinder
   */
  CylinderSDF(FLOAT_T radius, const Eigen::ParametrizedLine<FLOAT_T, 3> &axis) noexcept;

  /// Get the signed distance of a point from the cylinder
  /**
   * Get the signed distance of a point from the cylinder. The signed distance is the distance of the point from the
   * surface of the cylinder, with the sign indicating whether the point is inside(negative) or outside(positive) the
   * cylinder.
   *
   * @param point The point at which to calculate the signed distance
   * @return The signed distance of the point from the cylinder
   */
  [[nodiscard]] FLOAT_T signedDistance(const VectorF<3> &point) const noexcept;

  /// Get the normal of the cylinder at a point
  /**
   * Get the normal of the cylinder at a point. The normal is the direction from the point away from the axis.
   *
   * @param point The point at which to calculate the normal
   * @return The normal of the cylinder at the point
   */
  [[nodiscard]] Vector3F normal(const Vector3F &point) const noexcept;

  /// Get the bounding box of the cylinder
  /**
   * Get the bounding box of the cylinder. The bounding box is a minimal box that contains the cylinder.
   *
   * @return The bounding box of the cylinder
   */
  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, 3> boundingBox() const noexcept;
};
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