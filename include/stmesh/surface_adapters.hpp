#ifndef STMESH_SURFACE_ADAPTERS_HPP
#define STMESH_SURFACE_ADAPTERS_HPP

#include <array>
#include <concepts> // IWYU pragma: keep
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Geometry>

#include "edt.hpp" // IWYU pragma: keep
#include "marching_hypercubes.hpp"
#include "sdf.hpp"
#include "utility.hpp"

namespace stmesh {
/// A concept for a surface adapter
/**
 * A concept for a surface adapter. A surface adapter is a class that provides a way to interact with a surface. It
 * provides functions to find the closest point on the surface to a given point, to check if a point is inside the
 * surface, to check if the surface is intersected by a sphere, to get the bounding box of the surface, and to raycast
 * against the surface.
 *
 * @tparam T The type to check
 * @tparam D The dimension of the space
 */
template <typename T, unsigned D>
concept SurfaceAdapter = requires(const T t, HyperSphere<D> sphere, VectorF<D> vec,
                                  Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(D)> line, FLOAT_T num) {
  { t.closestPoint(vec) } -> std::convertible_to<VectorF<D>>;
  { t.inside(vec) } -> std::convertible_to<bool>;
  { t.intersectedBySphere(sphere) } -> std::convertible_to<bool>;
  { t.boundingBox() } -> std::convertible_to<Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)>>;
  { t.raycast(vec, vec, num) } -> std::convertible_to<std::optional<VectorF<D>>>;
  { t.raycast(line, num) } -> std::convertible_to<std::optional<VectorF<D>>>;
};

/// A concept for a surface adapter in 4 dimensions
/**
 * A concept for a surface adapter in 4 dimensions. This is a specialization of the SurfaceAdapter concept.
 *
 * @tparam T The type to check
 * @sa SurfaceAdapter
 */
template <typename T>
concept SurfaceAdapter4 = SurfaceAdapter<T, 4>;

/// A mixin used to provide the second form of the raycast function
/**
 * A mixin used to provide the second form of the raycast function. This mixin requires that the class has a raycast
 * function that takes a start point, a direction, and a maximum distance. It then provides a raycast function that
 * takes a parametrized line as the ray and a maximum distance.
 *
 * @tparam T The type to mix in
 * @tparam D The dimension of the space
 */
template <typename T, unsigned D> class LineRaycastMixin {
public:
  /// Raycast function that takes a parametrized line as the ray and a maximum distance
  /**
   * Raycast function that takes a parametrized line as the ray and a maximum distance. This function calls the raycast
   * function that takes a start point, a direction, and a maximum distance.
   *
   * @param ray The parametrized line to cast
   * @param max_distance The maximum distance to cast
   * @return The result of the raycast
   */
  [[nodiscard]] std::optional<VectorF<D>> raycast(const Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(D)> &ray,
                                                  FLOAT_T max_distance) const noexcept {
    return static_cast<const T *>(this)->raycast(ray.origin(), ray.direction(), max_distance);
  }
};

/// A surface adapter for an exact signed distance function
/**
 * A surface adapter for an exact signed distance function. The signed distance function is assumed to be exact, i.e.
 * it is not just an upper bound on the distance to the surface, but the exact distance to the surface. In order to
 * implement the raycast, the adapter uses a simple raymarching algorithm.
 *
 * @tparam T The exact signed distance function to adapt
 */
template <ExactSDF T> class SDFSurfaceAdapter : public LineRaycastMixin<SDFSurfaceAdapter<T>, T::kDimension> {
  using Base = LineRaycastMixin<SDFSurfaceAdapter<T>, T::kDimension>;

public:
  using Surface = T;                                                ///< The exact signed distance function
  static inline constexpr unsigned kDimension = T::kDimension;      ///< The dimension of the space
  constexpr static inline FLOAT_T kDt = static_cast<FLOAT_T>(1e-3); ///< The distance threshold for the raycast

private:
  T surface_;

public:
  /// Construct a new SDFSurfaceAdapter
  /**
   * Construct a new SDFSurfaceAdapter. This constructor forwards its arguments to the constructor of the exact signed
   * distance function.
   *
   * @tparam Args The types of the arguments
   * @param args The arguments to forward
   */
  template <typename... Args> explicit SDFSurfaceAdapter(Args &&...args) : surface_(std::forward<Args>(args)...) {}

  /// Find the closest point on the surface to a given point
  /**
   * Find the closest point on the surface to a given point. This function calls the signedDistance function of the
   * exact signed distance function and uses the normal function to find the closest point.
   *
   * @param point The point to find the closest point to
   * @return The closest point on the surface to the given point
   */
  [[nodiscard]] VectorF<kDimension> closestPoint(const VectorF<kDimension> &point) const noexcept {
    return point - surface_.signedDistance(point) * surface_.normal(point);
  }

  /// Check if a sphere intersects the surface
  /**
   * Check if a sphere intersects the surface. This function calls the signedDistance function of the exact signed
   * distance function and checks if the distance to the center of the sphere is less than the radius of the sphere.
   *
   * @param sphere The sphere to check
   * @return Whether the sphere intersects the surface
   */
  [[nodiscard]] bool intersectedBySphere(const HyperSphere<kDimension> &sphere) const noexcept {
    return sphere.radius() - surface_.distance(sphere.center()) > kDt;
  }

  /// Get the bounding box of the surface
  /**
   * Get the bounding box of the surface. This function calls the boundingBox function of the exact signed distance
   * function.
   *
   * @return The bounding box of the surface
   */
  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(kDimension)> boundingBox() const noexcept {
    return surface_.boundingBox();
  }

  /// Check if a point is inside the surface
  /**
   * Check if a point is inside the surface. This function calls the signedDistance function of the exact signed
   * distance function and checks if the distance to the point is negative.
   *
   * @param point The point to check
   * @return Whether the point is inside the surface
   */
  [[nodiscard]] bool inside(const VectorF<kDimension> &point) const noexcept {
    return surface_.signedDistance(point) < 0;
  }

  /// Raycast against the surface
  /**
   * Raycast against the surface. This function uses a simple raymarching algorithm to find the intersection of a ray
   * with the surface.
   *
   * @param start The start point of the ray
   * @param direction The direction of the ray
   * @param max_distance The maximum distance to cast
   * @return The result of the raycast
   */
  [[nodiscard]] std::optional<VectorF<kDimension>>
  raycast(const VectorF<kDimension> &start, const VectorF<kDimension> &direction, FLOAT_T max_distance) const noexcept {
    FLOAT_T total_distance = 0;
    while (total_distance < max_distance) {
      VectorF<kDimension> point = start + total_distance * direction;
      FLOAT_T distance = surface_.distance(point);
      if (distance < kDt)
        return point;
      total_distance += distance;
    }
    return std::nullopt;
  }

  using Base::raycast;
};

/// A surface adapter for a Euclidean distance transform
/**
 * A surface adapter for a Euclidean distance transform. This adapter uses a PLC determined from a Euclidean distance
 * transform to provide the functionality of a surface adapter. In order to determine the PLC, marching hypercubes is
 * used.
 * The adapter uses a simple raymarching algorithm to implement the raycast function.
 * It operates in a fixed dimension of 4.
 *
 * @tparam T The Euclidean distance transform to adapt
 */
template <EuclideanDistanceTransform<4> T> class EDTSurfaceAdapter : public LineRaycastMixin<EDTSurfaceAdapter<T>, 4> {
  using Base = LineRaycastMixin<EDTSurfaceAdapter<T>, 4>;

public:
  using Surface = T;                               ///< The Euclidean distance transform
  static inline constexpr unsigned kDimension = 4; ///< The dimension of the space

  constexpr static inline FLOAT_T kDt = static_cast<FLOAT_T>(1e-3); ///< The distance threshold for the raycast

private:
  std::shared_ptr<T> surface_;

public:
  /// Construct a new EDTSurfaceAdapter
  /**
   * Construct a new EDTSurfaceAdapter. This constructor takes a shared pointer to the Euclidean distance transform to
   * adapt.
   *
   * @param surface The Euclidean distance transform to adapt
   */
  explicit EDTSurfaceAdapter(std::shared_ptr<T> surface) : surface_(std::move(surface)) {}

  /// Find the closest point on the surface to a given point
  /**
   * Find the closest point on the surface to a given point. This function first finds the closest point on the surface
   * using the closestPoint function of the Euclidean distance transform. It then uses a raycast in the same direction
   * to find the closest point on the surface, as determined by the marching hypercubes algorithm.
   *
   * @param point The point to find the closest point to
   * @return The closest point on the surface to the given point
   */
  Vector4F closestPoint(const Vector4F &point, bool *point_inside = nullptr) const {
    const auto closest = surface_->closestAt(point);
    const auto line = Eigen::ParametrizedLine<FLOAT_T, 4>::Through(point, closest);
    // only modify the position if the ray hit from a positive direction
    const Vector4F min_corner =
        (line.direction().array() > FLOAT_T()).select(-Vector4F::Ones(), Vector4F::Zero()) + closest;
    std::array<bool, 1U << 4U> corner_values{};
    bool start_sign = surface_->signedDistanceAt(point) > 0;
    std::vector<Vector4F> tests{point - min_corner};
    for (size_t i = 0; i < 1U << 4U; i++) {
      Vector4F offset = Vector4F::Zero();
      for (Eigen::Index j = 0; j < 4; j++)
        if ((i & (1U << static_cast<unsigned>(j))) != 0U)
          offset(3 - j) = 1;
      // true corresponds to a binary 1 in the index, binary 1 in the index means that corner is inside
      // therefore check if the distance is positive, as that means inside
      Vector4F corner = min_corner + offset;
      corner_values.at(i) = surface_->signedDistanceAt(corner) > 0;
      if (corner_values.at(i) == start_sign)
        tests.emplace_back(offset);
    }
    bool *inside_ptr = nullptr;
    if (point_inside) {
      if ((point.array() >= min_corner.array()).all() && (point.array() <= min_corner.array() + 1.0).all())
        inside_ptr = point_inside;
      else
        *point_inside = surface_->signedDistanceAt(point) > 0;
    }
    for (const auto &test : tests) {
      if (auto intersection = surfaceRayIntersection(test, closest - min_corner, corner_values, inside_ptr);
          intersection)
        return *intersection + min_corner;
    }
    throw std::runtime_error("No intersection found");
  }

  /// Check if a sphere intersects the surface
  /**
   * Check if a sphere intersects the surface. This function uses the closestPoint function to find the closest point on
   * the surface to the center of the sphere and checks if the distance to the center of the sphere is less than the
   * radius of the sphere.
   *
   * @param sphere The sphere to check
   * @return Whether the sphere intersects the surface
   */
  [[nodiscard]] bool intersectedBySphere(const HyperSphere4 &sphere) const {
    return (closestPoint(sphere.center()) - sphere.center()).norm() < sphere.radius();
  }

  /// Get the bounding box of the surface
  /**
   * Get the bounding box of the surface. This function calls the boundingBox function of the Euclidean distance
   * transform.
   *
   * @return The bounding box of the surface
   */
  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, 4> boundingBox() const noexcept { return surface_->boundingBox(); }

  /// Check if a point is inside the surface
  /**
   * This function makes use of the same approach as the closestPoint function to determine if a point is inside the
   * surface.
   *
   * @param point The point to check
   * @return Whether the point is inside the surface
   */
  [[nodiscard]] bool inside(const Vector4F &point) const {
    bool point_inside{};
    closestPoint(point, &point_inside);
    return point_inside;
  }

  /// Raycast against the surface
  /**
   * Raycast against the surface. This function uses a simple raymarching algorithm to find the intersection of a ray
   * with the surface. To determine the distance, the fucntion uses the closestPoint function to find the closest point
   * on the surface to the point on the ray.
   *
   * @param start The start point of the ray
   * @param direction The direction of the ray
   * @param max_distance The maximum distance to cast
   * @return The result of the raycast
   */
  [[nodiscard]] std::optional<Vector4F> raycast(const Vector4F &start, const Vector4F &direction,
                                                FLOAT_T max_distance) const {
    FLOAT_T total_distance = 0;
    while (total_distance < max_distance) {
      VectorF<kDimension> point = start + total_distance * direction;
      bool point_inside{};
      FLOAT_T distance = (point - closestPoint(point, &point_inside)).norm();
      if (!point_inside)
        distance = -distance;
      if (distance < kDt)
        return point;
      total_distance += distance;
    }
    return std::nullopt;
  }

  using Base::raycast;
};
} // namespace stmesh
#endif
