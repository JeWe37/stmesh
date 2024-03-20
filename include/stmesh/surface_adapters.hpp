#ifndef STMESH_SURFACE_ADAPTERS_HPP
#define STMESH_SURFACE_ADAPTERS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts> // IWYU pragma: keep
#include <cstddef>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include <Eigen/Geometry>
#include <spdlog/spdlog.h>

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

  [[nodiscard]] static std::array<bool, 1U << 4U> surroundingCubesValues(const std::vector<bool> &insides,
                                                                         const std::span<const size_t, 4> strides,
                                                                         const std::span<const size_t, 4> offsets) {
    std::array<size_t, 4> stride_offsets{};
    stride_offsets[3] = 1;
    std::partial_sum(strides.rbegin(), strides.rbegin() + 3, stride_offsets.rbegin() + 1, std::multiplies<>());
    // offsets[0] = 1;
    // std::partial_sum(strides.begin() + 1, strides.end(), offsets.begin() + 1, std::multiplies<>());
    std::array<bool, 1U << 4U> result{};
    for (size_t i = 0; i < 1U << 4U; i++) {
      size_t index = 0;
      for (size_t j = 0; j < 4; j++)
        index += (offsets[j] + static_cast<size_t>((i & (1U << (3U - j))) != 0U)) * stride_offsets.at(j);
      result.at(i) = insides.at(index);
    }
    return result;
  }

  [[nodiscard]] std::vector<bool> getSigns(const Vector4F &point, const std::span<const size_t, 4> sizes) const {
    Vector4F spacing = surface_->spacing();
    Vector4F min_corner = point.cwiseProduct(spacing) - FLOAT_T(0.5) * spacing;
    std::vector<FLOAT_T> distances = surface_->signedDistanceAt(min_corner, sizes);
    std::vector<bool> signs;
    std::ranges::transform(distances, std::back_inserter(signs), [](FLOAT_T distance) { return distance > 0; });
    return signs;
  }

  template <typename F> bool allSteps(const F &f, const Vector4F &first, const Vector4F &second) const {
    const Array4F delta = (second - first).array().abs();
    Eigen::Index to_advance = delta.count();
    Eigen::Vector<size_t, 4> size = delta.cast<size_t>() + size_t{2};
    std::span<size_t, 4> sizes(size.data(), 4); // strides
    Vector4F smaller = first.cwiseMin(second);
    Vector4F larger = first.cwiseMax(second);
    std::vector<bool> signs = getSigns(smaller, sizes);
    for (size_t i = 0; i < 1U << static_cast<unsigned>(to_advance); ++i) {
      size_t j = 0;
      Vector4F combined;
      std::array<size_t, 4> offsets{};
      for (Eigen::Index k = 0; k < first.rows(); ++k) {
        if (delta[k] == FLOAT_T(1.0) && (i & (1U << j++))) {
          combined[k] = larger[k];
          offsets.at(static_cast<size_t>(k)) = 1;
        } else
          combined[k] = smaller[k];
      }
      if (combined == first)
        continue;
      std::array<bool, 1U << 4U> corner_values = surroundingCubesValues(signs, sizes, offsets);
      if (f(combined, corner_values))
        return true;
    }
    return false;
  }

  template <typename F> void iterateLine(const F &f, const Vector4F &start, const Vector4F &end) const {
    Vector4F direction = end.array() - start.array();
    Eigen::Index largest_index{};
    FLOAT_T largest = direction.array().abs().maxCoeff(&largest_index);
    FLOAT_T first_step =
        direction[largest_index] > FLOAT_T() ? std::ceil(largest) - largest : largest - std::floor(largest);
    Vector4F curr_pos = start.array().floor();
    constexpr std::array<size_t, 4> kSizes{2, 2, 2, 2};
    std::vector<bool> signs = getSigns(curr_pos, kSizes);
    std::array<bool, 1U << 4U> corner_values{};
    std::ranges::copy(signs, corner_values.begin());
    if (f(curr_pos, corner_values))
      return;
    // NOLINTNEXTLINE(cert-flp30-c)
    for (FLOAT_T acc = first_step / largest; acc <= FLOAT_T(1.0 + 1.0 / largest); acc += FLOAT_T(1.0) / largest) {
      Vector4F next_pos = (start + direction * acc).array().floor();
      next_pos[largest_index] = curr_pos[largest_index] + sgn(direction[largest_index]);
      if (allSteps(f, curr_pos, next_pos))
        return;
      curr_pos = next_pos;
    }
  }

  template <typename F> void spacingWrap(const F &f, const Vector4F &start, const Vector4F &end) const {
    const Vector4F spacing = surface_->spacing();
    Vector4F half = Vector4F::Constant(FLOAT_T(0.5));
    iterateLine(
        [&](const Vector4F &point, const std::array<bool, 1U << 4U> &corner_values) {
          return f(point.cwiseProduct(spacing), corner_values);
        },
        start.cwiseQuotient(spacing) + half, end.cwiseQuotient(spacing) + half);
  }

  std::optional<Vector4F> initialTest(const Vector4F &point, bool *point_inside) const {
    // logic: according to the marching cubes cell(with minimal corner min_corner) surrounding the point, the point must
    // lie on the same side of the surface as the next closest point(rounded_point) on the grid, as otherwise the
    // closest point(closest) and it would lie on the same side, meaning an intersection might not exist
    constexpr std::array<size_t, 4> kSizes{2, 2, 2, 2};
    Vector4F spacing = surface_->spacing();
    Array4F min_point = (point.array() / spacing.array() - FLOAT_T(0.5)).floor() + FLOAT_T(0.5);
    Vector4F min_corner = min_point * spacing.array();
    Vector4F rounded_point = ((point.array() / spacing.array()).floor() + FLOAT_T(0.5)) * spacing.array();
    std::vector<bool> signs = getSigns(min_point + FLOAT_T(0.5), kSizes);
    std::array<bool, 1U << 4U> corner_values{};
    std::ranges::copy(signs, corner_values.begin());
    std::optional<Vector4F> initial_result =
        surfaceRayIntersection((point - min_corner).cwiseQuotient(spacing),
                               (rounded_point - min_corner).cwiseQuotient(spacing), corner_values, point_inside);
    if (*point_inside != (surface_->signedDistanceAt(rounded_point) > 0))
      // should this not be the case, just find the intersection of the ray from the point to the rounded_point with the
      // surface, this should be close enough and should definitely exist
      return initial_result->cwiseProduct(spacing) + min_corner;
    return std::nullopt;
  }

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
  [[nodiscard]] Vector4F closestPoint(const Vector4F &point) const {
    Vector4F closest = surface_->closestAt(point);
    const Vector4F direction = closest - point;
    constexpr FLOAT_T kSkipMin = 4.0;
    FLOAT_T skip = std::max(FLOAT_T(0.0), direction.norm() - surface_->spacing().minCoeff() * kSkipMin);
    if (skip < FLOAT_T(0.5) * surface_->spacing().maxCoeff()) {
      bool initial_point_inside{};
      if (std::optional initial = initialTest(point, &initial_point_inside); initial)
        return *initial;
    }
    const Vector4F normalized_direction = direction.normalized();
    if (std::optional cast =
            raycast(point + normalized_direction * skip, normalized_direction,
                    direction.norm() + std::min(2.0 * surface_->spacing().maxCoeff() * kSkipMin - skip, FLOAT_T(0.0)));
        cast)
      return *cast;
    spdlog::warn("No intersection found for start point ({}), returning closest point", point.transpose());
    // kinda a hack, but extremely unlikely(impossible?) to happen
    return closest;
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
    initialTest(point, &point_inside);
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
    Vector4F spacing = surface_->spacing();
    std::optional<Vector4F> result;
    Vector4F end = start + max_distance * direction.normalized();
    spacingWrap(
        [&](const Vector4F &point, const std::array<bool, 1U << 4U> &corner_values) {
          Vector4F min_corner = point - spacing * 0.5;
          if (std::optional intersection =
                  surfaceRayIntersection((start - min_corner).cwiseQuotient(spacing),
                                         (end - min_corner).cwiseQuotient(spacing), corner_values);
              intersection) {
            result = *intersection + min_corner;
            return true;
          }
          return false;
        },
        start, end);
    return result;
  }

  using Base::raycast;
};
} // namespace stmesh
#endif
