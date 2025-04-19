#ifndef STMESH_RADIUS_SCHEMES_HPP
#define STMESH_RADIUS_SCHEMES_HPP

#include <filesystem>
#include <memory>
#include <type_traits>

#include "stmesh/edt.hpp"
#include "utility.hpp"

namespace stmesh::radius_schemes {
/// A concept for a radius scheme
/**
 * A concept for a radius scheme. This is a functional for the radius with a single argument of type Vector4F.
 *
 * @tparam F The type to check
 */
template <typename F>
concept RadiusScheme = std::is_invocable_r_v<FLOAT_T, F, Vector4F>;

/// A constant radius scheme
/**
 * A constant radius scheme. This scheme returns the same value for all points.
 */
class Constant {
  FLOAT_T value_;

public:
  /// Construct a constant radius scheme
  /**
   * Construct a constant radius scheme with the given value.
   *
   * @param value The value to return for all points
   */
  explicit Constant(const FLOAT_T value);

  /// Return the constant value for all points
  /**
   * Return the constant value for all points. Ignores the input point.
   *
   * @param vec The point to return the value for
   * @return The constant value
   */
  [[nodiscard]] FLOAT_T operator()(const Vector4F & /* unused */) const noexcept;
};

static_assert(RadiusScheme<Constant>);

namespace detail {
template <typename M, typename T>
requires std::is_invocable_r_v<FLOAT_T, M, FLOAT_T> && SignedDistance<T, 4>
class DistanceRadius {
protected:
  M mapper_;
  std::shared_ptr<T> signed_distance_;

public:
  explicit DistanceRadius(const std::shared_ptr<T> &signed_distance, M mapper = {})
      : mapper_(mapper), signed_distance_(signed_distance) {}
};
} // namespace detail

using Linear_t = decltype([](const FLOAT_T dist) { return dist; });

/// A boundary distance radius scheme
/**
 * A boundary distance radius scheme. This scheme returns the signed distance to the boundary for all points.
 * You can provide a custom mapper function to transform the distance.
 *
 * @tparam T The type from which to read the signed distance
 * @tparam M The type of the mapper function
 */
template <typename T, typename M = Linear_t> class BoundaryDistanceRadius : detail::DistanceRadius<M, T> {
  using Super = detail::DistanceRadius<M, T>;
  using Super::mapper_;
  using Super::signed_distance_;

public:
  /// Construct a boundary distance radius scheme
  /**
   * Construct a boundary distance radius scheme with the given EDT reader and mapper function.
   *
   * @param edt_reader The EDT reader to use for the distance
   * @param mapper The mapper function to transform the distance
   */
  using Super::Super;

  /// Return the signed distance to the boundary for the given point
  /**
   * Return the signed distance to the boundary for the given point. The distance is transformed by the mapper function.
   *
   * @param vec The point to return the distance for
   * @return The signed distance to the boundary
   */
  [[nodiscard]] FLOAT_T operator()(const Vector4F &vec) const noexcept {
    return mapper_(signed_distance_->signedDistance(vec));
  }
};

template <typename M, typename T> BoundaryDistanceRadius(const std::shared_ptr<T> &, M) -> BoundaryDistanceRadius<T, M>;

template <typename T> BoundaryDistanceRadius(const std::shared_ptr<T> &) -> BoundaryDistanceRadius<T>;

static_assert(RadiusScheme<BoundaryDistanceRadius<stmesh::EDTReader<4, false>>>);
static_assert(RadiusScheme<BoundaryDistanceRadius<stmesh::EDTReader<4, true>>>);

/// A local feature size radius scheme
/**
 * A local feature size radius scheme. This scheme returns the distance to the thinned image for all points.
 * You can provide a custom mapper function to transform the distance.
 *
 * @tparam M The type of the mapper function
 */
template <typename M = Linear_t> class LFSRadius : detail::DistanceRadius<M, stmesh::EDTReader<4, true>> {
  using Super = detail::DistanceRadius<M, stmesh::EDTReader<4, true>>;
  using Super::mapper_;
  using Super::signed_distance_;

public:
  /// Construct a local feature size radius scheme
  /**
   * Construct a local feature size radius scheme with the given EDT reader and mapper function.
   *
   * @param edt_reader The EDT reader to use for the distance
   * @param mapper The mapper function to transform the distance
   */
  using Super::Super;

  /// Return the distance to the thinned image for the given point
  /**
   * Return the distance to the thinned image for the given point. The distance is transformed by the mapper function.
   *
   * @param vec The point to return the distance for
   * @return The distance to the thinned image
   */
  [[nodiscard]] FLOAT_T operator()(const Vector4F &vec) const noexcept {
    return mapper_(signed_distance_->distanceToThinnedAt(vec));
  }
};

template <typename M> LFSRadius(const std::shared_ptr<stmesh::EDTReader<4, true>> &, M) -> LFSRadius<M>;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
LFSRadius(const std::shared_ptr<stmesh::EDTReader<4, true>> &) -> LFSRadius<>;
#pragma GCC diagnostic pop

static_assert(RadiusScheme<LFSRadius<>>);

/// An image radius scheme
/**
 * An image radius scheme. This scheme returns the value of the image at the given point.
 */
class ImageRadius {
  struct Impl;
  std::unique_ptr<Impl> pimpl_;

public:
  /// Construct an image radius scheme
  /**
   * Construct an image radius scheme with the given image path. The image must be in a format that can be read by ITK.
   *
   * @param path The path to the image to use
   */
  explicit ImageRadius(const std::filesystem::path &path);

  ~ImageRadius();
  ImageRadius(const ImageRadius &) = delete;
  ImageRadius &operator=(const ImageRadius &) = delete;
  ImageRadius(ImageRadius &&) noexcept;
  ImageRadius &operator=(ImageRadius &&) noexcept;

  /// Return the value of the image at the given point
  /**
   * Return the value of the image at the given point. The image is assumed to have an origin of 0.5 in all dimensions.
   *
   * @param vec The point to return the value for
   * @return The value of the image at the point
   */
  [[nodiscard]] FLOAT_T operator()(const Vector4F &vec) const noexcept;
};

static_assert(RadiusScheme<ImageRadius>);
} // namespace stmesh::radius_schemes

#endif
