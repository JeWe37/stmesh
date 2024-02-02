#ifndef STMESH_SURFACE_ADAPTERS_HPP
#define STMESH_SURFACE_ADAPTERS_HPP

#include <concepts> // IWYU pragma: keep
#include <optional>
#include <utility>

#include <Eigen/Geometry>

#include "sdf.hpp"
#include "utility.hpp"

namespace stmesh {
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

template <typename T>
concept SurfaceAdapter4 = SurfaceAdapter<T, 4>;

template <typename T, unsigned D> class LineRaycastMixin {
public:
  [[nodiscard]] std::optional<VectorF<D>> raycast(const Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(D)> &ray,
                                                  FLOAT_T max_distance) const noexcept {
    return static_cast<const T *>(this)->raycast(ray.origin(), ray.direction(), max_distance);
  }
};

template <ExactSDF T> class SDFSurfaceAdapter : public LineRaycastMixin<SDFSurfaceAdapter<T>, T::kDimension> {
  using Base = LineRaycastMixin<SDFSurfaceAdapter<T>, T::kDimension>;

public:
  using Surface = T;
  static inline constexpr unsigned kDimension = T::kDimension;
  // TODO: make this variable
  constexpr static inline FLOAT_T kDt = static_cast<FLOAT_T>(1e-3);

private:
  T surface_;

public:
  template <typename... Args> explicit SDFSurfaceAdapter(Args &&...args) : surface_(std::forward<Args>(args)...) {}

  [[nodiscard]] VectorF<kDimension> closestPoint(const VectorF<kDimension> &point) const noexcept {
    return point - surface_.signedDistance(point) * surface_.normal(point);
  }

  [[nodiscard]] bool intersectedBySphere(const HyperSphere<kDimension> &sphere) const noexcept {
    return sphere.radius() - surface_.distance(sphere.center()) > kDt;
  }

  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(kDimension)> boundingBox() const noexcept {
    return surface_.boundingBox();
  }

  [[nodiscard]] bool inside(const VectorF<kDimension> &point) const noexcept {
    return surface_.signedDistance(point) < 0;
  }

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

} // namespace stmesh
#endif