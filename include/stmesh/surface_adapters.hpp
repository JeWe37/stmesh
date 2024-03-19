#ifndef STMESH_SURFACE_ADAPTERS_HPP
#define STMESH_SURFACE_ADAPTERS_HPP

#include <boost/math/tools/toms748_solve.hpp>
#include <concepts> // IWYU pragma: keep
#include <cstddef>
#include <memory>
#include <optional>
#include <utility>

#include <Eigen/Geometry>

#include "edt.hpp" // IWYU pragma: keep
#include "marching_hypercubes.hpp"
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

template <EuclideanDistanceTransform<4> T> class EDTSurfaceAdapter : public LineRaycastMixin<EDTSurfaceAdapter<T>, 4> {
  using Base = LineRaycastMixin<EDTSurfaceAdapter<T>, 4>;

public:
  using surface_t = T;
  static inline constexpr unsigned kDimension = 4;

  constexpr static inline FLOAT_T kDt = static_cast<FLOAT_T>(1e-3);

private:
  std::shared_ptr<T> surface_;

public:
  template <typename... Args> explicit EDTSurfaceAdapter(std::shared_ptr<T> surface) : surface_(std::move(surface)) {}

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

  [[nodiscard]] bool intersectedBySphere(const HyperSphere4 &sphere) const {
    return (closestPoint(sphere.center()) - sphere.center()).norm() < sphere.radius();
  }

  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, 4> boundingBox() const noexcept { return surface_->boundingBox(); }

  [[nodiscard]] bool inside(const Vector4F &point) const {
    bool point_inside{};
    closestPoint(point, &point_inside);
    return point_inside;
  }

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
