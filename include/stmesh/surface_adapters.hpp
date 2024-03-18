#ifndef STMESH_SURFACE_ADAPTERS_HPP
#define STMESH_SURFACE_ADAPTERS_HPP

#include <boost/math/tools/toms748_solve.hpp>
#include <concepts> // IWYU pragma: keep
#include <optional>
#include <utility>

#include <Eigen/Geometry>
#include <boost/math/tools/roots.hpp>

#include "edt.hpp" // IWYU pragma: keep
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

template <unsigned D, EuclideanDistanceTransform<D> T>
class EDTSurfaceAdapter : public LineRaycastMixin<EDTSurfaceAdapter<D, T>, D> {
  using Base = LineRaycastMixin<EDTSurfaceAdapter<D, T>, D>;

public:
  using surface_t = T;
  static inline constexpr unsigned kDimension = D;

  constexpr static inline FLOAT_T kDt = static_cast<FLOAT_T>(1e-3);

private:
  T surface_;

  // NOLINTBEGIN(*-magic-numbers,*-member-init)
  [[nodiscard]] static FLOAT_T quadrilinearBinaryInterpolate(const std::array<bool, 16> &corner_values,
                                                             const Vector4F &pos) {
    std::array<FLOAT_T, 16> corner_weights;
    std::ranges::transform(corner_values, corner_weights.begin(), [](bool value) -> FLOAT_T { return value ? 1 : -1; });
    std::array<FLOAT_T, 8> edge_weights;
    for (size_t i = 0; i < 8; i++)
      edge_weights.at(i) = corner_weights.at(i + 8) * pos(3) + corner_weights.at(i) * (FLOAT_T(1) - pos(3));
    std::array<FLOAT_T, 4> face_weights;
    for (size_t i = 0; i < 4; i++)
      face_weights.at(i) = edge_weights.at(i + 4) * pos(2) + edge_weights.at(i) * (FLOAT_T(1) - pos(2));
    std::array<FLOAT_T, 2> volume_weights;
    for (size_t i = 0; i < 2; i++)
      volume_weights.at(i) = face_weights.at(i + 2) * pos(1) + face_weights.at(i) * (FLOAT_T(1) - pos(1));
    return volume_weights.at(1) * pos(0) + volume_weights.at(0) * (FLOAT_T(1) - pos(0));
  }
  // NOLINTEND(*-magic-numbers,*-member-init)

public:
  template <typename... Args> explicit EDTSurfaceAdapter(Args &&...args) : surface_(std::forward<Args>(args)...) {}

  VectorF<D> closestPoint(const VectorF<D> &point, bool *inside = nullptr) const {
    const auto closest = surface_.closestAt(point);
    const auto line = Eigen::ParametrizedLine<FLOAT_T, 4>::Through(point, closest);
    // only modify the position if the ray hit from a positive direction
    const Vector4F min_corner =
        (line.direction().array() > FLOAT_T()).select(-Vector4F::Ones(), Vector4F::Zero()) + closest;
    // NOLINTNEXTLINE(*-member-init)
    std::array<bool, 1U << 4U> corner_values;
    for (size_t i = 0; i < 1U << 4U; i++) {
      Vector4F offset = Vector4F::Zero();
      for (Eigen::Index j = 0; j < 4; j++)
        if ((i & (1U << static_cast<unsigned>(j))) != 0U)
          offset(j) = 1;
      corner_values.at(i) = surface_.signedDistanceAt(min_corner + offset) > 0;
    }
    const FLOAT_T initial_distance = (point - closest).norm();
    const auto point_side = [&](FLOAT_T x) -> FLOAT_T {
      return quadrilinearBinaryInterpolate(corner_values, line.pointAt(x) - min_corner);
    };
    // NOLINTNEXTLINE(*-magic-numbers)
    auto max_iter = static_cast<boost::uintmax_t>(100);
    const auto [lower, upper] =
        boost::math::tools::toms748_solve(point_side, initial_distance - FLOAT_T(2.0), initial_distance + FLOAT_T(2.0),
                                          boost::math::tools::eps_tolerance<FLOAT_T>(), max_iter);
    // NOLINTNEXTLINE(*-magic-numbers)
    if (max_iter > static_cast<boost::uintmax_t>(100))
      throw std::runtime_error("EDTSurfaceAdapter::closestPoint: max_iter exceeded");
    if (inside != nullptr)
      *inside = point_side(0) > 0;
    return line.pointAt((lower + upper) / 2);
  }

  [[nodiscard]] bool intersectedBySphere(const HyperSphere<D> &sphere) const {
    return (closestPoint(sphere.center()) - sphere.center()).norm() < sphere.radius();
  }

  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> boundingBox() const noexcept {
    return surface_.boundingBox();
  }

  [[nodiscard]] bool inside(const VectorF<D> &point) const {
    // NOLINTNEXTLINE(cppcoreguidelines-init-variables)
    bool point_inside;
    closestPoint(point, &point_inside);
    return point_inside;
  }

  [[nodiscard]] std::optional<VectorF<D>> raycast(const VectorF<D> &start, const VectorF<D> &direction,
                                                  FLOAT_T max_distance) const {
    FLOAT_T total_distance = 0;
    while (total_distance < max_distance) {
      VectorF<kDimension> point = start + total_distance * direction;
      // NOLINTNEXTLINE(cppcoreguidelines-init-variables)
      bool point_inside;
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