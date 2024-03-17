#ifndef STMESH_SDF_HPP
#define STMESH_SDF_HPP

#include <concepts> // IWYU pragma: keep

#include <Eigen/Geometry>

#include "sdf_mixins.hpp"
#include "utility.hpp"

namespace stmesh {
template <typename T, unsigned D>
concept SignedDistanceField = requires(T t, VectorF<D> vec) {
  { t.distance(vec) } -> std::convertible_to<FLOAT_T>;
  { t.signedDistance(vec) } -> std::convertible_to<FLOAT_T>;
  { t.normal(vec) } -> std::convertible_to<VectorF<D>>;
  { t.boundingBox() } -> std::convertible_to<Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)>>;
};

template <typename T>
requires SignedDistanceField<T, T::dimension>
struct ExactSDFTag {
  static constexpr bool value = false;
};

template <typename T>
concept ExactSDF = requires { ExactSDFTag<T>::value; };

template <unsigned D>
class HyperSphere : public stmesh::DistanceMixin<HyperSphere<D>, D>, public SamplingMixin<HyperSphere<D>, D> {
  FLOAT_T radius_;
  VectorF<D> center_;

public:
  static inline constexpr unsigned dimension = D;

  HyperSphere(FLOAT_T radius, const VectorF<D> &center) noexcept;

  [[nodiscard]] FLOAT_T radius() const noexcept;

  void scale(FLOAT_T factor) noexcept;

  [[nodiscard]] const VectorF<D> &center() const noexcept;

  [[nodiscard]] FLOAT_T signedDistance(const VectorF<D> &point) const noexcept;

  [[nodiscard]] VectorF<D> normal(const VectorF<D> &point) const noexcept;

  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> boundingBox() const noexcept;
};

template <unsigned D> struct ExactSDFTag<HyperSphere<D>> {
  static constexpr bool value = true;
};

using HyperSphere3 = HyperSphere<3>;
using HyperSphere4 = HyperSphere<4>;

template <unsigned D>
class HyperCube : public stmesh::DistanceMixin<HyperCube<D>, D>, public CentralDifferenceNormalMixin<HyperCube<D>, D> {
  VectorF<D> min_;
  VectorF<D> max_;

public:
  static inline constexpr unsigned dimension = D;

  HyperCube(const VectorF<D> &min, const VectorF<D> &max) noexcept;

  [[nodiscard]] FLOAT_T signedDistance(VectorF<D> point) const noexcept;

  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> boundingBox() const noexcept;
};

template <unsigned D> struct ExactSDFTag<HyperCube<D>> {
  static constexpr bool value = true;
};

using HyperCube4 = HyperCube<4>;
} // namespace stmesh
#endif