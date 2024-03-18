// NOLINTBEGIN(misc-include-cleaner)
#include <stmesh/sdf.hpp>

#include <algorithm>

namespace stmesh {
template <unsigned D>
HyperSphere<D>::HyperSphere(FLOAT_T radius, const VectorF<D> &center) noexcept : radius_(radius), center_(center) {}

template <unsigned D> [[nodiscard]] FLOAT_T HyperSphere<D>::radius() const noexcept { return radius_; }

template <unsigned D> void HyperSphere<D>::scale(FLOAT_T factor) noexcept { radius_ *= factor; }

template <unsigned D> [[nodiscard]] const VectorF<D> &HyperSphere<D>::center() const noexcept { return center_; }

template <unsigned D> [[nodiscard]] FLOAT_T HyperSphere<D>::signedDistance(const VectorF<D> &point) const noexcept {
  return (point - center_).norm() - radius_;
}

template <unsigned D> [[nodiscard]] VectorF<D> HyperSphere<D>::normal(const VectorF<D> &point) const noexcept {
  return (point - center_).normalized();
}

template <unsigned D>
[[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> HyperSphere<D>::boundingBox() const noexcept {
  return {center_ - radius_ * VectorF<D>::Ones(), center_ + radius_ * VectorF<D>::Ones()};
}

template <unsigned D> bool HyperSphere<D>::operator==(const HyperSphere<D> &other) const noexcept {
  return radius_ == other.radius_ && center_ == other.center_;
}

template class HyperSphere<3>;
template class HyperSphere<4>;

template <unsigned D>
HyperCube<D>::HyperCube(const VectorF<D> &min, const VectorF<D> &max) noexcept : min_(min), max_(max) {}

// https://iquilezles.org/articles/distfunctions/
template <unsigned D> [[nodiscard]] FLOAT_T HyperCube<D>::signedDistance(VectorF<D> point) const noexcept {
  const VectorF<D> center = (min_ + max_) / 2;
  point -= center;
  const VectorF<D> half_extent = (max_ - min_) / 2;
  const VectorF<D> q = point.cwiseAbs() - half_extent;
  const FLOAT_T max_dist = q.maxCoeff();
  return q.cwiseMax(VectorF<D>::Zero()).norm() + std::min(max_dist, FLOAT_T{0});
}

template <unsigned D>
[[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> HyperCube<D>::boundingBox() const noexcept {
  return {min_, max_};
}

template class HyperCube<4>;
} // namespace stmesh

// NOLINTEND(misc-include-cleaner)