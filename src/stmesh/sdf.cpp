#include "stmesh/sdf.hpp"

#include <algorithm>
#include <cmath>

#include <Eigen/Geometry>

#include "stmesh/utility.hpp"

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
  if (center_.isApprox(point))
    return VectorF<D>::Unit(0);
  return (point - center_).normalized();
}

template <unsigned D>
[[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> HyperSphere<D>::boundingBox() const noexcept {
  return {center_ - radius_ * VectorF<D>::Ones(), center_ + radius_ * VectorF<D>::Ones()};
}

template <unsigned D> bool HyperSphere<D>::operator==(const HyperSphere<D> &other) const noexcept {
  return radius_ == other.radius_ && center_ == other.center_;
}

template class HyperSphere<2>;
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

CylinderSDF::CylinderSDF(FLOAT_T radius, const Eigen::ParametrizedLine<FLOAT_T, 3> &axis) noexcept
    : radius_(radius), axis_(axis) {}

[[nodiscard]] FLOAT_T CylinderSDF::signedDistance(const VectorF<3> &point) const noexcept {
  // https://iquilezles.org/articles/distfunctions/
  const Vector3F ba = axis_.direction();
  const Vector3F pa = point - axis_.origin();
  const FLOAT_T baba = ba.squaredNorm();
  const FLOAT_T paba = ba.dot(pa);
  const FLOAT_T x = (pa * baba - ba * paba).norm() - radius_ * baba;
  const FLOAT_T y = std::abs(paba - baba * 0.5) - baba * 0.5;
  const FLOAT_T x2 = x * x;
  const FLOAT_T y2 = y * y * baba;
  const FLOAT_T d = (std::max(x, y) <= 0.0) ? -std::min(x2, y2) : (x > 0.0 ? x2 : 0.0) + (y > 0.0 ? y2 : 0.0);
  return (d > 0.0 ? 1.0 : -1.0) * std::sqrt(std::abs(d)) / baba;
}

[[nodiscard]] Vector3F CylinderSDF::normal(const Vector3F &point) const noexcept {
  // Extract axis components
  const Vector3F ba = axis_.direction();
  const Vector3F pa = point - axis_.origin();
  const FLOAT_T baba = ba.squaredNorm();
  const FLOAT_T paba = ba.dot(pa);

  // Compute perpendicular vector
  const Vector3F v = pa - (paba / baba) * ba;

  // Compute x and y as in SDF
  const FLOAT_T x = std::sqrt(baba) * (v.norm() - radius_);
  const FLOAT_T y = std::abs(paba - baba * 0.5) - baba * 0.5;

  // Determine normal based on dominant component
  if (x >= y)
    return v.normalized(); // Side normal
  else {
    const FLOAT_T sign = (paba > baba * 0.5) ? 1.0 : -1.0;
    return sign * ba.normalized(); // Cap normal
  }
}

[[nodiscard]] Eigen::AlignedBox<FLOAT_T, 3> CylinderSDF::boundingBox() const noexcept {
  const Vector3F end = axis_.pointAt(1.0);
  const Vector3F min = axis_.origin().cwiseMin(end);
  const Vector3F max = axis_.origin().cwiseMax(end);

  const Vector3F normal_dir = axis_.direction().normalized();
  Eigen::AlignedBox<FLOAT_T, 3> result;
  for (int i = 0; i < 3; ++i) {
    const Vector3F dirorth = Vector3F::Unit(i).cross(normal_dir);
    Array3F dirorthorth = dirorth.cross(normal_dir);
    dirorthorth *= dirorthorth.sign();
    result.min()[i] = min[i] - radius_ * dirorthorth[i];
    result.max()[i] = max[i] + radius_ * dirorthorth[i];
  }
  return result;
}
} // namespace stmesh
