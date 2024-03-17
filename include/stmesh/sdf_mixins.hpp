#ifndef STMESH_SDF_MIXINS_HPP
#define STMESH_SDF_MIXINS_HPP

#include <cstdlib>
#include <limits>
#include <random>

#include <Eigen/Geometry>

#include "utility.hpp"

namespace stmesh {
template <typename T, unsigned D> class CentralDifferenceNormalMixin {
  static inline constexpr FLOAT_T eps = sqrt(std::numeric_limits<FLOAT_T>::epsilon());

public:
  [[nodiscard]] VectorF<D> normal(const VectorF<D> &point) const noexcept {
    const T &self = static_cast<const T &>(*this);
    VectorF<D> normal;
    VectorF<D> delta = VectorF<D>::Zero();
    for (int i = 0; i < static_cast<int>(D); ++i) {
      delta[i] = eps;
      normal[i] = self.signedDistance(point + delta) - self.signedDistance(point - delta);
      delta[i] = 0;
    }
    return normal.normalized();
  }
};

template <typename T, unsigned D> class DistanceMixin {
public:
  [[nodiscard]] FLOAT_T distance(const VectorF<D> &point) const noexcept {
    return std::abs(static_cast<const T *>(this)->signedDistance(point));
  }
};

template <typename T, unsigned D> class SamplingMixin {
public:
  [[nodiscard]] VectorF<D> sample(std::uniform_random_bit_generator auto &gen) const noexcept {
    std::uniform_real_distribution<FLOAT_T> dist(0, 1);
    VectorF<D> unit_point;
    while (true) {
      for (int i = 0; i < static_cast<int>(D); ++i)
        unit_point[i] = dist(gen);
      Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> box = static_cast<const T *>(this)->boundingBox();
      VectorF<D> point = box.min() + unit_point.cwiseProduct(box.sizes());
      if (static_cast<const T *>(this)->signedDistance(point) < 0)
        return point;
    }
  }
};
} // namespace stmesh
#endif