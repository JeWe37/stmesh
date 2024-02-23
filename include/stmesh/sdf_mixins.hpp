#ifndef STMESH_SDF_MIXINS_HPP
#define STMESH_SDF_MIXINS_HPP

#include <cstdlib>
#include <limits>
#include <random>

#include <Eigen/Geometry>

#include "utility.hpp"

namespace stmesh {
/// A mixin for calculating the normal of a signed distance function
/**
 * A mixin for calculating the normal of a signed distance function. The normal is calculated using central differences.
 *
 * @tparam T The type of the class that the mixin is mixed into
 * @tparam D The dimension of the space
 */
template <typename T, unsigned D> class CentralDifferenceNormalMixin {
  static inline constexpr FLOAT_T kEps = sqrt(std::numeric_limits<FLOAT_T>::epsilon());

public:
  /// Calculate the normal of the signed distance function at a point
  /**
   * Calculate the normal of the signed distance function at a point. The normal is calculated using central
   * differences.
   *
   * @param point The point at which to calculate the normal
   * @return The normal of the signed distance function at the point
   */
  [[nodiscard]] VectorF<D> normal(const VectorF<D> &point) const noexcept {
    const T &self = static_cast<const T &>(*this);
    VectorF<D> normal;
    VectorF<D> delta = VectorF<D>::Zero();
    for (int i = 0; i < static_cast<int>(D); ++i) {
      delta[i] = kEps;
      normal[i] = self.signedDistance(point + delta) - self.signedDistance(point - delta);
      delta[i] = 0;
    }
    return normal.normalized();
  }
};

/// A mixin for adding the distance function to a signed distance function
/**
 * A mixin for adding the distance function to a signed distance function. The distance function is calculated as the
 * absolute value of the signed distance function in the subclass.
 *
 * @tparam T The type of the class that the mixin is mixed into
 * @tparam D The dimension of the space
 */
template <typename T, unsigned D> class DistanceMixin {
public:
  /// Calculate the distance of a point from the signed distance function
  /**
   * Calculate the distance of a point from the signed distance function. The distance is calculated as the absolute
   * value of the signed distance function in the subclass.
   *
   * @param point The point at which to calculate the distance
   * @return The distance of the point from the signed distance function
   */
  [[nodiscard]] FLOAT_T distance(const VectorF<D> &point) const noexcept {
    return std::abs(static_cast<const T *>(this)->signedDistance(point));
  }
};

/// A mixin for sampling a signed distance function
/**
 * A mixin for sampling a signed distance function. The sample function returns a random point that is inside the signed
 * distance function.
 *
 * @tparam T The type of the class that the mixin is mixed into
 * @tparam D The dimension of the space
 */
template <typename T, unsigned D> class SamplingMixin {
public:
  /// Sample a point that is inside the signed distance function
  /**
   * Sample a point that is inside the signed distance function. The point is sampled using rejection sampling in the
   * bounding box.
   *
   * @param gen The random number generator to use to sample the point
   * @return A point that is inside the signed distance function
   */
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