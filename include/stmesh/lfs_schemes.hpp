#ifndef STMESH_LFS_SCHEMES_HPP
#define STMESH_LFS_SCHEMES_HPP

#include <concepts> // IWYU pragma: keep
#include <memory>

#include "edt.hpp"
#include "utility.hpp"

namespace stmesh::lfs_schemes {
/// A concept for a local feature size approximation
/**
 * A concept for a local feature size approximation. This is a functional for the approximation with a maximum possible
 * value accessible using the max() method.
 *
 * @tparam T The type to check
 */
template <typename T>
concept LFSScheme = requires(T t, const Vector4F &vec) {
  { t(vec) } -> std::convertible_to<FLOAT_T>;
  { t.max() } -> std::convertible_to<FLOAT_T>;
};

/// A constant local feature size approximation
/**
 * A constant local feature size approximation. This approximation returns the same value for all points.
 */
class Constant {
  FLOAT_T value_;

public:
  /// Construct a constant local feature size approximation
  /**
   * Construct a constant local feature size approximation with the given value.
   *
   * @param value The value to return for all points
   */
  explicit Constant(const FLOAT_T value);

  /// Return the constant value for all points
  /**
   * Return the constant value for all points.
   *
   * @param vec The point to return the value for
   * @return The constant value
   */
  [[nodiscard]] FLOAT_T operator()(const stmesh::Vector4F & /* unused */) const noexcept;

  /// Return the constant value for all points
  /**
   * Return the constant value for all points.
   *
   * @return The constant value
   */
  [[nodiscard]] FLOAT_T max() const noexcept;
};

static_assert(LFSScheme<Constant>);

/// A local feature size approximation based on the distance to the thinned image
/**
 * A local feature size approximation based on the distance to the thinned image. This approximation returns the
 * distance to the thinned image multiplied by a constant value.
 *
 * @note This approximation requires an EDT reader to be passed in the constructor. It is thus not compatible with SDFs.
 */
class BinaryImageApproximation {
  FLOAT_T value_;
  std::shared_ptr<stmesh::EDTReader<4, true>> edt_reader_;

public:
  /// Construct a binary image approximation with the given value and EDT reader
  /**
   * Construct a binary image approximation with the given value and EDT reader.
   *
   * @param value The value to multiply the distance to the thinned image by
   * @param edt_reader The EDT reader to use for the distance to the thinned image
   */
  explicit BinaryImageApproximation(const FLOAT_T value, const std::shared_ptr<stmesh::EDTReader<4, true>> &edt_reader);

  /// Return the distance to the thinned image multiplied by the constant value
  /**
   * Return the distance to the thinned image multiplied by the constant value.
   *
   * @param vec The point to return the value for
   * @return The distance to the thinned image multiplied by the constant value
   */
  [[nodiscard]] FLOAT_T operator()(const stmesh::Vector4F &vec) const noexcept;

  /// Return the maximum value of the approximation
  /**
   * Return the maximum value of the approximation. This is the value multiplied by the maximum possible distance to the
   * thinned image.
   *
   * @return The maximum value of the approximation
   */
  [[nodiscard]] FLOAT_T max() const noexcept;
};

static_assert(LFSScheme<BinaryImageApproximation>);
} // namespace stmesh::lfs_schemes

#endif