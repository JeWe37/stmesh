#ifndef STMESH_EDT_READER_HPP
#define STMESH_EDT_READER_HPP

#include <array>
#include <concepts> // IWYU pragma: keep
#include <cstddef>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include <Eigen/Geometry>

#include "utility.hpp"

namespace stmesh {
/// A concept for an Euclidean distance transform
/**
 * A concept for an Euclidean distance transform. An Euclidean distance transform is a function that returns the signed
 * distance to the surface of a shape, the closest point on the surface of a shape, and the bounding box of the shape.
 *
 * @tparam T The type to check
 * @tparam D The dimension of the space
 */
template <typename T, unsigned D>
concept EuclideanDistanceTransform = requires(const T t, VectorF<D> vec, std::array<size_t, D> size) {
  { t.signedDistanceAt(vec) } -> std::convertible_to<FLOAT_T>;
  { t.signedDistanceAt(vec, size) } -> std::convertible_to<std::vector<FLOAT_T>>;
  { t.closestAt(vec) } -> std::convertible_to<VectorF<D>>;
  { t.boundingBox() } -> std::convertible_to<Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)>>;
  { t.spacing() } -> std::convertible_to<Vector4F>;
};

/// A reader for an Euclidean distance transform
/**
 * A reader for an Euclidean distance transform. This reader reads an Euclidean distance transform from a file and
 * provides functionality to query the Euclidean distance transform.
 * The file may be in any format, as long as it is supported by ITK. ITK is then used to compute the Euclidean distance
 * transform from the binary file. Outside values must be set to 0, the inside positive.
 * The centers of the voxels lie at integer coordinates + 0.5, so the bounding box lies at integer coordinates.
 *
 * @tparam D The dimension of the space
 */
template <unsigned D> class EDTReader {
  struct Impl;
  std::unique_ptr<Impl> pimpl_;

public:
  /// Construct an Euclidean distance transform reader
  /**
   * Construct an Euclidean distance transform reader. The reader reads the Euclidean distance transform from the file
   * and stores it in memory. The file may be in any format, as long as it is supported by ITK. ITK is then used to
   * compute the Euclidean distance transform from the binary file. Outside values must be set to 0, the inside
   * positive.
   *
   * @param filename The filename of the Euclidean distance transform
   */
  explicit EDTReader(const std::string &filename);

  // make pImpl unique_ptr deletable with incomplete type
  ~EDTReader();
  EDTReader(const EDTReader &) = delete;
  EDTReader &operator=(const EDTReader &) = delete;
  EDTReader(EDTReader &&) noexcept = default;
  EDTReader &operator=(EDTReader &&) noexcept;

  /// Get the spacing of the Euclidean distance transform
  /**
   * Get the spacing of the Euclidean distance transform. The spacing is the distance between the centers of the voxels
   * in the Euclidean distance transform.
   *
   * @return The spacing of the Euclidean distance transform
   */
  [[nodiscard]] Vector4F spacing() const noexcept;

  /// Get the signed distance at a point
  /**
   * Get the signed distance at a point. The signed distance is the distance to the surface of the shape. The sign of
   * the distance is positive if the point is outside the shape, and negative if the point is inside the shape.
   *
   * @param point The point to get the signed distance at
   * @return The signed distance at the point
   */
  [[nodiscard]] FLOAT_T signedDistanceAt(const VectorF<D> &point) const noexcept;

  /// Get the signed distance in a region
  /**
   * Get the signed distance in a region. The signed distance is the distance to the surface of the shape. The sign of
   * the distance is positive if the point is outside the shape, and negative if the point is inside the shape.
   *
   * @param point The point to get the signed distance at
   * @param size The size of the region to get the signed distance in
   * @return The signed distance in the region
   */
  [[nodiscard]] std::vector<FLOAT_T> signedDistanceAt(const VectorF<D> &point,
                                                      const std::span<const size_t, D> &size) const noexcept;

  /// An extension to use this as a BoundaryRegionManager
  /**
   * An extension to use this as a BoundaryRegionManager. This function returns the index of the original image at the
   * point.
   *
   * @param point The point to find the boundary region of
   * @return The index of the boundary region that the point is inside
   */
  [[nodiscard]] size_t findBoundaryRegion(const Vector4F &point) const noexcept;

  /// Get the closest point on the surface at a point
  /**
   * Get the closest point on the surface at a point. The closest point is the point on the surface that is closest to
   * the input point. The return value will always be the center of a voxel.
   *
   * @param point The point to get the closest point on the surface at
   * @return The closest point on the surface at the point
   */
  [[nodiscard]] VectorF<D> closestAt(const VectorF<D> &point) const noexcept;

  /// Get the bounding box of the Euclidean distance transform
  /**
   * Get the bounding box of the Euclidean distance transform. The bounding box is the bounding box of the positive
   * values in the Euclidean distance transform.
   *
   * @return The bounding box of the Euclidean distance transform
   */
  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> boundingBox() const noexcept;
};

using EDTReader4 = EDTReader<4>;
} // namespace stmesh
#endif