#ifndef STMESH_EDT_READER_HPP
#define STMESH_EDT_READER_HPP

#include <concepts> // IWYU pragma: keep
#include <string>
#include <vector>

#include <Eigen/Geometry>

#include "utility.hpp"

namespace stmesh {
template <typename T, unsigned D>
concept EuclideanDistanceTransform = requires(const T t, VectorF<D> vec) {
  { t.signedDistanceAt(vec) } -> std::convertible_to<FLOAT_T>;
  { t.closestAt(vec) } -> std::convertible_to<VectorF<D>>;
  { t.boundingBox() } -> std::convertible_to<Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)>>;
};

template <unsigned D> class EDTReader {
  std::vector<FLOAT_T> distance_map_;
  std::vector<Vector4F> vector_map_;
  std::array<size_t, D> sizes_;
  Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> bounding_box_;

  [[nodiscard]] Eigen::Vector<size_t, 4> clamp(const Vector4F &point) const noexcept;

  [[nodiscard]] size_t projection(const Vector4F &point) const noexcept;

public:
  explicit EDTReader(const std::string &filename);

  [[nodiscard]] FLOAT_T signedDistanceAt(const VectorF<D> &point) const noexcept;

  [[nodiscard]] VectorF<D> closestAt(const VectorF<D> &point) const noexcept;

  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> boundingBox() const noexcept;
};

using EDTReader4 = EDTReader<4>;
} // namespace stmesh
#endif