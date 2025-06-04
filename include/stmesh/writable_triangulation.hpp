#ifndef STMESH_WRITABLE_TRIANGULATION_HPP
#define STMESH_WRITABLE_TRIANGULATION_HPP

#include <concepts>
#include <cstddef>
#include <ranges>

#include "geometric_simplex.hpp"

namespace stmesh {
/// A cell in a triangulation that can be written
/**
 * A cell in a triangulation that can be written. This concept requires the cell to have a geometric simplex and a
 * method to check if a side is a surface side.
 */
template <typename T>
concept WritableCell = requires(T cell, size_t i) {
  { cell.geometricSimplex() } -> std::convertible_to<GeometricSimplex<4>>;
  { cell.isSurfaceSide(i) } -> std::same_as<bool>;
};

namespace detail {
struct DummyCell {
  [[nodiscard]] GeometricSimplex<4> geometricSimplex() const noexcept;
  [[nodiscard]] bool isSurfaceSide(size_t) const noexcept;
  // NOLINTNEXTLINE(hicpp-explicit-conversions)
  template <WritableCell T> operator T() const noexcept;
};

static_assert(WritableCell<DummyCell>);
} // namespace detail

/// A triangulation that can be written
/**
 * A triangulation that can be written. This concept requires the triangulation to have a bounding box and to be a
 * forward range of WritableCell.
 */
template <typename T>
concept WritableTriangulation = requires(T triangulation) {
  { triangulation.boundingBox() } -> std::same_as<Eigen::AlignedBox<double, 4>>;
} && std::ranges::forward_range<T> && WritableCell<std::ranges::range_value_t<T>>;
} // namespace stmesh
#endif