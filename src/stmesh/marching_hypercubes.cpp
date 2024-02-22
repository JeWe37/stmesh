#include "stmesh/marching_hypercubes.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include <CGAL/Dimension.h>
#include <CGAL/Epick_d.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <Eigen/Geometry>

#include "stmesh/table.hpp"
#include "stmesh/utility.hpp"

namespace stmesh {
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
std::optional<Vector4F> surfaceRayIntersection(const Vector4F &start, const Vector4F &end,
                                               const std::array<bool, (1U << 4U)> &corner_values, bool *inside) {
  using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<4>>;
  using Hyperplane = Kernel::Hyperplane_d;
  using Point = Kernel::Point_d;
  using Direction = Kernel::Direction_d;
  constexpr unsigned kHyperplaneStride = 5;
  constexpr FLOAT_T kEps = sqrt(std::numeric_limits<FLOAT_T>::epsilon());
  Kernel::Oriented_side_d side{};
  auto ray = Eigen::ParametrizedLine<FLOAT_T, 4>::Through(start, end);
  unsigned index = 0;
  for (unsigned i = 0; i < corner_values.size(); ++i) {
    if (corner_values.at(i))
      index |= 1U << (corner_values.size() - 1 - i);
  }
  // discard full and empty cubes
  if (index == 0 || index == (1U << corner_values.size()) - 1U) {
    if (inside != nullptr)
      *inside = index != 0;
    return std::nullopt;
  }
  auto hyperplane_from_index = [](unsigned idx) {
    const Vector4F normal{detail::kHyperplanes[idx], detail::kHyperplanes[idx + 1], detail::kHyperplanes[idx + 2],
                          detail::kHyperplanes[idx + 3]};
    const Direction direction(normal.begin(), normal.end());
    return std::make_pair(Hyperplane(direction, -detail::kHyperplanes[idx + 4]), normal);
  };
  std::vector<std::pair<std::vector<Hyperplane>, std::optional<std::pair<FLOAT_T, Point>>>> components;
  for (unsigned i = detail::kComponentIdxs[index]; i < detail::kComponentIdxs[index + 1]; ++i) {
    auto &component = components.emplace_back().first;
    std::vector<std::optional<std::pair<FLOAT_T, Point>>> intersections;
    // gather all hyperplanes and intersections
    for (unsigned j = detail::kHyperplaneIdxs[i]; j < detail::kHyperplaneIdxs[i + 1]; j += kHyperplaneStride) {
      const auto [hyperplane, normal] = hyperplane_from_index(j);
      component.emplace_back(hyperplane);

      const Eigen::Hyperplane<FLOAT_T, 4> plane(normal, detail::kHyperplanes[j + 4]);
      const FLOAT_T intersection_parameter = ray.intersectionParameter(plane);
      if (std::isfinite(intersection_parameter)) {
        Vector4F intersection_vector = ray.pointAt(intersection_parameter);
        if ((intersection_vector.array() <= FLOAT_T(1.0) + kEps * intersection_vector.array()).all() &&
            (intersection_vector.array() >= -kEps).all()) {
          const Point intersection{intersection_vector.begin(), intersection_vector.end()};
          intersections.emplace_back(std::pair{intersection_parameter, intersection});
          continue;
        }
      }
      intersections.emplace_back(std::nullopt);
    }
    // filter invalid intersections
    if (const auto it = std::ranges::find_if(
            intersections,
            [&, j = size_t{}](const std::optional<std::pair<FLOAT_T, Point>> &intersection) mutable {
              j++;
              return intersection && std::ranges::all_of(component, [&, k = size_t{}](const Hyperplane &plane) mutable {
                       return ++k == j || side(plane, intersection->second) == CGAL::Sign::ON_NEGATIVE_SIDE;
                     });
            });
        it != intersections.end())
      components.back().second = *it;
  }
  // check if inside
  if (inside != nullptr) {
    Point start_point(start.begin(), start.end());
    *inside = std::ranges::all_of(components, [&](const auto &component) {
      return std::ranges::all_of(
          component.first, [&](const auto &plane) { return side(plane, start_point) == CGAL::Sign::ON_NEGATIVE_SIDE; });
    });
  }
  // find first intersection
  if (const auto first_intersection = std::ranges::min_element(components, {},
                                                               [](const auto &component) {
                                                                 return component.second
                                                                            ? std::abs(component.second->first)
                                                                            : std::numeric_limits<FLOAT_T>::infinity();
                                                               });
      first_intersection != components.end() && first_intersection->second) {
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    const Point &point = first_intersection->second->second;
    return Vector4F{CGAL::to_double(point[0]), CGAL::to_double(point[1]), CGAL::to_double(point[2]),
                    CGAL::to_double(point[3])};
  }
  return std::nullopt;
}
} // namespace stmesh