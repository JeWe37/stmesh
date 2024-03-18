#ifndef STMESH_MESHING_ALGORITHM_HPP
#define STMESH_MESHING_ALGORITHM_HPP

#include <algorithm>
#include <concepts> // IWYU pragma: keep
#include <optional>
#include <random>
#include <tuple>
#include <variant>
#include <vector>

#include <Eigen/Geometry>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/iterator.hpp>
#include <spdlog/spdlog.h>

#include "meshing_cell.hpp"
#include "stmesh/geometric_simplex.hpp"
#include "stmesh/utility.hpp"
#include "surface_adapters.hpp" // IWYU pragma: keep
#include "triangulation.hpp"

namespace stmesh {
template <SurfaceAdapter4 Surface, std::uniform_random_bit_generator Random = std::mt19937_64> class MeshingAlgorithm {

  // may need map from full cell to fibonacci heap handle
  // need to check if those are ever invalidated
  boost::heap::fibonacci_heap<detail::Cell> queue_;

  detail::Triangulation triangulation_;
  Surface surface_;
  Random gen_;

  using Point = detail::Triangulation::BGPoint;
  using Box = detail::Triangulation::BGBox;
  using Tree = bg::index::rtree<std::tuple<Box, typename detail::Triangulation::FullCellHandle, HyperSphere4>,
                                // NOLINTNEXTLINE(*-magic-numbers)
                                bg::index::rstar<16>>;

  Tree tree_removal_;
  Tree tree_insertion_;

  FLOAT_T rho_bar_, tau_bar_, zeta_, b_, delta_;

  void updateHeap(const std::vector<typename detail::Triangulation::FullCellHandle> &inserted,
                  const std::vector<detail::CellHandle> &removed,
                  const std::set<typename detail::Triangulation::FullCellHandle> &dependent_full_cells) noexcept {
    for (const auto &full_cell : inserted)
      full_cell->data().extra_data.cell_handle =
          queue_.push({rulesSatisfied(full_cell), static_cast<int>(gen_()), full_cell});
    for (const auto &handle : removed) {
      std::visit([&](const auto &r) { r.remove(*this, (*handle).full_cell); }, (*handle).rule);
      queue_.erase(handle);
    }
    for (const auto &full_cell : dependent_full_cells) {
      std::visit([&](const auto &r) { r.remove(*this, full_cell); }, (*full_cell->data().extra_data.cell_handle).rule);
      queue_.erase(full_cell->data().extra_data.cell_handle);
      full_cell->data().extra_data.cell_handle =
          queue_.push({rulesSatisfied(full_cell), static_cast<int>(gen_()), full_cell});
    }
  }

  [[nodiscard]] static std::pair<std::vector<detail::CellHandle>,
                                 std::set<typename detail::Triangulation::FullCellHandle>>
  processRemovedCells(const std::vector<typename detail::Triangulation::FullCellHandle> &full_cells,
                      const std::vector<typename detail::Triangulation::FullCellHandle> &potentials) noexcept {
    std::vector<detail::CellHandle> handles;
    std::set<typename detail::Triangulation::FullCellHandle> dependent_full_cells(potentials.begin(), potentials.end());
    handles.reserve(full_cells.size());
    for (const auto &full_cell : full_cells) {
      handles.push_back(full_cell->data().extra_data.cell_handle);
      // NOLINTNEXTLINE(*-magic-numbers)
      for (unsigned i = 0; i < 5; ++i) {
        if ((full_cell->data().extra_data.dependents & static_cast<unsigned>(static_cast<unsigned char>(1) << i)) != 0U)
          dependent_full_cells.insert(full_cell->neighbor(static_cast<int>(i)));
      }
    }
    for (const auto &full_cell : full_cells)
      dependent_full_cells.erase(full_cell);
    return {handles, dependent_full_cells};
  }

  friend struct detail::Rule1;
  friend struct detail::Rule2;
  friend struct detail::Rule3;
  friend struct detail::Rule4;
  friend struct detail::Rule5;
  friend struct detail::Complete;

public:
  static void addDependency(const typename detail::Triangulation::FullCellHandle &dependency, int neighbor_index) {
    // if A(dependent) depends on B(dependency), then refresh A when B is removed
    // A depends on B if its rule is 5, which was confirmed using the mirror vertex of B
    // thus A can only depend on one B
    // B meanwhile can be the dependency of multiple A
    // thus viewing from the left, we have a one to many relationship, aka a map
    // however, viewing from the right, we have a many to one relationship, aka a multimap
    dependency->data().extra_data.dependents |=
        static_cast<unsigned char>(static_cast<unsigned char>(1) << static_cast<unsigned>(neighbor_index));
  }

  void addPointDependency(const HyperSphere4 &sphere, const typename detail::Triangulation::FullCellHandle dependency,
                          bool on_insert = true) {
    if (on_insert)
      tree_insertion_.insert({detail::Triangulation::boxFromAABB(sphere.boundingBox()), dependency, sphere});
    else
      tree_removal_.insert({detail::Triangulation::boxFromAABB(sphere.boundingBox()), dependency, sphere});
  }

  void removePointDependency(const HyperSphere4 &sphere,
                             const typename detail::Triangulation::FullCellHandle dependency, bool on_insert = true) {
    if (on_insert)
      tree_insertion_.remove({detail::Triangulation::boxFromAABB(sphere.boundingBox()), dependency, sphere});
    else
      tree_removal_.remove({detail::Triangulation::boxFromAABB(sphere.boundingBox()), dependency, sphere});
  }

  [[nodiscard]] std::vector<detail::Triangulation::FullCellHandle> potentialsInRadius(const Vector4F &point,
                                                                                      bool on_insert = true) const {
    Point p = detail::Triangulation::pointFromVector(point);
    auto query = bg::index::covers(p) &&
                 bg::index::satisfies([&](const auto &value) { return std::get<2>(value).signedDistance(point) < 0; });
    std::vector<detail::Triangulation::FullCellHandle> full_cells;
    auto output_iterator =
        boost::make_function_output_iterator([&](const auto &value) { full_cells.push_back(std::get<1>(value)); });
    if (on_insert)
      tree_insertion_.query(query, output_iterator);
    else
      tree_removal_.query(query, output_iterator);
    return full_cells;
  }

  // TODO: maybe move in surface rather than making a copy
  MeshingAlgorithm(const Surface &surface, FLOAT_T rho_bar, FLOAT_T tau_bar, FLOAT_T zeta, FLOAT_T b, FLOAT_T delta,
                   std::optional<typename Random::result_type> seed = std::nullopt)
      : triangulation_(calculateBoundingBox(surface, delta)), surface_(surface), gen_(), rho_bar_(rho_bar),
        tau_bar_(tau_bar), zeta_(zeta), b_(b), delta_(delta) {
    if (seed)
      gen_.seed(*seed);
    for (auto &full_cell : triangulation_) {
      typename detail::Triangulation::FullCellHandle full_cell_handle{&full_cell};
      full_cell_handle->data().extra_data.cell_handle =
          queue_.push({rulesSatisfied(full_cell_handle), static_cast<int>(gen_()), full_cell_handle});
    }
  }

  [[nodiscard]] static Eigen::AlignedBox<FLOAT_T, 4> calculateBoundingBox(const Surface &surface, FLOAT_T delta) {
    Eigen::AlignedBox<FLOAT_T, 4> bounding_box = surface.boundingBox();
    for (const Vector4F &corner : allCorners(bounding_box)) {
      Vector4F cfp = surface.closestPoint(corner);
      // TODO: fix if deltaSurface is ever changed
      FLOAT_T add_dist = std::max(FLOAT_T{0}, 2 * delta - (corner - cfp).norm());
      bounding_box.extend(corner + add_dist * (corner - cfp).normalized());
    }
    return bounding_box;
  }

  [[nodiscard]] const detail::Triangulation &triangulation() const noexcept { return triangulation_; }

  [[nodiscard]] const boost::heap::fibonacci_heap<detail::Cell> &queue() const noexcept { return queue_; }

private:
  template <typename Rule, typename... Rules>
  [[nodiscard]] detail::Rules rulesSatisfiedImpl(const typename detail::Triangulation::FullCellHandle &full_cell,
                                                 const detail::Rule &pre_rule) noexcept {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
    Rule rule(pre_rule);
#pragma GCC diagnostic pop
    if (rule.check(*this, full_cell))
      return rule;
    else if constexpr (sizeof...(Rules) > 0)
      return rulesSatisfiedImpl<Rules...>(full_cell, rule);
    // std::unreachable() not available in C++20
    __builtin_unreachable();
  }

public:
  [[nodiscard]] detail::Rules rulesSatisfied(const typename detail::Triangulation::FullCellHandle &full_cell) noexcept {
    return rulesSatisfiedImpl<detail::Rule1, detail::Rule2, detail::Rule3, detail::Rule4, detail::Rule5,
                              detail::Complete>(full_cell, {});
  }

  template <typename F = void (*)(void)> void triangulate(const F &callback = [] {}) {
    const detail::Rules *rule = nullptr;
    while (std::visit([](const auto &r) { return r.kIndex; }, *(rule = &queue_.top().rule)) !=
           detail::Complete::kIndex) {
      spdlog::info("Applying rule {}", std::visit([](const auto &r) { return r.kIndex; }, *rule) + 1);
      std::visit([&](auto &r) { r.apply(*this, queue_.top().full_cell); }, *rule);
      spdlog::info("Number of vertices: {}", triangulation_.vertexCount());
      callback();
    }
  }

  [[nodiscard]] std::optional<Vector4F>
  voronoiDual(const Eigen::Matrix<FLOAT_T, 4, 4> &vertices,
              std::pair<detail::Triangulation::FullCellHandle, int> *dependent_neighbor_info = nullptr) const {
    const auto [existing_side, optional_side] =
        triangulation_.facetMirrorVertices(triangulation_.facetFromVertices(vertices));
    const auto &[facet_mirror_vertex0, facet_side_mirror_index0, facet_side0] = existing_side;

    Vector4F mirror_vertex0 = triangulation_.pointToVec(facet_mirror_vertex0->point());

    GeometricSimplex<4, 4> simplex(vertices);
    Eigen::ParametrizedLine<FLOAT_T, 4> normal_ray = simplex.normalRay();
    Eigen::ParametrizedLine<FLOAT_T, 4> inverse_ray{normal_ray.origin(), -normal_ray.direction()};

    // we swap the rays if normal ray does not already point in the direction of mirror_vertex0, such that it does after
    if (normal_ray.direction().dot(mirror_vertex0 - normal_ray.origin()) < 0)
      std::swap(normal_ray, inverse_ray);

    const auto test =
        [&](const Eigen::ParametrizedLine<FLOAT_T, 4> &ray,
            const std::tuple<detail::Triangulation::VertexHandle, int, detail::Triangulation::FullCellHandle> &info)
        -> std::optional<Vector4F> {
      const auto &[vertex, index, full_cell] = info;
      std::optional<Vector4F> point = surface_.raycast(ray, triangulation_.boundingBox().diagonal().norm());
      if (point) {
        // issue: the insertion of a point can change the closest point, thus the face may no longer be restricted
        // however, this is only used in Rule5. It is not possible for a point to become further away due to
        // insertion this means it is sufficient to recheck the condition inside the application and do nothing if
        // no longer needed the issue is that rule5 also includes the removal of vertices. this means complete rules
        // could become invalid we need to store which (complete) rules depend on the existence of a point and
        // reevaluate them on removal NOLINTNEXTLINE(cppcoreguidelines-init-variables)
        Vector4F facet_mirror_vertex_vec = triangulation_.pointToVec(vertex->point());
        if ((*point - facet_mirror_vertex_vec).squaredNorm() > (*point - vertices.col(0)).squaredNorm()) {
          if (dependent_neighbor_info != nullptr)
            *dependent_neighbor_info = {full_cell, index};
          return point;
        }
      }
      return std::nullopt;
    };

    if (std::optional result = test(normal_ray, existing_side); result)
      return result;
    if (optional_side) {
      if (std::optional result = test(inverse_ray, *optional_side); result)
        return result;
    }
    return std::nullopt;
  }

  [[nodiscard]] HyperSphere4 surfaceBall(const Eigen::Matrix<FLOAT_T, 4, 4> &vertices) const {
    Vector4F voronoi_dual = voronoiDual(vertices).value();
    FLOAT_T radius = (vertices.col(0) - voronoi_dual).norm();
    return {radius, voronoi_dual};
  }

  [[nodiscard]] std::pair<Vector4F, FLOAT_T>
  // NOLINTNEXTLINE(*-magic-numbers)
  sampleFromPickingRegion(const Eigen::Matrix<FLOAT_T, 4, 5> &vertices) noexcept {
    HyperSphere<4> circumsphere = GeometricSimplex<4>(vertices).circumsphere();
    circumsphere.scale(zeta_);
    Vector4F sample;
    do {
      sample = circumsphere.sample(gen_);
    } while (!triangulation_.boundingBox().contains(sample));
    return {sample, circumsphere.radius()};
  }

  [[nodiscard]] std::pair<Vector4F, FLOAT_T> sampleFromPickingRegion(const Eigen::Matrix<FLOAT_T, 4, 4> &vertices) {
    HyperSphere<4> surface_ball = surfaceBall(vertices);
    surface_ball.scale(zeta_);
    Vector4F sample;
    do {
      sample = surface_.closestPoint(surface_ball.sample(gen_));
    } while (surface_ball.signedDistance(sample) > 0 || !triangulation_.boundingBox().contains(sample));
    return {sample, surface_ball.radius()};
  }

  // needs to take vertices by copy, as the cell calling this may be removed during iteration
  template <int N>
  detail::Triangulation::VertexHandle pickGoodPoint(const Eigen::Matrix<FLOAT_T, 4, N> vertices,
                                                    detail::Triangulation::FullCellHandle hint = {}) noexcept {
    // TODO: can be sped up, this is comically inefficient
    while (true) {
      auto [sample, radius] = sampleFromPickingRegion(vertices);
      detail::Triangulation::VertexHandle vertex = insert(sample, hint, N == 4);
      if (triangulation_.isGoodPoint(vertex, rho_bar_, tau_bar_, b_ * radius))
        return vertex;
      remove(vertex);
    }
  }

  [[nodiscard]] FLOAT_T deltaSurface(const Vector4F & /*unused*/) const noexcept { return delta_; }

  detail::Triangulation::VertexHandle insert(const Vector4F &point, detail::Triangulation::FullCellHandle hint = {},
                                             bool nonfree_vertex = false) noexcept {
    auto [removed, dependent_full_cells] =
        processRemovedCells(triangulation_.conflictZone(point, hint), potentialsInRadius(point));
    typename detail::Triangulation::VertexHandle vertex = triangulation_.insert(point, hint, nonfree_vertex);
    std::vector<typename detail::Triangulation::FullCellHandle> inserted = triangulation_.surroundingFullCells(vertex);
    updateHeap(inserted, removed, dependent_full_cells);
    return vertex;
  }

  void remove(const detail::Triangulation::VertexHandle &vertex) {
    Vector4F point = triangulation_.pointToVec(vertex->point());
    auto [removed, dependent_full_cells] =
        processRemovedCells(triangulation_.surroundingFullCells(vertex, false), potentialsInRadius(point, false));
    typename detail::Triangulation::FullCellHandle center_full_cell = triangulation_.remove(vertex);
    std::vector<typename detail::Triangulation::FullCellHandle> inserted =
        triangulation_.commitUncommitted(center_full_cell);
    updateHeap(inserted, removed, dependent_full_cells);
  }
};
} // namespace stmesh
#endif
