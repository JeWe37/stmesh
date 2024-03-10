#ifndef STMESH_MESHING_ALGORITHM_HPP
#define STMESH_MESHING_ALGORITHM_HPP

#include <algorithm>
#include <concepts> // IWYU pragma: keep
#include <limits>
#include <optional>
#include <random>
#include <set>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Geometry>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/parameters.hpp>
#include <boost/geometry/index/predicates.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/iterator.hpp>
#include <spdlog/spdlog.h>

#include "geometric_simplex.hpp"
#include "meshing_cell.hpp"
#include "sdf.hpp"
#include "surface_adapters.hpp" // IWYU pragma: keep
#include "triangulation.hpp"
#include "utility.hpp"

namespace stmesh {
/// The central class for the meshing algorithm
/**
 * The central class for the meshing algorithm. The meshing algorithm is a mesh generation algorithm that generates a
 * mesh of the interior of a surface. The meshing algorithm is based on the paper "4D Space-Time Delaunay Meshing for
 * Medical Images" by Panagiotis Foteinos and Nikos Chrisochoides. The meshing algorithm is based on the concept of
 * rules. The meshing algorithm is a template class, with the surface adapter and random number generator as template
 * parameters. The surface adapter is a class that provides the signed distance, normal, and bounding box of a surface.
 * The random number generator is a class that provides random numbers.
 * Fundamentally the algorithms goal is the production of a mesh that conforms well to the surface and contains only
 * good pentatopes, by some definition of good.
 * This is useful for space-time finite element algorithms, such as offered by XNS.
 *
 * @tparam Surface The surface adapter
 * @tparam Random The random number generator
 */
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

  FLOAT_T rho_bar_, tau_bar_, zeta_, b_, delta_, max_radius_;

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
  friend struct detail::Rule6;
  friend struct detail::Complete;

public:
  /// Add a dependency between two full cells
  /**
   * Adds a dependency between two full cells. This means that when the dependency is removed, the dependent will be
   * refreshed. This is used to ensure that for rules which depend on their surrounding cells, they will always be
   * correctly checked.
   *
   * @param dependency The full cell that the dependent depends on
   * @param neighbor_index The index of the neighbor that is the dependent
   */
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

  /// Add a depndency between a full cell and a region
  /**
   * Adds a dependency between a full cell and a region. This means that when a point is inserted into the sphere, the
   * full cell will be refreshed. This is used to ensure that for rules which depend on their surroundings more
   * generally. Such dependencies are tracked separately for full cells sensitive to removal and addition of vertices.
   *
   * @param sphere The sphere insertion into which the full cell is sensitive to
   * @param dependency The full cell that the sphere depends on
   * @param on_insert Whether the dependency is on insertion or removal
   */
  void addPointDependency(const HyperSphere4 &sphere, const typename detail::Triangulation::FullCellHandle dependency,
                          bool on_insert = true) {
    if (on_insert)
      tree_insertion_.insert({detail::Triangulation::boxFromAABB(sphere.boundingBox()), dependency, sphere});
    else
      tree_removal_.insert({detail::Triangulation::boxFromAABB(sphere.boundingBox()), dependency, sphere});
  }

  /// Remove a dependency between a full cell and a region
  /**
   * Removes a dependency between a full cell and a region. This means that when a point is inserted into the sphere,
   * the full cell will no longer be refreshed. This is used to ensure that for rules which depend on their surroundings
   * more generally. Such dependencies are tracked separately for full cells sensitive to removal and addition of
   * vertices.
   *
   * @param sphere The sphere insertion into which the full cell is sensitive to
   * @param dependency The full cell that the sphere depends on
   * @param on_insert Whether the dependency is on insertion or removal
   */
  void removePointDependency(const HyperSphere4 &sphere,
                             const typename detail::Triangulation::FullCellHandle dependency, bool on_insert = true) {
    if (on_insert)
      tree_insertion_.remove({detail::Triangulation::boxFromAABB(sphere.boundingBox()), dependency, sphere});
    else
      tree_removal_.remove({detail::Triangulation::boxFromAABB(sphere.boundingBox()), dependency, sphere});
  }

  /// Get the full cells that are sensitive to a point
  /**
   * Gets the full cells that are sensitive to a point. These are the full cells that have registered a region
   * dependendency, which this point lies inside of.
   *
   * @param point The point to check for
   * @param on_insert Whether to check for insertion or removal
   * @return The full cells that are sensitive to the point
   */
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

  /// Construct a meshing algorithm
  /**
   * Constructs a meshing algorithm. The meshing algorithm is constructed from a surface adapter, a rho bar, a tau bar,
   * a zeta, a b, and a delta. The rho bar is the maximum radius-edge ratio of a pentatope. The tau bar is the minimum
   * quality of a pentatope. The zeta is the factor by which the picking region is scaled. The b is the maximum radius
   * of small simplices. The delta is the maximum size of simplices. should maintain. The max radius is the maximum
   * radius of a full cell. The seed is an optional seed for the random number generator.
   *
   * @param surface The surface adapter
   * @param rho_bar The maximum radius-edge ratio of a pentatope
   * @param tau_bar The minimum quality of a pentatope
   * @param zeta The factor by which the picking region is scaled
   * @param b The maximum radius of small simplices
   * @param delta The maximum size of simplices
   * @param max_radius The maximum radius of a full cell
   * @param seed The seed for the random number generator
   */
  MeshingAlgorithm(const Surface &surface, FLOAT_T rho_bar, FLOAT_T tau_bar, FLOAT_T zeta, FLOAT_T b, FLOAT_T delta,
                   FLOAT_T max_radius, std::optional<typename Random::result_type> seed = std::nullopt)
      : triangulation_(calculateBoundingBox(surface, delta)), surface_(surface), gen_(), rho_bar_(rho_bar),
        tau_bar_(tau_bar), zeta_(zeta), b_(b), delta_(delta), max_radius_(max_radius) {
    if (seed)
      gen_.seed(*seed);
    for (auto &full_cell : triangulation_) {
      typename detail::Triangulation::FullCellHandle full_cell_handle{&full_cell};
      full_cell_handle->data().extra_data.cell_handle =
          queue_.push({rulesSatisfied(full_cell_handle), static_cast<int>(gen_()), full_cell_handle});
    }
  }

  /// Calculate the bounding box that is far enough from surface
  /**
   * Calculates the bounding box that is far enough from the surface. The bounding box is calculated by taking the
   * bounding box of the surface and extending it, such that none of its corners are closer to the surface than the
   * delta.
   *
   * @param surface The surface
   * @param delta The minimum distance from the surface
   * @return The bounding box that is far enough from the surface
   */
  [[nodiscard]] static Eigen::AlignedBox<FLOAT_T, 4> calculateBoundingBox(const Surface &surface, FLOAT_T delta) {
    Eigen::AlignedBox<FLOAT_T, 4> bounding_box = surface.boundingBox();
    bounding_box.min().array() -= 2 * delta;
    bounding_box.max().array() += 2 * delta;
    return bounding_box;
  }

  /// Get the underlying triangulation
  /**
   * Gets the underlying triangulation. The underlying triangulation is the triangulation that the meshing algorithm
   * operates on.
   *
   * @return The underlying triangulation
   */
  [[nodiscard]] const detail::Triangulation &triangulation() const noexcept { return triangulation_; }

  /// Get the central queue of the meshing algorithm
  /**
   * Gets the central queue of the meshing algorithm. The central queue is the queue that the meshing algorithm uses to
   * prioritize the application of rules.
   *
   * @return The central queue
   */
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
  /// Get the first rule that is not satisfied by a full cell
  /**
   * Gets the first rule that is not satisfied by a full cell. The first rule that is not satisfied is the first rule
   * that is not satisfied by the full cell, when the rules are checked in order. If all rules are satisfied, the
   * complete rule is returned. The return value is a variant of the rules and the complete rule.
   *
   * @param full_cell The full cell to check
   * @return The first rule that is not satisfied by the full cell
   */
  [[nodiscard]] detail::Rules rulesSatisfied(const typename detail::Triangulation::FullCellHandle &full_cell) noexcept {
    return rulesSatisfiedImpl<detail::Rule1, detail::Rule2, detail::Rule3, detail::Rule4, detail::Rule5, detail::Rule6,
                              detail::Complete>(full_cell, {});
  }

  /// Triangulate the mesh
  /**
   * Triangulates the mesh. The mesh is triangulated by applying rules until the complete rule is satisfied. The
   * triangulation is updated after each rule is applied. The callback is called after each rule is applied.
   *
   * @tparam F The type of the callback
   * @param callback The callback to call after each rule is applied
   */
  template <typename F = void (*)()> void triangulate(const F &callback = [] {}) {
    const detail::Rules *rule = nullptr;
    while (std::visit([](const auto &r) { return r.kIndex; }, *(rule = &queue_.top().rule)) !=
           detail::Complete::kIndex) {
      spdlog::debug("Applying rule {}", std::visit([](const auto &r) { return r.kIndex; }, *rule) + 1);
      std::visit([&](auto &r) { r.apply(*this, queue_.top().full_cell); }, *rule);
      spdlog::debug("Number of vertices: {}", triangulation_.vertexCount());
      callback();
    }
  }

  /// Get the voronoi dual of a face
  /**
   * Gets the voronoi dual of a face. The voronoi dual of a face is the point that lies on the intersection of the
   * surface and the normal ray of the face. The return value is an optional of the voronoi dual, as for some faces the
   * opposite vertices may lie closer than the surface. The dependent neighbor info is an optional output parameter that
   * is set to the neighbor which made this cell have a voronoi dual.
   *
   * @param vertices The vertices of the face
   * @param dependent_neighbor_info The dependent neighbor info
   */
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

  /// Get the surface ball of a face
  /**
   * Gets the surface ball of a face. The surface ball is the ball centered at the voronoi dual circumscribing the
   * veritices of the face.
   *
   * @param vertices The vertices of the face
   */
  [[nodiscard]] HyperSphere4 surfaceBall(const Eigen::Matrix<FLOAT_T, 4, 4> &vertices) const {
    Vector4F voronoi_dual = voronoiDual(vertices).value();
    FLOAT_T radius = (vertices.col(0) - voronoi_dual).norm();
    return {radius, voronoi_dual};
  }

  /// Sample from the picking region of a full cell
  /**
   * Samples from the picking region of a full cell. The picking region is the region from which points are picked to
   * check if they are good. It is the circumsphere of the full cell scaled by the zeta. The return value is a pair of
   * the sample and the radius of the circumsphere.
   *
   * @param vertices The vertices of the full cell
   * @return The sample and the radius of the circumsphere
   */
  [[nodiscard]] std::pair<Vector4F, FLOAT_T>
  sampleFromPickingRegion(const Eigen::Matrix<FLOAT_T, 4, 5> &vertices) noexcept {
    HyperSphere<4> circumsphere = GeometricSimplex<4>(vertices).circumsphere();
    circumsphere.scale(zeta_);
    Vector4F sample;
    do {
      sample = circumsphere.sample(gen_);
    } while (!triangulation_.boundingBox().contains(sample));
    return {sample, circumsphere.radius()};
  }

  /// Sample from the picking region of a facet
  /**
   * Samples from the picking region of a facet. The picking region is the region from which points are picked to check
   * if they are good. It is the surface ball of the facet scaled by the zeta. The return value is a pair of the sample
   * and the radius of the surface ball.
   *
   * @param vertices The vertices of the facet
   * @return The sample and the radius of the surface ball
   */
  [[nodiscard]] std::pair<Vector4F, FLOAT_T> sampleFromPickingRegion(const Eigen::Matrix<FLOAT_T, 4, 4> &vertices) {
    HyperSphere<4> surface_ball = surfaceBall(vertices);
    surface_ball.scale(zeta_);
    Vector4F sample;
    do {
      sample = surface_.closestPoint(surface_ball.sample(gen_));
    } while (surface_ball.signedDistance(sample) > 0 || !triangulation_.boundingBox().contains(sample));
    return {sample, surface_ball.radius()};
  }

  /// Pick a good point from a full cell or facet
  /**
   * Picks a good point from a full cell or facet. Returns the vertex handle of the inserted good point. To find it,
   * samples randomly from the picking region until a good point is found.
   *
   * @tparam N The dimension of the vertices
   * @param vertices The vertices of the full cell or facet
   * @param hint The hint for the full cell
   * @return The vertex handle of the good point
   */
  template <int N>
  // needs to take vertices by copy, as the cell calling this may be removed during iteration
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

  /// The delta of the surface
  /**
   * The delta of the surface. The delta of the surface is typically supposed to be the local feature size of the
   * surface, however as calculating the medial axis is difficult, it is instead set to a constant.
   *
   * @return The delta of the surface
   */
  [[nodiscard]] FLOAT_T deltaSurface(const Vector4F & /*unused*/) const noexcept { return delta_; }

  /// Insert a point into the triangulation
  /**
   * Inserts a point into the triangulation. The heap is updated with the inserted full cells and the removed cells.
   *
   * @param point The point to insert
   * @param hint The hint for the full cell
   * @param nonfree_vertex Whether the vertex is nonfree
   * @return The vertex handle of the inserted point
   */
  detail::Triangulation::VertexHandle insert(const Vector4F &point, detail::Triangulation::FullCellHandle hint = {},
                                             bool nonfree_vertex = false) noexcept {
    auto [removed, dependent_full_cells] =
        processRemovedCells(triangulation_.conflictZone(point, hint), potentialsInRadius(point));
    typename detail::Triangulation::VertexHandle vertex = triangulation_.insert(point, hint, nonfree_vertex);
    std::vector<typename detail::Triangulation::FullCellHandle> inserted = triangulation_.surroundingFullCells(vertex);
    updateHeap(inserted, removed, dependent_full_cells);
    return vertex;
  }

  /// Remove a point from the triangulation
  /**
   * Removes a point from the triangulation. The heap is updated with the inserted full cells and the removed cells.
   *
   * @param vertex The vertex to remove
   */
  void remove(const detail::Triangulation::VertexHandle &vertex) {
    Vector4F point = triangulation_.pointToVec(vertex->point());
    auto [removed, dependent_full_cells] =
        processRemovedCells(triangulation_.surroundingFullCells(vertex, false), potentialsInRadius(point, false));
    typename detail::Triangulation::FullCellHandle center_full_cell = triangulation_.remove(vertex);
    std::vector<typename detail::Triangulation::FullCellHandle> inserted =
        triangulation_.commitUncommitted(center_full_cell);
    updateHeap(inserted, removed, dependent_full_cells);
  }

  /// Min time of the triangulation
  /**
   * The minimum time of the triangulation. The minimum time of the triangulation is the minimum time of a vertex that
   * is part of a full cell whose circumcenter is inside the surface.
   *
   * @return The minimum time of the triangulation
   */
  [[nodiscard]] FLOAT_T minTime() const noexcept {
    FLOAT_T min_time = std::numeric_limits<FLOAT_T>::infinity();
    for (const auto &full_cell : triangulation_) {
      if (surface_.inside(
              triangulation_.fullCellSimplex(typename detail::Triangulation::FullCellConstHandle{&full_cell})
                  .circumsphere()
                  .center())) {
        for (int i = 0; i < 5; ++i)
          min_time = std::min(min_time, full_cell.vertex(i)->point()[3]);
      }
    }
    return min_time;
  }
};
} // namespace stmesh
#endif
