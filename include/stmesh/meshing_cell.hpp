#ifndef STMESH_MESHING_CELL_HPP
#define STMESH_MESHING_CELL_HPP
#include <algorithm>
#include <boost/heap/fibonacci_heap.hpp>
#include <compare>
#include <limits>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

#include "geometric_simplex.hpp"
#include "sdf.hpp"
#include "triangulation.hpp"
#include "utility.hpp"

namespace stmesh::detail {
// ordering for max-heap
struct Cell;

using CellHandle = typename boost::heap::fibonacci_heap<Cell>::handle_type;

struct ExtraData {
  CellHandle cell_handle{};
  unsigned char dependents{};
};

using Triangulation = Triangulation<ExtraData>;

struct Rule {
  [[nodiscard]] static unsigned priority() noexcept { return 0; }
  Vector4F z0;

  template <typename MeshingAlgorithm>
  static unsigned check(const MeshingAlgorithm & /*unused*/, Triangulation::FullCellHandle /*unused*/,
                        bool /* unused */ = false) {
    return 1;
  }

  template <typename MeshingAlgorithm>
  static void apply(MeshingAlgorithm & /*unused*/, Triangulation::FullCellHandle /*unused*/) {
    throw std::runtime_error("Invalid meshing rule");
  }

  template <typename MeshingAlgorithm>
  void remove(MeshingAlgorithm &meshing_algorithm, const Triangulation::FullCellHandle full_cell) const {
    // if a cell that has passed rule 0 is removed, it must be complete
    // the question: what is z, need to persist that somehow
    // or just recompute it, but that's expensive(surface.closestPoint might be slow)
    if (z0 != Vector4F::Constant(std::numeric_limits<FLOAT_T>::infinity())) {
      HyperSphere4 sphere{meshing_algorithm.deltaSurface(z0), z0};
      meshing_algorithm.removePointDependency(sphere, full_cell, false);
    }
  }
};

struct Rule1 : Rule {
  static constexpr inline int kIndex = 0;

  template <typename MeshingAlgorithm>
  unsigned check(MeshingAlgorithm &meshing_algorithm, Triangulation::FullCellHandle full_cell, bool dry_run = false) {
    const auto &surface = meshing_algorithm.surface_;
    const detail::Triangulation &triangulation = meshing_algorithm.triangulation_;
    GeometricSimplex<4> simplex = triangulation.fullCellSimplex(full_cell);
    HyperSphere<4> circumsphere = simplex.circumsphere();
    if (surface.intersectedBySphere(circumsphere)) {
      z0 = surface.closestPoint(circumsphere.center());
      // TODO: track potential points in kd tree
      // if potential point is within radius, make it a dependency and vice-versa
      // that assumes constant deltaSurface...
      // need dependency system that supports this...
      // might be easier to just recheck this before applying, this is the first rule anyway
      // this is going to cause issues when reordering to make rule3 first tho!
      // the bigger issue: vertex removal. that might make a simplex need rule1 again,
      // rather than completing it without doing anything

      // EXACT:
      //   if this condition is false, recheck when a vertex within deltaSurface(z) of (the potential) z is removed
      //   if this condition is true, recheck if a vertex within deltaSurface(z) of z is added

      // since deltaSurface is constant, this test can just be done via a nearest neighbor search with appropriate bound
      HyperSphere4 sphere{meshing_algorithm.deltaSurface(z0), z0};
      if (triangulation.verticesInRadius(z0, sphere.radius()).empty()) {
        // if we insert a point, this might start failing, making the cell complete
        if (!dry_run)
          meshing_algorithm.addPointDependency(sphere, full_cell);
        return 1;
      } else {
        // if we remove a point, this might start passing, making the cell incomplete
        if (!dry_run)
          meshing_algorithm.addPointDependency(sphere, full_cell, false);
      }
    } else
      z0 = Vector4F::Constant(std::numeric_limits<FLOAT_T>::infinity());
    return 0;
  }

  template <typename MeshingAlgorithm>
  void apply(MeshingAlgorithm &meshing_algorithm, Triangulation::FullCellHandle full_cell) const {
    meshing_algorithm.insert(z0, full_cell, true);
    // dependency is removed by subsequent removal of cell
  }

  template <typename MeshingAlgorithm>
  void remove(MeshingAlgorithm &meshing_algorithm, const Triangulation::FullCellHandle full_cell) const {
    // if a cell that is rule 0 is removed, it must be incomplete
    HyperSphere4 sphere{meshing_algorithm.deltaSurface(z0), z0};
    meshing_algorithm.removePointDependency(sphere, full_cell);
  }

  // issue: what about the removal of a cell that is complete for rule0?
  // basically all the cases of removal of something already passing rule0 are not covered
};

struct Rule2 : Rule {
  static constexpr inline int kIndex = 1;
  Vector4F z;
  bool z_outside;

  template <typename MeshingAlgorithm>
  unsigned check(MeshingAlgorithm &meshing_algorithm, Triangulation::FullCellHandle full_cell,
                 bool /*unused*/ = false) {
    const auto &surface = meshing_algorithm.surface_;
    const detail::Triangulation &triangulation = meshing_algorithm.triangulation_;
    GeometricSimplex<4> simplex = triangulation.fullCellSimplex(full_cell);
    HyperSphere<4> circumsphere = simplex.circumsphere();
    if (surface.intersectedBySphere(circumsphere)) {
      z = circumsphere.center();
      if (circumsphere.radius() >= 2 * meshing_algorithm.deltaSurface(surface.closestPoint(z))) {
        if ((z_outside = !triangulation.boundingBox().contains(z)))
          z = z.cwiseMax(triangulation.boundingBox().min()).cwiseMin(triangulation.boundingBox().max());
        return 1;
      }
    }
    return 0;
  }

  template <typename MeshingAlgorithm>
  void apply(MeshingAlgorithm &meshing_algorithm, Triangulation::FullCellHandle full_cell) const {
    meshing_algorithm.insert(z, full_cell, z_outside);
  }
};

struct Rule3 : Rule {
  static constexpr inline int kIndex = 2;
  Vector4F z;

  template <typename MeshingAlgorithm>
  unsigned check(MeshingAlgorithm &meshing_algorithm, Triangulation::FullCellHandle full_cell,
                 bool /*unused*/ = false) {
    const auto &surface = meshing_algorithm.surface_;
    const detail::Triangulation &triangulation = meshing_algorithm.triangulation_;
    GeometricSimplex<4> simplex = triangulation.fullCellSimplex(full_cell);
    z = simplex.circumsphere().center();
    return surface.inside(z) && simplex.radiusEdgeRatio() >= meshing_algorithm.rho_bar_;
  }

  template <typename MeshingAlgorithm>
  void apply(MeshingAlgorithm &meshing_algorithm, Triangulation::FullCellHandle full_cell) const {
    meshing_algorithm.insert(z, full_cell);
  }
};

struct Rule4 : Rule {
  static constexpr inline int kIndex = 3;
  Eigen::Matrix<FLOAT_T, 4, 5> vertices;
  unsigned prio;

  [[nodiscard]] unsigned priority() const noexcept { return prio; }

  template <typename MeshingAlgorithm>
  unsigned check(MeshingAlgorithm &meshing_algorithm, Triangulation::FullCellHandle full_cell,
                 bool /*unused*/ = false) {
    const auto &surface = meshing_algorithm.surface_;
    const detail::Triangulation &triangulation = meshing_algorithm.triangulation_;
    GeometricSimplex<4> simplex = triangulation.fullCellSimplex(full_cell);
    HyperSphere<4> circumsphere = simplex.circumsphere();
    if (unsigned dim = simplex.sliverSimplex(meshing_algorithm.rho_bar_, meshing_algorithm.tau_bar_);
        surface.inside(circumsphere.center()) && dim) {
      vertices = simplex.vertices();
      return prio = dim;
    }
    return 0;
  }

  template <typename MeshingAlgorithm>
  void apply(MeshingAlgorithm &meshing_algorithm, Triangulation::FullCellHandle full_cell) const {
    meshing_algorithm.pickGoodPoint(vertices, full_cell);
  }
};

struct Rule5 : Rule {
  static constexpr inline int kIndex = 4;
  Eigen::Matrix<FLOAT_T, 4, 4> vertices;

  template <typename MeshingAlgorithm>
  unsigned check(MeshingAlgorithm &meshing_algorithm, Triangulation::FullCellHandle full_cell, bool dry_run = false) {
    const auto &triangulation = meshing_algorithm.triangulation_;
    Eigen::Index idx{};
    auto nonfree_count = std::count_if(full_cell->vertices_begin(), full_cell->vertices_end(),
                                       [&, i = Eigen::Index{}](const auto &vertex) mutable {
                                         if (vertex->data().nonfree_vertex)
                                           idx = i;
                                         i++;
                                         return vertex->data().nonfree_vertex;
                                       });
    std::pair<Triangulation::FullCellHandle, int> dependent_neighbor_info;
    auto full_cell_simplex = triangulation.fullCellSimplex(full_cell);
    switch (nonfree_count) {
    case 0: {
      auto sub_simplices = full_cell_simplex.template subSimplices<4>();
      if (auto it = std::ranges::find_if(
              sub_simplices,
              [&](const auto &simplex) {
                return meshing_algorithm.voronoiDual(simplex.vertices(), &dependent_neighbor_info).has_value();
              });
          it != sub_simplices.end()) {
        if (!dry_run && dependent_neighbor_info.first != full_cell)
          std::apply(MeshingAlgorithm::addDependency, dependent_neighbor_info);
        vertices = it->vertices();
        return 1;
      }
      return 0;
    }
    case 1: {
      Eigen::Index j = 0;
      for (Eigen::Index i = 0; i < 5; ++i) {
        if (i != idx)
          vertices.col(j++) = full_cell_simplex.vertices().col(i);
      }
      if (meshing_algorithm.voronoiDual(vertices, &dependent_neighbor_info).has_value()) {
        if (!dry_run && dependent_neighbor_info.first != full_cell)
          std::apply(MeshingAlgorithm::addDependency, dependent_neighbor_info);
        return 1;
      }
      return 0;
    }
    default:
      return 0;
    }
  }

  template <typename MeshingAlgorithm>
  void apply(MeshingAlgorithm &meshing_algorithm, Triangulation::FullCellHandle full_cell) const {
    detail::Triangulation &triangulation = meshing_algorithm.triangulation_;
    typename detail::Triangulation::VertexHandle vertex = meshing_algorithm.pickGoodPoint(vertices, full_cell);
    const typename detail::Triangulation::Point &p = vertex->point();
    Vector4F point{static_cast<FLOAT_T>(p[0]), static_cast<FLOAT_T>(p[1]), static_cast<FLOAT_T>(p[2]),
                   static_cast<FLOAT_T>(p[3])};
    std::vector<typename detail::Triangulation::VertexHandle> close_vertices =
        triangulation.verticesInRadius(point, meshing_algorithm.deltaSurface(point));
    for (const auto &close_vertex : close_vertices) {
      if (close_vertex != vertex && !close_vertex->data().nonfree_vertex)
        meshing_algorithm.remove(close_vertex);
    }
  }
};

struct Complete : Rule {
  static constexpr inline int kIndex = 5;

  template <typename MeshingAlgorithm>
  static unsigned check(const MeshingAlgorithm &meshing_algorithm, Triangulation::FullCellHandle full_cell,
                        bool dry_run = false) {
    // TODO: optimally these rechecks only need to recheck rule5 i think
    if (!dry_run) {
      for (int i = 0; i < 5; ++i)
        meshing_algorithm.addDependency(full_cell->neighbor(i),
                                        meshing_algorithm.triangulation().mirrorIndex(full_cell, i));
    }
    return 1;
  }
};

using Rules = std::variant<Rule1, Rule2, Rule3, Rule4, Rule5, Complete>;

struct Cell {
  Rules rule{};
  int random{};
  Triangulation::FullCellHandle full_cell{};

  friend auto operator<=>(const Cell &lhs, const Cell &rhs) noexcept {
    int lhsindex = std::visit([](const auto &r) { return r.kIndex; }, lhs.rule);
    int rhsindex = std::visit([](const auto &r) { return r.kIndex; }, rhs.rule);
    // max heap => reverse order, lower rule index, higher priority
    if (std::strong_ordering cmp = rhsindex <=> lhsindex; cmp != std::strong_ordering::equal)
      return cmp;
    unsigned lhspriority = std::visit([](const auto &r) { return r.priority(); }, lhs.rule);
    unsigned rhspriority = std::visit([](const auto &r) { return r.priority(); }, rhs.rule);
    // max heap => reverse order, lower dimension, higher priority
    if (std::strong_ordering cmp = rhspriority <=> lhspriority; cmp != std::strong_ordering::equal)
      return cmp;
    return lhs.random <=> rhs.random;
  }
};
} // namespace stmesh::detail

#endif