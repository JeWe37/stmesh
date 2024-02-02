#ifndef STMESH_TRIANGULATION_HPP
#define STMESH_TRIANGULATION_HPP

#include <unordered_map>
#include <variant>
#include <vector>

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "geometric_simplex.hpp"
#include "utility.hpp"

namespace stmesh {
// for locking: only lock writes, reads are alaways safe while no writes going on, iirc that primitive exists in pthread
template <typename ExtraData = std::monostate> class Triangulation {
  Eigen::AlignedBox<FLOAT_T, 4> bounding_box_;

  /*struct FullCellData {
    std::vector<FullCellHandle> *deleted;
    FullCellHandle self;
    FullCellData() : deleted(nullptr) {}
    FullCellData(std::vector<FullCellHandle> *deleted, FullCellHandle self) : deleted(deleted), self(self) {}
    FullCellData(const FullCellData &other) = default;
    FullCellData(FullCellData &&other) noexcept : deleted(other.deleted), self(other.self) {
      new (&other) FullCellData();
    }
    FullCellData &operator=(const FullCellData &other) {
      new (this) FullCellData(other);
      return *this;
    }
    FullCellData &operator=(FullCellData &&other) noexcept {
      new (this) FullCellData(std::move(other));
      return *this;
    }
    ~FullCellData() {
      if (deleted != nullptr)
        deleted->emplace_back(self);
    }
  };*/

  struct FullCellData {
    bool committed;
    // maybe needs copy constructor that inits to false? tbd
    ExtraData extra_data;
  };

  struct VertexData {
    bool nonfree_vertex;
    // maybe needs copy constructor that inits to false? tbd
  };

public:
  using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<4>>;
  using Point = Kernel::Point_d;

  using TriangulationFullCell = CGAL::Triangulation_full_cell<Kernel, FullCellData>;
  using TriangulationVertex = CGAL::Triangulation_vertex<Kernel, VertexData>;
  using TriangulationDataStructure =
      CGAL::Triangulation_data_structure<Kernel::Dimension, TriangulationVertex, TriangulationFullCell>;
  using DelaunayTriangulation = CGAL::Delaunay_triangulation<Kernel, TriangulationDataStructure>;

  using VertexHandle = DelaunayTriangulation::Vertex_handle;
  using FullCell = DelaunayTriangulation::Full_cell;
  using FullCellHandle = DelaunayTriangulation::Full_cell_handle;
  using FullCellConstHandle = DelaunayTriangulation::Full_cell_const_handle;
  using Facet = DelaunayTriangulation::Facet;
  using Face = DelaunayTriangulation::Face;
  using LocateType = DelaunayTriangulation::Locate_type;

  using BGPoint = bg::model::point<FLOAT_T, 4, bg::cs::cartesian>;
  using BGBox = bg::model::box<BGPoint>;

private:
  DelaunayTriangulation triangulation_;

  using Tree = bg::index::rtree<BGPoint,
                                // NOLINTNEXTLINE(*-magic-numbers)
                                bg::index::rstar<16>>;

  Tree tree_;

  std::unordered_map<Vector4F, VertexHandle, Vector4FHash> vertex_handle_map_;

public:
  [[nodiscard]] static BGPoint pointFromVector(const Vector4F &vector) noexcept;

  [[nodiscard]] static Vector4F vectorFromPoint(const BGPoint &point) noexcept;

  [[nodiscard]] static BGBox boxFromAABB(const Eigen::AlignedBox<FLOAT_T, 4> &aabb) noexcept;

  [[nodiscard]] static Vector4F pointToVec(const Point &pt);

  explicit Triangulation(const Eigen::AlignedBox<FLOAT_T, 4> &bounding_box);

  // only deletes conflict zone and inserts surrounding
  VertexHandle insert(const Vector4F &point, FullCellHandle hint = {}, bool nonfree_vertex = false);

  // deletes surrounding, inserts is complicated, use commit
  FullCellHandle remove(VertexHandle vertex);

  [[nodiscard]] std::vector<FullCellHandle> conflictZone(const Vector4F &point, FullCellHandle hint = {}) const;

  std::vector<FullCellHandle> commitUncommitted(FullCellHandle start);

  std::vector<FullCellHandle> surroundingFullCells(VertexHandle vertex, bool commit = true);

  [[nodiscard]] int mirrorIndex(const FullCellHandle full_cell, int index) const;

  [[nodiscard]] Eigen::Matrix<FLOAT_T, 4, 4> facetVertices(Facet facet) const noexcept;

  [[nodiscard]] Facet facetFromVertices(const Eigen::Matrix<FLOAT_T, 4, 4> &vertices) const;

  [[nodiscard]] std::tuple<std::tuple<VertexHandle, int, FullCellHandle>,
                           std::optional<std::tuple<VertexHandle, int, FullCellHandle>>>
  facetMirrorVertices(Facet facet) const noexcept;

  [[nodiscard]] size_t vertexCount() const noexcept;

  [[nodiscard]] std::vector<VertexHandle> verticesInRadius(const Vector4F &point, FLOAT_T radius) const;

  [[nodiscard]] static GeometricSimplex<4> fullCellSimplex(FullCellConstHandle full_cell) noexcept;

  [[nodiscard]] bool isGoodPoint(VertexHandle vertex, FLOAT_T rho_bar, FLOAT_T tau_bar, FLOAT_T max_radius) const;

  [[nodiscard]] const Eigen::AlignedBox<FLOAT_T, 4> &boundingBox() const noexcept;

  [[nodiscard]] DelaunayTriangulation::Finite_full_cell_const_iterator cbegin() const;

  [[nodiscard]] DelaunayTriangulation::Finite_full_cell_const_iterator cend() const;

  [[nodiscard]] DelaunayTriangulation::Finite_full_cell_iterator begin();

  [[nodiscard]] DelaunayTriangulation::Finite_full_cell_iterator end();

  [[nodiscard]] DelaunayTriangulation::Finite_full_cell_const_iterator begin() const;

  [[nodiscard]] DelaunayTriangulation::Finite_full_cell_const_iterator end() const;
};
} // namespace stmesh
#endif