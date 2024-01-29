// NOLINTBEGIN(misc-include-cleaner)
#include <algorithm>
#include <stmesh/meshing_cell.hpp>
#include <stmesh/triangulation.hpp>

#include <array>
#include <iterator>
#include <tuple>

#include <boost/heap/fibonacci_heap.hpp>

namespace stmesh {
template <typename ExtraData> [[nodiscard]] Vector4F Triangulation<ExtraData>::pointToVec(const Point &pt) {
  return {static_cast<FLOAT_T>(pt[0]), static_cast<FLOAT_T>(pt[1]), static_cast<FLOAT_T>(pt[2]),
          static_cast<FLOAT_T>(pt[3])};
}

template <typename ExtraData>
Triangulation<ExtraData>::Triangulation(const Eigen::AlignedBox<FLOAT_T, 4> &bounding_box)
    : bounding_box_(bounding_box), triangulation_(4) {
  const std::array corners = allCorners(bounding_box);
  for (const auto &corner : corners)
    insert(corner, {}, true);
  for (auto &full_cell : *this)
    full_cell.data().committed = true;
}

// only deletes conflict zone and inserts surrounding
template <typename ExtraData>
auto Triangulation<ExtraData>::insert(const Vector4F &point, FullCellHandle hint, bool nonfree_vertex) -> VertexHandle {
  const Point p(point.begin(), point.end());
  // valid since iterators are not invalidated by insert
  // TODO: using hint from kdtree, should test with and without
  // then also wouldn't need closestPoint anymore
  (vertex_handle_map[point] = triangulation_.insert(p, hint))->data().nonfree_vertex = nonfree_vertex;
  tree_.insert(p);
  return vertex_handle_map.at(point);
}

// deletes surrounding, inserts is complicated, use commit
template <typename ExtraData> auto Triangulation<ExtraData>::remove(VertexHandle vertex) -> FullCellHandle {
  vertex_handle_map.erase(pointToVec(vertex->point()));
  tree_.remove(vertex->point());
  return triangulation_.remove(vertex);
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::conflictZone(const Vector4F &point, FullCellHandle hint) const
    -> std::vector<FullCellHandle> {
  std::vector<FullCellHandle> cells;
  const Point p(point.begin(), point.end());
  LocateType loc_type{};
  Face f(3);
  Facet ft;
  triangulation_.compute_conflict_zone(p, triangulation_.locate(p, loc_type, f, ft, hint), std::back_inserter(cells));
  std::vector<FullCellHandle> result;
  result.reserve(cells.size());
  std::copy_if(cells.begin(), cells.end(), std::back_inserter(result),
               [&](const auto &cell) { return !triangulation_.is_infinite(cell); });
  return result;
}

template <typename ExtraData>
auto Triangulation<ExtraData>::commitUncommitted(FullCellHandle start) -> std::vector<FullCellHandle> {
  std::vector<FullCellHandle> result;
  auto predicate = [&](const Facet &facet) {
    const FullCellHandle full_cell = facet.full_cell()->neighbor(facet.index_of_covertex());
    return !triangulation_.is_infinite(full_cell) && !full_cell->data().committed;
  };
  std::back_insert_iterator inserter(result);
  triangulation_.tds().gather_full_cells(start, predicate, inserter);
  for (const auto &full_cell : result)
    full_cell->data().committed = true;
  return result;
}

template <typename ExtraData>
auto Triangulation<ExtraData>::surroundingFullCells(VertexHandle vertex, bool commit) -> std::vector<FullCellHandle> {
  std::vector<FullCellHandle> cells;
  triangulation_.incident_full_cells(vertex, std::back_inserter(cells));
  std::vector<FullCellHandle> result;
  result.reserve(cells.size());
  std::copy_if(cells.begin(), cells.end(), std::back_inserter(result),
               [&](const auto &cell) { return !triangulation_.is_infinite(cell); });
  if (commit) {
    for (const auto &full_cell : result)
      full_cell->data().committed = true;
  }
  return result;
}

template <typename ExtraData>
[[nodiscard]] int Triangulation<ExtraData>::mirrorIndex(const FullCellHandle full_cell, int index) const {
  return triangulation_.tds().mirror_index(full_cell, index);
}

template <typename ExtraData>
[[nodiscard]] Eigen::Matrix<FLOAT_T, 4, 4> Triangulation<ExtraData>::facetVertices(Facet facet) const noexcept {
  Eigen::Matrix<FLOAT_T, 4, 4> vertices;
  const FullCellHandle cell = triangulation_.tds().full_cell(facet);
  int j = 0;
  // NOLINTNEXTLINE(*-magic-numbers)
  for (int i = 0; i < 5; ++i) {
    if (i == triangulation_.tds().index_of_covertex(facet))
      continue;
    const auto &vertex = triangulation_.tds().vertex(cell, i)->point();
    vertices.col(j++) = Vector4F{static_cast<FLOAT_T>(vertex[0]), static_cast<FLOAT_T>(vertex[1]),
                                 static_cast<FLOAT_T>(vertex[2]), static_cast<FLOAT_T>(vertex[3])};
  }
  return vertices;
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::facetFromVertices(const Eigen::Matrix<FLOAT_T, 4, 4> &vertices) const
    -> Facet {
  std::array<VertexHandle, 4> vertex_handles;
  auto *iter = vertex_handles.begin();
  for (const auto &vertex : vertices.colwise())
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    *iter++ = vertex_handle_map.at(vertex);
  std::vector<FullCellHandle> possible_full_cells;
  triangulation_.incident_full_cells(vertex_handles[0], std::back_inserter(possible_full_cells));
  for (const auto &full_cell : possible_full_cells) {
    Facet ret;
    bool done = true;
    // NOLINTNEXTLINE(*-magic-numbers)
    for (int i = 0; i < 5; ++i) {
      const Vector4F vertex = pointToVec(triangulation_.tds().vertex(full_cell, i)->point());
      if (!(vertices.array() == vertex.array().replicate<1, 4>()).colwise().all().any()) {
        if (ret.full_cell() == FullCellHandle())
          ret = {full_cell, i};
        else {
          done = false;
          break;
        }
      }
    }
    if (done)
      return ret;
  }
  throw std::runtime_error("Triangulation::facetFromVertices: no facet found");
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::facetMirrorVertices(Facet facet) const noexcept
    -> std::tuple<std::tuple<VertexHandle, int, FullCellHandle>,
                  std::optional<std::tuple<VertexHandle, int, FullCellHandle>>> {
  const FullCellHandle &full_cell = facet.full_cell();
  const int &covertexIndex = facet.index_of_covertex();
  const VertexHandle covertex = full_cell->vertex(covertexIndex);
  const FullCellHandle neighbor_full_cell = full_cell->neighbor(covertexIndex);
  const int mirror_vertex = triangulation_.tds().mirror_index(full_cell, covertexIndex);
  const VertexHandle neighbor_mirror_vertex = neighbor_full_cell->vertex(mirror_vertex);
  if (triangulation_.is_infinite(covertex))
    return {{neighbor_mirror_vertex, mirror_vertex, neighbor_full_cell}, std::nullopt};
  if (triangulation_.is_infinite(neighbor_mirror_vertex))
    return {{covertex, covertexIndex, full_cell}, std::nullopt};
  return {{covertex, covertexIndex, full_cell}, std::tuple{neighbor_mirror_vertex, mirror_vertex, neighbor_full_cell}};
}

template <typename ExtraData> [[nodiscard]] size_t Triangulation<ExtraData>::vertexCount() const noexcept {
  return triangulation_.number_of_vertices();
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::verticesInRadius(const Vector4F &point, FLOAT_T radius) const
    -> std::vector<VertexHandle> {
  const Point p(point.begin(), point.end());
  std::vector<Point> points;
  tree_.search(std::back_inserter(points), FuzzySphere(p, radius));
  std::vector<VertexHandle> result;
  result.reserve(points.size());
  std::ranges::transform(points, std::back_inserter(result),
                         [&](const auto &pt) { return vertex_handle_map.at(pointToVec(pt)); });
  return result;
}

template <typename ExtraData>
[[nodiscard]] GeometricSimplex<4> Triangulation<ExtraData>::fullCellSimplex(FullCellConstHandle full_cell) noexcept {
  // NOLINTBEGIN(*-magic-numbers)
  Eigen::Matrix<FLOAT_T, 4, 5> vertices;
  for (int i = 0; i < 5; ++i)
    // NOLINTEND(*-magic-numbers)
    vertices.col(i) = pointToVec(full_cell->vertex(i)->point());
  return GeometricSimplex<4>{vertices};
}

template <typename ExtraData>
[[nodiscard]] bool Triangulation<ExtraData>::isGoodPoint(VertexHandle vertex, FLOAT_T rho_bar, FLOAT_T tau_bar,
                                                         FLOAT_T max_radius) const {
  // only cells themselves could be bad or 3d facets, incident full cells and faces should be enough
  std::vector<FullCellHandle> incident_full_cells;
  triangulation_.incident_full_cells(vertex, std::back_inserter(incident_full_cells));
  if (std::ranges::any_of(incident_full_cells, [&](const auto &full_cell) {
        return !triangulation_.is_infinite(full_cell) &&
               fullCellSimplex(full_cell).smallSliverSimplex(rho_bar, tau_bar, max_radius);
      }))
    return false;
  std::vector<Face> incident_faces;
  triangulation_.incident_faces(vertex, 3, std::back_inserter(incident_faces));
  for (const auto &face : incident_faces) {
    if (triangulation_.is_infinite(face))
      continue;
    Eigen::Matrix<FLOAT_T, 4, 4> vertices;
    for (int i = 0; i < 4; ++i)
      vertices.col(i) = pointToVec(face.vertex(i)->point());
    const GeometricSimplex<4, 4> simplex(vertices);
    if (simplex.smallSliverSimplex(rho_bar, tau_bar, max_radius))
      return false;
  }
  return true;
}

template <typename ExtraData>
[[nodiscard]] const Eigen::AlignedBox<FLOAT_T, 4> &Triangulation<ExtraData>::boundingBox() const noexcept {
  return bounding_box_;
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::cbegin() const -> DelaunayTriangulation::Finite_full_cell_const_iterator {
  return triangulation_.finite_full_cells_begin();
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::cend() const -> DelaunayTriangulation::Finite_full_cell_const_iterator {
  return triangulation_.finite_full_cells_end();
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::begin() -> DelaunayTriangulation::Finite_full_cell_iterator {
  return triangulation_.finite_full_cells_begin();
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::end() -> DelaunayTriangulation::Finite_full_cell_iterator {
  return triangulation_.finite_full_cells_end();
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::begin() const -> DelaunayTriangulation::Finite_full_cell_const_iterator {
  return triangulation_.finite_full_cells_begin();
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::end() const -> DelaunayTriangulation::Finite_full_cell_const_iterator {
  return triangulation_.finite_full_cells_end();
}

template class Triangulation<detail::ExtraData>;
template class Triangulation<std::monostate>;
} // namespace stmesh

// NOLINTEND(misc-include-cleaner)