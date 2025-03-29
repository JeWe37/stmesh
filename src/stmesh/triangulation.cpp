#include "stmesh/triangulation.hpp"

#include <Eigen/src/Core/util/Constants.h>
#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <iterator>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <boost/geometry/index/predicates.hpp>
#include <boost/iterator/function_output_iterator.hpp>

#include "stmesh/geometric_simplex.hpp"
#include "stmesh/geometry_helpers.hpp"
#include "stmesh/meshing_cell.hpp"
#include "stmesh/mixd.hpp"
#include "stmesh/problem_types.hpp"
#include "stmesh/sdf.hpp"
#include "stmesh/utility.hpp"

namespace stmesh {
template <typename Cell>
detail::MixdCell<Cell>::MixdCell(const TriangulationFromMixdGeneric<CellType> *triangulation, size_t index)
    : triangulation_(triangulation), index_(index) {}

template <typename Cell> [[nodiscard]] GeometricSimplex<4> detail::MixdCell<Cell>::geometricSimplex() const noexcept {
  Eigen::Matrix<FLOAT_T, 4, 5> simplex;
  for (size_t i = 0; i < 5; ++i)
    simplex.col(static_cast<Eigen::Index>(i)) =
        triangulation_->mxyz_[static_cast<size_t>(triangulation_->mien_[index_].at(i)) - 1];
  return GeometricSimplex<4>(simplex);
}

template <typename Cell> [[nodiscard]] bool detail::MixdCell<Cell>::isSurfaceSide(size_t i) const noexcept {
  return triangulation_->mrng_[index_].at(mixd::kOmitted.at(i)) > 0;
}

template class detail::MixdCell<>;
template class detail::MixdCell<detail::MixdDataCell>;

[[nodiscard]] Eigen::Matrix<FLOAT_T, Eigen::Dynamic, 5> detail::MixdDataCell::data() const noexcept {
  Eigen::Matrix<FLOAT_T, Eigen::Dynamic, 5> ret(
      // NOLINTNEXTLINE(cppcoreguidelines-pro-type-static-cast-downcast)
      static_cast<const TriangulationFromMixdWithData *>(triangulation_)->problem_type_.entries(), 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < static_cast<size_t>(ret.rows()); ++j)
      ret(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) =
          // NOLINTNEXTLINE(cppcoreguidelines-pro-type-static-cast-downcast)
          static_cast<const TriangulationFromMixdWithData *>(triangulation_)
              ->data_.at(static_cast<size_t>(triangulation_->mien_[index_].at(i) - 1))[j];
  }
  return ret;
}

template <typename MixdCell>
auto TriangulationFromMixdGeneric<MixdCell>::iterator::dereference() const noexcept -> MixdCell {
  return {triangulation_, static_cast<size_t>(&*this->base_reference() - triangulation_->mien_.data())};
}

template <typename MixdCell>
TriangulationFromMixdGeneric<MixdCell>::iterator::iterator(const TriangulationFromMixdGeneric *triangulation,
                                                           const std::vector<std::array<int, 5>>::const_iterator &it)
    : iterator::iterator_adaptor_(it), triangulation_(triangulation) {}

template <typename MixdCell> TriangulationFromMixdGeneric<MixdCell>::iterator::iterator() : triangulation_(nullptr) {}

TriangulationFromMixdWithData::TriangulationFromMixdWithData(const std::filesystem::path &minf_file,
                                                             const ProblemType &problem_type,
                                                             const std::filesystem::path &data_file)
    : TriangulationFromMixdGeneric<detail::MixdDataCell>(minf_file),
      data_(mixd::readData(data_file, problem_type.entries())), problem_type_(problem_type) {}

template <typename MixdCell>
TriangulationFromMixdGeneric<MixdCell>::TriangulationFromMixdGeneric(const std::filesystem::path &minf_file) {
  auto [mxyz_file, mien_file, mrng_file] = mixd::readMinf(minf_file);
  mxyz_ = mixd::readMxyz(mxyz_file);
  mien_ = mixd::readIntMixd(mien_file);
  mrng_ = mixd::readIntMixd(mrng_file);
}

template <typename MixdCell>
[[nodiscard]] Eigen::AlignedBox<double, 4> TriangulationFromMixdGeneric<MixdCell>::boundingBox() const noexcept {
  Eigen::AlignedBox<double, 4> box(mxyz_[0], mxyz_[0]);
  for (const auto &vertex : mxyz_)
    box.extend(vertex.cast<double>());
  return box;
}

template <typename MixdCell>
[[nodiscard]] auto TriangulationFromMixdGeneric<MixdCell>::begin() const noexcept -> iterator {
  return {this, mien_.begin()};
}

template <typename MixdCell>
[[nodiscard]] auto TriangulationFromMixdGeneric<MixdCell>::end() const noexcept -> iterator {
  return {this, mien_.end()};
}

template <typename MixdCell>
[[nodiscard]] size_t TriangulationFromMixdGeneric<MixdCell>::findBoundaryRegion(const MixdCell &cell,
                                                                                size_t j) const noexcept {
  return static_cast<size_t>(mrng_[cell.index_].at(mixd::kOmitted.at(j)));
}

template class TriangulationFromMixdGeneric<>;
template class TriangulationFromMixdGeneric<detail::MixdDataCell>;

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
  const BGPoint bg_point = pointFromVector(point);
  // valid since iterators are not invalidated by insert
  (vertex_handle_map_[point] = triangulation_.insert(p, hint))->data().nonfree_vertex = nonfree_vertex;
  tree_.insert(bg_point);
  return vertex_handle_map_.at(point);
}

// deletes surrounding, inserts is complicated, use commit
template <typename ExtraData> auto Triangulation<ExtraData>::remove(VertexHandle vertex) -> FullCellHandle {
  const Vector4F point = pointToVec(vertex->point());
  vertex_handle_map_.erase(point);
  tree_.remove(pointFromVector(point));
  return triangulation_.remove(vertex);
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::conflictZone(const Vector4F &point, FullCellHandle hint) const
    -> std::vector<FullCellHandle> {
  std::vector<FullCellHandle> result;
  const Point p(point.begin(), point.end());
  LocateType loc_type{};
  Face f(3);
  Facet ft;
  triangulation_.compute_conflict_zone(p, triangulation_.locate(p, loc_type, f, ft, hint),
                                       boost::make_function_output_iterator([&](const auto &cell) {
                                         if (!triangulation_.is_infinite(cell))
                                           result.push_back(cell);
                                       }));
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
  std::vector<FullCellHandle> result;
  triangulation_.incident_full_cells(vertex, boost::make_function_output_iterator([&](const auto &cell) {
                                       if (!triangulation_.is_infinite(cell)) {
                                         result.push_back(cell);
                                         if (commit)
                                           cell->data().committed = true;
                                       }
                                     }));
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
    *iter++ = vertex_handle_map_.at(vertex);
  std::vector<FullCellHandle> possible_full_cells;
  triangulation_.incident_full_cells(vertex_handles[0], std::back_inserter(possible_full_cells));
  for (const auto &full_cell : possible_full_cells) {
    Facet ret;
    bool done = true;
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
[[nodiscard]] auto Triangulation<ExtraData>::facetMirrorVertices(Facet facet) const
    -> std::tuple<std::tuple<VertexHandle, int, FullCellHandle>,
                  std::optional<std::tuple<VertexHandle, int, FullCellHandle>>> {
  const FullCellHandle &full_cell = facet.full_cell();
  const int &covertex_index = facet.index_of_covertex();
  const VertexHandle covertex = full_cell->vertex(covertex_index);
  const FullCellHandle neighbor_full_cell = full_cell->neighbor(covertex_index);
  const int mirror_vertex = triangulation_.tds().mirror_index(full_cell, covertex_index);
  const VertexHandle neighbor_mirror_vertex = neighbor_full_cell->vertex(mirror_vertex);
  if (triangulation_.is_infinite(covertex))
    return {{neighbor_mirror_vertex, mirror_vertex, neighbor_full_cell}, std::nullopt};
  if (triangulation_.is_infinite(neighbor_mirror_vertex))
    return {{covertex, covertex_index, full_cell}, std::nullopt};
  return {{covertex, covertex_index, full_cell}, std::tuple{neighbor_mirror_vertex, mirror_vertex, neighbor_full_cell}};
}

template <typename ExtraData> [[nodiscard]] size_t Triangulation<ExtraData>::vertexCount() const noexcept {
  return triangulation_.number_of_vertices();
}

template <typename ExtraData>
[[nodiscard]] auto Triangulation<ExtraData>::verticesInRadius(const Vector4F &point, FLOAT_T radius) const
    -> std::vector<VertexHandle> {
  std::vector<VertexHandle> result;
  const HyperSphere4 sphere(radius, point);
  tree_.query(bg::index::intersects(boxFromAABB(sphere.boundingBox())),
              boost::make_function_output_iterator([&](const auto &pt) {
                const Vector4F p = vectorFromPoint(pt);
                if (sphere.signedDistance(p) < 0)
                  result.push_back(vertex_handle_map_.at(p));
              }));
  return result;
}

template <typename ExtraData>
[[nodiscard]] GeometricSimplex<4> Triangulation<ExtraData>::fullCellSimplex(FullCellConstHandle full_cell) noexcept {
  Eigen::Matrix<FLOAT_T, 4, 5> vertices;
  for (int i = 0; i < 5; ++i)
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
[[nodiscard]] bool Triangulation<ExtraData>::isInfinite(FullCellConstHandle full_cell) const {
  return triangulation_.is_infinite(full_cell);
}

template <typename ExtraData>
[[nodiscard]] const Eigen::AlignedBox<FLOAT_T, 4> &Triangulation<ExtraData>::boundingBox() const noexcept {
  return bounding_box_;
}

template <typename ExtraData> [[nodiscard]] auto Triangulation<ExtraData>::cbegin() const -> const_iterator {
  return triangulation_.finite_full_cells_begin();
}

template <typename ExtraData> [[nodiscard]] auto Triangulation<ExtraData>::cend() const -> const_iterator {
  return triangulation_.finite_full_cells_end();
}

template <typename ExtraData> [[nodiscard]] auto Triangulation<ExtraData>::begin() -> iterator {
  return triangulation_.finite_full_cells_begin();
}

template <typename ExtraData> [[nodiscard]] auto Triangulation<ExtraData>::end() -> iterator {
  return triangulation_.finite_full_cells_end();
}

template <typename ExtraData> [[nodiscard]] auto Triangulation<ExtraData>::begin() const -> const_iterator {
  return triangulation_.finite_full_cells_begin();
}

template <typename ExtraData> [[nodiscard]] auto Triangulation<ExtraData>::end() const -> const_iterator {
  return triangulation_.finite_full_cells_end();
}

template class Triangulation<detail::ExtraData>;
template class Triangulation<std::monostate>;
} // namespace stmesh
