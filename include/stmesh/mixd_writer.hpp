#ifndef STMESH_MIXD_WRITER_HPP
#define STMESH_MIXD_WRITER_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <numeric>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <spdlog/spdlog.h>

#include "boundary_region_manager.hpp" // IWYU pragma: keep
#include "geometric_simplex.hpp"
#include "sdf.hpp"
#include "surface_adapters.hpp" // IWYU pragma: keep
#include "triangulation.hpp"
#include "utility.hpp"

namespace stmesh {
namespace detail {
void writeMinf(const std::filesystem::path &minf_file, const std::filesystem::path &mxyz_file,
               const std::filesystem::path &mien_file, const std::filesystem::path &mrng_file, size_t number_elements,
               size_t number_nodes);

std::filesystem::path writeMxyz(std::filesystem::path file, const std::vector<Vector4F> &vertices,
                                stmesh::FLOAT_T scale, stmesh::FLOAT_T min_time);

std::filesystem::path writeIntMixd(std::filesystem::path file, std::string_view extension,
                                   const std::vector<std::array<size_t, 5>> &full_cell_vertex_ids);

bool positivePentatopeElementDet(const std::array<size_t, 5> &vertex_ids, const std::vector<Vector4F> &vertices);
} // namespace detail

/// Write a MIXD file from a triangulation
/**
 * Writes a MIXD file from a triangulation. The MIXD file format is a file format used by XNS to represent a mesh.
 * The format consists of:
 * - A .mxyz file containing the coordinates of the vertices of the mesh.
 * - A .mien file containing the vertex indices of the full cells of the mesh.
 * - A .mrng file containing the boundary region indices of the faces of the mesh.
 * - A .minf file containing metadata about the mesh.
 * Only full cells whose circumsphere is inside the surface are written to the MIXD file.
 * Boundary regions are determined by the boundary_region_manager.
 *
 * @param file The path to the .minf file to write. The other MIXD files will be placed alongside it.
 * @param surface The surface adapter to use to determine which full cells to write to the MIXD file.
 * @param triangulation The triangulation to write to the MIXD file.
 * @param boundary_region_manager The boundary region manager to use to determine the boundary regions of the faces of
 * the mesh.
 * @param scale The scale for the mxyz file. Defaults to 1.
 * @param min_time The minimum time of the triangulation for offsetting the mxyz file. Defaults to 0.
 * @tparam ExtraData The type of extra data stored in the triangulation.
 */
template <typename ExtraData>
void writeMixd([[maybe_unused]] const std::filesystem::path &file, [[maybe_unused]] const SurfaceAdapter4 auto &surface,
               const Triangulation<ExtraData> &triangulation, const BoundaryRegionManager auto &boundary_region_manager,
               stmesh::FLOAT_T scale = 1, stmesh::FLOAT_T min_time = 0) {
  std::unordered_map<Vector4F, size_t, Vector4FHash> vertex_map;
  std::vector<Vector4F> vertices;
  std::vector<std::array<size_t, 5>> full_cell_vertex_ids;
  std::vector<std::array<size_t, 5>> full_cell_face_ids;
  for (const auto &cell : triangulation) {
    GeometricSimplex<4> simplex =
        triangulation.fullCellSimplex(typename Triangulation<ExtraData>::FullCellConstHandle{&cell});
    HyperSphere4 circumsphere = simplex.circumsphere();
    if (surface.inside(circumsphere.center())) {
      std::array<size_t, 5> vertex_ids{};
      for (Eigen::Index i = 0; i < 5; ++i) {
        auto [it, inserted] = vertex_map.try_emplace(simplex.vertices().col(i), vertices.size());
        if (inserted)
          vertices.emplace_back(simplex.vertices().col(i));
        vertex_ids.at(static_cast<size_t>(i)) = it->second;
      }
      if (!detail::positivePentatopeElementDet(vertex_ids, vertices)) {
        spdlog::warn("Pentatope has negative determinant! Mesh will be tangled.");
        break;
      }
      std::ranges::transform(vertex_ids, vertex_ids.begin(), [&](size_t idx) { return idx + 1; });
      full_cell_vertex_ids.emplace_back(vertex_ids);
      std::array<size_t, 5> face_boundary_ids{};
      std::array<size_t, 5> vertex_boundary_ids{};
      for (size_t j = 0; j < 5; ++j)
        vertex_boundary_ids.at(j) =
            boundary_region_manager.findBoundaryRegion(simplex.vertices().col(static_cast<Eigen::Index>(j)));
      std::array<size_t, 5> vertex_boundary_idxs{};
      std::iota(vertex_boundary_idxs.begin(), vertex_boundary_idxs.end(), 0);
      std::partial_sort(
          vertex_boundary_idxs.begin(), vertex_boundary_idxs.begin() + 2, vertex_boundary_idxs.end(),
          [&](size_t lhs, size_t rhs) { return vertex_boundary_ids.at(lhs) > vertex_boundary_ids.at(rhs); });
      constexpr std::array<size_t, 5> kToOmit{4, 3, 2, 0, 1};
      // first two elements are now first two largest values
      for (size_t i = 0; i < 5; ++i) {
        size_t j = kToOmit.at(i);
        face_boundary_ids.at(i) =
            surface.inside(triangulation.fullCellSimplex(cell.neighbor(j)).circumsphere().center()) &&
                    !triangulation.isInfinite(cell.neighbor(j))
                ? 0
            : j != vertex_boundary_idxs[0] ? vertex_boundary_ids.at(vertex_boundary_idxs[0])
                                           : vertex_boundary_ids.at(vertex_boundary_idxs[1]);
      }
      full_cell_face_ids.emplace_back(face_boundary_ids);
    }
  }
  std::filesystem::path mxyz_file = detail::writeMxyz(file, vertices, scale, min_time);
  std::filesystem::path mien_file = detail::writeIntMixd(file, ".mien", full_cell_vertex_ids);
  std::filesystem::path mrng_file = detail::writeIntMixd(file, ".mrng", full_cell_face_ids);
  detail::writeMinf(file, mxyz_file, mien_file, mrng_file, full_cell_vertex_ids.size(), vertices.size());
}
} // namespace stmesh
#endif