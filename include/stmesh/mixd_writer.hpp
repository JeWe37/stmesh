#ifndef STMESH_MIXD_WRITER_HPP
#define STMESH_MIXD_WRITER_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <unordered_map>
#include <vector>

#include <spdlog/spdlog.h>

#include "boundary_region_manager.hpp" // IWYU pragma: keep
#include "geometric_simplex.hpp"
#include "mixd.hpp"
#include "sdf.hpp"
#include "surface_adapters.hpp" // IWYU pragma: keep
#include "utility.hpp"
#include "writable_triangulation.hpp"

namespace stmesh {
/// Write a MIXD file from a triangulation
/**
 * Writes a MIXD file from a triangulation. The MIXD file format is a file format used by XNS to represent a mesh.
 * The format consists of:
 * - A .mxyz file containing the coordinates of the vertices of the mesh.
 * - A .mien file containing the vertex indices of the full cells of the mesh.
 * - A .mrng file containing the boundary region indices of the faces of the mesh.
 * - A .minf file containing metadata about the mesh.
 * Which cells are to be written is specified by the triangulation.
 * Boundary regions are determined by the boundary_region_manager.
 *
 * @param file The path to the .minf file to write. The other MIXD files will be placed alongside it.
 * @param surface The surface adapter to use to determine which full cells to write to the MIXD file.
 * @param triangulation The writable triangulation to write to the MIXD file.
 * @param boundary_region_manager The boundary region manager to use to determine the boundary regions of the faces of
 * the mesh.
 * @param scale The scale for the mxyz file. Defaults to 1.
 * @param min_time The minimum time of the triangulation for offsetting the mxyz file. Defaults to 0.
 */
void writeMixd([[maybe_unused]] const std::filesystem::path &file, [[maybe_unused]] const SurfaceAdapter4 auto &surface,
               const WritableTriangulation auto &triangulation,
               const BoundaryRegionManager auto &boundary_region_manager, stmesh::FLOAT_T scale = 1,
               stmesh::FLOAT_T min_time = 0) {
  std::unordered_map<Vector4F, int, Vector4FHash> vertex_map;
  std::vector<Vector4F> vertices;
  std::vector<std::array<int, 5>> full_cell_vertex_ids;
  std::vector<std::array<int, 5>> full_cell_face_ids;
  for (const auto &cell : triangulation) {
    GeometricSimplex<4> simplex = cell.geometricSimplex();
    HyperSphere4 circumsphere = simplex.circumsphere();
    if (surface.inside(circumsphere.center())) {
      std::array<int, 5> vertex_ids{};
      for (Eigen::Index i = 0; i < 5; ++i) {
        auto [it, inserted] = vertex_map.try_emplace(simplex.vertices().col(i), vertices.size());
        if (inserted)
          vertices.emplace_back(simplex.vertices().col(i));
        vertex_ids.at(static_cast<size_t>(i)) = it->second;
      }
      if (!mixd::positivePentatopeElementDet(vertex_ids, vertices)) {
        spdlog::warn("Pentatope has negative determinant! Mesh will be tangled.");
        break;
      }
      std::ranges::transform(vertex_ids, vertex_ids.begin(), [&](size_t idx) { return idx + 1; });
      full_cell_vertex_ids.emplace_back(vertex_ids);
      std::array<int, 5> &face_boundary_ids = full_cell_face_ids.emplace_back();
      for (size_t i = 0; i < 5; ++i) {
        size_t j = mixd::kToOmit.at(i);
        face_boundary_ids.at(i) = !cell.isSurfaceSide(j) ? 0 : boundary_region_manager.findBoundaryRegion(cell, j);
      }
    }
  }
  std::filesystem::path mxyz_file = mixd::writeMxyz(file, vertices, scale, min_time);
  std::filesystem::path mien_file = mixd::writeIntMixd(file, ".mien", full_cell_vertex_ids);
  std::filesystem::path mrng_file = mixd::writeIntMixd(file, ".mrng", full_cell_face_ids);
  mixd::writeMinf(file, mxyz_file, mien_file, mrng_file, full_cell_vertex_ids.size(), vertices.size());
}
} // namespace stmesh
#endif