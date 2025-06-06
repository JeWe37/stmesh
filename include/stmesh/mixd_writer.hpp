#ifndef STMESH_MIXD_WRITER_HPP
#define STMESH_MIXD_WRITER_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <ranges>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/container_hash/hash.hpp>
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
 * @param write_ideal Whether to write the closest boundary positions of the vertices to a file with a .ideal.minf.
 * Defaults to false
 * @param write_neim Whether to write the node element index mapping. Defaults to false
 * @param compute_dual Whether to include the dual in the mrng file. Defaults to false
 */
void writeMixd([[maybe_unused]] const std::filesystem::path &file, [[maybe_unused]] const SurfaceAdapter4 auto &surface,
               const WritableTriangulation auto &triangulation,
               const BoundaryRegionManager auto &boundary_region_manager, stmesh::FLOAT_T scale = 1,
               stmesh::FLOAT_T min_time = 0, bool write_ideal = false, bool write_neim = false,
               bool compute_dual = false) {
  constexpr static FLOAT_T kNaN = std::numeric_limits<FLOAT_T>::quiet_NaN();

  std::unordered_map<Vector4F, int, Vector4FHash> vertex_map;
  std::vector<Vector4F> vertices;
  std::vector<Vector4F> ideal_vertices;
  std::set<size_t> boundary_vertices;

  std::vector<std::array<int, 5>> full_cell_vertex_ids;
  std::vector<std::array<int, 5>> full_cell_face_ids;
  std::vector<std::vector<int>> neim;
  using Face = std::array<int, 4>;
  std::unordered_map<Face, int, boost::hash<Face>> face_map;
  auto compute_missing = [](const Face &face, const std::array<int, 5> &cell) {
    int missing_val =
        static_cast<int>(static_cast<unsigned>(std::accumulate(face.begin(), face.end(), 0, std::bit_xor<>{})) ^
                         static_cast<unsigned>(std::accumulate(cell.begin(), cell.end(), 0, std::bit_xor<>{})));
    return mixd::kOmitted.at(static_cast<size_t>(std::distance(cell.begin(), std::ranges::find(cell, missing_val))));
  };
  auto add_face = [&](Face &face, const int cell) {
    std::ranges::sort(face);
    if (auto [it, inserted] = face_map.try_emplace(face, cell); !inserted) {
      size_t missing_this = compute_missing(face, full_cell_vertex_ids[static_cast<size_t>(cell - 1)]);
      size_t missing_other = compute_missing(face, full_cell_vertex_ids[static_cast<size_t>(it->second - 1)]);
      full_cell_face_ids[static_cast<size_t>(cell - 1)].at(missing_this) = -it->second;
      full_cell_face_ids[static_cast<size_t>(it->second - 1)].at(missing_other) = -cell;
      face_map.erase(it);
    }
  };
  std::array<Face, 5> faces{};
  for (const auto &cell : triangulation) {
    GeometricSimplex<4> simplex = cell.geometricSimplex();
    HyperSphere4 circumsphere = simplex.circumsphere();
    if (surface.inside(circumsphere.center())) {
      std::array<int, 5> vertex_ids{};
      for (Eigen::Index i = 0; i < 5; ++i) {
        auto [it, inserted] = vertex_map.try_emplace(simplex.vertices().col(i), vertices.size());
        if (inserted) {
          vertices.emplace_back(simplex.vertices().col(i));
          if (write_ideal)
            ideal_vertices.emplace_back(surface.closestPoint(vertices.back()));
        }
        vertex_ids.at(static_cast<size_t>(i)) = it->second;
      }
      if (!mixd::positivePentatopeElementDet(vertex_ids, vertices)) {
        spdlog::warn("Pentatope has negative determinant! Mesh will be tangled.");
        break;
      }
      std::array<int, 5> &face_boundary_ids = full_cell_face_ids.emplace_back();
      size_t num_faces = 0;
      for (size_t i = 0; i < 5; ++i) {
        if (cell.isSurfaceSide(i)) {
          for (size_t j = 0; j < 5; ++j) {
            if (j != i)
              boundary_vertices.insert(static_cast<size_t>(vertex_ids.at(j)));
          }
        }
        size_t j = mixd::kToOmit.at(i);
        if (cell.isSurfaceSide(j))
          face_boundary_ids.at(i) = boundary_region_manager.findBoundaryRegion(cell, j);
        else if (compute_dual) {
          Face &face = faces.at(num_faces++);
          size_t l = 0;
          for (size_t k = 0; k < 5; ++k) {
            if (k != j)
              face.at(l++) = vertex_ids.at(k) + 1;
          }
        } else
          face_boundary_ids.at(i) = 0;
      }
      if (write_neim) {
        for (int vertex_id : vertex_ids) {
          if (static_cast<size_t>(vertex_id) >= neim.size())
            neim.resize(static_cast<size_t>(vertex_id) + 1);
          neim[static_cast<size_t>(vertex_id)].emplace_back(full_cell_vertex_ids.size() + 1);
        }
      }
      std::ranges::transform(vertex_ids, vertex_ids.begin(), [&](size_t idx) { return idx + 1; });
      full_cell_vertex_ids.emplace_back(vertex_ids);
      if (compute_dual) {
        for (Face &face : faces | std::views::take(num_faces))
          add_face(face, full_cell_vertex_ids.size());
      }
    }
  }
  if (!face_map.empty())
    throw std::runtime_error("Surface is not correctly declared.");
  std::filesystem::path mxyz_file = mixd::writeMxyz(file, vertices, scale, min_time);
  if (write_ideal) {
    size_t start = 0;
    for (size_t boundary_vertex : boundary_vertices) {
      for (; start < boundary_vertex; ++start)
        ideal_vertices[start] = Vector4F::Constant(kNaN);
      start++;
    }
    std::filesystem::path ideal_file = file;
    ideal_file.replace_extension(".ideal.minf");
    mixd::writeMxyz(ideal_file, vertices, scale, min_time);
  }
  if (write_neim) {
    std::filesystem::path neim_file = file;
    neim_file.replace_extension(".neim");
    mixd::writeNeim(neim_file, neim);
  }
  std::filesystem::path mien_file = mixd::writeIntMixd(file, ".mien", full_cell_vertex_ids);
  std::filesystem::path mrng_file = mixd::writeIntMixd(file, ".mrng", full_cell_face_ids);
  mixd::writeMinf(file, mxyz_file, mien_file, mrng_file, full_cell_vertex_ids.size(), vertices.size());
}
} // namespace stmesh
#endif