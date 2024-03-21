#ifndef STMESH_VTK_WRITER_HPP
#define STMESH_VTK_WRITER_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <Eigen/Geometry>
#include <spdlog/spdlog.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkType.h>

#include "boundary_region_manager.hpp" // IWYU pragma: keep
#include "geometric_simplex.hpp"
#include "surface_adapters.hpp" // IWYU pragma: keep
#include "triangulation.hpp"
#include "utility.hpp"

namespace stmesh {
namespace detail {
class PointStorage {
public:
  std::unordered_map<Vector3F, vtkIdType, Vector3FHash> positions_;
  vtkNew<vtkPoints> points_;

  [[nodiscard]] vtkIdType insert(const Vector3F &point) noexcept;
};

class PolyhedraStorage {
public:
  CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> polyhedron_;
  std::vector<vtkIdType> point_ids_;
  size_t polyhedron_id_;

  PolyhedraStorage(const CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> &polyhedron,
                   PointStorage &points, size_t polyhedron_id);
};

void writeRaw(const std::filesystem::path &file, const std::vector<double> &values);

template <typename ExtraData>
std::pair<std::vector<std::vector<detail::PolyhedraStorage>>, std::vector<detail::PointStorage>>
computePolyhedra(const std::vector<Eigen::Hyperplane<FLOAT_T, 4>> &planes,
                 const stmesh::Triangulation<ExtraData> &triangulation, const SurfaceAdapter4 auto &surface,
                 FLOAT_T start_time, FLOAT_T dt, size_t n_positions) {
  std::vector<std::vector<detail::PolyhedraStorage>> polyhedra(n_positions);
  std::vector<detail::PointStorage> points(n_positions);
  size_t full_cell_id = 0;
  for (const auto &full_cell : triangulation) {
    std::array<Vector4F, 5> vertices;
    std::transform(full_cell.vertices_begin(), full_cell.vertices_end(), vertices.begin(),
                   [](const auto &vertex) -> Vector4F {
                     const auto &point = vertex->point();
                     return Vector4F{static_cast<FLOAT_T>(point[0]), static_cast<FLOAT_T>(point[1]),
                                     static_cast<FLOAT_T>(point[2]), static_cast<FLOAT_T>(point[3])};
                   });
    GeometricSimplex<4> simplex(vertices);
    if (!surface.inside(simplex.circumsphere().center()))
      continue;
    const auto [min_time, max_time] =
        std::ranges::minmax_element(vertices, {}, [](const Vector4F &vert) -> FLOAT_T { return vert[3]; });
    auto min_idx = static_cast<size_t>(std::max(FLOAT_T(), std::ceil(((*min_time)[3] - start_time) / dt)));
    auto max_idx =
        std::min(n_positions, static_cast<size_t>(std::max(FLOAT_T(), std::ceil(((*max_time)[3] - start_time) / dt))));
    for (size_t i = min_idx; i < max_idx; ++i)
      polyhedra[i].emplace_back(simplex.planeCut(planes[i]), points[i], full_cell_id);
    full_cell_id++;
  }
  return {std::move(polyhedra), std::move(points)};
}

template <typename ExtraData, BoundaryRegionManager BoundaryManager>
std::pair<std::vector<std::vector<detail::PolyhedraStorage>>, std::vector<detail::PointStorage>>
computePolygons(const std::vector<Eigen::Hyperplane<FLOAT_T, 4>> &planes,
                const stmesh::Triangulation<ExtraData> &triangulation, const BoundaryManager &boundary_region_manager,
                const SurfaceAdapter4 auto &surface, FLOAT_T start_time, FLOAT_T dt, size_t n_positions) {
  std::vector<std::vector<detail::PolyhedraStorage>> polygons(n_positions);
  std::vector<detail::PointStorage> points(n_positions);
  for (const auto &full_cell : triangulation) {
    typename stmesh::Triangulation<ExtraData>::FullCellConstHandle full_cell_handle{&full_cell};
    if (surface.inside(triangulation.fullCellSimplex(full_cell_handle).circumsphere().center())) {
      std::array<Vector4F, 5> vertices;
      std::transform(full_cell.vertices_begin(), full_cell.vertices_end(), vertices.begin(),
                     [](const auto &vertex) -> Vector4F {
                       const auto &point = vertex->point();
                       return Vector4F{static_cast<FLOAT_T>(point[0]), static_cast<FLOAT_T>(point[1]),
                                       static_cast<FLOAT_T>(point[2]), static_cast<FLOAT_T>(point[3])};
                     });
      for (size_t i = 0; i < 5; ++i) {
        auto neighbor = full_cell.neighbor(i);
        if (triangulation.isInfinite(neighbor) ||
            !surface.inside(triangulation.fullCellSimplex(neighbor).circumsphere().center())) {
          std::array<Vector4F, 4> face_vertices;
          std::copy_n(vertices.begin(), i, face_vertices.begin());
          std::copy_n(vertices.begin() + i + 1, 4 - i, face_vertices.begin() + i);
          const auto [min_time, max_time] =
              std::ranges::minmax_element(face_vertices, {}, [](const Vector4F &vert) -> FLOAT_T { return vert[3]; });
          auto min_idx = static_cast<size_t>(std::max(FLOAT_T(), std::ceil(((*min_time)[3] - start_time) / dt)));
          auto max_idx = std::min(
              n_positions, static_cast<size_t>(std::max(FLOAT_T(), std::ceil(((*max_time)[3] - start_time) / dt))));
          GeometricSimplex<4, 4> simplex(face_vertices);
          size_t boundary_id{};
          if constexpr (!std::is_same_v<BoundaryManager, NoopBoundaryManager>) {
            std::array<size_t, 4> vertex_boundary_ids{};
            for (size_t j = 0; j < 4; ++j)
              vertex_boundary_ids.at(j) = boundary_region_manager.findBoundaryRegion(face_vertices.at(j));
            boundary_id = *std::ranges::max_element(vertex_boundary_ids);
          }
          for (size_t j = min_idx; j < max_idx; ++j) {
            if (const auto polygon = simplex.planeCut(planes[j]); !polygon.empty())
              polygons[j].emplace_back(polygon, points[j], boundary_id);
            else {
              spdlog::warn("Empty polygon ignored at time {}", start_time + static_cast<FLOAT_T>(j) * dt);
              spdlog::debug("Face simplex vertices:\n{}", simplex.vertices());
            }
          }
        }
      }
    }
  }
  return {std::move(polygons), std::move(points)};
}

void writeVTUFile(const std::filesystem::path &directory, const std::string_view &name_format, FLOAT_T dt,
                  const std::vector<std::vector<detail::PolyhedraStorage>> &polyhedra,
                  const std::vector<detail::PointStorage> &points, const std::string_view &out_coord_format,
                  stmesh::FLOAT_T scale, stmesh::FLOAT_T min_time, size_t block_pos, size_t n_positions,
                  FLOAT_T start_time);

void writeVTPFile(const std::filesystem::path &directory, const std::string_view &name_format,
                  const std::vector<std::vector<detail::PolyhedraStorage>> &polygons,
                  const std::vector<detail::PointStorage> &points, size_t block_pos, size_t n_positions,
                  bool write_boundary_ids);
} // namespace detail

/// Write a VTK file from a triangulation
/**
 * Writes a VTK file from a triangulation. The VTK file format is a file format used by ParaView and other visualization
 * software. In this case, the VTK file format is used to represent a mesh, or unstructured grid in VTK terms, thus with
 * a .vtu file extension.
 * As we are dealing with a time-dependent mesh, we write a series of VTK files, one for each time step. The time step
 * is determined by the dt parameter.
 * Only full cells whose circumsphere is inside the surface are written to the VTK file.
 * In addition to the VTK file, the coordinates of the vertices of the mesh can be written to a separate file. This is
 * useful for using them as the input to the mesh_projector. Its output can then be added to the VTK files.
 * The number of blocks dictates how many slices are computed per thread at once. This is useful for large meshes where
 * the VTK files are too large to fit in memory. The more blocks, the less memory is used, but the longer the
 * computation time.
 * Additionally, a VTP file can be written. This file contains the boundary of the mesh. This is useful for visualizing
 * the boundary regions in ParaView. Boundary regions are determined by the boundary_region_manager.
 *
 * @param directory The directory to write the VTK files to.
 * @param name_format The format string for the VTK file names. This format string should contain a single {} which will
 * be replaced with the time step number.
 * @param dt The time step for the VTK output.
 * @param surface The surface adapter to use to determine which full cells to write to the VTK file.
 * @param triangulation The triangulation to write to the VTK file.
 * @param scale The scale for the output coord files. Defaults to 1.
 * @param min_time The minimum time of the triangulation for offsetting the output coord files. Defaults to 0.
 * @param vtp_name_format The format string for the VTP file names. This format string should contain a single {} which
 * will be replaced with the time step number. If empty, no VTP file will be written.
 * @param boundary_region_manager The boundary region manager to use to determine the boundary regions of the faces of
 * the mesh. Boundary will not be specified in VTP file.
 * @param out_coord_format The format string for the coordinates file name. This format string should contain a single
 * {} which will be replaced with the time step number. If empty, no coordinates file will be written.
 * @param blocks The number of blocks to write the VTK files in. Defaults to 1.
 * @tparam ExtraData The type of extra data stored in the triangulation.
 */
template <typename ExtraData, BoundaryRegionManager BoundaryManager = NoopBoundaryManager>
void writeVTU(const std::filesystem::path &directory, std::string_view name_format, FLOAT_T dt,
              [[maybe_unused]] const SurfaceAdapter4 auto &surface, const Triangulation<ExtraData> &triangulation,
              stmesh::FLOAT_T scale = 1, stmesh::FLOAT_T min_time = 0, const std::string_view &vtp_name_format = "",
              const BoundaryManager &boundary_region_manager = {}, const std::string_view &out_coord_format = "",
              size_t blocks = 1) {
  FLOAT_T time_step = triangulation.boundingBox().sizes()[3] / blocks;
  FLOAT_T dt_per_step = time_step / dt;
  std::vector<std::tuple<size_t, size_t, FLOAT_T>> blocks_info;
  size_t curr_pos = 0;
  FLOAT_T dt_pos = 0;
  for (size_t block = 0; block < blocks; ++block) {
    FLOAT_T start_time = triangulation.boundingBox().min()[3] + dt * static_cast<FLOAT_T>(curr_pos);
    dt_pos += dt_per_step;
    auto n_positions = static_cast<size_t>(dt_pos) + 1;
    dt_pos -= static_cast<FLOAT_T>(n_positions);
    blocks_info.emplace_back(curr_pos, n_positions, start_time);
    curr_pos += n_positions;
  }
// NOLINTNEXTLINE(openmp-use-default-none)
#pragma omp parallel for
  for (const auto &[block_pos, n_positions, start_time] : blocks_info) {
    std::vector<Eigen::Hyperplane<FLOAT_T, 4>> planes;
    for (size_t i = 0; i < n_positions; ++i) {
      FLOAT_T time = static_cast<FLOAT_T>(i) * dt + start_time;
      planes.emplace_back(Vector4F{0, 0, 0, 1}, -time);
    }
    {
      const auto [polyhedra, points] =
          detail::computePolyhedra(planes, triangulation, surface, start_time, dt, n_positions);
      detail::writeVTUFile(directory, name_format, dt, polyhedra, points, out_coord_format, scale, min_time, block_pos,
                           n_positions, start_time);
    }
    if (!vtp_name_format.empty()) {
      const auto [polygons, points] =
          detail::computePolygons(planes, triangulation, boundary_region_manager, surface, start_time, dt, n_positions);
      detail::writeVTPFile(directory, vtp_name_format, polygons, points, block_pos, n_positions,
                           !std::is_same_v<BoundaryManager, NoopBoundaryManager>);
    }
    spdlog::debug("Block {}-{} done", start_time, start_time + static_cast<FLOAT_T>(n_positions) * dt);
  }
}
} // namespace stmesh
#endif