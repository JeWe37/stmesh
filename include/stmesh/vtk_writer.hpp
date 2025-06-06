#ifndef STMESH_VTK_WRITER_HPP
#define STMESH_VTK_WRITER_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#include <CGAL/Polyhedron_3.h>
#pragma GCC diagnostic pop
#include <Eigen/Geometry>
#include <spdlog/spdlog.h>

#include "boundary_region_manager.hpp" // IWYU pragma: keep
#include "geometric_simplex.hpp"
#include "mesh_project.hpp"
#include "surface_adapters.hpp" // IWYU pragma: keep
#include "triangulation.hpp"
#include "utility.hpp"
#include "writable_triangulation.hpp"

namespace stmesh {
namespace detail {
using vtkIdType = long long;

class PolyhedraStorage;

class PointStorage {
  struct Impl;
  std::unique_ptr<Impl> pimpl_;

  friend void writeVTPFile(const std::filesystem::path &directory, const std::string_view &name_format,
                           const std::vector<std::vector<PolyhedraStorage>> &polygons,
                           const std::vector<PointStorage> &points, size_t block_pos, size_t n_positions,
                           bool write_boundary_ids);
  friend void writeVTUFile(const std::filesystem::path &directory, const std::string_view &name_format, FLOAT_T dt,
                           const std::vector<std::vector<PolyhedraStorage>> &polyhedra,
                           const std::vector<PointStorage> &points,
                           const Eigen::Transform<FLOAT_T, 4, Eigen::AffineCompact> &transformation,
                           const std::string_view &out_coord_format, stmesh::FLOAT_T scale, stmesh::FLOAT_T min_time,
                           size_t block_pos, size_t n_positions, FLOAT_T start_time, MeshProjector *projector);

public:
  PointStorage();
  ~PointStorage();
  PointStorage(const PointStorage &) = delete;
  PointStorage &operator=(const PointStorage &) = delete;
  PointStorage(PointStorage &&) noexcept = default;
  PointStorage &operator=(PointStorage &&) = default;

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

std::pair<std::vector<std::vector<detail::PolyhedraStorage>>, std::vector<detail::PointStorage>>
computePolyhedra(const std::vector<Eigen::Hyperplane<FLOAT_T, 4>> &planes,
                 const WritableTriangulation auto &triangulation,
                 const Eigen::Transform<FLOAT_T, 4, Eigen::AffineCompact> &transformation, FLOAT_T start_time,
                 FLOAT_T dt, size_t n_positions) {
  std::vector<std::vector<detail::PolyhedraStorage>> polyhedra(n_positions);
  std::vector<detail::PointStorage> points(n_positions);
  size_t full_cell_id = 0;
  for (const auto &cell : triangulation) {
    GeometricSimplex<4> simplex = cell.geometricSimplex();
    simplex.transform(transformation);
    const auto colwise = simplex.vertices().colwise(); // See https://gitlab.com/libeigen/eigen/-/issues/2882
    const auto [min_time, max_time] =
        std::ranges::minmax_element(colwise, {}, [](const Vector4F &vert) -> FLOAT_T { return vert[3]; });
    auto min_idx = static_cast<size_t>(std::max(FLOAT_T(), std::ceil(((*min_time)[3] - start_time) / dt)));
    auto max_idx =
        std::min(n_positions, static_cast<size_t>(std::max(FLOAT_T(), std::ceil(((*max_time)[3] - start_time) / dt))));
    for (size_t i = min_idx; i < max_idx; ++i) {
      auto polyhedron = simplex.planeCut(planes[i]);
      if (polyhedron.is_empty()) {
        spdlog::warn("Empty polyhedron ignored at time {}", start_time + static_cast<FLOAT_T>(i) * dt);
        continue;
      }
      polyhedra[i].emplace_back(polyhedron, points[i], full_cell_id);
    }
    full_cell_id++;
  }
  return {std::move(polyhedra), std::move(points)};
}

std::pair<std::vector<std::vector<detail::PolyhedraStorage>>, std::vector<detail::PointStorage>> computePolygons(
    const std::vector<Eigen::Hyperplane<FLOAT_T, 4>> &planes, const WritableTriangulation auto &triangulation,
    const Eigen::Transform<FLOAT_T, 4, Eigen::AffineCompact> &transformation,
    const BoundaryRegionManager auto &boundary_region_manager, FLOAT_T start_time, FLOAT_T dt, size_t n_positions) {
  std::vector<std::vector<detail::PolyhedraStorage>> polygons(n_positions);
  std::vector<detail::PointStorage> points(n_positions);
  for (const auto &cell : triangulation) {
    GeometricSimplex<4> simplex = cell.geometricSimplex();
    simplex.transform(transformation);
    for (size_t i = 0; i < 5; ++i) {
      if (cell.isSurfaceSide(i)) {
        Eigen::Matrix4<FLOAT_T> face_vertices;
        face_vertices << simplex.vertices().leftCols(i), simplex.vertices().rightCols(4 - i);
        const auto colwise = face_vertices.colwise(); // See https://gitlab.com/libeigen/eigen/-/issues/2882
        const auto [min_time, max_time] =
            std::ranges::minmax_element(colwise, {}, [](const Vector4F &vert) -> FLOAT_T { return vert[3]; });
        auto min_idx = static_cast<size_t>(std::max(FLOAT_T(), std::ceil(((*min_time)[3] - start_time) / dt)));
        auto max_idx = std::min(
            n_positions, static_cast<size_t>(std::max(FLOAT_T(), std::ceil(((*max_time)[3] - start_time) / dt))));
        size_t boundary_id = boundary_region_manager.findBoundaryRegion(cell, i);
        GeometricSimplex<4, 4> face_simplex(face_vertices);
        for (size_t j = min_idx; j < max_idx; ++j) {
          if (const auto polygon = face_simplex.planeCut(planes[j]); !polygon.empty())
            polygons[j].emplace_back(polygon, points[j], boundary_id);
          else {
            spdlog::warn("Empty polygon ignored at time {}", start_time + static_cast<FLOAT_T>(j) * dt);
            spdlog::debug("Face simplex vertices:\n{}", face_simplex.vertices());
          }
        }
      }
    }
  }
  return {std::move(polygons), std::move(points)};
}

void writeVTUFile(const std::filesystem::path &directory, const std::string_view &name_format, FLOAT_T dt,
                  const std::vector<std::vector<detail::PolyhedraStorage>> &polyhedra,
                  const std::vector<detail::PointStorage> &points,
                  const Eigen::Transform<FLOAT_T, 4, Eigen::AffineCompact> &transformation,
                  const std::string_view &out_coord_format, stmesh::FLOAT_T scale, stmesh::FLOAT_T min_time,
                  size_t block_pos, size_t n_positions, FLOAT_T start_time, MeshProjector *projector);

void writeVTPFile(const std::filesystem::path &directory, const std::string_view &name_format,
                  const std::vector<std::vector<detail::PolyhedraStorage>> &polygons,
                  const std::vector<detail::PointStorage> &points, size_t block_pos, size_t n_positions,
                  bool write_boundary_ids);
} // namespace detail

/// Add data to existing VTU files with precomputed coordinates
/**
 * Adds data to existing VTU files. The data is obtained by projecting the mesh data onto the points of the VTU files.
 * For determining the real positions of the points, the coordinates file is used.
 *
 * @param directory The directory containing the existing VTU files.
 * @param name_format The format string for the existing VTU file names. This format string should contain a single {}
 * which will be replaced with the time step number.
 * @param out_name_format The format string for the new VTU file names. This format string should contain a single {}
 * which will be replaced with the time step number.
 * @param steps The number of steps of the mesh.
 * @param triangulation The triangulation with data to use for the projector.
 * @param out_coord_format The format string for the coordinates file name. This format string should contain a single
 * {} which will be replaced with the time step number.
 */
void addVTUData(const std::filesystem::path &directory, std::string_view name_format, std::string_view out_name_format,
                size_t steps, const TriangulationFromMixdWithData &triangulation,
                const std::string_view &out_coord_format);

/// Add data to existing VTU files without precomputed coordinates
/**
 * Adds data to existing VTU files. The data is obtained by projecting the mesh data onto the points of the VTU files.
 * For determining the real positions of the points, they are computed based on the time step, transformation, scale,
 * and minimum time.
 *
 * @param directory The directory containing the existing VTU files.
 * @param name_format The format string for the existing VTU file names. This format string should contain a single {}
 * which will be replaced with the time step number.
 * @param out_name_format The format string for the new VTU file names. This format string should contain a single {}
 * which will be replaced with the time step number.
 * @param dt The time step for the VTU output.
 * @param triangulation The triangulation with data to use for the projector.
 * @param transformation The transformation to apply to the vertices of the mesh. Defaults to the identity matrix.
 * @param scale The scale for the output coord files. Defaults to 1.
 * @param min_time The minimum time of the triangulation for offsetting the output coord files. Defaults to 0.
 */
void addVTUData(const std::filesystem::path &directory, std::string_view name_format, std::string_view out_name_format,
                FLOAT_T dt, const TriangulationFromMixdWithData &triangulation,
                const Eigen::Transform<FLOAT_T, 4, Eigen::AffineCompact> &transformation =
                    Eigen::Transform<FLOAT_T, 4, Eigen::AffineCompact>::Identity(),
                stmesh::FLOAT_T scale = 1, stmesh::FLOAT_T min_time = 0);

/// Write a VTK file from a triangulation
/**
 * Writes a VTK file from a triangulation. The VTK file format is a file format used by ParaView and other visualization
 * software. In this case, the VTK file format is used to represent a mesh, or unstructured grid in VTK terms, thus with
 * a .vtu file extension.
 * As we are dealing with a time-dependent mesh, we write a series of VTK files, one for each time step. The time step
 * is determined by the dt parameter.
 * Which full cells are written to the VTK files is determined by the writable triangulation.
 * In addition to the VTK file, the coordinates of the vertices of the mesh can be written to a separate file. This is
 * useful for using them as the input to the mesh_projector. Its output can then be added to the VTK files.
 * The number of blocks dictates how many slices are computed per thread at once. This is useful for large meshes where
 * the VTK files are too large to fit in memory. The more blocks, the less memory is used, but the longer the
 * computation time.
 * Additionally, a VTP file can be written. This file contains the boundary of the mesh. This is useful for visualizing
 * the boundary regions in ParaView. Boundary regions are determined by the boundary_region_manager.
 * If the triangulation is a TriangulationFromMixdWithData, the mesh_projector is used to project the data onto the
 * points of the VTK file, which is then written to the VTK file.
 *
 * @param directory The directory to write the VTK files to.
 * @param name_format The format string for the VTK file names. This format string should contain a single {} which will
 * be replaced with the time step number.
 * @param dt The time step for the VTK output.
 * @param triangulation The writable triangulation to write to the VTK file.
 * @param transformation The transformation to apply to the vertices of the mesh. Defaults to the identity matrix.
 * @param scale The scale for the output coord files. Defaults to 1.
 * @param min_time The minimum time of the triangulation for offsetting the output coord files. Defaults to 0.
 * @param vtp_name_format The format string for the VTP file names. This format string should contain a single {} which
 * will be replaced with the time step number. If empty, no VTP file will be written.
 * @param boundary_region_manager The boundary region manager to use to determine the boundary regions of the faces of
 * the mesh. Boundary will not be specified in VTP file.
 * @param out_coord_format The format string for the coordinates file name. This format string should contain a single
 * {} which will be replaced with the time step number. If empty, no coordinates file will be written.
 * @param blocks The number of blocks to write the VTK files in. Defaults to 1.
 */
void writeVTU(const std::filesystem::path &directory, std::string_view name_format, FLOAT_T dt,
              const WritableTriangulation auto &triangulation,
              const Eigen::Transform<FLOAT_T, 4, Eigen::AffineCompact> &transformation =
                  Eigen::Transform<FLOAT_T, 4, Eigen::AffineCompact>::Identity(),
              stmesh::FLOAT_T scale = 1, stmesh::FLOAT_T min_time = 0, const std::string_view &vtp_name_format = "",
              const BoundaryRegionManager auto &boundary_region_manager = {},
              const std::string_view &out_coord_format = "", size_t blocks = 1) {
  spdlog::debug("Transforming with matrix:\n{}", transformation.matrix());
  Eigen::AlignedBox<FLOAT_T, 4> transformed_box(transformation * triangulation.boundingBox().min(),
                                                transformation * triangulation.boundingBox().max());
  FLOAT_T time_step = transformed_box.sizes()[3] / static_cast<FLOAT_T>(blocks);
  FLOAT_T dt_per_step = time_step / dt;
  std::vector<std::tuple<size_t, size_t, FLOAT_T>> blocks_info;
  size_t curr_pos = 0;
  FLOAT_T dt_pos = 0;
  for (size_t block = 0; block < blocks; ++block) {
    FLOAT_T start_time = transformed_box.min()[3] + dt * static_cast<FLOAT_T>(curr_pos);
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
          detail::computePolyhedra(planes, triangulation, transformation, start_time, dt, n_positions);
      if constexpr (std::is_same_v<std::decay_t<decltype(triangulation)>, TriangulationFromMixdWithData>) {
        MeshProjector projector(&triangulation);
        detail::writeVTUFile(directory, name_format, dt, polyhedra, points, transformation, out_coord_format, scale,
                             min_time, block_pos, n_positions, start_time, &projector);
      } else
        detail::writeVTUFile(directory, name_format, dt, polyhedra, points, transformation, out_coord_format, scale,
                             min_time, block_pos, n_positions, start_time, nullptr);
    }
    if (!vtp_name_format.empty()) {
      const auto [polygons, points] = detail::computePolygons(planes, triangulation, transformation,
                                                              boundary_region_manager, start_time, dt, n_positions);
      detail::writeVTPFile(directory, vtp_name_format, polygons, points, block_pos, n_positions,
                           !std::is_same_v<decltype(boundary_region_manager), NoopBoundaryManager>);
    }
    spdlog::debug("Block {}-{} done", start_time, start_time + static_cast<FLOAT_T>(n_positions) * dt);
  }
}
} // namespace stmesh
#endif