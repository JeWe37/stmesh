#ifndef STMESH_VTK_WRITER_HPP
#define STMESH_VTK_WRITER_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iterator>
#include <random>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <Eigen/Geometry>
#include <fmt/core.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

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
  std::vector<vtkIdType> pointIds_;
  size_t polyhedron_id_;

  PolyhedraStorage(const CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> &polyhedron,
                   PointStorage &points, size_t polyhedron_id);
};
} // namespace detail

/// Write a VTK file from a triangulation
/**
 * Writes a VTK file from a triangulation. The VTK file format is a file format used by ParaView and other visualization
 * software. In this case, the VTK file format is used to represent a mesh, or unstructured grid in VTK terms, thus with
 * a .vtu file extension.
 * As we are dealing with a time-dependent mesh, we write a series of VTK files, one for each time step. The time step
 * is determined by the dt parameter.
 * Only full cells whose circumsphere is inside the surface are written to the VTK file.
 * The number of blocks dictates how many slices are computed per thread at once. This is useful for large meshes where
 * the VTK files are too large to fit in memory. The more blocks, the less memory is used, but the longer the
 * computation time.
 *
 * @param directory The directory to write the VTK files to.
 * @param name_format The format string for the VTK file names. This format string should contain a single {} which will
 * be replaced with the time step number.
 * @param dt The time step for the VTK output.
 * @param surface The surface adapter to use to determine which full cells to write to the VTK file.
 * @param triangulation The triangulation to write to the VTK file.
 * @param blocks The number of blocks to write the VTK files in. Defaults to 1.
 * @tparam ExtraData The type of extra data stored in the triangulation.
 */
template <typename ExtraData>
void writeVTU([[maybe_unused]] const std::filesystem::path &directory, [[maybe_unused]] std::string_view name_format,
              FLOAT_T dt, [[maybe_unused]] const SurfaceAdapter4 auto &surface,
              const Triangulation<ExtraData> &triangulation, size_t blocks = 1) {
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
    std::vector<std::vector<detail::PolyhedraStorage>> polyhedra(n_positions);
    std::vector<detail::PointStorage> points(n_positions);
    size_t full_cell_id = 0;
    for (const auto &full_cell : triangulation) {
      std::vector<Vector4F> vertices;
      std::transform(full_cell.vertices_begin(), full_cell.vertices_end(), std::back_inserter(vertices),
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
      auto max_idx = std::min(n_positions,
                              static_cast<size_t>(std::max(FLOAT_T(), std::ceil(((*max_time)[3] - start_time) / dt))));
      for (size_t i = min_idx; i < max_idx; ++i)
        polyhedra[i].emplace_back(simplex.planeCut(planes[i]), points[i], full_cell_id);
      full_cell_id++;
    }

    for (size_t i = 0; i < n_positions; ++i) {
      vtkNew<vtkUnstructuredGrid> grid;
      grid->SetPoints(points[i].points_);
      vtkNew<vtkIntArray> id_array;
      id_array->SetName("Polyhedron ID");
      id_array->SetNumberOfTuples(static_cast<vtkIdType>(polyhedra[i].size()));
      vtkNew<vtkIntArray> random_array;
      random_array->SetName("Random");
      random_array->SetNumberOfTuples(static_cast<vtkIdType>(polyhedra[i].size()));
      for (const auto &polyhedron : polyhedra[i]) {
        int num_facets = 0;
        vtkNew<vtkIdList> faces;
        for (auto it = polyhedron.polyhedron_.facets_begin(); it != polyhedron.polyhedron_.facets_end(); ++it) {
          faces->InsertNextId(static_cast<vtkIdType>(it->facet_degree()));
          auto it_halfedge = it->facet_begin();
          for (size_t j = 0; j < it->facet_degree(); ++j) {
            auto pt = it_halfedge->vertex()->point();
            Vector3F point{static_cast<FLOAT_T>(pt.x()), static_cast<FLOAT_T>(pt.y()), static_cast<FLOAT_T>(pt.z())};
            faces->InsertNextId(points[i].positions_[point]);
            it_halfedge++;
          }
          num_facets++;
        }
        vtkIdType id = grid->InsertNextCell(VTK_POLYHEDRON, static_cast<vtkIdType>(polyhedron.pointIds_.size()),
                                            polyhedron.pointIds_.data(), num_facets, faces->GetPointer(0));
        id_array->InsertTuple1(id, static_cast<double>(polyhedron.polyhedron_id_));
        random_array->InsertTuple1(id, static_cast<double>(std::mt19937(polyhedron.polyhedron_id_)()));
      }
      grid->GetCellData()->AddArray(id_array);
      grid->GetCellData()->AddArray(random_array);
      vtkNew<vtkXMLUnstructuredGridWriter> writer;
      writer->SetFileName((directory / fmt::vformat(name_format, fmt::make_format_args(block_pos + i))).c_str());
      writer->SetInputData(grid);
      writer->SetDataModeToBinary();
      writer->Write();
    }
  }
}
} // namespace stmesh
#endif