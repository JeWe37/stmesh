#include "stmesh/vtk_writer.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <random>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <fmt/core.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "stmesh/mesh_project.hpp"
#include "stmesh/mixd.hpp"
#include "stmesh/problem_types.hpp"
#include "stmesh/triangulation.hpp"
#include "stmesh/utility.hpp"

namespace stmesh {
namespace detail {
static_assert(std::is_same_v<vtkIdType, ::vtkIdType>, "vtkIdType forward declaration is incorrect.");

struct PointStorage::Impl {
  std::unordered_map<Vector3F, vtkIdType, Vector3FHash> positions_;
  vtkNew<vtkPoints> points_;

  [[nodiscard]] vtkIdType insert(const Vector3F &point) noexcept {
    if (positions_.contains(point))
      return positions_[point];
    return positions_[point] = points_->InsertNextPoint(static_cast<double>(point[0]), static_cast<double>(point[1]),
                                                        static_cast<double>(point[2]));
  }
};

PointStorage::PointStorage() : pimpl_(new Impl) {}

PointStorage::~PointStorage() = default;

[[nodiscard]] vtkIdType PointStorage::insert(const Vector3F &point) noexcept { return pimpl_->insert(point); }

PolyhedraStorage::PolyhedraStorage(
    const CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> &polyhedron, PointStorage &points,
    size_t polyhedron_id)
    : polyhedron_(polyhedron), polyhedron_id_(polyhedron_id) {
  for (auto it_vert = polyhedron.vertices_begin(); it_vert != polyhedron.vertices_end(); ++it_vert) {
    const auto &vertex = it_vert->point();
    const Vector3F point{vertex.x(), vertex.y(), vertex.z()};
    point_ids_.push_back(points.insert(point));
  }
}

void writeRaw(const std::filesystem::path &file, const std::vector<double> &values) {
  std::ofstream out(file);
  for (const double value : values) {
    auto value_bytes = std::bit_cast<std::array<char, sizeof(FLOAT_T)>>(value);
    if constexpr (std::endian::native != std::endian::big)
      std::reverse(value_bytes.begin(), value_bytes.end());
    out.write(value_bytes.data(), sizeof(double));
  }
}

void addData(const std::vector<Vector4F> &poses, vtkUnstructuredGrid *grid, const MeshProjector &projector) {
  const ProblemType &problem_type = projector.problemType();
  std::vector<vtkNew<vtkDoubleArray>> data(problem_type.data_entries.size());
  for (size_t j = 0; j < data.size(); ++j) {
    data[j]->SetName(problem_type.data_entries[j].name.c_str());
    data[j]->SetNumberOfComponents(static_cast<int>(problem_type.data_entries[j].length));
  }

  for (const Vector4F &pos : poses)
    problem_type.forEach(projector.project(pos),
                         [&](size_t j, std::span<const FLOAT_T> value) { data[j]->InsertNextTuple(value.data()); });

  for (const auto &arr : data)
    grid->GetPointData()->AddArray(arr);
}

template <typename F>
void editVTU(const std::filesystem::path &in_vtu, const std::filesystem::path &out_vtu, const MeshProjector &projector,
             const F &point_callback) {
  const vtkNew<vtkXMLUnstructuredGridReader> reader;
  reader->SetFileName(in_vtu.c_str());
  reader->Update();
  const vtkSmartPointer<vtkUnstructuredGrid> grid = reader->GetOutput();

  addData(point_callback(grid), grid, projector);

  const vtkNew<vtkXMLUnstructuredGridWriter> writer;
  writer->SetFileName(out_vtu.c_str());
  writer->SetInputData(grid);
  writer->SetDataModeToBinary();
  writer->Write();
}

void writeVTUFile(const std::filesystem::path &directory, const std::string_view &name_format, FLOAT_T dt,
                  const std::vector<std::vector<detail::PolyhedraStorage>> &polyhedra,
                  const std::vector<detail::PointStorage> &points,
                  const Eigen::Transform<FLOAT_T, 4, Eigen::AffineCompact> &transformation,
                  const std::string_view &out_coord_format, stmesh::FLOAT_T scale, stmesh::FLOAT_T min_time,
                  size_t block_pos, size_t n_positions, FLOAT_T start_time, MeshProjector *projector) {
  for (size_t i = 0; i < n_positions; ++i) {
    std::vector<Vector4F> poses;
    if (!out_coord_format.empty() || projector != nullptr) {
      if (projector != nullptr)
        poses.resize(points[i].pimpl_->positions_.size());
      const FLOAT_T time = static_cast<FLOAT_T>(i) * dt + start_time;
      std::vector<double> raw_point_data(points[i].pimpl_->positions_.size() * 4);
      for (const vtkIdType id : points[i].pimpl_->positions_ | std::views::values) {
        points[i].pimpl_->points_->GetPoint(id, &raw_point_data[static_cast<size_t>(id * 4)]);
        raw_point_data[static_cast<size_t>(id * 4 + 3)] = time - min_time;
        Eigen::Map<Vector4F> mapped(&raw_point_data[static_cast<size_t>(id) * 4]);
        mapped = transformation.inverse().scale(scale) * mapped;
        if (projector != nullptr)
          poses[static_cast<size_t>(id)] = mapped;
      }
      if (!out_coord_format.empty())
        detail::writeRaw(directory / fmt::vformat(out_coord_format, fmt::make_format_args(block_pos + i)),
                         raw_point_data);
    }
    const vtkNew<vtkUnstructuredGrid> grid;
    grid->SetPoints(points[i].pimpl_->points_);

    const vtkNew<vtkIntArray> id_array;
    id_array->SetName("Polyhedron ID");
    id_array->SetNumberOfTuples(static_cast<vtkIdType>(polyhedra[i].size()));
    const vtkNew<vtkIntArray> random_array;
    random_array->SetName("Random");
    random_array->SetNumberOfTuples(static_cast<vtkIdType>(polyhedra[i].size()));
    for (const auto &polyhedron : polyhedra[i]) {
      int num_facets = 0;
      const vtkNew<vtkIdList> faces;
      for (auto it = polyhedron.polyhedron_.facets_begin(); it != polyhedron.polyhedron_.facets_end(); ++it) {
        faces->InsertNextId(static_cast<vtkIdType>(it->facet_degree()));
        auto it_halfedge = it->facet_begin();
        for (size_t j = 0; j < it->facet_degree(); ++j) {
          auto pt = it_halfedge->vertex()->point();
          const Vector3F point{static_cast<FLOAT_T>(pt.x()), static_cast<FLOAT_T>(pt.y()),
                               static_cast<FLOAT_T>(pt.z())};
          faces->InsertNextId(points[i].pimpl_->positions_.at(point));
          it_halfedge++;
        }
        num_facets++;
      }
      const vtkIdType id = grid->InsertNextCell(VTK_POLYHEDRON, static_cast<vtkIdType>(polyhedron.point_ids_.size()),
                                                polyhedron.point_ids_.data(), num_facets, faces->GetPointer(0));
      id_array->InsertTuple1(id, static_cast<double>(polyhedron.polyhedron_id_));
      random_array->InsertTuple1(
          id, static_cast<double>(std::mt19937(static_cast<uint_fast32_t>(polyhedron.polyhedron_id_))()));
    }
    grid->GetCellData()->AddArray(id_array);
    grid->GetCellData()->AddArray(random_array);
    if (projector != nullptr)
      addData(poses, grid, *projector);

    const vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName((directory / fmt::vformat(name_format, fmt::make_format_args(block_pos + i))).c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToBinary();
    writer->Write();
  }
}

void writeVTPFile(const std::filesystem::path &directory, const std::string_view &name_format,
                  const std::vector<std::vector<detail::PolyhedraStorage>> &polygons,
                  const std::vector<detail::PointStorage> &points, size_t block_pos, size_t n_positions,
                  bool write_boundary_ids) {
  for (size_t i = 0; i < n_positions; ++i) {
    const vtkNew<vtkPolyData> poly_data;
    const vtkNew<vtkCellArray> polys;
    poly_data->SetPoints(points[i].pimpl_->points_);
    const vtkNew<vtkIntArray> id_array;
    if (write_boundary_ids) {
      id_array->SetName("Boundary ID");
      id_array->SetNumberOfTuples(static_cast<vtkIdType>(polygons[i].size()));
    }
    for (const auto &polygon : polygons[i]) {
      auto it = polygon.polyhedron_.facets_begin();
      auto it_halfedge = it->facet_begin();
      std::vector<vtkIdType> face;
      for (size_t j = 0; j < it->facet_degree(); ++j) {
        auto pt = it_halfedge->vertex()->point();
        const Vector3F point{static_cast<FLOAT_T>(pt.x()), static_cast<FLOAT_T>(pt.y()), static_cast<FLOAT_T>(pt.z())};
        face.push_back(points[i].pimpl_->positions_.at(point));
        it_halfedge++;
      }
      const vtkIdType id = polys->InsertNextCell(static_cast<vtkIdType>(face.size()), face.data());
      if (write_boundary_ids)
        id_array->InsertTuple1(id, static_cast<double>(polygon.polyhedron_id_));
    }
    poly_data->SetPolys(polys);
    if (write_boundary_ids)
      poly_data->GetCellData()->AddArray(id_array);
    const vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetFileName((directory / fmt::vformat(name_format, fmt::make_format_args(block_pos + i))).c_str());
    writer->SetInputData(poly_data);
    writer->SetDataModeToBinary();
    writer->Write();
  }
}
} // namespace detail

void addVTUData(const std::filesystem::path &directory, std::string_view name_format, std::string_view out_name_format,
                size_t steps, const TriangulationFromMixdWithData &triangulation,
                const std::string_view &out_coord_format) {
  const MeshProjector projector(&triangulation);
  for (size_t step = 0; step < steps; ++step) {
    if (const std::filesystem::path in_vtu = directory / fmt::vformat(name_format, fmt::make_format_args(step));
        std::filesystem::status(in_vtu).type() != std::filesystem::file_type::regular)
      throw std::runtime_error(fmt::format("{} is not a regular file", in_vtu.c_str()));
  }

#pragma omp parallel for default(none)                                                                                 \
    shared(directory, name_format, out_name_format, steps, projector, triangulation, out_coord_format)
  for (size_t i = 0; i < steps; ++i) {
    detail::editVTU(directory / fmt::vformat(name_format, fmt::make_format_args(i)),
                    directory / fmt::vformat(out_name_format, fmt::make_format_args(i)), projector,
                    [&](const vtkSmartPointer<vtkUnstructuredGrid> & /*grid*/) {
                      return mixd::readMxyz(directory / fmt::vformat(out_coord_format, fmt::make_format_args(i)));
                    });
  }
}

void addVTUData(const std::filesystem::path &directory, std::string_view name_format, std::string_view out_name_format,
                FLOAT_T dt, const TriangulationFromMixdWithData &triangulation,
                const Eigen::Transform<FLOAT_T, 4, Eigen::AffineCompact> &transformation, stmesh::FLOAT_T scale,
                stmesh::FLOAT_T min_time) {
  const MeshProjector projector(&triangulation);
  Eigen::AlignedBox<FLOAT_T, 4> transformed_box(transformation * triangulation.boundingBox().min(),
                                                transformation * triangulation.boundingBox().max());
  const size_t steps = static_cast<size_t>(std::ceil(transformed_box.sizes()[3] / dt));
  for (size_t step = 0; step < steps; ++step) {
    if (const std::filesystem::path in_vtu = directory / fmt::vformat(name_format, fmt::make_format_args(step));
        std::filesystem::status(in_vtu).type() != std::filesystem::file_type::regular)
      throw std::runtime_error(fmt::format("{} is not a regular file", in_vtu.c_str()));
  }

#pragma omp parallel for default(none)                                                                                 \
    shared(directory, name_format, out_name_format, projector, triangulation, min_time, steps, transformed_box, scale, \
               dt, transformation, Eigen::Dynamic)
  for (size_t step = 0; step < steps; ++step) {
    detail::editVTU(directory / fmt::vformat(name_format, fmt::make_format_args(step)),
                    directory / fmt::vformat(out_name_format, fmt::make_format_args(step)), projector,
                    [&](const vtkSmartPointer<vtkUnstructuredGrid> &grid) {
                      std::vector<Vector4F> poses;
                      for (vtkIdType i = 0; i < grid->GetPoints()->GetNumberOfPoints(); ++i) {
                        std::array<FLOAT_T, 4> point{};
                        grid->GetPoints()->GetPoint(i, point.data());
                        point[3] = transformed_box.min()[3] + static_cast<FLOAT_T>(step) * dt - min_time;
                        Eigen::Map<Vector4F> mapped(point.data());
                        mapped = transformation.inverse().scale(scale) * mapped;
                        poses.emplace_back(mapped);
                      }
                      return poses;
                    });
  }
}
} // namespace stmesh
