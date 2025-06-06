#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "stmesh/mesh_project.hpp"
#include "stmesh/problem_types.hpp"
#include "stmesh/triangulation.hpp"
#include "stmesh/utility.hpp"

namespace nb = nanobind;
using namespace nb::literals;

class PyMeshProject {
  stmesh::TriangulationFromMixdWithData triangulation_;
  stmesh::MeshProjector projector_;

public:
  PyMeshProject(const std::string &minf_file, const std::string &data_file,
                const std::optional<stmesh::ProblemTypeEnum> &problem_type)
      : triangulation_(minf_file, data_file,
                       problem_type.has_value() ? std::optional{stmesh::kProblemTypeMap.at(*problem_type)}
                                                : std::nullopt),
        projector_(&triangulation_) {}

  PyMeshProject(const std::string &minf_file, const std::string &data_file,
                const std::vector<stmesh::DataEntry> &data_entries)
      : triangulation_(minf_file, data_file, stmesh::ProblemType{std::move(data_entries)}),
        projector_(&triangulation_) {}

  nb::dict project(const nb::DRef<Eigen::Matrix<stmesh::FLOAT_T, Eigen::Dynamic, 4>> &points) {
    std::vector<Eigen::MatrixX<stmesh::FLOAT_T>> data;
    for (const auto &data_entry : projector_.problemType().data_entries)
      data.emplace_back(points.rows(), data_entry.length);
    for (Eigen::Index i = 0; i < points.rows(); ++i) {
      projector_.problemType().forEach(projector_.project(points.row(i)),
                                       [&](size_t j, std::span<const stmesh::FLOAT_T> row_data) {
                                         std::ranges::copy(row_data, data[j].row(i).begin());
                                       });
    }
    nb::dict result;
    auto it = data.begin();
    for (const auto &data_entry : projector_.problemType().data_entries)
      result[data_entry.name.c_str()] = std::move(*it++);
    return result;
  }
};

NB_MODULE(py4dproject, m) {
  nb::enum_<stmesh::ProblemTypeEnum>(m, "ProblemType")
      .value("Solid", stmesh::ProblemTypeEnum::kSolidProblem)
      .value("Viscoelastic", stmesh::ProblemTypeEnum::kViscoelasticProblem)
      .value("AdvectionDiffusion", stmesh::ProblemTypeEnum::kAdvectionDiffusionProblem)
      .value("INS", stmesh::ProblemTypeEnum::kINSProblem)
      .value("CNS", stmesh::ProblemTypeEnum::kCNSProblem)
      .value("RVTCNS", stmesh::ProblemTypeEnum::kRVTCNSProblem)
      .value("EMUM", stmesh::ProblemTypeEnum::kEMUMProblem)
      .value("SolidUV", stmesh::ProblemTypeEnum::kSolidUVProblem);

  nb::class_<stmesh::DataEntry>(m, "DataEntry")
      .def(nb::init<const std::string &, size_t>())
      .def_rw("name", &stmesh::DataEntry::name)
      .def_rw("length", &stmesh::DataEntry::length);

  nb::class_<PyMeshProject>(m, "MeshProjector")
      .def(nb::init<const std::string &, const std::string &, const std::optional<stmesh::ProblemTypeEnum> &>(),
           "minf_file"_a, "data_file"_a, "problem_type"_a = nb::none())
      .def(nb::init<const std::string &, const std::string &, const std::vector<stmesh::DataEntry> &>(), "minf_file"_a,
           "data_file"_a, "data_entries"_a)
      .def("project", &PyMeshProject::project, "points"_a);
}
