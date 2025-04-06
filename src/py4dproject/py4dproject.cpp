#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <algorithm>
#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include "stmesh/mesh_project.hpp"
#include "stmesh/problem_types.hpp"
#include "stmesh/triangulation.hpp"
#include "stmesh/utility.hpp"

namespace nb = nanobind;

class PyMeshProject {
  stmesh::TriangulationFromMixdWithData triangulation_;
  stmesh::MeshProjector projector_;

public:
  PyMeshProject(const std::string &minf_file, const stmesh::ProblemTypeEnum &problem_type, const std::string &data_file)
      : triangulation_(minf_file, stmesh::kProblemTypeMap.at(problem_type), data_file), projector_(&triangulation_) {}

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
      result[data_entry.name] = std::move(*it++);
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
      .value("RVTCNS", stmesh::ProblemTypeEnum::kRVTCNSProblem);

  nb::class_<PyMeshProject>(m, "MeshProjector")
      .def(nb::init<const std::string &, const stmesh::ProblemTypeEnum &, const std::string &>())
      .def("project", &PyMeshProject::project, nb::arg("points"));
}
