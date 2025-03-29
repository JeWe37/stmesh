#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <vector>

#include "stmesh/geometric_simplex.hpp"
#include "stmesh/mesh_project.hpp"
#include "stmesh/mixd.hpp"
#include "stmesh/problem_types.hpp"
#include "stmesh/triangulation.hpp"
#include "stmesh/utility.hpp"

TEST_CASE("Test mesh projector", "[mesh_project]") {
  stmesh::TriangulationFromMixdWithData triangulation("data/sample_mixd/values.minf", stmesh::kINSProblem,
                                                      "data/sample_mixd/values.mxyz");
  stmesh::MeshProjector projector(&triangulation);

  std::array simplices{stmesh::GeometricSimplex<4>{
                           stmesh::Vector4F{0.0, 0.0, 0.0, 0.0},
                           stmesh::Vector4F{1.0, 0.0, 0.0, 0.0},
                           stmesh::Vector4F{0.0, 1.0, 0.0, 0.0},
                           stmesh::Vector4F{0.0, 0.0, 1.0, 0.0},
                           stmesh::Vector4F{0.3, 0.3, 0.3, 1.0},
                       },
                       stmesh::GeometricSimplex<4>{
                           stmesh::Vector4F{0.0, 0.0, 0.0, 0.0},
                           stmesh::Vector4F{1.0, 0.0, 0.0, 0.0},
                           stmesh::Vector4F{0.0, 1.0, 0.0, 0.0},
                           stmesh::Vector4F{0.0, 0.0, 1.0, 0.0},
                           stmesh::Vector4F{0.3, 0.3, 0.3, -1.0},
                       }};

  // As mxyz is used as data, result from projector should be exactly the point coordinates unless the point is outside
  constexpr static stmesh::FLOAT_T kEps = stmesh::sqrt(std::numeric_limits<stmesh::FLOAT_T>::epsilon());
  for (size_t i = 0; i < 100; ++i) {
    const stmesh::Vector4F point = stmesh::Vector4F::Random();
    const stmesh::Vector4F result = projector.project(point);
    const stmesh::Vector4F expected = [&] {
      for (const auto &simplex : simplices) {
        if (simplex.boundingBox().exteriorDistance(point) < kEps)
          return point;
      }
      return stmesh::Vector4F::Constant(-std::numeric_limits<stmesh::FLOAT_T>::quiet_NaN()).eval();
    }();
    INFO("Expected point: " << expected.transpose());
    INFO("Got point: " << result.transpose());
    if (expected.array().isNaN().all())
      REQUIRE(result.array().isNaN().all());
    else
      REQUIRE(result.isApprox(expected));
  }
}

TEST_CASE("Test mesh projector against old version", "[mesh_project]") {
  stmesh::TriangulationFromMixdWithData triangulation("data/projection_sample/artery.minf", stmesh::kINSProblem,
                                                      "data/projection_sample/INS_SST.out");
  stmesh::MeshProjector projector(&triangulation);

  std::vector<stmesh::Vector4F> test_positions = stmesh::mixd::readMxyz("data/projection_sample/stacked.dat");
  std::vector<std::vector<stmesh::FLOAT_T>> test_values =
      stmesh::mixd::readData("data/projection_sample/result.dat", 4);

  size_t failures = 0;

  for (size_t i = 0; i < test_positions.size(); ++i) {
    const stmesh::Vector4F expected = stmesh::Vector4F::Map(test_values[i].data());
    const stmesh::Vector4F result = projector.project(test_positions[i]);
    if (((result - expected).array() / result.norm()).abs().maxCoeff() >= 1e-5)
      failures++;
  }

  if (static_cast<stmesh::FLOAT_T>(failures) / static_cast<stmesh::FLOAT_T>(test_positions.size()) >= 0.1)
    FAIL("Too many failures: " << failures << " out of " << test_positions.size());
  WARN("Number of failures: " << failures << " out of " << test_positions.size());
}