#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <limits>
#include <ranges>
#include <set>
#include <span>
#include <vector>

#include <Eigen/Core>

#include "stmesh/geometric_simplex.hpp"
#include "stmesh/mixd.hpp"
#include "stmesh/problem_types.hpp"
#include "stmesh/sdf.hpp"
#include "stmesh/sdf_mixins.hpp"
#include "stmesh/surface_adapters.hpp"
#include "stmesh/triangulation.hpp"
#include "stmesh/utility.hpp"

using namespace Catch;

struct DummySDF : public stmesh::DistanceMixin<DummySDF, 4>, public stmesh::CentralDifferenceNormalMixin<DummySDF, 4> {
  constexpr static double kInf = std::numeric_limits<stmesh::FLOAT_T>::infinity();
  static constexpr unsigned kDimension = 4;

  [[nodiscard]] static stmesh::FLOAT_T signedDistance(const stmesh::Vector4F &point) noexcept { return 0.5 - point[3]; }

  [[nodiscard]] static Eigen::AlignedBox<stmesh::FLOAT_T, 4> boundingBox() noexcept {
    return {stmesh::Vector4F{-kInf, -kInf, -kInf, 0.5}, stmesh::Vector4F::Constant(kInf)};
  }
};

template <> struct stmesh::ExactSDFTag<DummySDF> {
  static constexpr bool kValue = false;
};

static_assert(stmesh::ExactSDF<DummySDF>);

TEST_CASE("Test triangulation store", "[triangulation]") {
  const stmesh::Vector4F min{-stmesh::FLOAT_T(2.0), -stmesh::FLOAT_T(2.0), -stmesh::FLOAT_T(2.0),
                             -stmesh::FLOAT_T(2.0)};
  const stmesh::Vector4F max{stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0)};
  const Eigen::AlignedBox<stmesh::FLOAT_T, 4> bounding_box{min, max};
  stmesh::Triangulation triangulation(bounding_box);
  SECTION("Correctly constructed") {
    REQUIRE(triangulation.boundingBox().min() == bounding_box.min());
    REQUIRE(triangulation.boundingBox().max() == bounding_box.max());
    for (const auto &full_cell : triangulation) {
      REQUIRE(full_cell.data().committed);
      for (const auto &vertex : std::span(full_cell.vertices_begin(), full_cell.vertices_end())) {
        const stmesh::Triangulation<>::Point pt = vertex->point();
        const stmesh::Vector4F point{static_cast<stmesh::FLOAT_T>(pt[0]), static_cast<stmesh::FLOAT_T>(pt[1]),
                                     static_cast<stmesh::FLOAT_T>(pt[2]), static_cast<stmesh::FLOAT_T>(pt[3])};
        REQUIRE(point.cwiseAbs() == max);
        REQUIRE(vertex->data().nonfree_vertex);
      }
    }
  }
  const std::vector<stmesh::Vector4F> vertices{
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)},
      {stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0)}};
  SECTION("Addition and removal") {
    std::vector<stmesh::Triangulation<>::VertexHandle> vertex_handles;
    for (const auto &vertex : vertices) {
      const std::vector<stmesh::Triangulation<>::FullCellHandle> conflict_zone = triangulation.conflictZone(vertex);
      std::vector<stmesh::Triangulation<>::FullCellConstHandle> const_conflict_zone;
      std::transform(conflict_zone.begin(), conflict_zone.end(), std::back_inserter(const_conflict_zone),
                     [](const auto &full_cell) { return stmesh::Triangulation<>::FullCellConstHandle{&*full_cell}; });

      std::set<stmesh::Triangulation<>::FullCellConstHandle> cells_before;
      std::transform(triangulation.begin(), triangulation.end(), std::inserter(cells_before, cells_before.begin()),
                     [](const auto &full_cell) { return stmesh::Triangulation<>::FullCellConstHandle{&full_cell}; });
      const stmesh::Triangulation<>::VertexHandle vertex_handle = triangulation.insert(vertex);
      vertex_handles.push_back(vertex_handle);
      REQUIRE_FALSE(vertex_handle->data().nonfree_vertex);
      std::set<stmesh::Triangulation<>::FullCellConstHandle> cells_after;
      std::transform(triangulation.begin(), triangulation.end(), std::inserter(cells_after, cells_after.begin()),
                     [](const auto &full_cell) { return stmesh::Triangulation<>::FullCellConstHandle{&full_cell}; });

      std::vector<stmesh::Triangulation<>::FullCellConstHandle> cells_before_not_after;
      std::set_difference(cells_before.begin(), cells_before.end(), cells_after.begin(), cells_after.end(),
                          std::back_inserter(cells_before_not_after));

      REQUIRE_THAT(const_conflict_zone, Matchers::UnorderedEquals(cells_before_not_after));

      std::vector<stmesh::Triangulation<>::FullCellConstHandle> cells_after_not_before;
      std::set_difference(cells_after.begin(), cells_after.end(), cells_before.begin(), cells_before.end(),
                          std::back_inserter(cells_after_not_before));
      const std::vector<stmesh::Triangulation<>::FullCellHandle> surrounding_full_cells =
          triangulation.surroundingFullCells(vertex_handle);
      std::vector<stmesh::Triangulation<>::FullCellConstHandle> const_surrounding_full_cells;
      std::transform(surrounding_full_cells.begin(), surrounding_full_cells.end(),
                     std::back_inserter(const_surrounding_full_cells),
                     [](const auto &full_cell) { return stmesh::Triangulation<>::FullCellConstHandle{&*full_cell}; });
      REQUIRE_THAT(const_surrounding_full_cells, Matchers::UnorderedEquals(cells_after_not_before));

      for (const auto &full_cell : triangulation)
        REQUIRE(full_cell.data().committed);
    }

    SECTION("Writable iteration works") {
      stmesh::SDFSurfaceAdapter<DummySDF> surface_adapter;
      auto writable_triangulation = triangulation.writableTriangulation(&surface_adapter);

      REQUIRE(writable_triangulation.boundingBox().min() == triangulation.boundingBox().min());
      REQUIRE(writable_triangulation.boundingBox().max() == triangulation.boundingBox().max());

      auto expected = std::count_if(triangulation.begin(), triangulation.end(), [&](const auto &full_cell) {
        return surface_adapter.inside(
            triangulation.fullCellSimplex(stmesh::Triangulation<>::FullCellConstHandle(&full_cell))
                .circumsphere()
                .center());
      });
      REQUIRE(std::distance(writable_triangulation.begin(), writable_triangulation.end()) == expected);

      for (const auto &cell : writable_triangulation)
        REQUIRE(surface_adapter.inside(cell.geometricSimplex().circumsphere().center()));

      // TODO: surface side check, not currently feasible
    }

    SECTION("Finds correct facet and mirror vertices") {
      std::vector<stmesh::Vector4F> simplex3(vertices.begin(), vertices.end() - 2);
      Eigen::Matrix<stmesh::FLOAT_T, 4, 4> simplex3_matrix;
      for (Eigen::Index i = 0; i < 4; i++)
        simplex3_matrix.col(i) = simplex3[static_cast<size_t>(i)];
      const stmesh::Triangulation<>::Facet facet = triangulation.facetFromVertices(simplex3_matrix);
      std::vector<stmesh::Vector4F> facet_vertices;
      const stmesh::Triangulation<>::FullCellHandle facet_full_cell = facet.full_cell();
      for (int i = 0; i < 5; i++) {
        if (i != facet.index_of_covertex()) {
          const stmesh::Triangulation<>::Point pt = facet_full_cell->vertex(i)->point();
          facet_vertices.emplace_back(static_cast<stmesh::FLOAT_T>(pt[0]), static_cast<stmesh::FLOAT_T>(pt[1]),
                                      static_cast<stmesh::FLOAT_T>(pt[2]), static_cast<stmesh::FLOAT_T>(pt[3]));
        }
      }
      REQUIRE_THAT(facet_vertices, Matchers::UnorderedEquals(simplex3));
      auto [existing_side, optional_side] = triangulation.facetMirrorVertices(facet);
      REQUIRE(optional_side.has_value());
      std::array<stmesh::Triangulation<>::VertexHandle, 2> mirror_vertices = {std::get<0>(existing_side),
                                                                              std::get<0>(*optional_side)};
      std::array<int, 2> vertex_indices = {std::get<1>(existing_side), std::get<1>(*optional_side)};
      std::array<stmesh::Triangulation<>::FullCellHandle, 2> full_cells = {std::get<2>(existing_side),
                                                                           std::get<2>(*optional_side)};
      for (size_t i = 0; i < 2; ++i)
        REQUIRE(full_cells.at(i)->vertex(vertex_indices.at(i)) == mirror_vertices.at(i));
      for (int i = 0; i < 5; ++i) {
        if (i != facet.index_of_covertex()) {
          for (size_t j = 0; j < 2; ++j)
            REQUIRE(full_cells.at(j)->has_vertex(facet.full_cell()->vertex(i)));
        }
      }
      REQUIRE(full_cells.at(0)->has_neighbor(full_cells.at(1)));
    }

    SECTION("Correct kd-tree queries") {
      const std::vector<stmesh::Triangulation<>::VertexHandle> simplex4_handles(vertex_handles.begin(),
                                                                                vertex_handles.end() - 1);

      REQUIRE_THAT(triangulation.verticesInRadius(
                       {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
                       stmesh::FLOAT_T(1.5)),
                   Matchers::UnorderedEquals(simplex4_handles));
      std::vector<stmesh::Triangulation<>::VertexHandle> just_origin = triangulation.verticesInRadius(
          {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
          stmesh::FLOAT_T(0.3));
      REQUIRE(just_origin.size() == 1);
      REQUIRE(just_origin[0] == vertex_handles[0]);
    }
    SECTION("Test good point check") {
      REQUIRE_FALSE(triangulation.isGoodPoint(vertex_handles[0], stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.2),
                                              stmesh::FLOAT_T(2.7)));
      REQUIRE_FALSE(triangulation.isGoodPoint(vertex_handles[1], stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.2),
                                              stmesh::FLOAT_T(1.0)));
      REQUIRE(triangulation.isGoodPoint(vertex_handles[0], stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.2),
                                        stmesh::FLOAT_T(0.5)));
      REQUIRE(triangulation.isGoodPoint(vertex_handles[1], stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.2),
                                        stmesh::FLOAT_T(0.5)));
    }

    for (const auto &vertex_handle : vertex_handles) {
      std::vector<stmesh::Triangulation<>::FullCellHandle> surrounding_full_cells =
          triangulation.surroundingFullCells(vertex_handle);
      std::vector<stmesh::Triangulation<>::FullCellConstHandle> const_surrounding_full_cells;
      std::transform(surrounding_full_cells.begin(), surrounding_full_cells.end(),
                     std::back_inserter(const_surrounding_full_cells),
                     [](const auto &full_cell) { return stmesh::Triangulation<>::FullCellConstHandle{&*full_cell}; });

      std::set<stmesh::Triangulation<>::FullCellConstHandle> cells_before;
      std::transform(triangulation.begin(), triangulation.end(), std::inserter(cells_before, cells_before.begin()),
                     [](const auto &full_cell) { return stmesh::Triangulation<>::FullCellConstHandle{&full_cell}; });
      const stmesh::Triangulation<>::FullCellHandle full_cell_handle = triangulation.remove(vertex_handle);
      std::set<stmesh::Triangulation<>::FullCellConstHandle> cells_after;
      std::transform(triangulation.begin(), triangulation.end(), std::inserter(cells_after, cells_after.begin()),
                     [](const auto &full_cell) { return stmesh::Triangulation<>::FullCellConstHandle{&full_cell}; });

      std::vector<stmesh::Triangulation<>::FullCellConstHandle> cells_before_not_after;
      std::set_difference(cells_before.begin(), cells_before.end(), cells_after.begin(), cells_after.end(),
                          std::back_inserter(cells_before_not_after));
      REQUIRE_THAT(const_surrounding_full_cells, Matchers::UnorderedEquals(cells_before_not_after));

      std::vector<stmesh::Triangulation<>::FullCellConstHandle> cells_after_not_before;
      std::set_difference(cells_after.begin(), cells_after.end(), cells_before.begin(), cells_before.end(),
                          std::back_inserter(cells_after_not_before));
      const std::vector<stmesh::Triangulation<>::FullCellHandle> uncommitted =
          triangulation.commitUncommitted(full_cell_handle);
      std::vector<stmesh::Triangulation<>::FullCellConstHandle> const_uncommitted;
      std::transform(uncommitted.begin(), uncommitted.end(), std::back_inserter(const_uncommitted),
                     [](const auto &full_cell) { return stmesh::Triangulation<>::FullCellConstHandle{&*full_cell}; });

      REQUIRE_THAT(const_uncommitted, Matchers::UnorderedEquals(cells_after_not_before));

      for (const auto &full_cell : triangulation)
        REQUIRE(full_cell.data().committed);
    }
  }
}

TEMPLATE_TEST_CASE("Test triangulation from MIXD", "[triangulation]", stmesh::TriangulationFromMixd,
                   stmesh::TriangulationFromMixdWithData) {
  std::vector<stmesh::Vector4F> mxyz{
      {0.0, 0.0, 0.0, 0.0}, {1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0}, {0.3, 0.3, 0.3, 1.0}, {0.3, 0.3, 0.3, -1.0},
  };
  std::vector<std::array<int, 5>> mien{
      {1, 2, 3, 4, 5},
      {1, 2, 3, 4, 6},
  };
  std::vector<std::array<int, 5>> mrng{
      {-1, 1, 1, 1, 1},
      {-1, 1, 1, 1, 1},
  };

  auto get_cell = [&](size_t i) {
    Eigen::Matrix<double, 4, 5> simplex;
    for (size_t j = 0; j < 5; ++j)
      simplex.col(static_cast<Eigen::Index>(j)) = mxyz[static_cast<size_t>(mien[i].at(j) - 1)];
    return stmesh::GeometricSimplex<4>(simplex);
  };

  TestType triangulation = [] {
    if constexpr (std::is_same_v<TestType, stmesh::TriangulationFromMixdWithData>)
      return TestType("data/sample_mixd/values.minf", "data/sample_mixd/values.mxyz", stmesh::kINSProblem);
    else
      return TestType("data/sample_mixd/values.minf");
  }();

  REQUIRE(triangulation.boundingBox().min() == stmesh::Vector4F{0.0, 0.0, 0.0, -1.0});
  REQUIRE(triangulation.boundingBox().max() == stmesh::Vector4F{1.0, 1.0, 1.0, 1.0});

  auto it = std::ranges::begin(triangulation);
  auto first_cell = *it++;
  auto second_cell = *it++;
  REQUIRE(it == std::ranges::end(triangulation));
  REQUIRE(first_cell.geometricSimplex() == get_cell(0));
  REQUIRE(second_cell.geometricSimplex() == get_cell(1));
  if constexpr (std::is_same_v<TestType, stmesh::TriangulationFromMixdWithData>) {
    REQUIRE(first_cell.geometricSimplex() == stmesh::GeometricSimplex<4>(first_cell.data()));
    REQUIRE(second_cell.geometricSimplex() == stmesh::GeometricSimplex<4>(second_cell.data()));
  }

  for (size_t j = 0; j < 5; ++j) {
    if (j == 4) {
      REQUIRE_FALSE(first_cell.isSurfaceSide(j));
      REQUIRE_FALSE(second_cell.isSurfaceSide(j));
    } else {
      REQUIRE(first_cell.isSurfaceSide(j));
      REQUIRE(second_cell.isSurfaceSide(j));
      REQUIRE(triangulation.findBoundaryRegion(first_cell, j) ==
              static_cast<size_t>(mrng[0].at(stmesh::mixd::kOmitted.at(j))));
      REQUIRE(triangulation.findBoundaryRegion(second_cell, j) ==
              static_cast<size_t>(mrng[1].at(stmesh::mixd::kOmitted.at(j))));
    }
  }
}
