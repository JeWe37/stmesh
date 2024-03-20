#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <tuple>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

#include "stmesh/geometric_simplex.hpp"
#include "stmesh/sdf.hpp"
#include "stmesh/utility.hpp"

using namespace Catch;

constexpr stmesh::FLOAT_T kEps = std::numeric_limits<stmesh::FLOAT_T>::epsilon();

TEST_CASE("Test 4-simplex in 4D", "[4simplex][geometric_simplex]") {
  const std::vector<stmesh::Vector4F> simplex_vertices{
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}};
  const stmesh::GeometricSimplex<4> simplex(simplex_vertices);
  SECTION("Correctly constructed") {
    REQUIRE(simplex.boundingBox().min() ==
            stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)});
    REQUIRE(simplex.boundingBox().max() ==
            stmesh::Vector4F{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0)});
  }
  SECTION("Basic geometric constructions") { REQUIRE(simplex.shortestEdgeLength() == stmesh::FLOAT_T(1.0)); }
  SECTION("Correct Caley-Menger calculations") {
    Eigen::Matrix<stmesh::FLOAT_T, 6, 6> matrix = Eigen::Matrix<stmesh::FLOAT_T, 6, 6>::Ones();
    matrix(Eigen::seq(2, Eigen::last), Eigen::seq(2, Eigen::last)) += Eigen::Matrix<stmesh::FLOAT_T, 4, 4>::Ones();
    matrix.diagonal() = Eigen::Matrix<stmesh::FLOAT_T, 6, 1>::Zero();
    REQUIRE(simplex.caleyMengerMatrix() == matrix);
    REQUIRE(simplex.caleyMengerDeterminant() == -stmesh::FLOAT_T(16.0));
    REQUIRE(simplex.content() == stmesh::FLOAT_T(1.0) / stmesh::FLOAT_T(24.0));
    SECTION("Correct circumcenter") {
      const stmesh::HyperSphere<4> circumsphere = simplex.circumsphere();
      for (const stmesh::Vector4F vertex : simplex.vertices().colwise())
        // Catch::Approx does not work well with stmesh::FLOAT_T(0.0)
        REQUIRE(circumsphere.distance(vertex) < std::sqrt(kEps));
    }
  }
  SECTION("Correct quality metrics") {
    REQUIRE(simplex.radiusEdgeRatio() == Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(simplex.quality() == Approx(stmesh::FLOAT_T(1.0) / stmesh::FLOAT_T(24.0)));
  }
  SECTION("Finds all sub simplices") {
    // clang-format off
    std::vector faces1{
        Eigen::Matrix<stmesh::FLOAT_T, 2, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 2, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 2, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 2, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 2, 4>{{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 2, 4>{{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 2, 4>{{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 2, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 2, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 2, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
      };
    std::vector faces2{
        Eigen::Matrix<stmesh::FLOAT_T, 3, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 3, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 3, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 3, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 3, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 3, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 3, 4>{{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 3, 4>{{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 3, 4>{{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 3, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
      };
    std::vector faces3{
        Eigen::Matrix<stmesh::FLOAT_T, 4, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 4, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 4, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 4, 4>{{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval(),
        Eigen::Matrix<stmesh::FLOAT_T, 4, 4>{{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}, {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}}.transpose().eval()
      };
    // clang-format on
    SECTION("3D faces") {
      const auto sub_simplices = simplex.subSimplices<4>();
      STATIC_REQUIRE(sub_simplices.size() == 5U);
      std::vector<Eigen::Matrix<stmesh::FLOAT_T, 4, 4>> vertices;
      std::ranges::transform(sub_simplices, std::back_inserter(vertices),
                             [](const auto &sub_simplex) { return sub_simplex.vertices(); });
      REQUIRE_THAT(vertices, Matchers::UnorderedRangeEquals(faces3));
    }
    SECTION("2D faces") {
      const auto sub_simplices = simplex.subSimplices<3>();
      STATIC_REQUIRE(sub_simplices.size() == 10U);
      std::vector<Eigen::Matrix<stmesh::FLOAT_T, 4, 3>> vertices;
      std::ranges::transform(sub_simplices, std::back_inserter(vertices),
                             [](const auto &sub_simplex) { return sub_simplex.vertices(); });
      REQUIRE_THAT(vertices, Matchers::UnorderedRangeEquals(faces2));
    }
    SECTION("1D faces") {
      const auto sub_simplices = simplex.subSimplices<2>();
      STATIC_REQUIRE(sub_simplices.size() == 10U);
      std::vector<Eigen::Matrix<stmesh::FLOAT_T, 4, 2>> vertices;
      std::ranges::transform(sub_simplices, std::back_inserter(vertices),
                             [](const auto &sub_simplex) { return sub_simplex.vertices(); });
      REQUIRE_THAT(vertices, Matchers::UnorderedRangeEquals(faces1));
    }
    SECTION("All at once") {
      const auto [sub_simplices1, sub_simplices2, sub_simplices3] = simplex.allSubSimplicies<4, 2>();
      STATIC_REQUIRE(std::tuple_size_v<decltype(sub_simplices3)> == 5U);
      STATIC_REQUIRE(std::tuple_size_v<decltype(sub_simplices2)> == 10U);
      STATIC_REQUIRE(std::tuple_size_v<decltype(sub_simplices1)> == 10U);
      std::vector<Eigen::Matrix<stmesh::FLOAT_T, 4, 4>> vertices3;
      std::ranges::transform(sub_simplices3, std::back_inserter(vertices3),
                             [](const auto &sub_simplex) { return sub_simplex.vertices(); });
      REQUIRE_THAT(vertices3, Matchers::UnorderedRangeEquals(faces3));
      std::vector<Eigen::Matrix<stmesh::FLOAT_T, 4, 3>> vertices2;
      std::ranges::transform(sub_simplices2, std::back_inserter(vertices2),
                             [](const auto &sub_simplex) { return sub_simplex.vertices(); });
      REQUIRE_THAT(vertices2, Matchers::UnorderedRangeEquals(faces2));
      std::vector<Eigen::Matrix<stmesh::FLOAT_T, 4, 2>> vertices1;
      std::ranges::transform(sub_simplices1, std::back_inserter(vertices1),
                             [](const auto &sub_simplex) { return sub_simplex.vertices(); });
      REQUIRE_THAT(vertices1, Matchers::UnorderedRangeEquals(faces1));
    }
  }
  SECTION("Correctly cuts by plane") {
    const stmesh::Vector4F normal{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0),
                                  stmesh::FLOAT_T(1.0)};
    const stmesh::Vector4F point{0.0, 0.0, 0.0, stmesh::FLOAT_T(0.5)};
    const Eigen::Hyperplane<stmesh::FLOAT_T, 4> plane{normal, point};
    const CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> polyhedron = simplex.planeCut(plane);
    std::vector<stmesh::Vector3F> vertices;
    std::transform(polyhedron.vertices_begin(), polyhedron.vertices_end(), std::back_inserter(vertices),
                   [](const auto &vertex) {
                     const auto &pt = vertex.point();
                     return stmesh::Vector3F{pt.x(), pt.y(), pt.z()};
                   });
    const std::vector<stmesh::Vector3F> expected_vertices{
        {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
        {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.5)},
        {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.0)},
        {stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}};
    REQUIRE_THAT(vertices, Matchers::UnorderedEquals(expected_vertices));
    REQUIRE(polyhedron.size_of_facets() == 4U);
    for (const auto &facet_handle : polyhedron.facet_handles())
      REQUIRE(facet_handle->facet_degree() == 3);
  }
}

TEST_CASE("4D sliver tests", "[sliver][geometric_simplex]") {
  const stmesh::GeometricSimplex<4> simplex(std::vector<stmesh::Vector4F>{
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.9), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}});
  const auto sub_simplices = simplex.subSimplices<4>();
  SECTION("Count well-shaped") {
    REQUIRE(std::ranges::count_if(sub_simplices, [](const auto &sub_simplex) {
              return sub_simplex.wellShaped(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.2));
            }) == 3);
    REQUIRE(std::ranges::count_if(sub_simplices, [](const auto &sub_simplex) {
              return sub_simplex.wellShaped(stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.2));
            }) == 0);
    REQUIRE(std::ranges::count_if(sub_simplices, [](const auto &sub_simplex) {
              return sub_simplex.wellShaped(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.1));
            }) == 5);
  }
  SECTION("Test sliver 4D") {
    REQUIRE(simplex.sliver(stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.1)));
    REQUIRE_FALSE(simplex.sliver(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.2)));
    REQUIRE_FALSE(simplex.sliver(stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.01)));
    REQUIRE_FALSE(simplex.sliver(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.1)));
  }
  SECTION("Test sliver 3D") {
    const stmesh::GeometricSimplex<4, 4> &sub_simplex = sub_simplices[4];
    REQUIRE(sub_simplex.sliver(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.3)));
    REQUIRE_FALSE(sub_simplex.sliver(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.1)));
    REQUIRE_FALSE(sub_simplex.sliver(stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.2)));
  }
  SECTION("Test sliver simplex") {
    REQUIRE(simplex.sliverSimplex(stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.1)) == 4U);
    REQUIRE(simplex.sliverSimplex(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.3)) == 3U);
    REQUIRE(simplex.sliverSimplex(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.01)) == 0U);
    REQUIRE(simplex.sliverSimplex(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.1)) == 0U);
  }
  SECTION("Test small sliver simplex") {
    REQUIRE(simplex.smallSliverSimplex(stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.1), stmesh::FLOAT_T(2.0)));
    REQUIRE(simplex.smallSliverSimplex(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.3), stmesh::FLOAT_T(2.0)));
    REQUIRE_FALSE(simplex.smallSliverSimplex(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.01), stmesh::FLOAT_T(2.0)));
    REQUIRE_FALSE(simplex.smallSliverSimplex(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.1), stmesh::FLOAT_T(2.0)));
    REQUIRE_FALSE(simplex.smallSliverSimplex(stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.1), stmesh::FLOAT_T(0.5)));
    REQUIRE_FALSE(simplex.smallSliverSimplex(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.2), stmesh::FLOAT_T(0.5)));
    REQUIRE_FALSE(simplex.smallSliverSimplex(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.01), stmesh::FLOAT_T(0.5)));
    REQUIRE_FALSE(simplex.smallSliverSimplex(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.1), stmesh::FLOAT_T(0.5)));
  }
}

TEST_CASE("Test 3-simplex in 4D", "[3simplex][geometric_simplex]") {
  const std::vector<stmesh::Vector4F> simplex_vertices{
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
      {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}};
  const stmesh::GeometricSimplex<4, 4> simplex(simplex_vertices);
  SECTION("Correctly constructed") {
    REQUIRE(simplex.boundingBox().min() ==
            stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)});
    REQUIRE(simplex.boundingBox().max() ==
            stmesh::Vector4F{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)});
  }
  SECTION("Basic geometric constructions") { REQUIRE(simplex.shortestEdgeLength() == stmesh::FLOAT_T(1.0)); }
  SECTION("Correct Caley-Menger calculations") {
    Eigen::Matrix<stmesh::FLOAT_T, 5, 5> matrix = Eigen::Matrix<stmesh::FLOAT_T, 5, 5>::Ones();
    matrix(Eigen::seq(2, Eigen::last), Eigen::seq(2, Eigen::last)) += Eigen::Matrix<stmesh::FLOAT_T, 3, 3>::Ones();
    matrix.diagonal() = Eigen::Matrix<stmesh::FLOAT_T, 5, 1>::Zero();
    REQUIRE(simplex.caleyMengerMatrix() == matrix);
    REQUIRE(simplex.caleyMengerDeterminant() == stmesh::FLOAT_T(8.0));
    REQUIRE(simplex.content() == stmesh::FLOAT_T(1.0) / stmesh::FLOAT_T(6.0));
    SECTION("Correct circumcenter") {
      const stmesh::HyperSphere<4> circumsphere = simplex.circumsphere();
      for (const stmesh::Vector4F vertex : simplex.vertices().colwise())
        // Catch::Approx does not work well with stmesh::FLOAT_T(0.0)
        REQUIRE(circumsphere.distance(vertex) < std::sqrt(kEps));
    }
  }
  SECTION("Correct quality metrics") {
    REQUIRE(simplex.radiusEdgeRatio() == Approx(std::sqrt(0.75)));
    REQUIRE(simplex.quality() == Approx(stmesh::FLOAT_T(1.0) / stmesh::FLOAT_T(6.0)));
  }
  SECTION("Check correct normal ray") {
    const Eigen::ParametrizedLine<stmesh::FLOAT_T, 4> ray = simplex.normalRay();
    REQUIRE((ray.origin() -
             stmesh::Vector4F{stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.0)})
                .cwiseAbs()
                .maxCoeff() < std::sqrt(kEps));
    REQUIRE(ray.direction() ==
            stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)});
  }
}