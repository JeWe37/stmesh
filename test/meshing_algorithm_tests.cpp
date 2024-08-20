#include <catch2/catch_approx.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <limits>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

#include <boost/heap/fibonacci_heap.hpp>

#include "stmesh/geometric_simplex.hpp"
#include "stmesh/lfs_schemes.hpp"
#include "stmesh/meshing_algorithm.hpp"
#include "stmesh/meshing_cell.hpp"
#include "stmesh/radius_schemes.hpp"
#include "stmesh/sdf.hpp"
#include "stmesh/surface_adapters.hpp"
#include "stmesh/triangulation.hpp"
#include "stmesh/utility.hpp"

using namespace Catch;

constexpr stmesh::FLOAT_T kEps = std::numeric_limits<stmesh::FLOAT_T>::epsilon();

void verifyMeshingAlgorithm(auto &meshing_algorithm, bool disable_rule6 = false) {
  const stmesh::detail::Triangulation &triangulation = meshing_algorithm.triangulation();
  std::array<stmesh::detail::Rules, 7> rule_list{
      stmesh::detail::Rule1{}, stmesh::detail::Rule2{}, stmesh::detail::Rule3{},   stmesh::detail::Rule4{},
      stmesh::detail::Rule5{}, stmesh::detail::Rule6{}, stmesh::detail::Complete{}};
  using Heap = boost::heap::fibonacci_heap<stmesh::detail::Cell>;
  const Heap &queue = meshing_algorithm.queue();
  std::vector<Heap::handle_type> handles;
  handles.reserve(queue.size());
  for (auto it = queue.begin(); it != queue.end(); ++it)
    handles.emplace_back(queue.s_handle_from_iterator(it));
  size_t full_cell_count = 0;
  for (const auto &full_cell : triangulation) {
    REQUIRE(full_cell.data().committed);
    auto it = std::ranges::find(handles, full_cell.data().extra_data.cell_handle);
    REQUIRE(it != handles.end());
    stmesh::detail::Triangulation::FullCellHandle full_cell_handle{
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
        const_cast<stmesh::detail::Triangulation::FullCell *>(&full_cell)};
    REQUIRE((**it).full_cell == full_cell_handle);
    const auto rule_index = static_cast<size_t>(std::visit([](const auto &r) { return r.kIndex; }, (**it).rule));
    const unsigned rule_priority = std::visit([](const auto &r) { return r.priority(); }, (**it).rule);
    if (rule_index != stmesh::detail::Complete::kIndex || disable_rule6) {
      if (!std::visit([&](auto &r) { return r.check(meshing_algorithm, full_cell_handle, true); },
                      rule_list.at(rule_index)))
        WARN("Issue!");
      REQUIRE(std::visit([&](auto &r) { return r.check(meshing_algorithm, full_cell_handle, true); },
                         rule_list.at(rule_index)));
      REQUIRE(std::visit([&](auto &r) { return r.priority(); }, rule_list.at(rule_index)) == rule_priority);
    }
    for (size_t i = 0; i < rule_index - (disable_rule6 && rule_index >= stmesh::detail::Rule6::kIndex); ++i) {
      if (!std::visit([&](auto &r) { return !r.check(meshing_algorithm, full_cell_handle, true); }, rule_list.at(i)))
        WARN("Issue!");
      REQUIRE(
          std::visit([&](auto &r) { return !r.check(meshing_algorithm, full_cell_handle, true); }, rule_list.at(i)));
    }
    full_cell_count++;
  }
  REQUIRE(handles.size() == full_cell_count);
}

TEST_CASE("Test meshing algorithm base functionality", "[meshing_base][meshing_algorithm]") {
  const stmesh::SDFSurfaceAdapter<stmesh::HyperSphere4> sdf_surface_adapter(
      stmesh::FLOAT_T(1.0),
      stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)});
  stmesh::MeshingAlgorithm meshing_algorithm(
      sdf_surface_adapter, stmesh::FLOAT_T(20.0), stmesh::FLOAT_T(0.001), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(5.0),
      stmesh::lfs_schemes::Constant(stmesh::FLOAT_T(2.0)), stmesh::radius_schemes::Constant(stmesh::FLOAT_T(0.5)));
  SECTION("Correctly constructed") {
    REQUIRE(meshing_algorithm.triangulation().boundingBox().min() == -stmesh::Vector4F::Constant(5.0));
    REQUIRE(meshing_algorithm.triangulation().boundingBox().max() == stmesh::Vector4F::Constant(5.0));
    REQUIRE(meshing_algorithm.deltaSurface(stmesh::Vector4F::Zero()) == stmesh::FLOAT_T(2.0));
    verifyMeshingAlgorithm(meshing_algorithm);
  }
  SECTION("Test insertion") {
    const Eigen::Matrix<stmesh::FLOAT_T, 4, 4> points{
        {stmesh::FLOAT_T(0.3), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
        {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.3), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
        {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.3), stmesh::FLOAT_T(0.0)},
        {stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.3)}};
    std::vector<stmesh::detail::Triangulation::VertexHandle> vertex_handles;
    for (const auto &point : points.colwise()) {
      vertex_handles.push_back(meshing_algorithm.insert(point));
      verifyMeshingAlgorithm(meshing_algorithm);
    }
    SECTION("Test voronoi dual") {
      std::pair<stmesh::detail::Triangulation::FullCellHandle, int> dependent_neighbor_info;
      std::optional voronoi_dual = meshing_algorithm.voronoiDual(points, &dependent_neighbor_info);
      REQUIRE(voronoi_dual.has_value());
      REQUIRE((*voronoi_dual -
               stmesh::Vector4F{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0)}
                   .normalized())
                  .cwiseAbs()
                  .maxCoeff() < 4.0 * kEps);
      const auto &[full_cell, index_of_covertex] = dependent_neighbor_info;
      for (const auto &vertex_handle : vertex_handles) {
        int i{};
        REQUIRE(full_cell->has_vertex(vertex_handle, i));
        REQUIRE(i != index_of_covertex);
      }

      SECTION("Other direction") {
        const stmesh::Vector4F additional_point = {-stmesh::FLOAT_T(0.3), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0),
                                                   stmesh::FLOAT_T(0.0)};
        Eigen::Matrix<stmesh::FLOAT_T, 4, 4> points2;
        points2 << additional_point, points.rightCols<3>();
        meshing_algorithm.insert(additional_point);
        verifyMeshingAlgorithm(meshing_algorithm);
        std::optional voronoi_dual2 = meshing_algorithm.voronoiDual(points2);
        REQUIRE(voronoi_dual2.has_value());
        REQUIRE((*voronoi_dual2 - stmesh::Vector4F{-stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0),
                                                   stmesh::FLOAT_T(1.0)}
                                      .normalized())
                    .cwiseAbs()
                    .maxCoeff() < 4.0 * kEps);
      }

      SECTION("Test surface ball") {
        stmesh::HyperSphere4 surface_ball = meshing_algorithm.surfaceBall(points);
        REQUIRE(surface_ball.radius() ==
                Approx((points.col(0) - stmesh::Vector4F{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0),
                                                         stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0)}
                                            .normalized())
                           .norm()));
        REQUIRE(surface_ball.center() == *voronoi_dual);
        SECTION("Test picking region") {
          surface_ball.scale(stmesh::FLOAT_T(0.5));
          SECTION("Test 3-simplex picking region sampling") {
            for (int i = 0; i < 100; ++i) {
              const auto [point, radius] = meshing_algorithm.sampleFromPickingRegion(points);
              REQUIRE(surface_ball.signedDistance(point) < stmesh::FLOAT_T(0.0));
              REQUIRE(point.norm() == Approx(stmesh::FLOAT_T(1.0)));
              REQUIRE(radius == surface_ball.radius());
            }
          }
          SECTION("Pick good point from 3-simplex") {
            const auto vertex_handle = meshing_algorithm.pickGoodPoint(points);
            const auto pt = vertex_handle->point();
            const stmesh::Vector4F point{static_cast<stmesh::FLOAT_T>(pt[0]), static_cast<stmesh::FLOAT_T>(pt[1]),
                                         static_cast<stmesh::FLOAT_T>(pt[2]), static_cast<stmesh::FLOAT_T>(pt[3])};
            REQUIRE(surface_ball.signedDistance(point) < stmesh::FLOAT_T(0.0));
            REQUIRE(point.norm() == Approx(stmesh::FLOAT_T(1.0)));
            REQUIRE(meshing_algorithm.triangulation().isGoodPoint(vertex_handle, stmesh::FLOAT_T(20.0),
                                                                  stmesh::FLOAT_T(0.001),
                                                                  stmesh::FLOAT_T(5.0) * surface_ball.radius()));
            verifyMeshingAlgorithm(meshing_algorithm);
          }
        }
      }
    }

    Eigen::Matrix<stmesh::FLOAT_T, 4, 5> points2;
    const stmesh::Vector4F zero_seven{stmesh::FLOAT_T(0.7), stmesh::FLOAT_T(0.7), stmesh::FLOAT_T(0.7),
                                      stmesh::FLOAT_T(0.7)};
    points2 << points, zero_seven;
    meshing_algorithm.insert(zero_seven);
    verifyMeshingAlgorithm(meshing_algorithm);
    meshing_algorithm.insert(-zero_seven);
    verifyMeshingAlgorithm(meshing_algorithm);

    REQUIRE_FALSE(meshing_algorithm.voronoiDual(points).has_value());
    SECTION("Test 4-simplex picking region sampling") {
      stmesh::HyperSphere4 picking_region = stmesh::GeometricSimplex<4>(points2).circumsphere();
      picking_region.scale(stmesh::FLOAT_T(0.5));

      for (int i = 0; i < 100; ++i) {
        const auto [point, radius] = meshing_algorithm.sampleFromPickingRegion(points2);
        REQUIRE(picking_region.signedDistance(point) < stmesh::FLOAT_T(0.0));
        REQUIRE(radius == picking_region.radius());
      }
      SECTION("Pick good point from 4-simplex") {
        const auto vertex_handle = meshing_algorithm.pickGoodPoint(points2);
        const auto pt = vertex_handle->point();
        const stmesh::Vector4F point{static_cast<stmesh::FLOAT_T>(pt[0]), static_cast<stmesh::FLOAT_T>(pt[1]),
                                     static_cast<stmesh::FLOAT_T>(pt[2]), static_cast<stmesh::FLOAT_T>(pt[3])};
        REQUIRE(picking_region.signedDistance(point) < stmesh::FLOAT_T(0.0));
        REQUIRE(meshing_algorithm.triangulation().isGoodPoint(vertex_handle, stmesh::FLOAT_T(20.0),
                                                              stmesh::FLOAT_T(0.001),
                                                              stmesh::FLOAT_T(5.0) * picking_region.radius()));
        verifyMeshingAlgorithm(meshing_algorithm);
      }
    }
    for (const auto &vertex_handle : vertex_handles) {
      meshing_algorithm.remove(vertex_handle);
      verifyMeshingAlgorithm(meshing_algorithm);
    }
  }
}

template <typename... Ts> struct Overload : Ts... {
  using Ts::operator()...;
};
template <class... Ts> Overload(Ts...) -> Overload<Ts...>;

TEST_CASE("Test meshing rules", "[meshing_rules][meshing_algorithm]") {
  const stmesh::SDFSurfaceAdapter<stmesh::HyperSphere4> sdf_surface_adapter(
      stmesh::FLOAT_T(40.0),
      stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)});
  stmesh::MeshingAlgorithm meshing_algorithm(
      sdf_surface_adapter, stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(5.0),
      stmesh::lfs_schemes::Constant(stmesh::FLOAT_T(5.0)), stmesh::radius_schemes::Constant(stmesh::FLOAT_T(30.0)));
  for (const auto &point :
       {stmesh::Vector4F{stmesh::FLOAT_T(28.4), stmesh::FLOAT_T(11.6), stmesh::FLOAT_T(-32), stmesh::FLOAT_T(39.2)},
        stmesh::Vector4F{stmesh::FLOAT_T(-37.2), stmesh::FLOAT_T(-22), stmesh::FLOAT_T(-32.8), stmesh::FLOAT_T(30.8)},
        stmesh::Vector4F{stmesh::FLOAT_T(25.6), stmesh::FLOAT_T(32.8), stmesh::FLOAT_T(-35.6), stmesh::FLOAT_T(-13.6)},
        stmesh::Vector4F{stmesh::FLOAT_T(33.2), stmesh::FLOAT_T(-14), stmesh::FLOAT_T(-30.4), stmesh::FLOAT_T(-8)},
        stmesh::Vector4F{stmesh::FLOAT_T(12.4), stmesh::FLOAT_T(1.2), stmesh::FLOAT_T(27.6), stmesh::FLOAT_T(38.4)},
        stmesh::Vector4F{stmesh::FLOAT_T(-10.8), stmesh::FLOAT_T(-23.2), stmesh::FLOAT_T(16), stmesh::FLOAT_T(-4.4)},
        stmesh::Vector4F{stmesh::FLOAT_T(-17.6), stmesh::FLOAT_T(8), stmesh::FLOAT_T(-0.4), stmesh::FLOAT_T(-28)},
        stmesh::Vector4F{stmesh::FLOAT_T(4.4), stmesh::FLOAT_T(2.8), stmesh::FLOAT_T(-26), stmesh::FLOAT_T(-7.2)},
        stmesh::Vector4F{stmesh::FLOAT_T(12.8), stmesh::FLOAT_T(36), stmesh::FLOAT_T(-20), stmesh::FLOAT_T(18.4)},
        stmesh::Vector4F{stmesh::FLOAT_T(-16.8), stmesh::FLOAT_T(4), stmesh::FLOAT_T(-32), stmesh::FLOAT_T(22.4)},
        stmesh::Vector4F{stmesh::FLOAT_T(33.6), stmesh::FLOAT_T(4.8), stmesh::FLOAT_T(-34.4), stmesh::FLOAT_T(4)},
        stmesh::Vector4F{stmesh::FLOAT_T(39.6), stmesh::FLOAT_T(12.8), stmesh::FLOAT_T(-14), stmesh::FLOAT_T(-19.6)},
        stmesh::Vector4F{stmesh::FLOAT_T(19.2), stmesh::FLOAT_T(6.4), stmesh::FLOAT_T(22), stmesh::FLOAT_T(16.4)},
        stmesh::Vector4F{stmesh::FLOAT_T(-29.2), stmesh::FLOAT_T(38.4), stmesh::FLOAT_T(-27.6), stmesh::FLOAT_T(10.8)},
        stmesh::Vector4F{stmesh::FLOAT_T(-24.4), stmesh::FLOAT_T(-40), stmesh::FLOAT_T(19.2), stmesh::FLOAT_T(-40)},
        stmesh::Vector4F{stmesh::FLOAT_T(-15.6), stmesh::FLOAT_T(19.6), stmesh::FLOAT_T(10.8), stmesh::FLOAT_T(13.6)},
        stmesh::Vector4F{stmesh::FLOAT_T(-11.2), stmesh::FLOAT_T(-4), stmesh::FLOAT_T(29.6), stmesh::FLOAT_T(-7.6)},
        stmesh::Vector4F{stmesh::FLOAT_T(30), stmesh::FLOAT_T(-20), stmesh::FLOAT_T(10.4), stmesh::FLOAT_T(0)},
        stmesh::Vector4F{stmesh::FLOAT_T(-20.8), stmesh::FLOAT_T(30.4), stmesh::FLOAT_T(4.8), stmesh::FLOAT_T(-18.4)},
        stmesh::Vector4F{stmesh::FLOAT_T(-26.8), stmesh::FLOAT_T(9.6), stmesh::FLOAT_T(-24.8), stmesh::FLOAT_T(-6.8)},
        stmesh::Vector4F{stmesh::FLOAT_T(-20.8), stmesh::FLOAT_T(8.4), stmesh::FLOAT_T(-11.2), stmesh::FLOAT_T(28.8)},
        stmesh::Vector4F{stmesh::FLOAT_T(-27.2), stmesh::FLOAT_T(-34), stmesh::FLOAT_T(16), stmesh::FLOAT_T(18)},
        stmesh::Vector4F{stmesh::FLOAT_T(-30.4), stmesh::FLOAT_T(-6.4), stmesh::FLOAT_T(-16), stmesh::FLOAT_T(-12.4)},
        stmesh::Vector4F{stmesh::FLOAT_T(28), stmesh::FLOAT_T(-23.6), stmesh::FLOAT_T(-10.8), stmesh::FLOAT_T(-20.8)},
        stmesh::Vector4F{stmesh::FLOAT_T(14), stmesh::FLOAT_T(34.4), stmesh::FLOAT_T(-1.2), stmesh::FLOAT_T(-37.2)},
        stmesh::Vector4F{stmesh::FLOAT_T(-18.4), stmesh::FLOAT_T(38.8), stmesh::FLOAT_T(-28.8), stmesh::FLOAT_T(-7.6)},
        stmesh::Vector4F{stmesh::FLOAT_T(34.8), stmesh::FLOAT_T(19.2), stmesh::FLOAT_T(-15.6), stmesh::FLOAT_T(28.8)},
        stmesh::Vector4F{stmesh::FLOAT_T(3.2), stmesh::FLOAT_T(24.8), stmesh::FLOAT_T(-17.2), stmesh::FLOAT_T(-4.4)},
        stmesh::Vector4F{stmesh::FLOAT_T(34.8), stmesh::FLOAT_T(18.8), stmesh::FLOAT_T(-32.8), stmesh::FLOAT_T(-26.4)},
        stmesh::Vector4F{stmesh::FLOAT_T(-23.2), stmesh::FLOAT_T(-36), stmesh::FLOAT_T(-6.4), stmesh::FLOAT_T(-3.6)},
        stmesh::Vector4F{stmesh::FLOAT_T(-40), stmesh::FLOAT_T(-10), stmesh::FLOAT_T(-34), stmesh::FLOAT_T(34)}})
    meshing_algorithm.insert(point);
  for (const auto &point :
       {stmesh::Vector4F{stmesh::FLOAT_T(38.4), stmesh::FLOAT_T(15.2), stmesh::FLOAT_T(-39.2), stmesh::FLOAT_T(17.6)},
        stmesh::Vector4F{stmesh::FLOAT_T(-10.4), stmesh::FLOAT_T(12), stmesh::FLOAT_T(30.4), stmesh::FLOAT_T(33.2)},
        stmesh::Vector4F{stmesh::FLOAT_T(-39.2), stmesh::FLOAT_T(-29.6), stmesh::FLOAT_T(-15.2),
                         stmesh::FLOAT_T(-38.8)}})
    meshing_algorithm.insert(sdf_surface_adapter.closestPoint(point), {}, true);
  SECTION("Check rules") {
    verifyMeshingAlgorithm(meshing_algorithm);
    const stmesh::detail::Triangulation &triangulation = meshing_algorithm.triangulation();
    bool checked_rule1{}, checked_rule2{}, checked_rule3{}, checked_rule4{}, checked_rule5{}, checked_rule6{};
    for (const auto &cell : meshing_algorithm.queue())
      std::visit(Overload{
                     [&](const stmesh::detail::Rule1 &r) {
                       const stmesh::GeometricSimplex<4> simplex = triangulation.fullCellSimplex(cell.full_cell);
                       const stmesh::HyperSphere4 circumsphere = simplex.circumsphere();
                       REQUIRE(sdf_surface_adapter.intersectedBySphere(circumsphere));
                       const stmesh::Vector4F z = sdf_surface_adapter.closestPoint(circumsphere.center());
                       REQUIRE(z == r.z0);
                       REQUIRE(std::ranges::all_of(triangulation.verticesInRadius(z, meshing_algorithm.deltaSurface(z)),
                                                   [](const auto &vertex) { return !vertex->data().nonfree_vertex; }));
                       checked_rule1 = true;
                     },
                     [&](const stmesh::detail::Rule2 &r) {
                       const stmesh::GeometricSimplex<4> simplex = triangulation.fullCellSimplex(cell.full_cell);
                       const stmesh::HyperSphere4 circumsphere = simplex.circumsphere();
                       REQUIRE(sdf_surface_adapter.intersectedBySphere(circumsphere));
                       stmesh::Vector4F z = circumsphere.center();
                       REQUIRE_FALSE(triangulation.boundingBox().contains(z) == r.z_outside);
                       if (r.z_outside)
                         REQUIRE(
                             r.z ==
                             z.cwiseMax(triangulation.boundingBox().min()).cwiseMin(triangulation.boundingBox().max()));
                       else
                         REQUIRE(z == r.z);
                       REQUIRE(2 * meshing_algorithm.deltaSurface(sdf_surface_adapter.closestPoint(z)) <=
                               circumsphere.radius());
                       checked_rule2 = true;
                     },
                     [&](const stmesh::detail::Rule3 &r) {
                       const stmesh::GeometricSimplex<4> simplex = triangulation.fullCellSimplex(cell.full_cell);
                       const stmesh::HyperSphere4 circumsphere = simplex.circumsphere();
                       const stmesh::Vector4F &z = circumsphere.center();
                       REQUIRE(sdf_surface_adapter.inside(z));
                       REQUIRE(circumsphere.radius() > stmesh::FLOAT_T(30.0));
                       REQUIRE(r.z == z);
                       checked_rule3 = true;
                     },
                     [&](const stmesh::detail::Rule4 &r) {
                       const stmesh::GeometricSimplex<4> simplex = triangulation.fullCellSimplex(cell.full_cell);
                       const stmesh::HyperSphere4 circumsphere = simplex.circumsphere();
                       const stmesh::Vector4F &z = circumsphere.center();
                       REQUIRE(sdf_surface_adapter.inside(z));
                       REQUIRE(simplex.radiusEdgeRatio() >= stmesh::FLOAT_T(1.0));
                       REQUIRE(r.z == z);
                       checked_rule4 = true;
                     },
                     [&](const stmesh::detail::Rule5 &r) {
                       const stmesh::GeometricSimplex<4> simplex = triangulation.fullCellSimplex(cell.full_cell);
                       const stmesh::HyperSphere4 circumsphere = simplex.circumsphere();
                       REQUIRE(sdf_surface_adapter.inside(circumsphere.center()));
                       REQUIRE(simplex.sliverSimplex(stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.5)));
                       REQUIRE(r.vertices == simplex.vertices());
                       checked_rule5 = true;
                     },
                     [&](const stmesh::detail::Rule6 &r) {
                       const stmesh::GeometricSimplex<4> simplex = triangulation.fullCellSimplex(cell.full_cell);
                       auto sub_simplices = simplex.subSimplices<4>();
                       REQUIRE(std::ranges::find_if(sub_simplices, [&](const auto &sub_simplex) {
                                 return sub_simplex.vertices() == r.vertices;
                               }) != sub_simplices.end());
                       REQUIRE(meshing_algorithm.voronoiDual(r.vertices).has_value());
                       auto facet = meshing_algorithm.triangulation().facetFromVertices(r.vertices);
                       REQUIRE(std::any_of(facet.full_cell()->vertices_begin(), facet.full_cell()->vertices_end(),
                                           [&, i = 0](const auto &vertex) mutable {
                                             return i++ != facet.index_of_covertex() && !vertex->data().nonfree_vertex;
                                           }));
                       checked_rule6 = true;
                     },
                     [&](const stmesh::detail::Complete & /* unused */) {},
                 },
                 cell.rule);
    REQUIRE(checked_rule1);
    REQUIRE(checked_rule2);
    REQUIRE(checked_rule3);
    REQUIRE(checked_rule4);
    REQUIRE(checked_rule5);
    REQUIRE(checked_rule6);
  }
  const size_t size = meshing_algorithm.queue().size();
  auto i = GENERATE_COPY(Catch::Generators::range(static_cast<size_t>(0), size));
  DYNAMIC_SECTION("Test applying rule: " << i) {
    auto it = meshing_algorithm.queue().begin();
    std::advance(it, i);
    std::visit(
        [&]<typename Rule>(const Rule &r) {
          // picking region based rules cannot be tested, picking may not terminate for unbounded aspect ratio
          // similarly, if Complete were applied it would throw an error, so also ignore that
          if constexpr (Rule::kIndex < stmesh::detail::Rule5::kIndex)
            r.apply(meshing_algorithm, it->full_cell);
        },
        it->rule);
    verifyMeshingAlgorithm(meshing_algorithm);
  }
}

TEST_CASE("End to end sphere meshing", "[sphere_meshing][meshing_algorithm]") {
  const stmesh::SDFSurfaceAdapter<stmesh::HyperSphere4> sdf_surface_adapter(
      stmesh::FLOAT_T(0.2),
      stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)});
  stmesh::MeshingAlgorithm meshing_algorithm(
      sdf_surface_adapter, stmesh::FLOAT_T(20.0), stmesh::FLOAT_T(1e-6), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(5.0),
      stmesh::lfs_schemes::Constant(stmesh::FLOAT_T(0.15)), stmesh::radius_schemes::Constant(stmesh::FLOAT_T(0.12)),
      std::nullopt, true);
  meshing_algorithm.triangulate([&] { verifyMeshingAlgorithm(meshing_algorithm, true); });
}
