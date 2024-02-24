#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <memory>

#include "stmesh/edt.hpp"
#include "stmesh/sdf.hpp"
#include "stmesh/surface_adapters.hpp"
#include "stmesh/utility.hpp"

constexpr stmesh::FLOAT_T kEps = std::numeric_limits<stmesh::FLOAT_T>::epsilon();

TEST_CASE("Test sdf adapter with sphere sdf", "[sdf_adapter][surface_adapter]") {
  const auto radius = stmesh::FLOAT_T(1.0);
  const stmesh::Vector3F center{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)};
  const stmesh::SDFSurfaceAdapter<stmesh::HyperSphere3> adapter(radius, center);
  SECTION("Correctly constructed") {
    REQUIRE(adapter.boundingBox().min() ==
            stmesh::Vector3F{-stmesh::FLOAT_T(1.0), -stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)});
    REQUIRE(adapter.boundingBox().max() ==
            stmesh::Vector3F{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(2.0)});
  }
  SECTION("Correctly determines if inside") {
    REQUIRE_FALSE(adapter.inside({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}));
    REQUIRE(adapter.inside({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}));
    REQUIRE_FALSE(adapter.inside({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0)}));
    REQUIRE_FALSE(adapter.inside({stmesh::FLOAT_T(0.0), -stmesh::FLOAT_T(3.0), stmesh::FLOAT_T(0.0)}));
  }
  SECTION("Correctly computes closest point") {
    REQUIRE(adapter.closestPoint({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(1.0)}) ==
            stmesh::Vector3F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0)});
    REQUIRE(adapter.closestPoint({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(3.0)}) ==
            stmesh::Vector3F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0)});
    REQUIRE(adapter.closestPoint({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0)}) ==
            stmesh::Vector3F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0)});
    REQUIRE((adapter.closestPoint({stmesh::FLOAT_T(0.0), -stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(3.0)}) - center -
             stmesh::Vector3F{stmesh::FLOAT_T(0.0), -1.0, 1.0}.normalized())
                .cwiseAbs()
                .maxCoeff() < kEps);
  }
  SECTION("Correctly checks sphere intersection") {
    const auto radius2 = stmesh::FLOAT_T(0.25);
    const stmesh::Vector3F center2{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)};
    const stmesh::HyperSphere3 sphere2{radius2, center2};
    REQUIRE(adapter.intersectedBySphere(sphere2));
    const auto radius3 = stmesh::FLOAT_T(0.25);
    const stmesh::Vector3F center3{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), -stmesh::FLOAT_T(1.0)};
    const stmesh::HyperSphere3 sphere3{radius3, center3};
    REQUIRE_FALSE(adapter.intersectedBySphere(sphere3));
  }
  SECTION("Correctly computes raycast") {
    const auto max_distance = stmesh::FLOAT_T(10.0);
    const stmesh::Vector3F direction1{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)};
    const stmesh::Vector3F point1{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)};
    REQUIRE(*adapter.raycast(point1, direction1, max_distance) ==
            stmesh::Vector3F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)});
    const stmesh::Vector3F direction2{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), -stmesh::FLOAT_T(1.0)};
    const stmesh::Vector3F point2{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(3.0)};
    REQUIRE(*adapter.raycast(point2, direction2, max_distance) ==
            stmesh::Vector3F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0)});
    const stmesh::Vector3F direction3{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)};
    const stmesh::Vector3F point3{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)};
    REQUIRE(*adapter.raycast(point3, direction3, max_distance) ==
            stmesh::Vector3F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0)});
  }
}

TEST_CASE("Test edt adapter with sphere edt", "[edt_adapter][surface_adapter]") {
  const auto edt_reader = std::make_shared<stmesh::EDTReader<4>>("data/sphere.mha");
  const stmesh::EDTSurfaceAdapter adapter(edt_reader);
  {
    INFO("Correctly constructed");
    REQUIRE(adapter.boundingBox().min() == stmesh::Vector4F::Constant(stmesh::FLOAT_T(3.0)));
    REQUIRE(adapter.boundingBox().max() == stmesh::Vector4F::Constant(stmesh::FLOAT_T(42.0)));
  }
  {
    INFO("Correctly computes closest point");
    REQUIRE(
        adapter.closestPoint(
            {stmesh::FLOAT_T(12.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5)}) ==
        stmesh::Vector4F{stmesh::FLOAT_T(3.0), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5)});
    REQUIRE(
        adapter.closestPoint(
            {stmesh::FLOAT_T(1.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5)}) ==
        stmesh::Vector4F{stmesh::FLOAT_T(3.0), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5)});
    REQUIRE(
        adapter.closestPoint(
            {stmesh::FLOAT_T(2.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5)}) ==
        stmesh::Vector4F{stmesh::FLOAT_T(3.0), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5)});
    REQUIRE(
        (adapter.closestPoint(stmesh::Vector4F::Constant(18.0)) - stmesh::Vector4F::Constant(stmesh::FLOAT_T(12.625)))
            .cwiseAbs()
            .maxCoeff() < std::sqrt(kEps));
  }
  {
    INFO("Correctly determines if inside");
    REQUIRE(
        adapter.inside({stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5)}));
    REQUIRE(
        adapter.inside({stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(3.01)}));
    REQUIRE_FALSE(
        adapter.inside({stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(2.99)}));
    REQUIRE_FALSE(
        adapter.inside({stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(1.5)}));
    REQUIRE_FALSE(
        adapter.inside({stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(43.0), stmesh::FLOAT_T(22.5)}));
  }
  {
    INFO("Correctly checks sphere intersection");
    const auto radius2 = stmesh::FLOAT_T(5.0);
    const stmesh::Vector4F center2{stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5),
                                   stmesh::FLOAT_T(6.5)};
    const stmesh::HyperSphere4 sphere2{radius2, center2};
    REQUIRE(adapter.intersectedBySphere(sphere2));
    const auto radius3 = stmesh::FLOAT_T(5.0);
    const stmesh::Vector4F center3{stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5),
                                   stmesh::FLOAT_T(22.5)};
    const stmesh::HyperSphere4 sphere3{radius3, center3};
    REQUIRE_FALSE(adapter.intersectedBySphere(sphere3));
  }
  {
    INFO("Correctly computes raycast");
    const auto max_distance = stmesh::FLOAT_T(100.0);
    const stmesh::Vector4F direction1{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0),
                                      stmesh::FLOAT_T(1.0)};
    const stmesh::Vector4F point1{stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5),
                                  stmesh::FLOAT_T(22.5)};
    REQUIRE(
        *adapter.raycast(point1, direction1, max_distance) ==
        stmesh::Vector4F{stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(42.0)});
    const stmesh::Vector4F direction2{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0),
                                      stmesh::FLOAT_T(-1.0)};
    const stmesh::Vector4F point2{stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5),
                                  stmesh::FLOAT_T(22.5)};
    REQUIRE(
        *adapter.raycast(point2, direction2, max_distance) ==
        stmesh::Vector4F{stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(3.0)});
    const stmesh::Vector4F point3{stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5),
                                  stmesh::FLOAT_T(22.5)};
    const stmesh::Vector4F direction3{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0),
                                      stmesh::FLOAT_T(0.0)};
    REQUIRE(
        *adapter.raycast(point3, direction3, max_distance) ==
        stmesh::Vector4F{stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(42.0), stmesh::FLOAT_T(22.5)});
  }
}