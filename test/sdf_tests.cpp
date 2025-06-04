#include <catch2/catch_approx.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <unordered_set>
#include <vector>

#include "stmesh/sdf.hpp"
#include "stmesh/utility.hpp"

using namespace Catch;

constexpr stmesh::FLOAT_T kEps = std::numeric_limits<stmesh::FLOAT_T>::epsilon();

TEST_CASE("Test hypersphere functionality in 3D", "[hypersphere][sdf]") {
  const auto radius = stmesh::FLOAT_T(1.0);
  const stmesh::Vector3F center{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)};
  stmesh::HyperSphere3 sphere{radius, center};
  SECTION("Correctly constructed") {
    REQUIRE(sphere.radius() == radius);
    REQUIRE(sphere.center() == center);
    REQUIRE(sphere.boundingBox().min() ==
            stmesh::Vector3F{-stmesh::FLOAT_T(1.0), -stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)});
    REQUIRE(sphere.boundingBox().max() ==
            stmesh::Vector3F{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(2.0)});
  }
  SECTION("Scales correctly") {
    sphere.scale(stmesh::FLOAT_T(2.0));
    REQUIRE(sphere.radius() == stmesh::FLOAT_T(2.0));
    REQUIRE(sphere.center() == center);
    REQUIRE(sphere.boundingBox().min() ==
            stmesh::Vector3F{-stmesh::FLOAT_T(2.0), -stmesh::FLOAT_T(2.0), -stmesh::FLOAT_T(1.0)});
    REQUIRE(sphere.boundingBox().max() ==
            stmesh::Vector3F{stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(3.0)});
    sphere.scale(stmesh::FLOAT_T(0.5));
    REQUIRE(sphere.radius() == stmesh::FLOAT_T(1.0));
    REQUIRE(sphere.center() == center);
    REQUIRE(sphere.boundingBox().min() ==
            stmesh::Vector3F{-stmesh::FLOAT_T(1.0), -stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)});
    REQUIRE(sphere.boundingBox().max() ==
            stmesh::Vector3F{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(2.0)});
  }
  SECTION("Computes distances correctly") {
    REQUIRE(sphere.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}) ==
            Approx(stmesh::FLOAT_T(0.0)));
    REQUIRE(sphere.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}) ==
            Approx(stmesh::FLOAT_T(0.0)));
    REQUIRE(sphere.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}) ==
            Approx(-stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}) ==
            Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(3.0)}) ==
            Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(3.0)}) ==
            Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(1.0)}) ==
            Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(1.0)}) ==
            Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(4.0), stmesh::FLOAT_T(4.0)}) ==
            Approx(stmesh::FLOAT_T(4.0)));
    REQUIRE(sphere.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(4.0), stmesh::FLOAT_T(4.0)}) ==
            Approx(stmesh::FLOAT_T(4.0)));
  }
  SECTION("Computes normals correctly") {
    REQUIRE(sphere.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.5)}) ==
            stmesh::Vector3F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), -stmesh::FLOAT_T(1.0)});
    REQUIRE(sphere.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}) ==
            stmesh::Vector3F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), -stmesh::FLOAT_T(1.0)});
    REQUIRE(sphere.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0)}) ==
            stmesh::Vector3F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)});
    REQUIRE(sphere.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(1.0)}) ==
            stmesh::Vector3F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)});
    REQUIRE(sphere.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(3.0)}) ==
            stmesh::Vector3F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0)}.normalized());
  }
  SECTION("All samples lie within circle and all are different") {
    // NOLINTNEXTLINE(cert-msc51-cpp,cert-msc32-c)
    std::mt19937_64 gen{0};
    std::unordered_set<stmesh::Vector3F, stmesh::Vector3FHash> samples;
    for (size_t i = 0; i < 1000; ++i) {
      const auto sample = sphere.sample(gen);
      REQUIRE((sample - center).norm() <= radius);
      REQUIRE_FALSE(samples.contains(sample));
      samples.insert(sample);
    }
  }
}

TEST_CASE("Test hypersphere functionality in 4D", "[hypersphere][sdf]") {
  const auto radius = stmesh::FLOAT_T(1.0);
  const stmesh::Vector4F center{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)};
  stmesh::HyperSphere4 sphere{radius, center};
  SECTION("Correctly constructed") {
    REQUIRE(sphere.radius() == radius);
    REQUIRE(sphere.center() == center);
    REQUIRE(sphere.boundingBox().min() == stmesh::Vector4F{-stmesh::FLOAT_T(1.0), -stmesh::FLOAT_T(1.0),
                                                           -stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)});
    REQUIRE(sphere.boundingBox().max() ==
            stmesh::Vector4F{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(2.0)});
  }
  SECTION("Scales correctly") {
    sphere.scale(stmesh::FLOAT_T(2.0));
    REQUIRE(sphere.radius() == stmesh::FLOAT_T(2.0));
    REQUIRE(sphere.center() == center);
    REQUIRE(sphere.boundingBox().min() == stmesh::Vector4F{-stmesh::FLOAT_T(2.0), -stmesh::FLOAT_T(2.0),
                                                           -stmesh::FLOAT_T(2.0), -stmesh::FLOAT_T(1.0)});
    REQUIRE(sphere.boundingBox().max() ==
            stmesh::Vector4F{stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(3.0)});
    sphere.scale(stmesh::FLOAT_T(0.5));
    REQUIRE(sphere.radius() == stmesh::FLOAT_T(1.0));
    REQUIRE(sphere.center() == center);
    REQUIRE(sphere.boundingBox().min() == stmesh::Vector4F{-stmesh::FLOAT_T(1.0), -stmesh::FLOAT_T(1.0),
                                                           -stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)});
    REQUIRE(sphere.boundingBox().max() ==
            stmesh::Vector4F{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(2.0)});
  }
  SECTION("Computes distances correctly") {
    REQUIRE(sphere.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0),
                                   stmesh::FLOAT_T(0.0)}) == Approx(stmesh::FLOAT_T(0.0)));
    REQUIRE(sphere.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}) ==
            Approx(stmesh::FLOAT_T(0.0)));
    REQUIRE(sphere.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0),
                                   stmesh::FLOAT_T(1.0)}) == Approx(-stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}) ==
            Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0),
                                   stmesh::FLOAT_T(3.0)}) == Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(3.0)}) ==
            Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0),
                                   stmesh::FLOAT_T(1.0)}) == Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(1.0)}) ==
            Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(sphere.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(4.0), stmesh::FLOAT_T(0.0),
                                   stmesh::FLOAT_T(4.0)}) == Approx(stmesh::FLOAT_T(4.0)));
    REQUIRE(sphere.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(4.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(4.0)}) ==
            Approx(stmesh::FLOAT_T(4.0)));
  }
  SECTION("Computes normals correctly") {
    REQUIRE(sphere.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.5)}) ==
            stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), -stmesh::FLOAT_T(1.0)});
    REQUIRE(sphere.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}) ==
            stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), -stmesh::FLOAT_T(1.0)});
    REQUIRE(sphere.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0)}) ==
            stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)});
    REQUIRE(sphere.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(1.0)}) ==
            stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)});
    REQUIRE(sphere.normal({stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(3.0)}) ==
            stmesh::Vector4F{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0)}
                .normalized());
  }
  SECTION("All samples lie within circle and all are different") {
    // NOLINTNEXTLINE(cert-msc51-cpp,cert-msc32-c)
    std::mt19937_64 gen{0};
    std::unordered_set<stmesh::Vector4F, stmesh::Vector4FHash> samples;
    for (size_t i = 0; i < 1000; ++i) {
      const auto sample = sphere.sample(gen);
      REQUIRE((sample - center).norm() <= radius);
      REQUIRE_FALSE(samples.contains(sample));
      samples.insert(sample);
    }
  }
}

TEST_CASE("Test hypercube functionality in 4D", "[hypercube][sdf]") {
  const stmesh::Vector4F min{-stmesh::FLOAT_T(1.0), -stmesh::FLOAT_T(1.0), -stmesh::FLOAT_T(1.0),
                             -stmesh::FLOAT_T(1.0)};
  const stmesh::Vector4F max{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0)};
  const stmesh::HyperCube4 cube{min, max};
  SECTION("Correctly constructed") {
    REQUIRE(cube.boundingBox().min() == min);
    REQUIRE(cube.boundingBox().max() == max);
  }
  SECTION("Computes distances correctly") {
    REQUIRE(cube.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0),
                                 stmesh::FLOAT_T(0.0)}) == Approx(-stmesh::FLOAT_T(1.0)));
    REQUIRE(cube.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}) ==
            Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(cube.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0),
                                 stmesh::FLOAT_T(1.0)}) == Approx(stmesh::FLOAT_T(0.0)));
    REQUIRE(cube.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)}) ==
            Approx(stmesh::FLOAT_T(0.0)));
    REQUIRE(cube.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0),
                                 stmesh::FLOAT_T(3.0)}) == Approx(stmesh::FLOAT_T(2.0)));
    REQUIRE(cube.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(3.0)}) ==
            Approx(stmesh::FLOAT_T(2.0)));
    REQUIRE(cube.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0),
                                 stmesh::FLOAT_T(1.0)}) == Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(cube.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(1.0)}) ==
            Approx(stmesh::FLOAT_T(1.0)));
    REQUIRE(cube.signedDistance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(5.0), stmesh::FLOAT_T(0.0),
                                 stmesh::FLOAT_T(4.0)}) == Approx(stmesh::FLOAT_T(5.0)));
    REQUIRE(cube.distance({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(5.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(4.0)}) ==
            Approx(stmesh::FLOAT_T(5.0)));
  }
  SECTION("Computes normals correctly") {
    REQUIRE(cube.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.5)}) ==
            stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)});
    REQUIRE(cube.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0)}) ==
            stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0)});
    REQUIRE(cube.normal({-stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)}) ==
            stmesh::Vector4F{-stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)});
    REQUIRE(cube.normal({stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(1.0)}) ==
            stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)});
    REQUIRE((cube.normal({stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(0.0)}) -
             stmesh::Vector4F{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0)}
                 .normalized())
                .maxCoeff() < kEps);
  }
}

TEST_CASE("Test CylinderSDF functionality", "[cylinder][sdf]") {
  const stmesh::FLOAT_T radius = 1.0;
  const Eigen::ParametrizedLine<stmesh::FLOAT_T, 3> axis(stmesh::Vector3F{0.0, 0.0, 0.0},
                                                         stmesh::Vector3F{0.0, 0.0, 1.0});
  stmesh::CylinderSDF cylinder(radius, axis);

  SECTION("Correctly constructed") {
    REQUIRE(cylinder.boundingBox().min() == stmesh::Vector3F{-1.0, -1.0, 0.0});
    REQUIRE(cylinder.boundingBox().max() == stmesh::Vector3F{1.0, 1.0, 1.0});
  }

  SECTION("Computes signed distances correctly") {
    // Inside the cylinder
    REQUIRE(cylinder.signedDistance({0.0, 0.0, 0.5}) == Approx(-0.5));
    // On the side
    REQUIRE(cylinder.signedDistance({1.0, 0.0, 0.5}) == Approx(0.0));
    // Outside the side
    REQUIRE(cylinder.signedDistance({2.0, 0.0, 0.5}) == Approx(1.0));
    // On the bottom cap
    REQUIRE(cylinder.signedDistance({0.0, 0.0, 0.0}) == Approx(0.0));
    // On the top cap
    REQUIRE(cylinder.signedDistance({0.0, 0.0, 1.0}) == Approx(0.0));
    // Below the bottom cap
    REQUIRE(cylinder.signedDistance({0.0, 0.0, -1.0}) == Approx(1.0));
    // Above the top cap
    REQUIRE(cylinder.signedDistance({0.0, 0.0, 2.0}) == Approx(1.0));
    // On the edge (bottom)
    REQUIRE(cylinder.signedDistance({1.0, 0.0, 0.0}) == Approx(0.0));
    // On the edge (top)
    REQUIRE(cylinder.signedDistance({1.0, 0.0, 1.0}) == Approx(0.0));
    // Outside both radially and axially
    REQUIRE(cylinder.signedDistance({2.0, 0.0, 2.0}) == Approx(std::sqrt(2.0)));
  }

  SECTION("Computes distances correctly") {
    // Inside the cylinder (negative signed distance, so distance is positive)
    REQUIRE(cylinder.distance({0.0, 0.0, 0.5}) == Approx(0.5));
    // On the side
    REQUIRE(cylinder.distance({1.0, 0.0, 0.5}) == Approx(0.0));
    // Outside the side
    REQUIRE(cylinder.distance({2.0, 0.0, 0.5}) == Approx(1.0));
    // On the bottom cap
    REQUIRE(cylinder.distance({0.0, 0.0, 0.0}) == Approx(0.0));
    // On the top cap
    REQUIRE(cylinder.distance({0.0, 0.0, 1.0}) == Approx(0.0));
    // Below the bottom cap
    REQUIRE(cylinder.distance({0.0, 0.0, -1.0}) == Approx(1.0));
    // Above the top cap
    REQUIRE(cylinder.distance({0.0, 0.0, 2.0}) == Approx(1.0));
    // On the edge (bottom)
    REQUIRE(cylinder.distance({1.0, 0.0, 0.0}) == Approx(0.0));
    // On the edge (top)
    REQUIRE(cylinder.distance({1.0, 0.0, 1.0}) == Approx(0.0));
    // Outside both radially and axially
    REQUIRE(cylinder.distance({2.0, 0.0, 2.0}) == Approx(std::sqrt(2.0)));
  }

  SECTION("Computes normals correctly") {
    // On the side
    REQUIRE(cylinder.normal({1.0, 0.0, 0.5}) == stmesh::Vector3F{1.0, 0.0, 0.0});
    // On the bottom cap
    REQUIRE(cylinder.normal({0.0, 0.0, 0.0}) == stmesh::Vector3F{0.0, 0.0, -1.0});
    // On the top cap
    REQUIRE(cylinder.normal({0.0, 0.0, 1.0}) == stmesh::Vector3F{0.0, 0.0, 1.0});
    // On the edge (bottom, prefers side normal)
    REQUIRE(cylinder.normal({1.0, 0.0, 0.0}) == stmesh::Vector3F{1.0, 0.0, 0.0});
    // On the edge (top, prefers side normal)
    REQUIRE(cylinder.normal({1.0, 0.0, 1.0}) == stmesh::Vector3F{1.0, 0.0, 0.0});
  }
}

TEST_CASE("Test ExtrudedSDF with HyperSphere<2> against CylinderSDF", "[extruded][sdf]") {
  // Define the base 2D circle
  const stmesh::FLOAT_T radius = 1.0;
  stmesh::HyperSphere<2> circle(1.0, {0.0, 0.0});

  // Define the extrusion axis (z-axis, length 1)
  const Eigen::ParametrizedLine<stmesh::FLOAT_T, 3> axis(stmesh::Vector3F{0.0, 0.0, 0.0}, // Origin
                                                         stmesh::Vector3F{0.0, 0.0, 2.0}  // Direction
  );

  // Create the extruded SDF and cylinder SDF
  stmesh::ExtrudedSDF<stmesh::HyperSphere<2>> extruded_circle(circle, axis);
  stmesh::CylinderSDF cylinder(radius, axis);

  // Test signed distances
  SECTION("Signed distances match") {
    for (const auto &point : std::vector<stmesh::Vector3F>{
             {0.0, 0.0, 1.0},  // Inside the cylinder
             {1.0, 0.0, 1.0},  // On the curved surface
             {2.0, 0.0, 1.0},  // Outside radially
             {0.0, 0.0, 0.0},  // On the bottom cap
             {0.0, 0.0, 2.0},  // On the top cap
             {0.0, 0.0, -1.0}, // Below the cylinder
             {0.0, 0.0, 3.0},  // Above the cylinder
             {1.0, 0.0, 0.0},  // On bottom edge
             {1.0, 0.0, 2.0},  // On top edge
             {2.0, 0.0, 3.0}   // Outside both radially and axially
         }) {
      INFO("Point: " << point.transpose());
      stmesh::FLOAT_T extruded_distance = extruded_circle.signedDistance(point);
      stmesh::FLOAT_T cylinder_distance = cylinder.signedDistance(point);
      REQUIRE(extruded_distance == Approx(cylinder_distance).epsilon(1e-6));
    }
  }

  // Test normal vectors
  SECTION("Normals match") {
    for (const auto &point : std::vector<stmesh::Vector3F>{
             {1.0, 0.0, 1.0},  // On the curved surface
             {0.0, 0.0, 0.0},  // On the bottom cap
             {0.0, 0.0, 2.0},  // On the top cap
             {1.01, 0.0, 0.0}, // On bottom edge
             {1.01, 0.0, 2.0}  // On top edge
         }) {
      INFO("Point: " << point.transpose());
      stmesh::Vector3F extruded_normal = extruded_circle.normal(point);
      stmesh::Vector3F cylinder_normal = cylinder.normal(point);
      REQUIRE((extruded_normal - cylinder_normal).norm() < 1e-6);
    }
  }

  SECTION("Random normals and distances match") {
    // NOLINTNEXTLINE(cert-msc51-cpp,cert-msc32-c)
    std::mt19937_64 gen{0};
    std::unordered_set<stmesh::Vector3F, stmesh::Vector3FHash> samples;
    for (size_t i = 0; i < 1000; ++i) {
      const auto sample = cylinder.sample(gen);
      INFO("Sample: " << sample.transpose());
      stmesh::FLOAT_T extruded_distance = extruded_circle.signedDistance(sample);
      stmesh::FLOAT_T cylinder_distance = cylinder.signedDistance(sample);
      INFO("Extruded distance: " << extruded_distance);
      INFO("Cylinder distance: " << cylinder_distance);
      REQUIRE(extruded_distance == Approx(cylinder_distance).epsilon(1e-6));
      stmesh::Vector3F extruded_normal = extruded_circle.normal(sample);
      stmesh::Vector3F cylinder_normal = cylinder.normal(sample);
      INFO("Extruded normal: " << extruded_normal.transpose());
      INFO("Cylinder normal: " << cylinder_normal.transpose());
      REQUIRE((extruded_normal - cylinder_normal).norm() < 1e-6);
    }
  }

  // Test bounding boxes
  SECTION("Bounding boxes match") {
    REQUIRE(extruded_circle.boundingBox().min() == stmesh::Vector3F{-1.0, -1.0, 0.0});
    REQUIRE(extruded_circle.boundingBox().max() == stmesh::Vector3F{1.0, 1.0, 2.0});
  }
}