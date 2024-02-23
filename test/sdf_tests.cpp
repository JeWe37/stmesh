#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <limits>
#include <random>
#include <unordered_set>

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