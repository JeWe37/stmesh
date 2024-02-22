// NOLINTBEGIN(*-magic-numbers,misc-include-cleaner,readability-function-cognitive-complexity)
#include "stmesh/edt.hpp"
#include "stmesh/geometric_simplex.hpp"
#include "stmesh/marching_hypercubes.hpp"
#include "stmesh/meshing_algorithm.hpp"
#include "stmesh/meshing_cell.hpp"
#include "stmesh/sdf.hpp"
#include "stmesh/surface_adapters.hpp"
#include "stmesh/triangulation.hpp"
#include "stmesh/utility.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <Eigen/src/Geometry/Hyperplane.h>
#include <Eigen/src/Geometry/ParametrizedLine.h>
#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <array>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <memory>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <random>
#include <stmesh/stmesh.hpp>

using namespace Catch;

constexpr stmesh::FLOAT_T kEps = std::numeric_limits<stmesh::FLOAT_T>::epsilon();

TEST_CASE("Factorial is computed", "[factorial][utility]") {
  REQUIRE(stmesh::factorial(0U) == 1U);
  REQUIRE(stmesh::factorial(1U) == 1U);
  REQUIRE(stmesh::factorial(2U) == 2U);
  REQUIRE(stmesh::factorial(3U) == 6U);
  REQUIRE(stmesh::factorial(10U) == 3628800U);
}

TEST_CASE("Sqrt is computed", "[sqrt][utility]") {
  REQUIRE(std::abs(stmesh::sqrt(kEps * kEps) - kEps) < kEps);
  REQUIRE(std::abs(stmesh::sqrt(kEps) * stmesh::sqrt(kEps) - kEps) < kEps);
  REQUIRE(stmesh::sqrt(stmesh::FLOAT_T(0.0)) == stmesh::FLOAT_T(0.0));
  REQUIRE(stmesh::sqrt(stmesh::FLOAT_T(1.0)) == stmesh::FLOAT_T(1.0));
  REQUIRE(stmesh::sqrt(stmesh::FLOAT_T(4.0)) == stmesh::FLOAT_T(2.0));
  REQUIRE(stmesh::sqrt(stmesh::FLOAT_T(16.0)) == stmesh::FLOAT_T(4.0));
  REQUIRE(std::abs(stmesh::sqrt(stmesh::FLOAT_T(2.0)) * stmesh::sqrt(stmesh::FLOAT_T(2.0)) - stmesh::FLOAT_T(2.0)) <
          stmesh::FLOAT_T(3.0) * kEps);
  REQUIRE(std::abs(stmesh::sqrt(stmesh::FLOAT_T(3.0)) * stmesh::sqrt(stmesh::FLOAT_T(3.0)) - stmesh::FLOAT_T(3.0)) <
          stmesh::FLOAT_T(3.0) * kEps);
  REQUIRE(std::abs(stmesh::sqrt(stmesh::FLOAT_T(5.0)) * stmesh::sqrt(stmesh::FLOAT_T(5.0)) - stmesh::FLOAT_T(5.0)) <
          stmesh::FLOAT_T(5.0) * kEps);
}

TEST_CASE("nChoosek is computed", "[nChoosek][utility]") {
  REQUIRE(stmesh::nChoosek(0U, 0U) == 1U);
  REQUIRE(stmesh::nChoosek(1U, 0U) == 1U);
  REQUIRE(stmesh::nChoosek(1U, 1U) == 1U);
  REQUIRE(stmesh::nChoosek(2U, 0U) == 1U);
  REQUIRE(stmesh::nChoosek(2U, 1U) == 2U);
  REQUIRE(stmesh::nChoosek(2U, 2U) == 1U);
  REQUIRE(stmesh::nChoosek(4U, 0U) == 1U);
  REQUIRE(stmesh::nChoosek(4U, 1U) == 4U);
  REQUIRE(stmesh::nChoosek(4U, 2U) == 6U);
  REQUIRE(stmesh::nChoosek(4U, 3U) == 4U);
  REQUIRE(stmesh::nChoosek(4U, 4U) == 1U);
  REQUIRE(stmesh::nChoosek(6U, 0U) == 1U);
  REQUIRE(stmesh::nChoosek(6U, 1U) == 6U);
  REQUIRE(stmesh::nChoosek(6U, 2U) == 15U);
  REQUIRE(stmesh::nChoosek(6U, 3U) == 20U);
  REQUIRE(stmesh::nChoosek(6U, 4U) == 15U);
  REQUIRE(stmesh::nChoosek(6U, 5U) == 6U);
  REQUIRE(stmesh::nChoosek(6U, 6U) == 1U);
  REQUIRE(stmesh::nChoosek(10U, 0U) == 1U);
  REQUIRE(stmesh::nChoosek(10U, 1U) == 10U);
  REQUIRE(stmesh::nChoosek(10U, 2U) == 45U);
  REQUIRE(stmesh::nChoosek(10U, 5U) == 252U);
  REQUIRE(stmesh::nChoosek(10U, 8U) == 45U);
  REQUIRE(stmesh::nChoosek(10U, 9U) == 10U);
  REQUIRE(stmesh::nChoosek(10U, 10U) == 1U);
}

TEST_CASE("All corners are found", "[allCorners][utility]") {
  stmesh::Vector4F min{-stmesh::FLOAT_T(1.0), -stmesh::FLOAT_T(1.0), -stmesh::FLOAT_T(1.0), -stmesh::FLOAT_T(1.0)};
  stmesh::Vector4F max{stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(1.0)};
  const Eigen::AlignedBox<stmesh::FLOAT_T, 4> box{min, max};
  const auto corners = stmesh::allCorners<4U>(box);
  STATIC_REQUIRE(corners.size() == 16U);
  REQUIRE_THAT(corners, Matchers::UnorderedRangeEquals(std::vector<stmesh::Vector4F>{min,
                                                                                     {min[0], min[1], min[2], max[3]},
                                                                                     {min[0], min[1], max[2], min[3]},
                                                                                     {min[0], min[1], max[2], max[3]},
                                                                                     {min[0], max[1], min[2], min[3]},
                                                                                     {min[0], max[1], min[2], max[3]},
                                                                                     {min[0], max[1], max[2], min[3]},
                                                                                     {min[0], max[1], max[2], max[3]},
                                                                                     {max[0], min[1], min[2], min[3]},
                                                                                     {max[0], min[1], min[2], max[3]},
                                                                                     {max[0], min[1], max[2], min[3]},
                                                                                     {max[0], min[1], max[2], max[3]},
                                                                                     {max[0], max[1], min[2], min[3]},
                                                                                     {max[0], max[1], min[2], max[3]},
                                                                                     {max[0], max[1], max[2], min[3]},
                                                                                     max}));
}

TEST_CASE("Non-square full rank matrix kernel is computed correctly", "[kernel][utility]") {
  Eigen::Matrix<stmesh::FLOAT_T, 4, 3> matrix;
  matrix << stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0),
      stmesh::FLOAT_T(1.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(1.0),
      stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0);
  const Eigen::Matrix<stmesh::FLOAT_T, 4, 1> result = stmesh::kernel<4U, 3U>(matrix);
  STATIC_REQUIRE(result.rows() == 4U);
  STATIC_REQUIRE(result.cols() == 1U);
  REQUIRE(result == Eigen::Matrix<stmesh::FLOAT_T, 4, 1>::Unit(3));
}

TEST_CASE("Different hashes for different vectors 3D", "[hash][utility]") {
  const std::vector<stmesh::Vector3F> vectors = {
      {1.0, 1.0, 1.0}, {4.0, 1.0, 2.0}, {4.0, 2.0, 1.0}, {1.0, 3.0, 2.0}, {2.0, 1.0, 1.0}};
  for (size_t i = 0; i < vectors.size(); ++i) {
    for (size_t j = i + 1; j < vectors.size(); ++j)
      REQUIRE(stmesh::Vector3FHash{}(vectors[i]) != stmesh::Vector3FHash{}(vectors[j]));
  }
}

TEST_CASE("Different hashes for different vectors 4D", "[hash][utility]") {
  const std::vector<stmesh::Vector4F> vectors = {
      {1.0, 1.0, 1.0, 1.0}, {4.0, 1.0, 2.0, 3.0}, {4.0, 2.0, 1.0, 3.0}, {1.0, 3.0, 2.0, 4.0}, {2.0, 1.0, 1.0, 1.0}};
  for (size_t i = 0; i < vectors.size(); ++i) {
    for (size_t j = i + 1; j < vectors.size(); ++j)
      REQUIRE(stmesh::Vector4FHash{}(vectors[i]) != stmesh::Vector4FHash{}(vectors[j]));
  }
}

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

// NOLINTBEGIN(bugprone-unchecked-optional-access)
TEST_CASE("Test marching hypercubes", "[marching_hypercubes]") {
  {
    const stmesh::Vector4F start{stmesh::FLOAT_T(0.1), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.5),
                                 stmesh::FLOAT_T(0.5)};
    const stmesh::Vector4F end{stmesh::FLOAT_T(0.8), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.5)};
    bool inside = true;
    std::array<bool, 16> corner_values{};
    std::fill_n(corner_values.begin() + 8, 8, true);
    REQUIRE(stmesh::surfaceRayIntersection(start, end, corner_values, &inside).value() ==
            stmesh::Vector4F::Constant(stmesh::FLOAT_T(0.5)));
    REQUIRE_FALSE(inside);
  }
  {
    const stmesh::Vector4F start{stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.5),
                                 stmesh::FLOAT_T(0.1)};
    const stmesh::Vector4F end{stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(0.8)};
    bool inside = true;
    std::array<bool, 16> corner_values{};
    for (size_t i = 0; i < corner_values.size(); i += 2)
      corner_values.at(i) = true;
    REQUIRE(stmesh::surfaceRayIntersection(start, end, corner_values, &inside).value() ==
            stmesh::Vector4F::Constant(stmesh::FLOAT_T(0.5)));
    REQUIRE(inside);
  }
  {
    const auto start1 = stmesh::Vector4F::Constant(stmesh::FLOAT_T(0.4));
    const auto start2 = stmesh::Vector4F::Constant(stmesh::FLOAT_T(0.6));
    const auto end1 = stmesh::Vector4F::Constant(stmesh::FLOAT_T(0.0));
    const auto end2 = stmesh::Vector4F::Constant(stmesh::FLOAT_T(1.0));
    bool inside = true;
    std::array<bool, 16> corner_values{};
    std::fill_n(corner_values.begin() + 1, 14, true);
    REQUIRE(stmesh::surfaceRayIntersection(start1, end1, corner_values, &inside).value() ==
            stmesh::Vector4F::Constant(stmesh::FLOAT_T(0.125)));
    REQUIRE(inside);
    REQUIRE(stmesh::surfaceRayIntersection(end1, start1, corner_values, &inside).value() ==
            stmesh::Vector4F::Constant(stmesh::FLOAT_T(0.125)));
    REQUIRE_FALSE(inside);
    REQUIRE(stmesh::surfaceRayIntersection(start2, end2, corner_values, &inside).value() ==
            stmesh::Vector4F::Constant(stmesh::FLOAT_T(0.875)));
    REQUIRE(inside);
    REQUIRE(stmesh::surfaceRayIntersection(end2, start2, corner_values, &inside).value() ==
            stmesh::Vector4F::Constant(stmesh::FLOAT_T(0.875)));
    REQUIRE_FALSE(inside);
  }
  {
    const auto start = stmesh::Vector4F::Zero();
    const auto end = stmesh::Vector4F::Ones();
    bool inside = true;
    std::array<bool, 16> corner_values{};
    std::fill_n(corner_values.begin(), 3, true);
    const auto result = stmesh::Vector4F::Constant(stmesh::FLOAT_T(0.25));
    REQUIRE(stmesh::surfaceRayIntersection(start, end, corner_values, &inside).value() == result);
    REQUIRE(inside);
    REQUIRE(stmesh::surfaceRayIntersection(end, start, corner_values, &inside).value() == result);
    REQUIRE_FALSE(inside);
  }
}
// NOLINTEND(bugprone-unchecked-optional-access)

TEST_CASE("Test euclidean distance transform with sphere", "[sphere_edt][edt]") {
  const stmesh::EDTReader<4> reader("data/sphere.mha");
  REQUIRE(reader.boundingBox().min() == stmesh::Vector4F::Ones() * stmesh::FLOAT_T(2.0));
  REQUIRE(reader.boundingBox().max() == stmesh::Vector4F::Ones() * stmesh::FLOAT_T(42.0));
  REQUIRE(reader.signedDistanceAt(
              {stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0)}) -
              stmesh::FLOAT_T(20.0) <
          stmesh::FLOAT_T(std::sqrt(2.0)));
  REQUIRE(reader.signedDistanceAt(
              {stmesh::FLOAT_T(12.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0)}) -
              stmesh::FLOAT_T(10.0) <
          stmesh::FLOAT_T(std::sqrt(2.0)));
  REQUIRE(reader.signedDistanceAt(
              {stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0)}) -
              stmesh::FLOAT_T(-20.0) <
          stmesh::FLOAT_T(std::sqrt(2.0)));
  REQUIRE(reader.closestAt({stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0)}) ==
          stmesh::Vector4F{stmesh::FLOAT_T(12.5), stmesh::FLOAT_T(12.5), stmesh::FLOAT_T(12.5), stmesh::FLOAT_T(13.5)});
  REQUIRE(
      reader.closestAt({stmesh::FLOAT_T(12.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0)}) ==
      stmesh::Vector4F{stmesh::FLOAT_T(2.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5)});
  REQUIRE(reader.findBoundaryRegion({stmesh::FLOAT_T(12.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0),
                                     stmesh::FLOAT_T(22.0)}) == size_t{1});
  REQUIRE(reader.findBoundaryRegion({stmesh::FLOAT_T(32.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0),
                                     stmesh::FLOAT_T(22.0)}) == size_t{2});
}

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
    // NOLINTBEGIN(bugprone-unchecked-optional-access)
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
    // NOLINTEND(bugprone-unchecked-optional-access)
  }
}

TEST_CASE("Test edt adapter with sphere edt", "[edt_adapter][surface_adapter]") {
  const auto edt_reader = std::make_shared<stmesh::EDTReader<4>>("data/sphere.mha");
  const stmesh::EDTSurfaceAdapter adapter(edt_reader);
  {
    INFO("Correctly constructed");
    REQUIRE(adapter.boundingBox().min() ==
            stmesh::Vector4F{stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0)});
    REQUIRE(adapter.boundingBox().max() == stmesh::Vector4F{stmesh::FLOAT_T(42.0), stmesh::FLOAT_T(42.0),
                                                            stmesh::FLOAT_T(42.0), stmesh::FLOAT_T(42.0)});
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
    // NOLINTBEGIN(bugprone-unchecked-optional-access)
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
    // NOLINTEND(bugprone-unchecked-optional-access)
  }
}

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
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      for (const auto *it = full_cell.vertices_begin(); it != full_cell.vertices_end(); ++it) {
        const stmesh::Triangulation<>::Point pt = (*it)->point();
        const stmesh::Vector4F vertex{static_cast<stmesh::FLOAT_T>(pt[0]), static_cast<stmesh::FLOAT_T>(pt[1]),
                                      static_cast<stmesh::FLOAT_T>(pt[2]), static_cast<stmesh::FLOAT_T>(pt[3])};
        REQUIRE(vertex.cwiseAbs() == max);
        REQUIRE((*it)->data().nonfree_vertex);
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
      // NOLINTBEGIN(bugprone-unchecked-optional-access)
      std::array<stmesh::Triangulation<>::VertexHandle, 2> mirror_vertices = {std::get<0>(existing_side),
                                                                              std::get<0>(*optional_side)};
      std::array<int, 2> vertex_indices = {std::get<1>(existing_side), std::get<1>(*optional_side)};
      std::array<stmesh::Triangulation<>::FullCellHandle, 2> full_cells = {std::get<2>(existing_side),
                                                                           std::get<2>(*optional_side)};
      // NOLINTEND(bugprone-unchecked-optional-access)
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

void verifyMeshingAlgorithm(auto &meshing_algorithm) {
  const stmesh::detail::Triangulation &triangulation = meshing_algorithm.triangulation();
  std::array<stmesh::detail::Rules, 6> rule_list{stmesh::detail::Rule1{}, stmesh::detail::Rule2{},
                                                 stmesh::detail::Rule3{}, stmesh::detail::Rule4{},
                                                 stmesh::detail::Rule5{}, stmesh::detail::Complete{}};
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
    if (rule_index != stmesh::detail::Complete::kIndex) {
      if (!std::visit([&](auto &r) { return r.check(meshing_algorithm, full_cell_handle, true); },
                      rule_list.at(rule_index)))
        WARN("Issue!");
      REQUIRE(std::visit([&](auto &r) { return r.check(meshing_algorithm, full_cell_handle, true); },
                         rule_list.at(rule_index)));
      REQUIRE(std::visit([&](auto &r) { return r.priority(); }, rule_list.at(rule_index)) == rule_priority);
    }
    for (size_t i = 0; i < rule_index; ++i) {
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
  stmesh::MeshingAlgorithm meshing_algorithm(sdf_surface_adapter, stmesh::FLOAT_T(20.0), stmesh::FLOAT_T(0.001),
                                             stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(5.0), stmesh::FLOAT_T(2.0));
  SECTION("Correctly constructed") {
    stmesh::Vector4F two_five{stmesh::FLOAT_T(2.5), stmesh::FLOAT_T(2.5), stmesh::FLOAT_T(2.5), stmesh::FLOAT_T(2.5)};
    REQUIRE(meshing_algorithm.triangulation().boundingBox().min() == -two_five);
    REQUIRE(meshing_algorithm.triangulation().boundingBox().max() == two_five);
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
      // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
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
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
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
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
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
      stmesh::FLOAT_T(0.2),
      stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)});
  stmesh::MeshingAlgorithm meshing_algorithm(sdf_surface_adapter, stmesh::FLOAT_T(20.0), stmesh::FLOAT_T(0.01),
                                             stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(5.0), stmesh::FLOAT_T(0.15));
  for (const auto &point : {
           stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
           stmesh::Vector4F{stmesh::FLOAT_T(0.05), stmesh::FLOAT_T(0.1), stmesh::FLOAT_T(0.1), stmesh::FLOAT_T(0.05)},
           stmesh::Vector4F{stmesh::FLOAT_T(0.05), stmesh::FLOAT_T(0.1), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.05)},
           stmesh::Vector4F{stmesh::FLOAT_T(0.2), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
           stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.2), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)},
           stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.2), stmesh::FLOAT_T(0.0)},
           stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.2)},
           stmesh::Vector4F{-stmesh::FLOAT_T(0.1), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.1)},
           stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), -stmesh::FLOAT_T(0.1), stmesh::FLOAT_T(0.0)},
           stmesh::Vector4F{stmesh::FLOAT_T(0.002), stmesh::FLOAT_T(0.002), stmesh::FLOAT_T(0.002),
                            stmesh::FLOAT_T(0.199)},
       })
    meshing_algorithm.insert(point);
  SECTION("Check rules") {
    verifyMeshingAlgorithm(meshing_algorithm);
    const stmesh::detail::Triangulation &triangulation = meshing_algorithm.triangulation();
    for (const auto &cell : meshing_algorithm.queue())
      std::visit(Overload{
                     [&](const stmesh::detail::Rule1 &r) {
                       const stmesh::GeometricSimplex<4> simplex = triangulation.fullCellSimplex(cell.full_cell);
                       const stmesh::HyperSphere4 circumsphere = simplex.circumsphere();
                       REQUIRE(sdf_surface_adapter.intersectedBySphere(circumsphere));
                       const stmesh::Vector4F z = sdf_surface_adapter.closestPoint(circumsphere.center());
                       REQUIRE(z == r.z0);
                       REQUIRE(triangulation.verticesInRadius(z, meshing_algorithm.deltaSurface(z)).empty());
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
                     },
                     [&](const stmesh::detail::Rule3 &r) {
                       const stmesh::GeometricSimplex<4> simplex = triangulation.fullCellSimplex(cell.full_cell);
                       const stmesh::HyperSphere4 circumsphere = simplex.circumsphere();
                       const stmesh::Vector4F &z = circumsphere.center();
                       REQUIRE(sdf_surface_adapter.inside(z));
                       REQUIRE(simplex.radiusEdgeRatio() >= stmesh::FLOAT_T(20.0));
                       REQUIRE(r.z == z);
                     },
                     [&](const stmesh::detail::Rule4 &r) {
                       const stmesh::GeometricSimplex<4> simplex = triangulation.fullCellSimplex(cell.full_cell);
                       const stmesh::HyperSphere4 circumsphere = simplex.circumsphere();
                       REQUIRE(sdf_surface_adapter.inside(circumsphere.center()));
                       REQUIRE(simplex.sliverSimplex(20.0, stmesh::FLOAT_T(0.01)));
                       REQUIRE(r.vertices == simplex.vertices());
                     },
                     [&](const stmesh::detail::Rule5 &r) {
                       const stmesh::GeometricSimplex<4> simplex = triangulation.fullCellSimplex(cell.full_cell);
                       auto sub_simplices = simplex.subSimplices<4>();
                       REQUIRE(std::ranges::find_if(sub_simplices, [&](const auto &sub_simplex) {
                                 return sub_simplex.vertices() == r.vertices;
                               }) != sub_simplices.end());
                       REQUIRE(meshing_algorithm.voronoiDual(r.vertices).has_value());
                     },
                     [&](const stmesh::detail::Complete & /* unused */) {},
                 },
                 cell.rule);
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
          if constexpr (Rule::kIndex < stmesh::detail::Rule4::kIndex)
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
  stmesh::MeshingAlgorithm meshing_algorithm(sdf_surface_adapter, stmesh::FLOAT_T(20.0), stmesh::FLOAT_T(1e-6),
                                             stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(5.0), stmesh::FLOAT_T(0.15));
  meshing_algorithm.triangulate([&] { verifyMeshingAlgorithm(meshing_algorithm); });
}
// NOLINTEND(*-magic-numbers,misc-include-cleaner,readability-function-cognitive-complexity)
