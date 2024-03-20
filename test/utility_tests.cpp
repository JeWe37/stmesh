#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "stmesh/utility.hpp"

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