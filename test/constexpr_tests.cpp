// NOLINTBEGIN(*-magic-numbers,misc-include-cleaner)
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>

#include <stmesh/stmesh.hpp>

TEST_CASE("Factorial is computed with constexpr", "[factorial][utility]") {
  STATIC_REQUIRE(stmesh::factorial(0U) == 1U);
  STATIC_REQUIRE(stmesh::factorial(1U) == 1U);
  STATIC_REQUIRE(stmesh::factorial(2U) == 2U);
  STATIC_REQUIRE(stmesh::factorial(3U) == 6U);
  STATIC_REQUIRE(stmesh::factorial(10U) == 3628800U);
}

TEST_CASE("Sqrt is computed with constexpr", "[sqrt][utility]") {
  constexpr stmesh::FLOAT_T kEps = std::numeric_limits<stmesh::FLOAT_T>::epsilon();
  STATIC_REQUIRE(std::abs(stmesh::sqrt(kEps * kEps) - kEps) < kEps);
  STATIC_REQUIRE(std::abs(stmesh::sqrt(kEps) * stmesh::sqrt(kEps) - kEps) < kEps);
  STATIC_REQUIRE(stmesh::sqrt(stmesh::FLOAT_T(0.0)) == stmesh::FLOAT_T(0.0));
  STATIC_REQUIRE(stmesh::sqrt(stmesh::FLOAT_T(1.0)) == stmesh::FLOAT_T(1.0));
  STATIC_REQUIRE(stmesh::sqrt(stmesh::FLOAT_T(4.0)) == stmesh::FLOAT_T(2.0));
  STATIC_REQUIRE(stmesh::sqrt(stmesh::FLOAT_T(16.0)) == stmesh::FLOAT_T(4.0));
  STATIC_REQUIRE(std::abs(stmesh::sqrt(stmesh::FLOAT_T(2.0)) * stmesh::sqrt(stmesh::FLOAT_T(2.0)) -
                          stmesh::FLOAT_T(2.0)) < stmesh::FLOAT_T(3.0) * kEps);
  STATIC_REQUIRE(std::abs(stmesh::sqrt(stmesh::FLOAT_T(3.0)) * stmesh::sqrt(stmesh::FLOAT_T(3.0)) -
                          stmesh::FLOAT_T(3.0)) < stmesh::FLOAT_T(3.0) * kEps);
  STATIC_REQUIRE(std::abs(stmesh::sqrt(stmesh::FLOAT_T(5.0)) * stmesh::sqrt(stmesh::FLOAT_T(5.0)) -
                          stmesh::FLOAT_T(5.0)) < stmesh::FLOAT_T(5.0) * kEps);
}

TEST_CASE("nChoosek is computed with constexpr", "[nChoosek][utility]") {
  STATIC_REQUIRE(stmesh::nChoosek(0U, 0U) == 1U);
  STATIC_REQUIRE(stmesh::nChoosek(1U, 0U) == 1U);
  STATIC_REQUIRE(stmesh::nChoosek(1U, 1U) == 1U);
  STATIC_REQUIRE(stmesh::nChoosek(2U, 0U) == 1U);
  STATIC_REQUIRE(stmesh::nChoosek(2U, 1U) == 2U);
  STATIC_REQUIRE(stmesh::nChoosek(2U, 2U) == 1U);
  STATIC_REQUIRE(stmesh::nChoosek(4U, 0U) == 1U);
  STATIC_REQUIRE(stmesh::nChoosek(4U, 1U) == 4U);
  STATIC_REQUIRE(stmesh::nChoosek(4U, 2U) == 6U);
  STATIC_REQUIRE(stmesh::nChoosek(4U, 3U) == 4U);
  STATIC_REQUIRE(stmesh::nChoosek(4U, 4U) == 1U);
  STATIC_REQUIRE(stmesh::nChoosek(6U, 0U) == 1U);
  STATIC_REQUIRE(stmesh::nChoosek(6U, 1U) == 6U);
  STATIC_REQUIRE(stmesh::nChoosek(6U, 2U) == 15U);
  STATIC_REQUIRE(stmesh::nChoosek(6U, 3U) == 20U);
  STATIC_REQUIRE(stmesh::nChoosek(6U, 4U) == 15U);
  STATIC_REQUIRE(stmesh::nChoosek(6U, 5U) == 6U);
  STATIC_REQUIRE(stmesh::nChoosek(6U, 6U) == 1U);
  STATIC_REQUIRE(stmesh::nChoosek(10U, 0U) == 1U);
  STATIC_REQUIRE(stmesh::nChoosek(10U, 1U) == 10U);
  STATIC_REQUIRE(stmesh::nChoosek(10U, 2U) == 45U);
  STATIC_REQUIRE(stmesh::nChoosek(10U, 5U) == 252U);
  STATIC_REQUIRE(stmesh::nChoosek(10U, 8U) == 45U);
  STATIC_REQUIRE(stmesh::nChoosek(10U, 9U) == 10U);
  STATIC_REQUIRE(stmesh::nChoosek(10U, 10U) == 1U);
}
// NOLINTEND(*-magic-numbers,misc-include-cleaner)
