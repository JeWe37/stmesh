#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cstddef>

#include "stmesh/marching_hypercubes.hpp"
#include "stmesh/utility.hpp"

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
      corner_values[i] = true;
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