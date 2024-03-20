#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>

#include "stmesh/edt.hpp"
#include "stmesh/utility.hpp"

TEST_CASE("Test euclidean distance transform with sphere", "[sphere_edt][edt]") {
  const stmesh::EDTReader<4> reader("data/sphere.mha");
  REQUIRE(reader.boundingBox().min() == stmesh::Vector4F::Constant(stmesh::FLOAT_T(2.5)));
  REQUIRE(reader.boundingBox().max() == stmesh::Vector4F::Constant(stmesh::FLOAT_T(42.5)));
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