#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <memory>

#include "stmesh/edt.hpp"
#include "stmesh/lfs_schemes.hpp"
#include "stmesh/utility.hpp"

TEST_CASE("Test local feature size constant scheme", "[lfs_schemes][constant_lfs]") {
  stmesh::lfs_schemes::Constant constant(1.0);
  REQUIRE(constant(stmesh::Vector4F::Zero()) == 1.0);
  REQUIRE(constant.max() == 1.0);
}

TEST_CASE("Test local feature size binary image approximation scheme", "[lfs_schemes][binary_image_approximation]") {
  const auto reader = std::make_shared<stmesh::EDTReader<4, true>>("data/sphere.mha");
  const stmesh::lfs_schemes::BinaryImageApproximation binary_image_approximation(stmesh::FLOAT_T(2.0), reader);
  const stmesh::FLOAT_T radius = reader->boundingBox().sizes().mean() / 2.0 - 1.0;
  const stmesh::Vector4F surface_point = reader->closestAt(reader->boundingBox().center());
  REQUIRE(std::abs(binary_image_approximation(surface_point) - stmesh::FLOAT_T(2.0) * radius) <= stmesh::FLOAT_T(2.0));
}
