#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <memory>

#include "stmesh/edt.hpp"
#include "stmesh/radius_schemes.hpp"
#include "stmesh/utility.hpp"

using namespace Catch;

TEST_CASE("Constant Radius Scheme", "[radius_schemes][constant]") {
  stmesh::radius_schemes::Constant constant_scheme(5.0);

  SECTION("Returns constant value regardless of input") {
    stmesh::Vector4F vec1{1.0, 2.0, 3.0, 4.0};
    stmesh::Vector4F vec2{-1.0, -2.0, -3.0, -4.0};

    REQUIRE(constant_scheme(vec1) == 5.0);
    REQUIRE(constant_scheme(vec2) == 5.0);
  }
}

TEST_CASE("BoundaryDistanceRadius Scheme", "[radius_schemes][boundary]") {
  const auto reader = std::make_shared<stmesh::EDTReader<4, false>>("data/sphere.mha");
  stmesh::Vector4F vec{0.5, 0.5, 0.5, 0.0};

  stmesh::radius_schemes::BoundaryDistanceRadius boundary_scheme(reader);

  stmesh::FLOAT_T result = boundary_scheme(vec);
  stmesh::FLOAT_T expected = reader->signedDistanceAt(vec);
  REQUIRE(result == Approx(expected));

  auto double_mapper = [](stmesh::FLOAT_T dist) { return dist * 2.0; };
  stmesh::radius_schemes::BoundaryDistanceRadius custom_boundary_scheme(reader, double_mapper);

  stmesh::FLOAT_T result_custom = custom_boundary_scheme(vec);
  stmesh::FLOAT_T expected_custom = reader->signedDistanceAt(vec) * 2.0;
  REQUIRE(result_custom == Approx(expected_custom));
}

TEST_CASE("LFSRadius Scheme", "[radius_schemes][lfs]") {
  const auto reader = std::make_shared<stmesh::EDTReader<4, true>>("data/sphere.mha");
  stmesh::Vector4F vec{0.5, 0.5, 0.5, 0.0};

  stmesh::radius_schemes::LFSRadius lfs_scheme(reader);

  stmesh::FLOAT_T result = lfs_scheme(vec);
  stmesh::FLOAT_T expected = reader->distanceToThinnedAt(vec);
  REQUIRE(result == Approx(expected));

  auto half_mapper = [](stmesh::FLOAT_T dist) { return dist * 0.5; };
  stmesh::radius_schemes::LFSRadius custom_lfs_scheme(reader, half_mapper);

  stmesh::FLOAT_T result_custom = custom_lfs_scheme(vec);
  stmesh::FLOAT_T expected_custom = reader->distanceToThinnedAt(vec) * 0.5;
  REQUIRE(result_custom == Approx(expected_custom));
}

TEST_CASE("ImageRadius Scheme", "[radius_schemes][image]") {
  const std::filesystem::path image_path = "data/gradient.mha";

  SECTION("Reads image correctly") {
    stmesh::radius_schemes::ImageRadius image_radius(image_path);

    stmesh::Vector4F vec{10.0, 10.0, 10.0, 5.0};
    stmesh::FLOAT_T result = image_radius(vec);

    REQUIRE(result == Approx(10.0));
  }

  SECTION("Clamps values to image size") {
    stmesh::radius_schemes::ImageRadius image_radius(image_path);

    stmesh::Vector4F vec{-1.0, -1.0, -1.0, -1.0};
    stmesh::FLOAT_T result = image_radius(vec);

    REQUIRE(result == Approx(0.0));

    stmesh::Vector4F vec2{20.0, 20.0, 20.0, 20.0};
    result = image_radius(vec2);

    REQUIRE(result == Approx(19.0));
  }
}
