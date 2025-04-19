#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>

#include "stmesh/edt.hpp"
#include "stmesh/geometric_simplex.hpp"
#include "stmesh/utility.hpp"
#include "stmesh/writable_triangulation.hpp"

struct DummyWritableCell {
  stmesh::GeometricSimplex<4> simplex;

  [[nodiscard]] stmesh::GeometricSimplex<4> geometricSimplex() const noexcept { return simplex; }
  [[nodiscard]] static bool isSurfaceSide(size_t) noexcept { return false; }
};

static_assert(stmesh::WritableCell<DummyWritableCell>);

TEST_CASE("Test euclidean distance transform with sphere", "[sphere_edt][edt]") {
  const stmesh::EDTReader<4, true> reader("data/sphere.mha");
  REQUIRE(reader.boundingBox().min() == stmesh::Vector4F::Constant(stmesh::FLOAT_T(2.5)));
  REQUIRE(reader.boundingBox().max() == stmesh::Vector4F::Constant(stmesh::FLOAT_T(42.5)));
  REQUIRE(reader.signedDistance(
              {stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0)}) -
              stmesh::FLOAT_T(-20.0) <
          stmesh::FLOAT_T(std::sqrt(2.0)));
  REQUIRE(reader.signedDistance(
              {stmesh::FLOAT_T(12.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0)}) -
              stmesh::FLOAT_T(-10.0) <
          stmesh::FLOAT_T(std::sqrt(2.0)));
  REQUIRE(
      reader.signedDistance({stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0)}) -
          stmesh::FLOAT_T(20.0) <
      stmesh::FLOAT_T(std::sqrt(2.0)));
  REQUIRE(reader.closestAt({stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0), stmesh::FLOAT_T(2.0)}) ==
          stmesh::Vector4F{stmesh::FLOAT_T(12.5), stmesh::FLOAT_T(12.5), stmesh::FLOAT_T(12.5), stmesh::FLOAT_T(13.5)});
  REQUIRE(
      reader.closestAt({stmesh::FLOAT_T(12.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0)}) ==
      stmesh::Vector4F{stmesh::FLOAT_T(2.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5), stmesh::FLOAT_T(22.5)});
  stmesh::Vector4F region1vec{stmesh::FLOAT_T(12.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0),
                              stmesh::FLOAT_T(22.0)};
  stmesh::Vector4F region2vec{stmesh::FLOAT_T(32.0), stmesh::FLOAT_T(22.0), stmesh::FLOAT_T(22.0),
                              stmesh::FLOAT_T(22.0)};
  REQUIRE(reader.findBoundaryRegion(DummyWritableCell({region1vec, region1vec, region1vec, region1vec, region1vec}),
                                    0) == size_t{1});
  REQUIRE(reader.findBoundaryRegion(DummyWritableCell({region2vec, region2vec, region2vec, region2vec, region2vec}),
                                    0) == size_t{2});
  REQUIRE(reader.findBoundaryRegion(DummyWritableCell({region1vec, region2vec, region1vec, region1vec, region1vec}),
                                    0) == size_t{2});
  REQUIRE(reader.findBoundaryRegion(DummyWritableCell({region2vec, region1vec, region1vec, region1vec, region1vec}),
                                    0) == size_t{1});

  const stmesh::FLOAT_T radius = reader.boundingBox().sizes().mean() / 2.0 - 1.0;
  for (int i = 0; i < 100; ++i) {
    const stmesh::Vector4F random_in_bounding_box =
        stmesh::Vector4F::Random().cwiseProduct(reader.boundingBox().sizes()) + reader.boundingBox().min();
    const stmesh::Vector4F random_on_sphere = reader.closestAt(random_in_bounding_box);
    REQUIRE(std::abs(reader.distanceToThinnedAt(random_on_sphere) - radius) <= stmesh::FLOAT_T(2.0));
  }
}
