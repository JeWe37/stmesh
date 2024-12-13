#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <array>
#include <filesystem>
#include <vector>

#include "stmesh/mixd.hpp"
#include "stmesh/utility.hpp"

TEST_CASE("Read-write test", "[mixd]") {
  std::vector<stmesh::Vector4F> mxyz{
      {0.0, 0.0, 0.0, 0.0}, {1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0}, {0.3, 0.3, 0.3, 1.0}, {0.3, 0.3, 0.3, -1.0},
  };
  std::vector<std::array<int, 5>> mien{
      {1, 2, 3, 4, 5},
      {1, 2, 3, 4, 6},
  };
  std::vector<std::array<int, 5>> mrng{
      {-1, 1, 1, 1, 1},
      {-1, 1, 1, 1, 1},
  };

  auto base_dir = std::filesystem::temp_directory_path() / "stmesh_test";
  std::filesystem::remove_all(base_dir);
  std::filesystem::create_directory(base_dir);

  auto minf_path = base_dir / "values.minf";
  stmesh::mixd::writeMxyz(minf_path, mxyz, 1.0, 0.0);
  stmesh::mixd::writeIntMixd(minf_path, ".mien", mien);
  stmesh::mixd::writeIntMixd(minf_path, ".mrng", mrng);

  stmesh::mixd::writeMinf(minf_path, base_dir / "values.mxyz", base_dir / "values.mien", base_dir / "values.mrng", 2,
                          6);

  REQUIRE_THAT(stmesh::mixd::readMxyz(base_dir / "values.mxyz"), Catch::Matchers::Equals(mxyz));
  REQUIRE_THAT(stmesh::mixd::readIntMixd(base_dir / "values.mien"), Catch::Matchers::Equals(mien));
  REQUIRE_THAT(stmesh::mixd::readIntMixd(base_dir / "values.mrng"), Catch::Matchers::Equals(mrng));

  auto [mxyz_file, mien_file, mrng_file] = stmesh::mixd::readMinf(minf_path);
  REQUIRE(mxyz_file == base_dir / "values.mxyz");
  REQUIRE(mien_file == base_dir / "values.mien");
  REQUIRE(mrng_file == base_dir / "values.mrng");
}