#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstddef>
#include <set>
#include <vector>

#include "stmesh/voxel_complex.hpp"

TEST_CASE("VoxelComplex Face operations", "[voxelcomplex][face]") {
  SECTION("Face setPos and getPos") {
    stmesh::VoxelComplex<3>::Face face;
    face.setPos(0, 1);
    face.setPos(1, -1);
    face.setPos(2, 0);

    REQUIRE(face.getPos(0) == 1);
    REQUIRE(face.getPos(1) == -1);
    REQUIRE(face.getPos(2) == 0);
  }

  SECTION("Face getSet and getUnset") {
    stmesh::VoxelComplex<3>::Face face;
    face.setPos(0, 1);
    face.setPos(1, -1);
    face.setPos(2, 0);

    REQUIRE(face.getSet() == 0b00000101);
    REQUIRE(face.getUnset() == 0b00010000);
  }

  SECTION("Face getDim") {
    stmesh::VoxelComplex<3>::Face face;
    REQUIRE(face.getDim() == 3);

    face.setPos(0, 1);
    REQUIRE(face.getDim() == 2);

    face.setPos(1, -1);
    REQUIRE(face.getDim() == 1);

    face.setPos(2, 1);
    REQUIRE(face.getDim() == 0);
  }
}

TEST_CASE("VoxelComplex construction and basic operations", "[voxelcomplex]") {
  SECTION("Construction from dimensions") {
    stmesh::VoxelComplex<3> voxel_complex({4, 5, 6});
    REQUIRE(voxel_complex.index({0, 0, 0}) == 0);
    REQUIRE(voxel_complex.index({3, 4, 5}) == 4 * 5 * 6 - 1);
  }

  SECTION("Construction from vector") {
    std::vector<std::vector<std::vector<bool>>> data = {{{true, false}, {false, true}}, {{false, true}, {true, false}}};
    stmesh::VoxelComplex<3> voxel_complex(data);
    REQUIRE(voxel_complex[{0, 0, 0}] == true);
    REQUIRE(voxel_complex[{1, 0, 0}] == false);
    REQUIRE(voxel_complex[{0, 1, 0}] == false);
    REQUIRE(voxel_complex[{1, 1, 0}] == true);
    REQUIRE(voxel_complex[{0, 0, 1}] == false);
    REQUIRE(voxel_complex[{1, 0, 1}] == true);
    REQUIRE(voxel_complex[{0, 1, 1}] == true);
    REQUIRE(voxel_complex[{1, 1, 1}] == false);
  }
}

TEST_CASE("VoxelComplex offsetByFace", "[voxelcomplex]") {
  stmesh::VoxelComplex<3> voxel_complex({4, 5, 6});
  stmesh::VoxelComplex<3>::Face face;

  SECTION("Offset in positive direction") {
    face.setPos(0, 1);
    face.setPos(1, 0);
    face.setPos(2, 0);
    REQUIRE(voxel_complex.offsetByFace(0, face) == 1);
  }

  SECTION("Offset in negative direction") {
    face.setPos(0, 0);
    face.setPos(1, -1);
    face.setPos(2, 0);
    REQUIRE(voxel_complex.offsetByFace(5, face) == 1);
  }

  SECTION("Offset in multiple directions") {
    face.setPos(0, 1);
    face.setPos(1, -1);
    face.setPos(2, 1);
    REQUIRE(voxel_complex.offsetByFace(25, face) == 25 + 1 - 4 + 20);
  }
}

TEST_CASE("VoxelComplex index and operator[]", "[voxelcomplex]") {
  std::vector<std::vector<std::vector<bool>>> data = {{{true, false}, {false, true}}, {{false, true}, {true, false}}};
  stmesh::VoxelComplex<3> voxel_complex(data);

  SECTION("Index calculation") {
    REQUIRE(voxel_complex.index({0, 0, 0}) == 0);
    REQUIRE(voxel_complex.index({1, 0, 0}) == 1);
    REQUIRE(voxel_complex.index({0, 1, 0}) == 2);
    REQUIRE(voxel_complex.index({1, 1, 1}) == 7);
  }

  SECTION("Operator[] access") {
    REQUIRE(voxel_complex[{0, 0, 0}] == true);
    REQUIRE(voxel_complex[{1, 0, 0}] == false);
    REQUIRE(voxel_complex[{0, 1, 0}] == false);
    REQUIRE(voxel_complex[{1, 1, 0}] == true);
    REQUIRE(voxel_complex[{0, 0, 1}] == false);
    REQUIRE(voxel_complex[{1, 0, 1}] == true);
    REQUIRE(voxel_complex[{0, 1, 1}] == true);
    REQUIRE(voxel_complex[{1, 1, 1}] == false);
  }
}

TEST_CASE("Core works correctly", "[voxelcomplex][core]") {
  std::vector<std::vector<std::vector<bool>>> data = {
      {{true, false, false}, {false, false, false}, {false, true, true}},
      {{false, false, false}, {false, true, true}, {false, true, false}},
      {{false, true, false}, {false, false, false}, {false, true, false}}};
  stmesh::VoxelComplex<3> voxel_complex(data);
  auto core = voxel_complex.core(voxel_complex.index({1, 1, 1}));
  stmesh::VoxelComplex<3>::Face face;
  for (const int i : {-1, 0, 1}) {
    face.setPos(0, i);
    for (const int j : {-1, 0, 1}) {
      face.setPos(1, j);
      for (const int k : {-1, 0, 1}) {
        face.setPos(2, k);
        std::set<size_t> unique_x_pos;
        std::set<size_t> unique_y_pos;
        std::set<size_t> unique_z_pos;
        for (const size_t x : {static_cast<size_t>(1), static_cast<size_t>(1 + i)}) {
          for (const size_t y : {static_cast<size_t>(1), static_cast<size_t>(1 + j)}) {
            for (const size_t z : {static_cast<size_t>(1), static_cast<size_t>(1 + k)}) {
              if (voxel_complex[{x, y, z}]) {
                unique_x_pos.insert(x);
                unique_y_pos.insert(y);
                unique_z_pos.insert(z);
              }
            }
          }
        }
        unique_x_pos.erase(1);
        unique_y_pos.erase(1);
        unique_z_pos.erase(1);
        REQUIRE(core[face.pos] == unique_x_pos.size() + unique_y_pos.size() + unique_z_pos.size());
      }
    }
  }
}

TEST_CASE("End to end thinning", "[voxelcomplex][thinning]") {
  std::vector<std::vector<std::vector<bool>>> data(3, std::vector<std::vector<bool>>(5, std::vector<bool>(7, false)));
  data[1][2][3] = true; // center
  data[1][1][3] = true; // right
  data[1][3][3] = true; // left
  data[1][2][1] = true; // top1
  data[1][2][2] = true; // top2
  data[1][2][4] = true; // bottom1
  data[1][2][5] = true; // bottom2
  stmesh::VoxelComplex<3> voxel_complex(data);

  SECTION("Thinning twice") {
    REQUIRE(voxel_complex.thinningStep());
    std::set<std::array<size_t, 3>> expected = {{3, 2, 1}, {2, 2, 1}, {4, 2, 1}};
    voxel_complex.iterateVoxels(
        [&](const std::array<size_t, 3> &coords) { REQUIRE(expected.contains(coords) == voxel_complex[coords]); });
    REQUIRE(voxel_complex.thinningStep());
    voxel_complex.iterateVoxels([&](const std::array<size_t, 3> &coords) {
      REQUIRE((coords[0] == 3 && coords[1] == 2 && coords[2] == 1) == voxel_complex[coords]);
    });
    REQUIRE_FALSE(voxel_complex.thinningStep());
    voxel_complex.iterateVoxels([&](const std::array<size_t, 3> &coords) {
      REQUIRE((coords[0] == 3 && coords[1] == 2 && coords[2] == 1) == voxel_complex[coords]);
    });
  }

  SECTION("Thinning with fixing") {
    voxel_complex.fixOneNeighbor();
    std::set<std::array<size_t, 3>> fixed = {{1, 2, 1}, {5, 2, 1}};
    voxel_complex.iterateVoxels([&](const std::array<size_t, 3> &coords) {
      REQUIRE(voxel_complex.fixed()[voxel_complex.index(coords)] == fixed.contains(coords));
    });
    REQUIRE(voxel_complex.thinningStep());
    std::set<std::array<size_t, 3>> expected = {{1, 2, 1}, {2, 2, 1}, {3, 2, 1}, {4, 2, 1}, {5, 2, 1}};
    voxel_complex.iterateVoxels(
        [&](const std::array<size_t, 3> &coords) { REQUIRE(expected.contains(coords) == voxel_complex[coords]); });
    REQUIRE_FALSE(voxel_complex.thinningStep());
    voxel_complex.iterateVoxels(
        [&](const std::array<size_t, 3> &coords) { REQUIRE(expected.contains(coords) == voxel_complex[coords]); });
  }
}

TEST_CASE("Multithreaded sphere thinning", "[voxelcomplex][sphere_thin]") {
  constexpr ptrdiff_t kSize = 50;
  std::vector data(kSize, std::vector(kSize, std::vector(kSize, std::vector<bool>(kSize, false))));

  for (ptrdiff_t i = 0; i < kSize; ++i) {
    for (ptrdiff_t j = 0; j < kSize; ++j) {
      for (ptrdiff_t k = 0; k < kSize; ++k) {
        for (ptrdiff_t l = 0; l < kSize; ++l) {
          data[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)][static_cast<size_t>(l)] =
              (i - kSize / 2) * (i - kSize / 2) + (j - kSize / 2) * (j - kSize / 2) +
                  (k - kSize / 2) * (k - kSize / 2) + (l - kSize / 2) * (l - kSize / 2) <
              20 * 20;
        }
      }
    }
  }

  stmesh::VoxelComplex complex(data);

  while (complex.thinningStep(4))
    ;

  int count = 0;
  complex.iterateVoxels([&](const auto &coords) {
    if (complex[coords])
      count++;
  });
  REQUIRE(count == 1);
}
