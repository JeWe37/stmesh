#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "stmesh/bitset.hpp"

TEST_CASE("Bitset construction and basic operations", "[bitset][construction_ops]") {
  SECTION("Default constructor") {
    stmesh::Bitset b;
    REQUIRE(b.size() == 0);
  }

  SECTION("Size constructor") {
    stmesh::Bitset b(100);
    REQUIRE(b.size() == 100);
    for (size_t i = 0; i < 100; ++i)
      REQUIRE(b[i] == false);
  }

  SECTION("Vector<bool> constructor") {
    std::vector<bool> v = {true, false, true, true, false};
    stmesh::Bitset b(v);
    REQUIRE(b.size() == 5);
    REQUIRE(b[0] == true);
    REQUIRE(b[1] == false);
    REQUIRE(b[2] == true);
    REQUIRE(b[3] == true);
    REQUIRE(b[4] == false);
  }

  SECTION("Copy constructor") {
    stmesh::Bitset b1(100);
    b1.set(50);
    stmesh::Bitset b2(b1);
    REQUIRE(b2.size() == 100);
    REQUIRE(b2[50] == true);
    for (size_t i = 0; i < 100; ++i) {
      if (i != 50)
        REQUIRE(b2[i] == false);
    }
  }

  SECTION("Move constructor") {
    stmesh::Bitset b1(100);
    b1.set(50);
    stmesh::Bitset b2(std::move(b1));
    REQUIRE(b2.size() == 100);
    REQUIRE(b2[50] == true);
    for (size_t i = 0; i < 100; ++i) {
      if (i != 50)
        REQUIRE(b2[i] == false);
    }
  }
}

TEST_CASE("Bitset set operations", "[bitset][set_ops]") {
  SECTION("Set single bit") {
    stmesh::Bitset b(100);
    b.set(50);
    REQUIRE(b[50] == true);
    for (size_t i = 0; i < 100; ++i) {
      if (i != 50)
        REQUIRE(b[i] == false);
    }
  }

  SECTION("Set range of bits") {
    stmesh::Bitset b(100);
    b.setRange(25, 75);
    for (size_t i = 0; i < 100; ++i)
      REQUIRE(b[i] == (i >= 25 && i <= 75));
  }

  SECTION("Set from vector<bool>") {
    stmesh::Bitset b(10);
    std::vector<bool> v = {true, false, true, true, false, true, false, true, true, true};
    b.setFrom(v);
    for (size_t i = 0; i < 10; ++i)
      REQUIRE(b[i] == v[i]);
  }
}

TEST_CASE("Bitset edge cases", "[bitset][edge_cases]") {
  SECTION("Out of range access") {
    stmesh::Bitset b(100);
    REQUIRE_THROWS_AS(b[100], std::out_of_range);
  }

  SECTION("Set range out of bounds") {
    stmesh::Bitset b(100);
    REQUIRE_THROWS_AS(b.setRange(90, 110), std::out_of_range);
  }

  SECTION("Set range with swapped indices") {
    stmesh::Bitset b(100);
    b.setRange(75, 25); // This should work, just swap internally
    for (size_t i = 0; i < 100; ++i)
      REQUIRE(b[i] == (i >= 25 && i <= 75));
  }
}

TEST_CASE("Bitset assignment operators", "[bitset][assignment]") {
  SECTION("Copy assignment") {
    stmesh::Bitset b1(100);
    b1.set(50);
    stmesh::Bitset b2;
    b2 = b1;
    REQUIRE(b2.size() == 100);
    REQUIRE(b2[50] == true);
    for (size_t i = 0; i < 100; ++i) {
      if (i != 50)
        REQUIRE(b2[i] == false);
    }
  }

  SECTION("Move assignment") {
    stmesh::Bitset b1(100);
    b1.set(50);
    stmesh::Bitset b2;
    b2 = std::move(b1);
    REQUIRE(b2.size() == 100);
    REQUIRE(b2[50] == true);
    for (size_t i = 0; i < 100; ++i) {
      if (i != 50)
        REQUIRE(b2[i] == false);
    }
  }
}
