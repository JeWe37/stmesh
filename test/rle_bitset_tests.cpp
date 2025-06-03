#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <forward_list>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <vector>

#include "stmesh/rle_bitset.hpp"

TEST_CASE("RleBitset construction and basic operations", "[rlebitset][construction]") {
  SECTION("Default constructor") {
    stmesh::RleBitset b;
    REQUIRE(b.size() == 0);
  }

  SECTION("Size constructor") {
    stmesh::RleBitset b(100);
    REQUIRE(b.size() == 100);
    for (size_t i = 0; i < 100; ++i)
      REQUIRE(b[i] == false);
  }

  SECTION("Vector<bool> constructor") {
    std::vector<bool> v = {true, false, true, true, false, true, true, true};
    stmesh::RleBitset b(v);
    REQUIRE(b.size() == 8);
    REQUIRE(b[0] == true);
    REQUIRE(b[1] == false);
    REQUIRE(b[2] == true);
    REQUIRE(b[3] == true);
    REQUIRE(b[4] == false);
    REQUIRE(b[5] == true);
    REQUIRE(b[6] == true);
    REQUIRE(b[7] == true);
  }
}

TEST_CASE("RleBitset set operations", "[rlebitset]") {
  SECTION("Set single bit") {
    stmesh::RleBitset b(100);
    b.registerThreads(1);
    b.set(50, 0);
    b.commit(0);
    b.unregisterThreads();
    REQUIRE(b[50] == true);
    for (size_t i = 0; i < 100; ++i) {
      if (i != 50)
        REQUIRE(b[i] == false);
    }
  }

  SECTION("Set range of bits") {
    stmesh::RleBitset b(100);
    std::forward_list<stmesh::RleBitset::Run> runs;
    b.setRange(10, 20, runs.before_begin(), runs);
    REQUIRE(std::distance(runs.begin(), runs.end()) == 1);
    b.setRange(15, 25, runs.before_begin(), runs);
    REQUIRE(std::distance(runs.begin(), runs.end()) == 1);
    b.setRange(25, 50, runs.before_begin(), runs);
    REQUIRE(std::distance(runs.begin(), runs.end()) == 1);
    b.setRange(75, 99, runs.before_begin(), runs);
    REQUIRE(std::distance(runs.begin(), runs.end()) == 2);
    for (size_t i = 0; i < 100; ++i)
      REQUIRE(b[i] == ((i >= 10 && i <= 50) || i >= 75));
  }
}

TEST_CASE("RleBitset threading support", "[rlebitset][threading]") {
  SECTION("Multi-threaded set operation") {
    stmesh::RleBitset b(1000);
    b.registerThreads(4);

    auto setOperation = [&](size_t start, size_t end, size_t thread_id) {
      for (size_t i = start; i < end; ++i) {
        if (i % 2 == 0)
          b.set(i, thread_id);
      }
      b.commit(thread_id);
    };

    std::vector<std::thread> threads;
    for (size_t i = 0; i < 4; ++i)
      threads.emplace_back(setOperation, i * 250, (i + 1) * 250, i);

    for (auto &thread : threads)
      thread.join();

    b.unregisterThreads();

    for (size_t i = 0; i < 1000; ++i)
      REQUIRE(b[i] == (i % 2 == 0));
  }
}

TEST_CASE("RleBitset iterateSet", "[rlebitset][iterate_set]") {
  SECTION("Single-threaded iteration") {
    stmesh::RleBitset b(100);
    b.registerThreads(1);
    for (size_t i = 0; i < 100; i += 2)
      b.set(i, 0);
    b.commit(0);
    b.unregisterThreads();

    std::vector<size_t> setIndices;
    b.iterateSet([&](size_t idx, size_t) { setIndices.push_back(idx); });

    REQUIRE(setIndices.size() == 50);
    for (size_t i = 0; i < 50; ++i)
      REQUIRE(setIndices[i] == i * 2);
  }

  SECTION("Multi-threaded iteration") {
    stmesh::RleBitset b(1000);
    b.registerThreads(1);
    for (size_t i = 0; i < 1000; i += 3)
      b.set(i, 0);
    b.commit(0);
    b.unregisterThreads();

    std::vector<std::vector<size_t>> threadSetIndices(4);
    b.iterateSet([&](size_t idx, size_t thread_id) { threadSetIndices[thread_id].push_back(idx); },
                 [](size_t) {}, // initialize
                 [](size_t) {}, // finalize
                 4              // n_threads
    );

    const size_t totalSetIndices =
        std::accumulate(threadSetIndices.begin(), threadSetIndices.end(), size_t{},
                        [](size_t total, const std::vector<size_t> &indices) { return total + indices.size(); });

    REQUIRE(totalSetIndices == 334); // 1000 / 3 rounded up

    // Verify that each thread's indices are in order and divisible by 3
    for (const auto &indices : threadSetIndices) {
      REQUIRE(std::is_sorted(indices.begin(), indices.end()));
      for (size_t idx : indices)
        REQUIRE(idx % 3 == 0);
    }
  }
}

TEST_CASE("RleBitset edge cases", "[rlebitset][edge_cases]") {
  SECTION("Set out of range") {
    stmesh::RleBitset b(100);
    b.registerThreads(1);
    REQUIRE_THROWS_AS(b.set(100, 0), std::out_of_range);
    b.unregisterThreads();
  }

  SECTION("Set range out of bounds") {
    stmesh::RleBitset b(100);
    std::forward_list<stmesh::RleBitset::Run> runs;
    REQUIRE_THROWS_AS(b.setRange(90, 110, runs.before_begin(), runs), std::out_of_range);
  }
}
