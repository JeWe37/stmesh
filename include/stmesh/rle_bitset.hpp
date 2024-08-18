#ifndef STMESH_RLE_BITSET_HPP
#define STMESH_RLE_BITSET_HPP
#include <algorithm>
#include <bits/ranges_algo.h>
#include <concepts>
#include <cstddef>
#include <forward_list>
#include <iterator>
#include <thread>
#include <vector>

#include "bitset.hpp"

namespace stmesh {
class RleBitset : public Bitset {
public:
  struct Run {
    size_t start;
    size_t length;
    [[nodiscard]] bool operator<=>(const Run &other) const = default;
  };

  using hint_type = std::forward_list<Run>::iterator;

private:
  struct CurrentUpdate {
    hint_type hint;
    std::forward_list<Run> thread_runs;
    size_t start{};
    size_t length{};
  };

  template <std::forward_iterator T> [[nodiscard]] static T locateRun(T hint, size_t idx);

public:
  RleBitset();
  RleBitset(const RleBitset &other);
  RleBitset(RleBitset &&other) noexcept;
  RleBitset &operator=(const RleBitset &other);
  RleBitset &operator=(RleBitset &&other) noexcept;
  ~RleBitset();

  explicit RleBitset(size_t size);

  explicit RleBitset(const std::vector<bool> &data);

  // cppcheck-suppress duplInheritedMember
  void setFrom(const std::vector<bool> &data);

  template <std::invocable<size_t, size_t> F, std::invocable<size_t> E = decltype([](size_t /*unused*/) {}),
            std::invocable<size_t> B = decltype([](size_t /*unused*/) {})>
  constexpr void iterateSet(const F &func, const B &initialize = {}, const E &finalize = {},
                            size_t n_threads = 1) const {
    std::vector<std::jthread> threads;
    for (size_t i = 0; i < n_threads; ++i) {
      auto task = [&, i] {
        const size_t start = size_ / n_threads * i;
        auto it = std::ranges::lower_bound(runs_, start, {}, [](const Run &run) { return run.start + run.length; });
        initialize(i);
        for (; it != runs_.end(); it++) {
          for (size_t j = std::max(it->start, start); j < it->start + it->length; ++j) {
            if (j >= size_ / n_threads * (i + 1)) {
              finalize(i);
              return;
            }
            func(j, i);
          }
        }
        finalize(i);
      };
      if (n_threads == 1)
        task();
      else
        threads.emplace_back(task);
    }
  }

  void commit(size_t thread_id);

  void registerThreads(size_t n_threads);

  void restartIteration(size_t thread_id);

  void finalizeThread(size_t thread_id);

  void unregisterThreads();

  void set(size_t idx, size_t thread_id);

  hint_type setRange(size_t idx_start, size_t idx_end, hint_type hint, std::forward_list<Run> &runs);

  [[nodiscard]] bool operator==(const RleBitset &other) const;

private:
  std::vector<CurrentUpdate> current_updates_;
  std::vector<Run> runs_;
};
} // namespace stmesh

#endif
