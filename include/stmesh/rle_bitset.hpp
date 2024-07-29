#ifndef STMESH_RLE_BITSET_HPP
#define STMESH_RLE_BITSET_HPP
#include <algorithm>
#include <atomic>
#include <bits/ranges_algo.h>
#include <concepts>
#include <cstddef>
#include <forward_list>
#include <iterator>
#include <thread>
#include <vector>

namespace stmesh {
class RleBitset {
  struct Run {
    size_t start;
    size_t length;
  };

  using hint_type = std::forward_list<Run>::iterator;

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

  void setFrom(const std::vector<bool> &data);

  template <std::invocable<size_t, size_t> F, std::invocable<size_t> E, std::invocable<size_t> B>
  constexpr void iterateSet(const F &func, const B &initialize, const E &finalize, size_t n_threads) const {
    std::vector<std::jthread> threads;
    for (size_t i = 0; i < n_threads; ++i)
      threads.emplace_back([&, i] {
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
      });
  }

  [[nodiscard]] bool operator[](size_t idx) const;

  void commit(size_t thread_id);

  void registerThreads(size_t n_threads);

  void restartIteration(size_t thread_id);

  void finalizeThread(size_t thread_id);

  void unregisterThreads();

  void set(size_t idx, size_t thread_id);

  hint_type setRange(size_t idx_start, size_t idx_end, hint_type hint, std::forward_list<Run> &runs);

  [[nodiscard]] size_t size() const noexcept;

private:
  std::vector<CurrentUpdate> current_updates_;
  size_t size_;
  std::vector<std::atomic_uint64_t> data_;
  std::vector<Run> runs_;
};
} // namespace stmesh

#endif
