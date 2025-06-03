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
/// A run-length encoded bitset
/**
 * A run-length encoded bitset that allows for efficient iteration over the set bits. The bitset is stored as a list of
 * runs of set bits. The runs are stored in a forward list, which allows for efficient insertion and deletion of runs.
 * For efficient access by index, the same data is also stored as a classical bitset.
 *
 * Methods are provided for allowing the registration of threads for setting bits in parallel.
 */
class RleBitset : public Bitset {
public:
  /// A run of set bits
  /**
   * A run of set bits. A run is defined by a start index and a length. The run starts at the start index and has the
   * given length.
   */
  struct Run {
    size_t start;
    size_t length;

    /// Compare two runs
    /**
     * Compare two runs. The runs are compared by their start index, then by their length.
     */
    [[nodiscard]] bool operator<=>(const Run &other) const = default;
  };

  using hint_type = std::forward_list<Run>::iterator; ///< The type of the hint for the setRange function

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

  /// Construct a bitset of a given size
  /**
   * Construct a bitset of a given size. The bitset is initialized to all bits being unset.
   *
   * @param size The size of the bitset
   */
  explicit RleBitset(size_t size);

  /// Construct a bitset from a vector of booleans
  /**
   * Construct a bitset from a vector of booleans. The bitset is initialized to the values of the vector of booleans.
   *
   * @param data The vector of booleans to initialize the bitset from
   */
  explicit RleBitset(const std::vector<bool> &data);

  /// Set the bitset from a vector of booleans
  /**
   * Set the bitset from a vector of booleans. The bitset is set to the values of the vector of booleans.
   *
   * @param data The vector of booleans to set the bitset from
   */
  // cppcheck-suppress duplInheritedMember
  void setFrom(const std::vector<bool> &data);

  /// Iterate over the set bits
  /**
   * Iterate over the set bits. The function is called for each set bit with the index of the bit as the first argument
   * and the thread id as the second argument. The function is called in the order of the set bits.
   *
   * @param func The function to call for each set bit
   * @param initialize The function to call before the iteration of a thread starts, defaults to a no-op
   * @param finalize The function to call after the iteration of a thread ends, defaults to a no-op
   * @param n_threads The number of threads to use for the iteration
   * @tparam F The type of the function to call for each set bit
   * @tparam E The type of the function to call before the iteration of a thread starts
   * @tparam B The type of the function to call after the iteration of a thread ends
   */
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

  /// Commits the current update for a thread
  /**
   * Commits the current update for a thread. The current update is the range of bits that have been set since the last
   * commit. The range is set in the bitset and the runs are updated accordingly. This function should be called for
   * every thread before unregistering the threads.
   *
   * @param thread_id The id of the thread to commit the update for
   * @sa registerThreads, restartIteration, finalizeThread
   */
  void commit(size_t thread_id);

  /// Registers threads for iteration
  /**
   * Registers threads for iteration. This function should be called before starting the iteration over the set bits.
   * Only thread_ids from 0 to n_threads - 1 are valid, for any functions that require them.
   *
   * @param n_threads The number of threads to register
   */
  void registerThreads(size_t n_threads);

  /// Restarts the iteration for a thread
  /**
   * Restarts the iteration for a thread. This function should be called after a thread has finished iterating over the
   * set bits and before starting the iteration again. Without it, bits can only be set in order.
   *
   * @param thread_id The id of the thread to restart the iteration for
   */
  void restartIteration(size_t thread_id);

  /// Unregisters threads after iteration
  /**
   * Unregisters threads after iteration. This function should be called after all threads have finished iterating over
   * the set bits. It copies the runs from the threads to the main list of runs.
   */
  void unregisterThreads();

  /// Set a bit
  /**
   * Set a bit in the bitset. The bit at the given index is set to true. To be used only with threads that were
   * previously registered.
   *
   * @param idx The index of the bit to set
   * @param thread_id The id of the thread to set the bit for
   * @sa registerThreads, commit, restartIteration, finalizeThread
   */
  void set(size_t idx, size_t thread_id);

  /// Set a range of bits
  /**
   * Set a range of bits in the bitset. The bits in the range from idx_start to idx_end are set to true. This only
   * operates on the runs given to it, and is unrelated to threading.
   *
   * @param idx_start The start index of the range to set
   * @param idx_end The end index of the range to set
   * @param hint The hint for the setRange function, an iterator into runs
   * @param runs The list of runs to set the range in
   */
  hint_type setRange(size_t idx_start, size_t idx_end, hint_type hint, std::forward_list<Run> &runs);

  /// Compares two RleBitsets
  /**
   * Compares two RleBitsets. The RleBitsets are compared by their runs, therefore no threads should be registered.
   *
   * @param other The RleBitset to compare to
   * @return True if the RleBitsets are equal, false otherwise
   */
  [[nodiscard]] bool operator==(const RleBitset &other) const;

private:
  std::vector<CurrentUpdate> current_updates_;
  std::vector<Run> runs_;
};
} // namespace stmesh

#endif
