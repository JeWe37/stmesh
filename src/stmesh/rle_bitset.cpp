#include "stmesh/rle_bitset.hpp"

#include <algorithm>
#include <cstddef>
#include <forward_list>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>

#include "stmesh/bitset.hpp"

namespace stmesh {
static constexpr size_t k64Mask = 0x3F;
static constexpr size_t k64Bits = 6;

RleBitset::RleBitset() = default;

// cppcheck-suppress missingMemberCopy
RleBitset::RleBitset(const RleBitset &other) : Bitset(other), runs_(other.runs_) {}

RleBitset::RleBitset(RleBitset &&other) noexcept
    : current_updates_(std::move(other.current_updates_)), runs_(std::move(other.runs_)) {}

RleBitset &RleBitset::operator=(const RleBitset &other) {
  if (this == &other)
    return *this;
  Bitset::operator=(other);
  current_updates_ = std::vector<CurrentUpdate>();
  runs_ = other.runs_;
  return *this;
}

RleBitset &RleBitset::operator=(RleBitset &&other) noexcept {
  current_updates_ = std::move(other.current_updates_);
  size_ = other.size_;
  data_ = std::move(other.data_);
  runs_ = std::move(other.runs_);
  return *this;
}

RleBitset::~RleBitset() = default;

RleBitset::RleBitset(size_t size) : Bitset(size) { std::ranges::fill(data_, 0UL); }

RleBitset::RleBitset(const std::vector<bool> &data) : RleBitset(data.size()) { setFrom(data); }

// cppcheck-suppress duplInheritedMember
void RleBitset::setFrom(const std::vector<bool> &data) {
  size_t start = 0;
  size_t length = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    if (data[i]) {
      data_[i >> k64Bits] |= 1UL << (i & k64Mask);
      if (length == 0)
        start = i;
      length++;
    } else {
      if (length != 0)
        runs_.emplace_back(start, length);
      length = 0;
    }
  }
  if (length != 0)
    runs_.emplace_back(start, length);
}

template <std::forward_iterator T> [[nodiscard]] T RleBitset::locateRun(T hint, size_t idx) {
  auto it_old = hint;
  // actually valid to use end from runs_ here
  for (auto it = ++hint; it != std::forward_list<Run>().end(); it++) {
    if (it->start + it->length > idx)
      return it_old;
    it_old = it;
  }
  return it_old;
}

void RleBitset::commit(size_t thread_id) {
  CurrentUpdate &current_update = current_updates_[thread_id];
  if (current_update.length == 0)
    return;
  current_update.hint = setRange(current_update.start, current_update.start + current_update.length - 1,
                                 current_update.hint, current_update.thread_runs);
  current_update.length = 0;
}

void RleBitset::registerThreads(size_t n_threads) {
  current_updates_ = std::vector<CurrentUpdate>(n_threads);
  size_t thread_id = 0;
  for (CurrentUpdate &current_update : current_updates_) {
    const size_t initial_pos = size_ / current_updates_.size() * thread_id++;
    current_update.hint = current_update.thread_runs.before_begin();
    current_update.start = initial_pos;
    current_update.length = 0;
  }
}

void RleBitset::restartIteration(size_t thread_id) {
  current_updates_[thread_id].hint = current_updates_[thread_id].thread_runs.before_begin();
}

void RleBitset::unregisterThreads() {
  for (const CurrentUpdate &current_update : current_updates_)
    std::ranges::copy(current_update.thread_runs, std::back_inserter(runs_));
  current_updates_.clear();
}

void RleBitset::set(size_t idx, size_t thread_id) {
  CurrentUpdate &current_update = current_updates_[thread_id];
  if (idx >= size_)
    throw std::out_of_range("Index out of range");
  if (current_update.length == 0 || current_update.start + current_update.length != idx) {
    if (current_update.length != 0)
      commit(thread_id);
    current_update.hint = locateRun(current_update.hint, idx);
    current_update.start = idx;
    current_update.length = 1;
  } else
    ++current_update.length;
}

auto RleBitset::setRange(size_t idx_start, size_t idx_end, hint_type hint, std::forward_list<Run> &runs) -> hint_type {
  if (idx_start >= size_ || idx_end >= size_)
    throw std::out_of_range("Index out of range");
  if (idx_start > idx_end)
    std::swap(idx_start, idx_end);

  Bitset::setRange(idx_start, idx_end);

  // Update the run-length encoding
  auto it = locateRun(hint, idx_start);
  auto next = std::next(it);

  // Handle the case where we're setting to 1
  if (next == runs.end() || next->start > idx_end) {
    // Insert a new run
    if (it != runs.before_begin() && it->start + it->length >= idx_start)
      // Extend the previous run
      it->length = std::max(it->length, idx_end - it->start + 1);
    else
      // Insert a new run
      it = runs.insert_after(it, Run{idx_start, idx_end - idx_start + 1});
  } else {
    // Extend the existing run
    const size_t new_start = std::min(next->start, idx_start);
    const size_t new_end = std::max(next->start + next->length - 1, idx_end);
    next->start = new_start;
    next->length = new_end - new_start + 1;

    // Merge with subsequent runs if necessary
    auto current = next;
    ++next;
    while (next != runs.end() && next->start <= current->start + current->length) {
      current->length = std::max(current->length, next->start + next->length - current->start);
      next = runs.erase_after(current);
    }
  }

  return it;
}

[[nodiscard]] bool RleBitset::operator==(const RleBitset &other) const { return runs_ == other.runs_; }
} // namespace stmesh
