#include "stmesh/bitset.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace stmesh {
static constexpr size_t k64Mask = 0x3F;
static constexpr size_t k64Bits = 6;

Bitset::Bitset() : size_(0) {}

// cppcheck-suppress missingMemberCopy
Bitset::Bitset(const Bitset &other) : size_(other.size_), data_(other.data_.size()) {
  for (size_t i = 0; i < other.data_.size(); ++i)
    data_[i] = other.data_[i].load();
}

Bitset::Bitset(Bitset &&other) noexcept : size_(other.size_), data_(std::move(other.data_)) {}

Bitset &Bitset::operator=(const Bitset &other) {
  if (this == &other)
    return *this;
  size_ = other.size_;
  data_ = std::vector<std::atomic_uint64_t>(other.data_.size());
  for (size_t i = 0; i < other.data_.size(); ++i)
    data_[i] = other.data_[i].load();
  return *this;
}

Bitset &Bitset::operator=(Bitset &&other) noexcept {
  size_ = other.size_;
  data_ = std::move(other.data_);
  return *this;
}

Bitset::~Bitset() = default;

Bitset::Bitset(size_t size) : size_(size), data_((size >> k64Bits) + 1UL) { std::ranges::fill(data_, 0UL); }

Bitset::Bitset(const std::vector<bool> &data) : Bitset(data.size()) { setFrom(data); }

void Bitset::setFrom(const std::vector<bool> &data) {
  for (size_t i = 0; i < data.size(); ++i) {
    if (data[i])
      data_[i >> k64Bits] |= 1UL << (i & k64Mask);
  }
}

[[nodiscard]] bool Bitset::operator[](size_t idx) const {
  if (idx >= size_)
    throw std::out_of_range("Index out of range");
  return static_cast<bool>(data_[idx >> k64Bits] & (1UL << (idx & k64Mask)));
}

void Bitset::set(size_t idx) { data_[idx >> k64Bits] |= 1UL << (idx & k64Mask); }

void Bitset::setRange(size_t idx_start, size_t idx_end) {
  if (idx_start >= size_ || idx_end >= size_)
    throw std::out_of_range("Index out of range");
  if (idx_start > idx_end)
    std::swap(idx_start, idx_end);

  // Set the bits in the underlying data structure
  const size_t start_word = idx_start >> 6ULL;
  const size_t end_word = idx_end >> 6ULL;
  const uint64_t start_mask = ~0ULL << (idx_start & 63ULL);
  const uint64_t end_mask = ~0ULL >> (63ULL - (idx_end & 63ULL));

  if (start_word == end_word) {
    const uint64_t mask = start_mask & end_mask;
    data_[start_word] |= mask;
  } else {
    data_[start_word] |= start_mask;
    for (size_t i = start_word + 1; i < end_word; ++i)
      data_[i] = ~0ULL;
    data_[end_word] |= end_mask;
  }
}

[[nodiscard]] size_t Bitset::size() const noexcept { return size_; }
} // namespace stmesh
