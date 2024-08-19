#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <vector>

#include <fuzzer/FuzzedDataProvider.h>

#include "stmesh/rle_bitset.hpp"

// NOLINTNEXTLINE(performance-unnecessary-value-param)
bool verifyBitset(std::vector<bool> bits, const stmesh::RleBitset &rle_bitset) {
  if (bits.size() != rle_bitset.size())
    return false;
  // NOLINTNEXTLINE(misc-const-correctness)
  for (size_t i = 0; i < bits.size(); ++i)
    if (bits[i] != rle_bitset[i])
      return false;
  const ptrdiff_t last_idx = -1;
  // NOLINTNEXTLINE(misc-const-correctness)
  bool out_of_order = false;
  rle_bitset.iterateSet([&](size_t idx, size_t /*unused*/) {
    if (static_cast<ptrdiff_t>(idx) <= last_idx)
      out_of_order = true;
    bits[idx] = !bits[idx];
  });
  return !out_of_order && !std::ranges::any_of(bits, std::identity{});
}

// Fuzzer that attempts to invoke undefined behavior for signed integer overflow
// cppcheck-suppress unusedFunction symbolName=LLVMFuzzerTestOneInput
extern "C" int LLVMFuzzerTestOneInput(const uint8_t data[], size_t size) { // NOLINT(readability-identifier-naming)
  FuzzedDataProvider data_provider(data, size);

  const size_t num_bits = data_provider.ConsumeIntegralInRange<unsigned>(0U, 3000U);
  std::vector<bool> bits;
  bits.reserve(num_bits);
  std::generate_n(std::back_inserter(bits), num_bits, [&] { return data_provider.ConsumeBool(); });

  stmesh::RleBitset rle_bitset(num_bits);
  rle_bitset.registerThreads(1);
  for (size_t i = 0; i < num_bits; ++i) {
    if (bits[i])
      rle_bitset.set(i, 0);
  }
  rle_bitset.commit(0);

  if (bits.empty())
    return 0;

  rle_bitset.restartIteration(0);
  auto modification_start = data_provider.ConsumeIntegralInRange<size_t>(0, bits.size() - 1);
  auto modification_end = data_provider.ConsumeIntegralInRange<size_t>(modification_start, bits.size() - 1);
  for (size_t i = modification_start; i < modification_end; ++i) {
    rle_bitset.set(i, 0);
    bits[i] = true;
  }
  rle_bitset.commit(0);
  rle_bitset.unregisterThreads();

  if (!verifyBitset(bits, rle_bitset))
    __builtin_trap();

  return 0;
}
