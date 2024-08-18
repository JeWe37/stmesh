#ifndef STMESH_BITSET_HPP
#define STMESH_BITSET_HPP
#include <atomic>
#include <cstddef>
#include <vector>

namespace stmesh {
class Bitset {
public:
  Bitset();
  Bitset(const Bitset &other);
  Bitset(Bitset &&other) noexcept;
  Bitset &operator=(const Bitset &other);
  Bitset &operator=(Bitset &&other) noexcept;
  ~Bitset();

  explicit Bitset(size_t size);

  explicit Bitset(const std::vector<bool> &data);

  void setFrom(const std::vector<bool> &data);

  [[nodiscard]] bool operator[](size_t idx) const;

  void set(size_t idx);

  void setRange(size_t idx_start, size_t idx_end);

  [[nodiscard]] size_t size() const noexcept;

protected:
  size_t size_;
  std::vector<std::atomic_uint64_t> data_;
};
} // namespace stmesh

#endif
