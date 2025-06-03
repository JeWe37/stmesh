#ifndef STMESH_BITSET_HPP
#define STMESH_BITSET_HPP
#include <atomic>
#include <cstddef>
#include <vector>

namespace stmesh {
/// A simple bitset class
/**
 * A simple bitset class that allows setting and getting bits. Setting bits is thread-safe thanks to the use of
 * atomic variables.
 */
class Bitset {
public:
  Bitset();
  Bitset(const Bitset &other);
  Bitset(Bitset &&other) noexcept;
  Bitset &operator=(const Bitset &other);
  Bitset &operator=(Bitset &&other) noexcept;
  ~Bitset();

  /// Constructor
  /**
   * Constructor that initializes the bitset to a given size.
   *
   * @param size The size of the bitset
   */
  explicit Bitset(size_t size);

  /// Constructor from a vector of booleans
  /**
   * Constructor that initializes the bitset from a vector of booleans.
   *
   * @param data The vector of booleans to initialize the bitset from
   */
  explicit Bitset(const std::vector<bool> &data);

  /// Set the bitset from a vector of booleans
  /**
   * Set the bitset from a vector of booleans.
   *
   * @param data The vector of booleans to set the bitset from
   */
  void setFrom(const std::vector<bool> &data);

  /// Read a bit
  /**
   * Read a bit from the bitset.
   *
   * @param idx The index of the bit to read
   * @return The value of the bit at the given index
   */
  [[nodiscard]] bool operator[](size_t idx) const;

  /// Set a bit
  /**
   * Set a bit in the bitset.
   *
   * @param idx The index of the bit to set
   */
  void set(size_t idx);

  /// Set a range of bits
  /**
   * Set a range of bits in the bitset.
   *
   * @param idx_start The start index of the range to set
   * @param idx_end The end index of the range to set
   */
  void setRange(size_t idx_start, size_t idx_end);

  /// Get the size of the bitset
  /**
   * Get the size of the bitset.
   *
   * @return The size of the bitset
   */
  [[nodiscard]] size_t size() const noexcept;

protected:
  size_t size_;                            ///< The size of the bitset
  std::vector<std::atomic_uint64_t> data_; ///< The data contained in the bitset
};
} // namespace stmesh

#endif
