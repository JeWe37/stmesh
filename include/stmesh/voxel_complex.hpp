#ifndef STMESH_VOXEL_COMPLEX_HPP
#define STMESH_VOXEL_COMPLEX_HPP

#include <array>
#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

#include "rle_bitset.hpp"
#include "utility.hpp"

// 2h

namespace stmesh {
template <unsigned D> struct BoolVectorDim {
  using type = std::vector<typename BoolVectorDim<D - 1>::type>;
};

template <> struct BoolVectorDim<1> {
  using type = std::vector<bool>;
};

template <unsigned D> using BoolVectorDimType = BoolVectorDim<D>::type;

template <typename T> struct GetBoolVectorDim {
  static constexpr unsigned kValue = GetBoolVectorDim<typename T::value_type>::kValue + 1;
};

template <> struct GetBoolVectorDim<std::vector<bool>> {
  static constexpr unsigned kValue = 1;
};

template <typename T> constexpr unsigned kGetBoolVectorDimV = GetBoolVectorDim<T>::kValue;

template <unsigned D> class VoxelComplex {
  static_assert(D <= 4, "VoxelComplex only supports dimensions up to 4");

public:
  struct Face {
    constexpr static unsigned char kDimMask = (1U << (D << 1U)) - 1U;
    constexpr static size_t kMaxFaces = std::numeric_limits<unsigned char>::max() + 1;

    // A face is represented as the coefficients of the vectors spanning its (oriented) dual space
    // This data is compressed into a byte of 4 2-bit values, each representing -1, 0 or 1
    // -1 becomes 1(0b01), 0 becomes 2(0b10), 1 becomes 3(0b11)
    // NOLINTNEXTLINE(*-magic-numbers)
    unsigned char pos{0b10101010U & kDimMask};

    constexpr void setPos(unsigned i, int val);

    [[nodiscard]] constexpr unsigned char getSet() const;

    [[nodiscard]] constexpr unsigned char getUnset() const;

    [[nodiscard]] constexpr int getPos(unsigned i) const;

    [[nodiscard]] constexpr unsigned char getDim() const;
  };

  template <typename T> [[nodiscard]] static auto getDims(const T &data);

  template <typename T> static void flatten(const T &data, std::vector<bool> &result);

  explicit VoxelComplex(const BoolVectorDimType<D> &data);

  explicit VoxelComplex(const std::array<size_t, D> &dims);

  [[nodiscard]] size_t offsetByFace(size_t idx, const Face face) const;

  [[nodiscard]] size_t index(const std::array<size_t, D> &coords) const;

  [[nodiscard]] bool operator[](const std::array<size_t, D> &coords) const;

  // potentially it might be faster to use the precomputed dim faces, but really this can be optimized out fully
  // usually
  template <std::invocable<const Face> F>
  static consteval void iterateSubfacesRecursive(Face face, const F &func, unsigned d = 0, int func_dim = -1,
                                                 unsigned i = 0) {
    const unsigned dim = face.getDim();
    if (dim == d || i >= D) {
      // only look for faces down to dimension d
      if (func_dim == -1 || static_cast<int>(dim) == func_dim)
        // run for faces of all dimensions or if func_dim is set, only for faces of that dimension
        func(face);
      return;
    }
    if (face.getPos(i) != 0)
      // ensures these are actually subfaces
      iterateSubfacesRecursive(face, func, d, func_dim, i + 1);
    else {
      // iterates over possibilities for remaining free coords
      for (const int val : {-1, 0, 1}) {
        face.setPos(i, val);
        iterateSubfacesRecursive(face, func, d, func_dim, i + 1);
      }
    }
  }

private:
  template <unsigned Dim, typename F, int TargetDim = -1> static void iterateSubfaces(Face face, const F &func);

public:
  [[nodiscard]] std::array<unsigned char, Face::kMaxFaces> core(const size_t idx) const;

  // for each face, we store whether by setting the i-th coordinate to zero we reach an existing face
  // initially, this is always the case
  [[nodiscard]] consteval static std::array<unsigned char, Face::kMaxFaces> initialFaces();

  template <unsigned Dim> struct DimSubfaceStore {
    std::array<size_t, Dim + 1> dims{};
    std::array<Face, expN(Dim, 3)> faces;

    [[nodiscard]] constexpr std::span<const Face> asSpan(unsigned dim) const;
  };

  template <unsigned Dim> [[nodiscard]] consteval static DimSubfaceStore<Dim> dimFaces(const Face face);

  template <int Dim = D> [[nodiscard]] consteval static auto allSortedDims();

  constexpr static auto kAllSortedDims = allSortedDims();

  template <unsigned Dim>
  static constexpr bool essential(const Face face, const std::array<unsigned char, Face::kMaxFaces> &onto);

  template <unsigned Dim>
  static bool collapseFace(const Face face, const std::array<unsigned char, Face::kMaxFaces> &onto);

  // shortcut: if a voxel has other voxels in all directions, it is trivially critical => handles bulk
  [[nodiscard]] bool checkShortcut(size_t idx) const;

  template <std::invocable<std::array<size_t, D>> F, unsigned Dim = D - 1>
  void iterateVoxels(const F &func, std::array<size_t, D> coords = {}) const {
    // parallelize outermost loop, i.e. Dim == D - 1
    // might need an overload for that case
    for (size_t i = 0; i < dims_.at(Dim); ++i) {
      coords.at(Dim) = i;
      if constexpr (Dim == 0)
        func(coords);
      else
        iterateVoxels<F, Dim - 1>(func, coords);
    }
  }

  [[nodiscard]] constexpr bool checkContainingVoxels(size_t idx, const Face face, const RleBitset &removed,
                                                     unsigned i = 0);

  // idea: to deduplicate tests across the clique, set the entire clique in new_table
  // we can then just check for removed && !new_table at the start
  // we cant remove from removed as another face might still need it in checkContainingVoxels
  // current issue: we cannot immediately test for new_table changes, as those have yet to be committed
  // would need to add ability to check for changes in the current update
  // void setContainingVoxels(size_t idx, const Face face, RleBitset &removed, RleBitset::hint_type &hint, unsigned i
  // = 0,
  //                          size_t initial_idx = std::numeric_limits<size_t>::max()) {
  //   if (i == D) {
  //     if (idx == initial_idx || idx == initial_idx + 1)
  //       hint = removed.set(idx, hint, true);
  //     else
  //       removed.setImmediate(idx, hint, true);
  //     return;
  //   }
  //   initial_idx = initial_idx == std::numeric_limits<size_t>::max() ? idx : initial_idx;
  //   setContainingVoxels(idx, face, removed, i + 1, initial_idx);
  //   setContainingVoxels(
  //       static_cast<size_t>(static_cast<ptrdiff_t>(idx) + face.getPos(i) *
  //       static_cast<ptrdiff_t>(projection_.at(i))), face, removed, i + 1, initial_idx);
  // }

  void thinningStep(size_t n_threads = 1);

private:
  std::array<size_t, D> projection_;
  std::array<size_t, D> dims_;
  RleBitset table_;
};

template <size_t D> VoxelComplex(std::array<size_t, D>) -> VoxelComplex<static_cast<unsigned>(D)>;

template <typename T> VoxelComplex(const std::vector<T> &) -> VoxelComplex<kGetBoolVectorDimV<std::vector<T>>>;
} // namespace stmesh

#endif