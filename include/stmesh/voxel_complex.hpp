#ifndef STMESH_VOXEL_COMPLEX_HPP
#define STMESH_VOXEL_COMPLEX_HPP

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

#include "bitset.hpp"
#include "rle_bitset.hpp"
#include "utility.hpp"

// 4h
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

    constexpr void setPos(unsigned i, int val) {
      const unsigned shift = i << 1U;
      pos &= static_cast<unsigned char>(~(3U << shift));
      pos |= static_cast<unsigned char>(static_cast<unsigned>(val + 2) << shift);
    }

    [[nodiscard]] constexpr unsigned char getSet() const {
      // NOLINTNEXTLINE(*-magic-numbers)
      return pos & 0b01010101U;
    }

    [[nodiscard]] constexpr unsigned char getUnset() const {
      // NOLINTNEXTLINE(*-magic-numbers)
      return static_cast<unsigned char>(~pos) & (0b01010101U & kDimMask);
    }

    [[nodiscard]] constexpr int getPos(unsigned i) const {
      return static_cast<int>((static_cast<unsigned>(pos) >> (i << 1U)) & 3U) - 2;
    }

    [[nodiscard]] constexpr unsigned char getDim() const {
      return static_cast<unsigned char>(std::popcount(getUnset()));
    }

    [[nodiscard]] constexpr auto operator<=>(const Face &) const = default;
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
  [[nodiscard]] consteval static std::array<unsigned char, Face::kMaxFaces> initialFaces() {
    std::array<unsigned char, Face::kMaxFaces> result{};
    iterateSubfacesRecursive({}, [&](const Face face) constexpr { result.at(face.pos) = face.getSet(); });
    return result;
  }

  template <unsigned Dim> struct DimSubfaceStore {
    std::array<size_t, Dim + 1> dims{};
    std::array<Face, expN(Dim, 3)> faces;

    [[nodiscard]] constexpr auto operator<=>(const DimSubfaceStore &) const = default;

    [[nodiscard]] constexpr std::span<const Face> asSpan(unsigned dim) const {
      if (dim == Dim)
        return std::span(faces.begin() + dims.at(Dim), faces.end());
      return std::span(faces.begin() + dims.at(dim), faces.begin() + dims.at(dim + 1));
    }
  };

  template <unsigned Dim> [[nodiscard]] consteval static DimSubfaceStore<Dim> dimFaces(const Face face) {
    std::array<std::vector<Face>, Dim + 1> result{};
    iterateSubfacesRecursive(face,
                             [&](const Face subface) constexpr { result.at(subface.getDim()).push_back(subface); });
    // flatten into array and store indices of dims
    std::array<size_t, Dim + 1> dims{};
    std::array<Face, expN(Dim, 3)> faces{};
    for (size_t i = 0; i < Dim + 1; ++i) {
      std::ranges::copy(result.at(i), faces.begin() + dims.at(i));
      if (i < Dim)
        dims.at(i + 1) = result.at(i).size() + dims.at(i);
    }
    return {dims, faces};
  }

  template <int Dim = D> [[nodiscard]] consteval static auto allSortedDims() {
    if constexpr (Dim == -1)
      return std::tuple<>();
    else {
      // NOLINTNEXTLINE(bugprone-misplaced-widening-cast)
      constexpr auto kSize = static_cast<std::size_t>((1U << static_cast<unsigned>(D - Dim)) * nChoosek(D, Dim));
      std::array<std::pair<unsigned char, unsigned short>, kSize> result_map{};
      std::array<DimSubfaceStore<static_cast<unsigned>(Dim)>, kSize> result_values{};
      size_t i = 0;
      iterateSubfacesRecursive(
          {},
          [&](const Face face) consteval {
            result_map.at(i) = {face.pos, i};
            result_values.at(i) = dimFaces<Dim>(face);
            i++;
          },
          Dim, Dim);
      return std::tuple_cat(allSortedDims<Dim - 1>(), std::make_tuple(std::pair{result_map, result_values}));
    }
  }

  constexpr static auto kAllSortedDims = allSortedDims();

  template <unsigned Dim>
  static constexpr bool essential(const Face face, const std::array<unsigned char, Face::kMaxFaces> &onto) {
    return onto.at(face.pos) >= D - Dim;
  }

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

  [[nodiscard]] constexpr unsigned checkContainingVoxels(size_t idx, const Face face, const RleBitset &removed,
                                                         unsigned i = 0);

  bool thinningStep(size_t n_threads = 1);

  void fixOneNeighbor(size_t n_threads = 1);

  [[nodiscard]] const Bitset &fixed() const;

private:
  std::array<size_t, D> projection_;
  std::array<size_t, D> dims_;
  RleBitset table_;
  Bitset fixed_;
};

template <size_t D> VoxelComplex(std::array<size_t, D>) -> VoxelComplex<static_cast<unsigned>(D)>;

template <typename T> VoxelComplex(const std::vector<T> &) -> VoxelComplex<kGetBoolVectorDimV<std::vector<T>>>;
} // namespace stmesh

#endif