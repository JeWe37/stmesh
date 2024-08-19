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
/// Has a member with a D-dimensional vector of bools
/**
 * Computes the type of a D-dimensional vector of bools via recursion.
 */
template <unsigned D> struct BoolVectorDim {
  using type = std::vector<typename BoolVectorDim<D - 1>::type>;
};

template <> struct BoolVectorDim<1> {
  using type = std::vector<bool>;
};

template <unsigned D> using BoolVectorDimType = BoolVectorDim<D>::type; ///< The type of a D-dimensional vector of bools

/// Gets the dimension of a D-dimensional vector of bools
/**
 * Computes the dimension of a D-dimensional vector of bools via recursion.
 */
template <typename T> struct GetBoolVectorDim {
  static constexpr unsigned kValue = GetBoolVectorDim<typename T::value_type>::kValue + 1;
};

template <> struct GetBoolVectorDim<std::vector<bool>> {
  static constexpr unsigned kValue = 1;
};

template <typename T>
constexpr unsigned kGetBoolVectorDimV =
    GetBoolVectorDim<T>::kValue; ///< The dimension of a D-dimensional vector of bools

/// A D-dimensional voxel complex
/**
 * This class provides methods for the storage and thinning of a D-dimensional voxel complex. The voxel complex is
 * stored as a D-dimensional RleBitset. This format permits efficient thinning via critical thinning.
 *
 * @tparam D The dimension(<= 4) of the voxel complex
 */
template <unsigned D> class VoxelComplex {
  static_assert(D <= 4, "VoxelComplex only supports dimensions up to 4");

public:
  /// A face of a D-dimensional voxel
  /**
   * A face is represented as the coefficients of the vectors spanning its (oriented) dual space. This data is
   * compressed into a byte of 4 2-bit values, each representing -1, 0 or 1. -1 becomes 1(0b01), 0 becomes 2(0b10),
   * 1 becomes 3(0b11). This permits the storage of a face in a single byte and the efficient computation of
   * several important properties.
   */
  struct Face {
    constexpr static unsigned char kDimMask = (1U << (D << 1U)) - 1U; ///< A mask for the dimension bits
    constexpr static size_t kMaxFaces = std::numeric_limits<unsigned char>::max() + 1; ///< The maximum number of faces

    // NOLINTNEXTLINE(*-magic-numbers)
    unsigned char pos{0b10101010U & kDimMask}; ///< The face stored in the encoded form

    /// Sets the i-th coordinate of the face
    /**
     * Can be set to -1, 0 or 1.
     *
     * @param i The index of the coordinate
     * @param val The value of the coordinate
     */
    constexpr void setPos(unsigned i, int val) {
      const unsigned shift = i << 1U;
      pos &= static_cast<unsigned char>(~(3U << shift));
      pos |= static_cast<unsigned char>(static_cast<unsigned>(val + 2) << shift);
    }

    /// Gets the set coordinates of the face
    /**
     * The set coordinates are the coordinates that are not zero.
     *
     * @return The set coordinates
     */
    [[nodiscard]] constexpr unsigned char getSet() const {
      // NOLINTNEXTLINE(*-magic-numbers)
      return pos & 0b01010101U;
    }

    /// Gets the unset coordinates of the face
    /**
     * The unset coordinates are the coordinates that are zero.
     *
     * @return The unset coordinates
     */
    [[nodiscard]] constexpr unsigned char getUnset() const {
      // NOLINTNEXTLINE(*-magic-numbers)
      return static_cast<unsigned char>(~pos) & (0b01010101U & kDimMask);
    }

    /// Gets the i-th coordinate of the face
    /**
     * Can be -1, 0 or 1.
     *
     * @param i The index of the coordinate
     * @return The i-th coordinate
     */
    [[nodiscard]] constexpr int getPos(unsigned i) const {
      return static_cast<int>((static_cast<unsigned>(pos) >> (i << 1U)) & 3U) - 2;
    }

    /// Gets the dimension of the face
    /**
     * The dimension is the number of unset coordinates.
     *
     * @return The dimension
     */
    [[nodiscard]] constexpr unsigned char getDim() const {
      return static_cast<unsigned char>(std::popcount(getUnset()));
    }

    /// Allows comparison of faces
    /**
     * Faces are compared by their encoded form.
     *
     * @param other The other face
     * @return The comparison result
     */
    [[nodiscard]] constexpr auto operator<=>(const Face &) const = default;
  };

private:
  template <typename T> [[nodiscard]] static auto getDims(const T &data);

  template <typename T> static void flatten(const T &data, std::vector<bool> &result);

public:
  /// Constructs a D-dimensional voxel complex
  /**
   * The voxel complex is initialized with the given data.
   *
   * @param data The data of the voxel complex
   */
  explicit VoxelComplex(const BoolVectorDimType<D> &data);

  /// Constructs a D-dimensional voxel complex
  /**
   * The voxel complex is created with the given sizes and initialized to being empty.
   *
   * @param dims The dimensions of the voxel complex
   */
  explicit VoxelComplex(const std::array<size_t, D> &dims);

  /// Gets the neighbor of a voxel given by an index
  /**
   * The neighbor is given by the index of the voxel and the face of the neighbor. The neighbor selected is the neighbor
   * lying across the face.
   *
   * @param idx The index of the voxel
   * @param face The face shared with the neighbor
   */
  [[nodiscard]] size_t offsetByFace(size_t idx, const Face face) const;

  /// Gets the index of a face based on its coordinates
  /**
   * The index is the index of the face in the voxel complex. Uses the projection to map the coordinates to the index.
   *
   * @param coords The coordinates of the face
   * @return The index of the face
   */
  [[nodiscard]] size_t index(const std::array<size_t, D> &coords) const;

  /// Gets the value of a voxel given by a position
  /**
   * The value is true if the voxel is set and false if the voxel is unset. The index is computed from the position.
   *
   * @param idx The index of the voxel
   * @return The value of the voxel
   */
  [[nodiscard]] bool operator[](const std::array<size_t, D> &coords) const;

  /// Recursively iterates over all subfaces of a face
  /**
   * The function is called for each subface of the face. The function is called with the subface as an argument.
   * Using d and func_dim, the iteration can be limited to faces of a certain dimension. By default, the function
   * is called for all subfaces.
   *
   * @tparam F The type of the function
   * @param face The face
   * @param func The function
   * @param d The dimension to iterate down to
   * @param func_dim The dimension to run the function for
   * @param i The current coordinate
   */
  template <std::invocable<const Face> F>
  static consteval void iterateSubfacesRecursive(Face face, const F &func, unsigned d = 0, int func_dim = -1,
                                                 unsigned i = 0) {
    // potentially it might be faster to use the precomputed dim faces, but really this can be optimized out fully
    // usually
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
  /// Gets the core of a voxel
  /**
   * The core is represented by an array of the number of unique coordinate directions in which there are other voxels,
   * that lie within the clique of a given face. The array can be indexed using the pos of the face.
   *
   * @param idx The index of the voxel
   */
  [[nodiscard]] std::array<unsigned char, Face::kMaxFaces> core(const size_t idx) const;

  /// Initial superface array
  /**
   * Gets the initial superface array. This array indicates which positions, if reset, would lead to a superface  of the
   * current facethat is still present.
   */
  [[nodiscard]] consteval static std::array<unsigned char, Face::kMaxFaces> initialFaces() {
    std::array<unsigned char, Face::kMaxFaces> result{};
    iterateSubfacesRecursive({}, [&](const Face face) constexpr { result.at(face.pos) = face.getSet(); });
    return result;
  }

  /// Store for the subfaces of a face
  /**
   * Stores the subfaces of a face in a struct. Stores faces sorted by dimension.
   */
  template <unsigned Dim> struct DimSubfaceStore {
    std::array<size_t, Dim + 1> dims{};   ///< Boundaries between the dimensions in faces
    std::array<Face, expN(Dim, 3)> faces; ///< Faces sorted by dimension

    /// Compares two DimSubfaceStores
    /**
     * Compares two DimSubfaceStores by their dims and faces.
     *
     * @param other The other DimSubfaceStore
     * @return The comparison result
     */
    [[nodiscard]] constexpr auto operator<=>(const DimSubfaceStore &) const = default;

    /// Gets a span of the faces of a certain dimension
    /**
     * Gets a span of the faces of a certain dimension. The span is a view into the faces array.
     *
     * @param dim The dimension of the faces
     * @return The span of the faces
     */
    [[nodiscard]] constexpr std::span<const Face> asSpan(unsigned dim) const {
      if (dim == Dim)
        return std::span(faces.begin() + dims.at(Dim), faces.end());
      return std::span(faces.begin() + dims.at(dim), faces.begin() + dims.at(dim + 1));
    }
  };

  /// Gets the subfaces of a face
  /**
   * Gets the subfaces of a face in a DimSubfaceStore. The faces are sorted by dimension.
   *
   * @param face The face
   * @return The DimSubfaceStore of the subfaces
   */
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

  /// Gets the subfaces of the faces for all dimensions
  /**
   * Gets the subfaces of the faces for all dimensions. The faces are stored in a tuple for each dimension. In each,
   * They contain a map that can be used with mph to determin the index of a face inside of the faces array, containing
   * DimSubfaceStores.
   *
   * @return The tuple of pairs of mph maps and arrays of DimSubfaceStores
   */
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

  constexpr static auto kAllSortedDims = allSortedDims(); ///< Result of allSortedDims

  /// Checks whether a face is essential
  /**
   * Whether a face is essential depends on the dimension we are collapsing in.
   *
   * @tparam Dim The dimension we are collapsing in
   * @param face The face to check
   * @param onto The onto array from the core
   * @return Whether the face is essential
   */
  template <unsigned Dim>
  static constexpr bool essential(const Face face, const std::array<unsigned char, Face::kMaxFaces> &onto) {
    return onto.at(face.pos) >= D - Dim;
  }

  /// Checks for collapse onto a subcomplex
  /**
   * Checks whether the voxel complex formed by the faces of face can be collapsed onto its subcomplex given
   * by the essential faces of the onto subcomplex.
   *
   * @tparam Dim The dimension we are collapsing in
   * @param face The face to check
   * @param onto The onto array from the core
   * @return Whether the face can be collapsed
   */
  template <unsigned Dim>
  static bool collapseFace(const Face face, const std::array<unsigned char, Face::kMaxFaces> &onto);

  /// Checks whether we can immediately keep a voxel
  /**
   * Checks whether we can immediately keep a voxel. If a voxel has other voxels in all directions, it is trivially
   * critical and can be kept.
   *
   * @param idx The index of the voxel
   * @return Whether the voxel can be kept
   */
  [[nodiscard]] bool checkShortcut(size_t idx) const;

  /// Iterates over all voxels
  /**
   * Iterates over all voxels and calls the function for each voxel.
   *
   * @tparam F The type of the function
   * @param func The function
   * @param coords The coordinates of the voxel, not needed for external calls
   */
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

  /// Counts the number of set voxels in a clique
  /**
   * Counts the number of set voxels that are set in removed in a clique. The clique is given by the face.
   *
   * @param idx The index of the voxel
   * @param face The face of the clique
   * @param removed The voxel complex to count in
   * @return The number of set voxels
   */
  [[nodiscard]] constexpr unsigned checkContainingVoxels(size_t idx, const Face face, const RleBitset &removed,
                                                         unsigned i = 0);

  /// Performs one step of thinning
  /**
   * Performs one step of thinning. The function returns true if the voxel complex was changed.
   *
   * To do the general thinning, critical thinning is used.
   *
   * @param n_threads The number of threads to use
   * @return Whether the voxel complex was changed
   */
  [[nodiscard]] bool thinningStep(size_t n_threads = 1);

  /// Fixes voxels with exactly one neighbor
  /**
   * Fixes voxels with exactly one neighbor. This is useful to keep additional branches of the skeleton. Such voxels
   * will not be removed by the thinning.
   *
   * @param n_threads The number of threads to use
   */
  void fixOneNeighbor(size_t n_threads = 1);

  /// Gets the set of fixed voxels
  /**
   * Gets the set of fixed voxels. These voxels will not be removed by the thinning.
   *
   * @return The set of fixed voxels
   */
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