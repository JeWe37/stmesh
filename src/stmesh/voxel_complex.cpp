#include "stmesh/voxel_complex.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <mph>

#include "stmesh/rle_bitset.hpp"
#include "stmesh/utility.hpp"

namespace stmesh {
template <unsigned D> constexpr void VoxelComplex<D>::Face::setPos(unsigned i, int val) {
  const unsigned shift = i << 1U;
  pos &= static_cast<unsigned char>(~(3U << shift));
  pos |= static_cast<unsigned char>(static_cast<unsigned>(val + 2) << shift);
}

template <unsigned D> [[nodiscard]] constexpr unsigned char VoxelComplex<D>::Face::getSet() const {
  // NOLINTNEXTLINE(*-magic-numbers)
  return pos & 0b01010101U;
}

template <unsigned D> [[nodiscard]] constexpr unsigned char VoxelComplex<D>::Face::getUnset() const {
  // NOLINTNEXTLINE(*-magic-numbers)
  return static_cast<unsigned char>(~pos) & (0b01010101U & kDimMask);
}

template <unsigned D> [[nodiscard]] constexpr int VoxelComplex<D>::Face::getPos(unsigned i) const {
  return static_cast<int>((static_cast<unsigned>(pos) >> (i << 1U)) & 3U) - 2;
}

template <unsigned D> [[nodiscard]] constexpr unsigned char VoxelComplex<D>::Face::getDim() const {
  return static_cast<unsigned char>(std::popcount(getUnset()));
}

template <unsigned D> template <typename T> [[nodiscard]] auto VoxelComplex<D>::getDims(const T &data) {
  static constexpr unsigned kDim = kGetBoolVectorDimV<T>;
  if constexpr (kDim == 1)
    return std::array{data.size()};
  else {
    std::array<size_t, kDim> result{};
    std::ranges::copy(getDims(data[0]), result.begin());
    result[kDim - 1] = data.size();
    return result;
  }
}

template <unsigned D> template <typename T> void VoxelComplex<D>::flatten(const T &data, std::vector<bool> &result) {
  static constexpr unsigned kDim = kGetBoolVectorDimV<T>;
  for (const auto &vec : data) {
    if constexpr (kDim == 2)
      std::ranges::copy(vec, std::back_inserter(result));
    else
      flatten(vec, result);
  }
}

template <unsigned D> VoxelComplex<D>::VoxelComplex(const BoolVectorDimType<D> &data) : VoxelComplex(getDims(data)) {
  std::vector<bool> flattend;
  flatten(data, flattend);
  table_.setFrom(flattend);
}

template <unsigned D>
VoxelComplex<D>::VoxelComplex(const std::array<size_t, D> &dims)
    : dims_(dims), table_(std::accumulate(dims.begin(), dims.end(), size_t{1}, std::multiplies<>())) {
  projection_[0] = 1;
  std::partial_sum(dims.begin(), dims.begin() + (D - 1), projection_.begin() + 1, std::multiplies<>());
}

template <unsigned D> [[nodiscard]] size_t VoxelComplex<D>::offsetByFace(size_t idx, const Face face) const {
  for (unsigned i = 0; i < D; ++i)
    idx += static_cast<size_t>(face.getPos(i)) * projection_.at(i);
  return idx;
}

template <unsigned D> [[nodiscard]] size_t VoxelComplex<D>::index(const std::array<size_t, D> &coords) const {
  size_t index = 0;
  for (size_t i = 0; i < D; ++i)
    index += coords.at(i) * projection_.at(i);
  return index;
}

template <unsigned D> [[nodiscard]] bool VoxelComplex<D>::operator[](const std::array<size_t, D> &coords) const {
  return table_[index(coords)];
}

template <unsigned D>
template <unsigned Dim, typename F, int TargetDim>
void VoxelComplex<D>::iterateSubfaces(Face face, const F &func) {
  static constexpr auto kDimSortedDims = std::get<Dim>(kAllSortedDims);
#ifdef _CLANGD // the mph call is atrociously slow
  auto &dim_subface_store = kDimSortedDims.second.at(0);
#else
  auto &dim_subface_store = kDimSortedDims.second.at(mph::lookup<kDimSortedDims.first>(face.pos));
#endif
  if constexpr (TargetDim == -1)
    [&]<unsigned... Vals>(std::integer_sequence<unsigned, Vals...>) {
      (
          [&] {
            // NOLINTNEXTLINE(clang-analyzer-core.CallAndMessage)
            for (const Face &subface : dim_subface_store.asSpan(Vals))
              func(subface, std::integral_constant<unsigned, Vals>{});
          }(),
          ...);
    }(std::make_integer_sequence<unsigned, Dim>());
  else if constexpr (TargetDim == Dim)
    func(face);
  else {
    // NOLINTNEXTLINE(clang-analyzer-core.CallAndMessage)
    for (const Face &subface : dim_subface_store.asSpan(static_cast<unsigned>(TargetDim)))
      func(subface);
  }
}

template <unsigned D>
[[nodiscard]] auto VoxelComplex<D>::core(const size_t idx) const -> std::array<unsigned char, Face::kMaxFaces> {
  std::array<unsigned char, Face::kMaxFaces> result{};
  iterateSubfaces<D>({}, [&]<unsigned SubDim>(const Face &face, std::integral_constant<unsigned, SubDim>) {
    // offsetByFace is unnecessraily expensive
    if (table_[offsetByFace(idx, face)]) {
      const unsigned char set = face.getSet();
      result.at(face.pos) |= set;
      iterateSubfaces<SubDim>(face, [&](const Face &subface, auto /*unused*/) { result.at(subface.pos) |= set; });
    }
  });
  std::ranges::transform(result, result.begin(), [](unsigned char val) { return std::popcount(val); });
  // the core of a (N-d)-face consists of all faces contained in it with result >= d
  // an (N-d)-face is essential if it has result > D - d
  return result;
}

// for each face, we store whether by setting the i-th coordinate to zero we reach an existing face
// initially, this is always the case
template <unsigned D>
[[nodiscard]] consteval auto VoxelComplex<D>::initialFaces() -> std::array<unsigned char, Face::kMaxFaces> {
  std::array<unsigned char, Face::kMaxFaces> result{};
  iterateSubfacesRecursive({}, [&](const Face face) constexpr { result.at(face.pos) = face.getSet(); });
  return result;
}

template <unsigned D>
template <unsigned Dim>
[[nodiscard]] constexpr auto VoxelComplex<D>::DimSubfaceStore<Dim>::asSpan(unsigned dim) const
    -> std::span<const Face> {
  if (dim == Dim)
    return std::span(faces.begin() + dims.at(Dim), faces.end());
  return std::span(faces.begin() + dims.at(dim), faces.begin() + dims.at(dim + 1));
}

template <unsigned D>
template <unsigned Dim>
[[nodiscard]] consteval auto VoxelComplex<D>::dimFaces(const Face face) -> DimSubfaceStore<Dim> {
  std::array<std::vector<Face>, Dim + 1> result{};
  iterateSubfacesRecursive(face, [&](const Face subface) constexpr { result.at(subface.getDim()).push_back(subface); });
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

template <unsigned D> template <int Dim> [[nodiscard]] consteval auto VoxelComplex<D>::allSortedDims() {
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

template <unsigned D>
template <unsigned Dim>
constexpr bool VoxelComplex<D>::essential(const Face face, const std::array<unsigned char, Face::kMaxFaces> &onto) {
  return onto.at(face.pos) >= D - Dim;
}

template <unsigned D>
template <unsigned Dim>
bool VoxelComplex<D>::collapseFace(const Face face, const std::array<unsigned char, Face::kMaxFaces> &onto) {
  constexpr unsigned kRemoved = 0xFF;

  static constexpr auto kDimSortedDims = std::get<Dim>(kAllSortedDims);
#ifdef _CLANGD // the mph call is atrociously slow
  auto &dim_subface_store = kDimSortedDims.second.at(0);
#else
  auto &dim_subface_store = kDimSortedDims.second.at(mph::lookup<kDimSortedDims.first>(face.pos));
#endif
  auto faces = initialFaces();

  // zero those directions that are already in face
  // kinda a mess to do this every time
  if constexpr (Dim != D) {
    unsigned char mask = ~face.getSet();
    faces.at(face.pos) &= mask;
    iterateSubfaces<Dim>(face, [&](const Face subface, auto /*unused*/) { faces.at(subface.pos) &= mask; });
  }

  auto remove_face = [&](Face subface) {
    // remove faces by setting them to kRemoved
    faces.at(subface.pos) = kRemoved;
    // subfaces of removed faces are now no longer contained in thenit as it has been removed
    for (unsigned i = 0; i < D; ++i) {
      // this is an annoying condition, maybe we should just use the lookup table and do this basically statically
      if (subface.getPos(i) == 0) {
        Face face_cpy = subface;
        face_cpy.setPos(i, 1);
        faces.at(face_cpy.pos) &= static_cast<unsigned char>(~(1U << (i << 1U)));
        face_cpy.setPos(i, -1);
        faces.at(face_cpy.pos) &= static_cast<unsigned char>(~(1U << (i << 1U)));
      }
    }
  };

  // assumes popcount is 1, removes both this face and the face containing it
  auto try_collapse = [&](Face subface) {
    // isn't essential, can be collapsed
    // non-essential implies unique containing face is not essential
    if (essential<Dim - 1>(subface, onto) || std::popcount(faces.at(subface.pos)) != 1)
      return false;
    // get the index of the coordinate that the parent to be removed does not have set
    const unsigned i = static_cast<unsigned>(std::countr_zero(faces.at(subface.pos))) >> 1U;
    Face superface = subface;
    superface.setPos(i, 0);
    // for this to be removable, the superface must not be contained in any other faces
    // this also implcitly means it must not have been removed before
    if (faces.at(superface.pos) != 0U)
      return false;
    remove_face(subface);
    remove_face(superface);
    return true;
  };

  // remain means: not essential and not removed aka has the upper bits set in faces
  auto remaining = [&](const Face subface) {
    return !essential<Dim - 1>(subface, onto) && faces.at(subface.pos) < 1U << 4U;
  };

  // collapse dimension by dimension, breadth first, for as long as a dimension still contains non-essential faces.
  // as we always collapse faces of dims d-1 and d and we detect collapsibility based on d-1 faces, start with
  // dimension n-1 if we iterate over a dimension, it still contains non-essential faces and we didn't collapse
  // anything, we fail the collapse if we never failed and reached dim 0, we succeed

  // in d == Dim - 1 exactly one collapse must happen
  // count number of essential faces that were not removed per dimension.
  // this must be equal to the number of collapses in the next dimension
  // at the end, the essential faces remaining number must be 0
  unsigned prior_dim_remaining = 1;
  for (int d = static_cast<int>(Dim) - 1; d >= 0; --d) {
    bool collapsed{};
    unsigned current_dim_remaining = 0;
    do {
      current_dim_remaining = 0;
      collapsed = false;
      // not totally optimal, removed faces get retried even if efficiently skipped
      for (const Face dim_face : dim_subface_store.asSpan(static_cast<unsigned>(d))) {
        if (try_collapse(dim_face)) {
          prior_dim_remaining--;
          collapsed = true;
        } else if (remaining(dim_face))
          ++current_dim_remaining;
      }
    } while (collapsed); // at least one face has been collapsed this iteration
    // we can fail if there remain any faces in the dimension above the current one
    if (static_cast<bool>(prior_dim_remaining))
      return false;
    prior_dim_remaining = current_dim_remaining;
  }
  // no faces should remain in the final dimension
  return !prior_dim_remaining;
}

template <unsigned D> [[nodiscard]] bool VoxelComplex<D>::checkShortcut(size_t idx) const {
  for (size_t i = 0; i < D; ++i) {
    if (!table_[idx + projection_.at(i)] || !table_[idx - projection_.at(i)])
      return false;
  }
  return true;
}

template <unsigned D>
[[nodiscard]] constexpr bool VoxelComplex<D>::checkContainingVoxels(size_t idx, const Face face,
                                                                    const RleBitset &removed, unsigned i) {
  if (i == D) // leaf
    return removed[idx];
  // should see if we need to manually optimize face.getPos(i) == 0 case
  // check the case of using this for the change and not using it
  return checkContainingVoxels(static_cast<size_t>(static_cast<ptrdiff_t>(idx) +
                                                   face.getPos(i) * static_cast<ptrdiff_t>(projection_.at(i))),
                               face, removed, i + 1) &&
         checkContainingVoxels(idx, face, removed, i + 1);
}

template <unsigned D> void VoxelComplex<D>::thinningStep(size_t n_threads) {
  RleBitset new_table(table_.size());
  new_table.registerThreads(n_threads);
  table_.iterateSet(
      [&](size_t idx, size_t thread_id) {
        if (checkShortcut(idx))
          new_table.set(idx, thread_id);
      },
      [&](size_t /*unused*/) { /* just registered the threads, already at start */ },
      [&](size_t thread_id) { new_table.finalizeThread(thread_id); }, n_threads);

  // it might not look it, but this is actually a loop(at compile time, for kD=D..0)
  [&]<unsigned... Vals>(std::integer_sequence<unsigned, Vals...>) {
    (
        [&] {
          constexpr unsigned kD = D - Vals;

          RleBitset removed(table_.size());
          removed.registerThreads(n_threads);
          table_.iterateSet(
              [&](size_t idx, size_t thread_id) {
                if (!new_table[idx])
                  removed.set(idx, thread_id);
              },
              [&](size_t /*unused*/) { /* just registered the threads, already at start */ },
              [&](size_t thread_id) { removed.finalizeThread(thread_id); }, n_threads);
          removed.unregisterThreads();

          removed.iterateSet(
              [&](size_t idx, size_t thread_id) constexpr {
                // issue: we recompute the core a lot which might cost performance
                auto onto = core(idx);
                auto voxel_test = [&](const Face &face) {
                  if (essential<kD>(face, onto) && checkContainingVoxels(idx, face, removed) &&
                      !collapseFace<kD>(face, onto)) // essential in this dim + not collapsible == critical
                    new_table.set(idx, thread_id);
                };
                iterateSubfaces<D, decltype(voxel_test), kD>({}, voxel_test);
              },
              [&](size_t thread_id) { new_table.restartIteration(thread_id); },
              [&](size_t thread_id) { new_table.finalizeThread(thread_id); }, n_threads);
        }(),
        ...);
  }(std::make_integer_sequence<unsigned, D + 1>());
  new_table.unregisterThreads();
  table_ = std::move(new_table);
}

template class VoxelComplex<2U>;
template class VoxelComplex<3U>;
template class VoxelComplex<4U>;
} // namespace stmesh
