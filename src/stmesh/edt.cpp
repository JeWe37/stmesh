#include "stmesh/edt.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>
#include <span>
#include <string>
#include <thread>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <itkAddImageFilter.h>
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImportImageFilter.h>
#ifndef __clang_analyzer__
#include <itkImageRegionConstIterator.h>
#endif
#include <itkIndex.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkOffset.h>
#include <itkSubtractImageFilter.h>

#include "stmesh/utility.hpp"
#include "stmesh/voxel_complex.hpp"

namespace stmesh {

template <unsigned D, bool LFS> struct EDTReader<D, LFS>::Impl {
  using Index = itk::Index<D>;

  using Image = itk::Image<unsigned char, D>;
  using FloatImage = itk::Image<FLOAT_T, D>;
  using OffsetImage = itk::Image<itk::Offset<D>, D>;

  using BinaryThresholdImageFilter = itk::BinaryThresholdImageFilter<Image, Image>;
  using DanielssonDistanceMapImageFilter = itk::DanielssonDistanceMapImageFilter<Image, FloatImage>;
  using ImportFilter = itk::ImportImageFilter<unsigned char, D>;
  using AddImageFilter = itk::AddImageFilter<OffsetImage>;
  using SubtractImageFilter = itk::SubtractImageFilter<FloatImage>;

  typename BinaryThresholdImageFilter::Pointer threshold_image_filter;
  typename DanielssonDistanceMapImageFilter::Pointer thinned_distance_map_image_filter;
  typename DanielssonDistanceMapImageFilter::Pointer distance_map_image_filter;
  typename DanielssonDistanceMapImageFilter::Pointer not_distance_map_image_filter;
  typename ImportFilter::Pointer thinned_import_filter;
  typename AddImageFilter::Pointer add_image_filter;
  typename SubtractImageFilter::Pointer subtract_image_filter;

  VoxelComplex<D> voxel_complex;

  Image::Pointer read_image;
  FloatImage::Pointer distance_map;
  FloatImage::Pointer thinned_distance_map;
  OffsetImage::Pointer vector_map;
  Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> bounding_box;

  Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)>
  boundingBoxFromImage(const itk::Image<unsigned char, D> *image) const {
    using BinaryImageToLabelMapFilter = itk::BinaryImageToLabelMapFilter<Image>;
    const typename BinaryImageToLabelMapFilter::Pointer binary_image_to_label_map_filter =
        BinaryImageToLabelMapFilter::New();
    binary_image_to_label_map_filter->SetInputForegroundValue(0);
    binary_image_to_label_map_filter->SetInput(image);
    binary_image_to_label_map_filter->Update();

    using LabelMapToLabelImageFilter =
        itk::LabelMapToLabelImageFilter<typename BinaryImageToLabelMapFilter::OutputImageType, Image>;
    const typename LabelMapToLabelImageFilter::Pointer label_map_to_label_image_filter =
        LabelMapToLabelImageFilter::New();
    label_map_to_label_image_filter->SetInput(binary_image_to_label_map_filter->GetOutput());
    label_map_to_label_image_filter->Update();

    using LabelStatisticsImageFilter = itk::LabelStatisticsImageFilter<Image, Image>;
    const typename LabelStatisticsImageFilter::Pointer label_statistics_image_filter =
        LabelStatisticsImageFilter::New();
    label_statistics_image_filter->SetLabelInput(label_map_to_label_image_filter->GetOutput());
    label_statistics_image_filter->SetInput(image);
    label_statistics_image_filter->Update();

    Index min;
    Index max;
    for (size_t i = 0; i < D; ++i) {
      min[static_cast<Index::size_type>(i)] = label_statistics_image_filter->GetBoundingBox(1)[i * 2];
      max[static_cast<Index::size_type>(i)] = label_statistics_image_filter->GetBoundingBox(1)[i * 2 + 1];
    }
    return {indexToVector(min) - 1.0 * spacing(), indexToVector(max) + 1.0 * spacing()};
  }

  template <size_t I = 0> auto vectorFromImage(const Image *image, Index index = {}) const {
    if constexpr (I == D)
      return static_cast<bool>(image->GetPixel(index));
    else {
      std::vector<decltype(vectorFromImage<I + 1>(image, index))> result;
      result.reserve(image->GetLargestPossibleRegion().GetSize()[I]);
      for (typename Index::value_type i = 0;
           i < static_cast<Index::value_type>(image->GetLargestPossibleRegion().GetSize()[I]); ++i) {
        index[I] = i;
        result.push_back(vectorFromImage<I + 1>(image, index));
      }
      return result;
    }
  }

  explicit Impl(const std::string &filename, size_t n_threads = 1)
      : threshold_image_filter(BinaryThresholdImageFilter::New()),
        thinned_distance_map_image_filter(DanielssonDistanceMapImageFilter::New()),
        distance_map_image_filter(DanielssonDistanceMapImageFilter::New()),
        not_distance_map_image_filter(DanielssonDistanceMapImageFilter::New()),
        thinned_import_filter(ImportFilter::New()), add_image_filter(AddImageFilter::New()),
        subtract_image_filter(SubtractImageFilter::New()), read_image(itk::ReadImage<Image>(filename)),
        distance_map(subtract_image_filter->GetOutput()),
        thinned_distance_map(thinned_distance_map_image_filter->GetOutput()),
        vector_map(add_image_filter->GetOutput()) {
    const stmesh::Vector4F origin = spacing() * 0.5;
    read_image->SetOrigin(origin.data());

    threshold_image_filter->SetLowerThreshold(1);
    threshold_image_filter->SetInsideValue(0);
    threshold_image_filter->SetOutsideValue(1);
    threshold_image_filter->SetInput(read_image);
    threshold_image_filter->ReleaseDataFlagOn();

    // NOLINTNEXTLINE(misc-const-correctness)
    std::thread th;
    if constexpr (LFS) {
      voxel_complex = VoxelComplex<D>(vectorFromImage(read_image));

      th = std::thread([&]() {
        do {
          voxel_complex.fixOneNeighbor(n_threads);
        } while (voxel_complex.thinningStep(n_threads));
      });
    }

    distance_map_image_filter->SetInput(read_image);
    distance_map_image_filter->ReleaseDataFlagOn();

    not_distance_map_image_filter->SetInput(threshold_image_filter->GetOutput());
    not_distance_map_image_filter->ReleaseDataFlagOn();

    add_image_filter->SetInput1(not_distance_map_image_filter->GetVectorDistanceMap());
    add_image_filter->SetInput2(distance_map_image_filter->GetVectorDistanceMap());
    add_image_filter->ReleaseDataFlagOn();
    add_image_filter->Update();

    subtract_image_filter->SetInput1(distance_map_image_filter->GetDistanceMap());
    subtract_image_filter->SetInput2(not_distance_map_image_filter->GetDistanceMap());
    subtract_image_filter->ReleaseDataFlagOn();
    subtract_image_filter->Update();

    // NOLINTNEXTLINE(cppcoreguidelines-prefer-member-initializer)
    bounding_box = boundingBoxFromImage(threshold_image_filter->GetOutput());

    if constexpr (LFS) {
      th.join();

      const size_t total_size = read_image->GetLargestPossibleRegion().GetNumberOfPixels();
      // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
      auto *thinned_data = new unsigned char[total_size]{};
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      voxel_complex.table().iterateSet([&](size_t idx, size_t /*unused*/) { thinned_data[idx] = 1; }, {}, {},
                                       n_threads);

      thinned_import_filter->SetOrigin(read_image->GetOrigin());
      thinned_import_filter->SetSpacing(read_image->GetSpacing());
      thinned_import_filter->SetRegion(read_image->GetLargestPossibleRegion());
      thinned_import_filter->SetImportPointer(thinned_data, total_size, true);

      thinned_distance_map_image_filter->SetInput(thinned_import_filter->GetOutput());
      thinned_distance_map_image_filter->ReleaseDataFlagOn();
      thinned_distance_map_image_filter->Update();
    }
  }

  [[nodiscard]] Index vectorToIndex(const Vector4F &vec) const noexcept {
    using Point = Image::PointType;
    using Size = Image::SizeType;

    Size size = read_image->GetLargestPossibleRegion().GetSize();

    const Point point(vec.data());
    Index result = read_image->TransformPhysicalPointToIndex(point);
    for (size_t i = 0; i < D; ++i)
      result[static_cast<Index::size_type>(i)] =
          std::clamp(result[static_cast<Index::size_type>(i)], typename Index::IndexValueType(0),
                     typename Index::IndexValueType(size[static_cast<Size::size_type>(i)] - 1));
    return result;
  }

  [[nodiscard]] VectorF<D> indexToVector(const Index &index) const noexcept {
    using Point = Image::PointType;

    Point point;
    read_image->TransformIndexToPhysicalPoint(index, point);
    Vector4F result;
    std::ranges::copy(point, result.begin());
    return result;
  }

  [[nodiscard]] VectorF<D> spacing() const noexcept {
    using Spacing = Image::SpacingType;

    Spacing spacing = read_image->GetSpacing();
    Vector4F result;
    std::ranges::copy(spacing, result.begin());
    return result;
  }

  template <size_t N = 0>
  void getMultiIndex(Index min_corner, const std::span<const size_t, D> &size, auto &&it, auto image) const {
    if constexpr (N == D)
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      *it++ = image->GetPixel(min_corner);
    else {
      for (size_t i = 0; i < size[N]; ++i) {
        getMultiIndex<N + 1>(min_corner, size, it, image);
        if (min_corner[N] < static_cast<Index::size_type>(read_image->GetLargestPossibleRegion().GetSize()[N] - 1))
          min_corner[N]++;
      }
    }
  }
};

template <unsigned D, bool LFS>
EDTReader<D, LFS>::EDTReader(const std::string &filename, size_t n_threads) : pimpl_(new Impl(filename, n_threads)) {}

template <unsigned D, bool LFS> EDTReader<D, LFS>::~EDTReader() = default;

template <unsigned D, bool LFS> Vector4F EDTReader<D, LFS>::spacing() const noexcept { return pimpl_->spacing(); }

template <unsigned D, bool LFS> FLOAT_T EDTReader<D, LFS>::signedDistance(const VectorF<D> &point) const noexcept {
  return pimpl_->distance_map->GetPixel(pimpl_->vectorToIndex(point));
}

template <unsigned D, bool LFS>
[[nodiscard]] std::vector<FLOAT_T>
EDTReader<D, LFS>::signedDistance(const VectorF<D> &point, const std::span<const size_t, D> &size) const noexcept {
  const typename Impl::Index min_corner = pimpl_->vectorToIndex(point);
  std::vector<FLOAT_T> result;
  result.reserve(std::accumulate(size.begin(), size.end(), size_t{1}, std::multiplies<>()));
  pimpl_->getMultiIndex(min_corner, size, std::back_insert_iterator(result), pimpl_->distance_map);
  return result;
}

template <unsigned D, bool LFS> size_t EDTReader<D, LFS>::pointBoundaryRegion(const Vector4F &point) const noexcept {
  const typename Impl::Index index = pimpl_->vectorToIndex(point - 0.5 * spacing());
  std::array<unsigned char, 1U << 4U> distances{};
  std::array<size_t, D> size{};
  std::ranges::fill(size, 2U);
  pimpl_->getMultiIndex(index, size, distances.begin(), pimpl_->read_image);
  return *std::ranges::max_element(distances);
}

template <unsigned D, bool LFS> VectorF<D> EDTReader<D, LFS>::closestAt(const VectorF<D> &point) const noexcept {
  typename Impl::Index index = pimpl_->vectorToIndex(point);
  index += pimpl_->vector_map->GetPixel(index);
  return pimpl_->indexToVector(index);
}

template <unsigned D, bool LFS>
FLOAT_T EDTReader<D, LFS>::distanceToThinnedAt(const VectorF<D> &point) const noexcept
requires LFS
{
  return pimpl_->thinned_distance_map->GetPixel(pimpl_->vectorToIndex(point));
}

template <unsigned D, bool LFS>
Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> EDTReader<D, LFS>::boundingBox() const noexcept {
  return pimpl_->bounding_box;
}

template class EDTReader<4, false>;
template class EDTReader<4, true>;
} // namespace stmesh
