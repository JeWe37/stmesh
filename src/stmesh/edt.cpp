// NOLINTBEGIN(misc-include-cleaner)
#include <stmesh/edt.hpp>

#include <functional>
#include <numeric>
#include <type_traits>

#include <itkAddImageFilter.h>
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkBinaryNotImageFilter.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkOffset.h>
#include <itkSmartPointer.h>
#include <itkSubtractImageFilter.h>

namespace stmesh {

namespace {
template <typename TImage, typename F>
std::vector<std::remove_cvref_t<std::invoke_result_t<F, typename TImage::PixelType>>>
deepCopy(TImage *input, const F &conversion) noexcept {
  std::vector<std::remove_cvref_t<std::invoke_result_t<F, typename TImage::PixelType>>> output;

#ifndef __clang_analyzer__
  for (itk::ImageRegionConstIterator<TImage> input_iterator(input, input->GetLargestPossibleRegion());
       !input_iterator.IsAtEnd(); ++input_iterator)
    output.push_back(conversion(input_iterator.Get()));
#else
  (void)input;
  (void)conversion;
#endif

  return output;
}

template <unsigned D>
Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> boundingBoxFromImage(const itk::Image<unsigned char, D> *image) {
  using Image = itk::Image<unsigned char, D>;
  using BinaryImageToLabelMapFilter = itk::BinaryImageToLabelMapFilter<Image>;
  const typename BinaryImageToLabelMapFilter::Pointer binary_image_to_label_map_filter =
      BinaryImageToLabelMapFilter::New();
  binary_image_to_label_map_filter->SetInputForegroundValue(1);
  binary_image_to_label_map_filter->SetInput(image);
  binary_image_to_label_map_filter->Update();

  using LabelMapToLabelImageFilter =
      itk::LabelMapToLabelImageFilter<typename BinaryImageToLabelMapFilter::OutputImageType, Image>;
  const typename LabelMapToLabelImageFilter::Pointer label_map_to_label_image_filter =
      LabelMapToLabelImageFilter::New();
  label_map_to_label_image_filter->SetInput(binary_image_to_label_map_filter->GetOutput());
  label_map_to_label_image_filter->Update();

  using LabelStatisticsImageFilter = itk::LabelStatisticsImageFilter<Image, Image>;
  const typename LabelStatisticsImageFilter::Pointer label_statistics_image_filter = LabelStatisticsImageFilter::New();
  label_statistics_image_filter->SetLabelInput(label_map_to_label_image_filter->GetOutput());
  label_statistics_image_filter->SetInput(image);
  label_statistics_image_filter->Update();

  Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> result;
  for (size_t i = 0; i < D; ++i) {
    result.min()[static_cast<Eigen::Index>(i)] = FLOAT_T(label_statistics_image_filter->GetBoundingBox(1)[i * 2] - 1);
    result.max()[static_cast<Eigen::Index>(i)] =
        FLOAT_T(label_statistics_image_filter->GetBoundingBox(1)[i * 2 + 1] + 1);
  }
  return result;
}
} // namespace

template <unsigned D> Eigen::Vector<size_t, 4> EDTReader<D>::clamp(const Vector4F &point) const noexcept {
  return point.cwiseMax(Vector4F::Zero())
      .template cast<size_t>()
      .cwiseMin(Eigen::Map<const Eigen::Vector<size_t, 4>>(sizes_.data()) -
                Eigen::Vector<size_t, 4>::Constant(size_t(1)));
}

template <unsigned D> size_t EDTReader<D>::projection(const Vector4F &point) const noexcept {
  Eigen::Vector<size_t, static_cast<int>(D)> proj;
  proj[0] = size_t{1};
  std::partial_sum(sizes_.begin(), sizes_.end() - 1, proj.begin() + 1, std::multiplies<>());
  return clamp(point).dot(proj);
}

template <unsigned D> EDTReader<D>::EDTReader(const std::string &filename) {
  using Image = itk::Image<unsigned char, D>;
  using FloatImage = itk::Image<FLOAT_T, D>;
  using OffsetImage = itk::Image<itk::Offset<D>, D>;

  auto read = itk::ReadImage<Image>(filename);
  std::ranges::copy(read->GetLargestPossibleRegion().GetSize(), sizes_.begin());
  bounding_box_ = boundingBoxFromImage(read.GetPointer());

  using DanielssonDistanceMapImageFilter = itk::DanielssonDistanceMapImageFilter<Image, FloatImage>;
  const typename DanielssonDistanceMapImageFilter::Pointer distance_map_image_filter =
      DanielssonDistanceMapImageFilter::New();
  distance_map_image_filter->SetInput(read);

  using BinaryNotImageFilter = itk::BinaryNotImageFilter<Image>;
  const typename BinaryNotImageFilter::Pointer not_image_filter = BinaryNotImageFilter::New();
  not_image_filter->SetBackgroundValue(0);
  not_image_filter->SetForegroundValue(1);
  not_image_filter->SetInput(read);

  const typename DanielssonDistanceMapImageFilter::Pointer not_distance_map_image_filter =
      DanielssonDistanceMapImageFilter::New();
  not_distance_map_image_filter->SetInput(not_image_filter->GetOutput());

  using AddImageFilter = itk::AddImageFilter<OffsetImage>;
  const typename AddImageFilter::Pointer add_image_filter = AddImageFilter::New();
  add_image_filter->SetInput1(not_distance_map_image_filter->GetVectorDistanceMap());
  add_image_filter->SetInput2(distance_map_image_filter->GetVectorDistanceMap());
  add_image_filter->Update();

  using SubtractImageFilter = itk::SubtractImageFilter<FloatImage>;
  const typename SubtractImageFilter::Pointer subtract_image_filter = SubtractImageFilter::New();
  subtract_image_filter->SetInput1(not_distance_map_image_filter->GetDistanceMap());
  subtract_image_filter->SetInput2(distance_map_image_filter->GetDistanceMap());
  subtract_image_filter->Update();

  distance_map_ = deepCopy(subtract_image_filter->GetOutput(), std::identity());
  vector_map_ = deepCopy(add_image_filter->GetOutput(), [](const itk::Offset<D> &offset) {
    VectorF<D> result;
    for (unsigned i = 0; i < D; ++i)
      result[i] = static_cast<FLOAT_T>(offset[i]);
    return result;
  });
}

template <unsigned D> FLOAT_T EDTReader<D>::signedDistanceAt(const VectorF<D> &point) const noexcept {
  return distance_map_[projection(point)];
}

template <unsigned D> VectorF<D> EDTReader<D>::closestAt(const VectorF<D> &point) const noexcept {
  // NOLINTNEXTLINE(*-magic-numbers)
  return (clamp(point).array().template cast<FLOAT_T>() + FLOAT_T(0.5)).matrix() + vector_map_[projection(point)];
}

template <unsigned D> Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> EDTReader<D>::boundingBox() const noexcept {
  return bounding_box_;
}

template class EDTReader<4>;
} // namespace stmesh
  // NOLINTEND(misc-include-cleaner)