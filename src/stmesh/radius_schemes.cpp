#include "stmesh/radius_schemes.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>

#include <itkImage.h>
#include <itkImageFileReader.h>

#include "stmesh/utility.hpp"

namespace stmesh::radius_schemes {
Constant::Constant(const FLOAT_T value) : value_(value) {}

FLOAT_T Constant::operator()(const Vector4F & /* unused */) const noexcept { return value_; }

struct ImageRadius::Impl {
  using Image = itk::Image<float, 4>;
  using Index = Image::IndexType;
  Image::Pointer read_image;

  explicit Impl(const std::filesystem::path &path) : read_image(itk::ReadImage<Image>(path)) {
    auto spacing = read_image->GetSpacing();
    std::array<float, 4> origin = {static_cast<float>(spacing[0] * 0.5), static_cast<float>(spacing[1] * 0.5),
                                   static_cast<float>(spacing[2] * 0.5), static_cast<float>(spacing[3] * 0.5)};
    read_image->SetOrigin(origin.data());
  }

  [[nodiscard]] Index vectorToIndex(const Vector4F &vec) const noexcept {
    using Point = Image::PointType;
    using Size = Image::SizeType;

    Size size = read_image->GetLargestPossibleRegion().GetSize();

    const Point point(vec.data());
    Index result = read_image->TransformPhysicalPointToIndex(point);
    for (size_t i = 0; i < 4; ++i)
      result[static_cast<Index::size_type>(i)] =
          std::clamp(result[static_cast<Index::size_type>(i)], typename Index::IndexValueType(0),
                     typename Index::IndexValueType(size[static_cast<Size::size_type>(i)] - 1));
    return result;
  }
};

ImageRadius::~ImageRadius() = default;
ImageRadius::ImageRadius(ImageRadius &&) noexcept = default;
ImageRadius &ImageRadius::operator=(ImageRadius &&) noexcept = default;

ImageRadius::ImageRadius(const std::filesystem::path &path) : pimpl_(new Impl(path)) {}

FLOAT_T ImageRadius::operator()(const Vector4F &vec) const noexcept {
  return static_cast<FLOAT_T>(pimpl_->read_image->GetPixel(pimpl_->vectorToIndex(vec)));
}
} // namespace stmesh::radius_schemes
