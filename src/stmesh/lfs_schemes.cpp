#include "stmesh/lfs_schemes.hpp"

#include <memory>

#include "stmesh/edt.hpp"
#include "stmesh/utility.hpp"

namespace stmesh::lfs_schemes {
Constant::Constant(const FLOAT_T value) : value_(value) {}

FLOAT_T Constant::operator()(const stmesh::Vector4F & /* unused */) const noexcept { return value_; }

FLOAT_T Constant::max() const noexcept { return value_; }

BinaryImageApproximation::BinaryImageApproximation(const FLOAT_T value,
                                                   const std::shared_ptr<EDTReader<4, true>> &edt_reader)
    : value_(value), edt_reader_(edt_reader) {}

FLOAT_T BinaryImageApproximation::operator()(const stmesh::Vector4F &vec) const noexcept {
  return value_ * edt_reader_->distanceToThinnedAt(vec);
}

FLOAT_T BinaryImageApproximation::max() const noexcept {
  return value_ * edt_reader_->boundingBox().sizes().maxCoeff() / 2.0;
}
} // namespace stmesh::lfs_schemes
