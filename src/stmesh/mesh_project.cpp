#include "stmesh/mesh_project.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <boost/geometry/index/inserter.hpp>
#include <boost/geometry/index/predicates.hpp>

#include "stmesh/geometric_simplex.hpp"
#include "stmesh/geometry_helpers.hpp"
#include "stmesh/problem_types.hpp"
#include "stmesh/triangulation.hpp"
#include "stmesh/utility.hpp"

namespace stmesh {
MeshProjector::MeshProjector(const TriangulationFromMixdWithData *triangulation) : triangulation_(triangulation) {
  std::transform(triangulation_->begin(), triangulation_->end(), bg::index::inserter(tree_),
                 [i = size_t{}](const auto &cell) mutable -> Tree::value_type {
                   return {boxFromAABB(cell.geometricSimplex().boundingBox()), i++};
                 });
}

Eigen::Vector<FLOAT_T, Eigen::Dynamic> MeshProjector::project(const Vector4F &point) const {
  constexpr static FLOAT_T kEps = sqrt(std::numeric_limits<FLOAT_T>::epsilon());
  constexpr static FLOAT_T kInf = std::numeric_limits<FLOAT_T>::infinity();
  const BGBox bg_box = boxFromAABB({point - Vector4F::Constant(kEps), point + Vector4F::Constant(kEps)});
  std::vector<std::pair<BGBox, size_t>> result;
  tree_.query(bg::index::intersects(bg_box), std::back_inserter(result));
  VectorF<5> best_coords = VectorF<5>::Constant(-kInf);
  Eigen::Matrix<FLOAT_T, Eigen::Dynamic, 5> best_data =
      Eigen::Matrix<FLOAT_T, Eigen::Dynamic, 5>::Zero(static_cast<Eigen::Index>(problemType().entries()), 5);
  for (const auto &[box, idx] : result) {
    auto cell = triangulation_->begin() + static_cast<TriangulationFromMixdWithData::iterator::difference_type>(idx);
    const VectorF<5> barycentric = cell->geometricSimplex().barycentricCoordinates(point);
    if (barycentric.minCoeff() > best_coords.minCoeff()) {
      best_coords = barycentric;
      best_data = cell->data();
      if (best_coords.minCoeff() >= 0.0)
        break;
    }
  }
  return best_data * best_coords;
}

const ProblemType &MeshProjector::problemType() const noexcept { return triangulation_->problem_type_; }
} // namespace stmesh