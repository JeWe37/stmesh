#include "stmesh/vtk_writer.hpp"

#include <cstddef>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <vtkType.h>

#include "stmesh/utility.hpp"

namespace stmesh::detail {
[[nodiscard]] vtkIdType PointStorage::insert(const Vector3F &point) noexcept {
  if (positions_.contains(point))
    return positions_[point];
  return positions_[point] = points_->InsertNextPoint(static_cast<double>(point[0]), static_cast<double>(point[1]),
                                                      static_cast<double>(point[2]));
}

PolyhedraStorage::PolyhedraStorage(
    const CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> &polyhedron, PointStorage &points,
    size_t polyhedron_id)
    : polyhedron_(polyhedron), polyhedron_id_(polyhedron_id) {
  for (auto it_vert = polyhedron.vertices_begin(); it_vert != polyhedron.vertices_end(); ++it_vert) {
    const auto &vertex = it_vert->point();
    const Vector3F point{vertex.x(), vertex.y(), vertex.z()};
    pointIds_.push_back(points.insert(point));
  }
}
} // namespace stmesh::detail
