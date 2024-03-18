#ifndef STMESH_STATISTICS_COLLECTOR_HPP
#define STMESH_STATISTICS_COLLECTOR_HPP

#include <filesystem>
#include <fstream>
#include <iterator>

#include <Eigen/Geometry>
#include <fmt/core.h>
#include <fmt/ostream.h>

#include "geometric_simplex.hpp"
#include "surface_adapters.hpp" // IWYU pragma: keep
#include "triangulation.hpp"

namespace stmesh {
template <typename ExtraData>
void writeStatistics([[maybe_unused]] const std::filesystem::path &file,
                     [[maybe_unused]] const SurfaceAdapter4 auto &surface,
                     const Triangulation<ExtraData> &triangulation) {
  std::ofstream out(file);
  size_t full_cell_id = 0;
  out << "full_cell_id,content,shortest_edge_length,circumsphere_radius,radius_edge_ratio,quality\n";
  std::ostream_iterator<char> out_it(out);
  for (const auto &cell : triangulation) {
    GeometricSimplex<4> simplex =
        triangulation.fullCellSimplex(typename Triangulation<ExtraData>::FullCellConstHandle{&cell});
    HyperSphere4 circumsphere = simplex.circumsphere();
    if (surface.inside(circumsphere.center()))
      fmt::format_to(out_it, "{},{},{},{},{},{}\n", full_cell_id++, simplex.content(), simplex.shortestEdgeLength(),
                     circumsphere.radius(), simplex.radiusEdgeRatio(), simplex.quality());
  }
}
} // namespace stmesh
#endif