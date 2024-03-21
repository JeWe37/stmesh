#ifndef STMESH_STATISTICS_COLLECTOR_HPP
#define STMESH_STATISTICS_COLLECTOR_HPP

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>

#include <Eigen/Geometry>
#include <fmt/core.h>
#include <fmt/ostream.h>

#include "geometric_simplex.hpp"
#include "sdf.hpp"
#include "surface_adapters.hpp" // IWYU pragma: keep
#include "triangulation.hpp"

namespace stmesh {
/// Write statistics to a file from a triangulation
/**
 * Writes statistics to a file from a triangulation. The statistics written are:
 * - The full cell id
 * - The content(volume) of the full cell
 * - The shortest edge length of the full cell
 * - The radius of the circumsphere of the full cell
 * - The ratio of the radius of the circumsphere to the shortest edge length
 * - The quality of the full cell
 * Only full cells whose circumsphere is inside the surface are written to the statistics file.
 * The file is written in CSV format.
 *
 * @param file The path to the file to write the statistics to
 * @param surface The surface adapter to use to determine which full cells to write to the statistics file
 * @param triangulation The triangulation to write the statistics from
 * @tparam ExtraData The type of extra data stored in the triangulation
 */
template <typename ExtraData>
void writeStatistics([[maybe_unused]] const std::filesystem::path &file,
                     [[maybe_unused]] const SurfaceAdapter4 auto &surface,
                     const Triangulation<ExtraData> &triangulation) {
  std::ofstream out(file);
  size_t full_cell_id = 0;
  out << "full_cell_id,content,shortest_edge_length,circumsphere_radius,radius_edge_ratio,quality,metric1,metric2,"
         "metric3\n";
  std::ostream_iterator<char> out_it(out);
  for (const auto &cell : triangulation) {
    GeometricSimplex<4> simplex =
        triangulation.fullCellSimplex(typename Triangulation<ExtraData>::FullCellConstHandle{&cell});
    HyperSphere4 circumsphere = simplex.circumsphere();
    if (surface.inside(circumsphere.center()))
      fmt::format_to(out_it, "{},{},{},{},{},{},{},{},{}\n", full_cell_id++, simplex.content(),
                     simplex.shortestEdgeLength(), circumsphere.radius(), simplex.radiusEdgeRatio(), simplex.quality(),
                     simplex.metric1(), simplex.metric2(), simplex.metric3());
  }
}
} // namespace stmesh
#endif