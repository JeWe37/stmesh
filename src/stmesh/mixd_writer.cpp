#include "stmesh/mixd_writer.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string_view>
#include <vector>

#include "stmesh/utility.hpp"

namespace stmesh::detail {
void writeMinf(const std::filesystem::path &minf_file, const std::filesystem::path &mxyz_file,
               const std::filesystem::path &mien_file, const std::filesystem::path &mrng_file, size_t number_elements,
               size_t number_nodes) {
  std::ofstream out(minf_file);
  out << "ne " << number_elements << "\n";
  out << "nn " << number_nodes << "\n";
  out << "nen " << 5 << "\n";
  out << "nef " << 5 << "\n";
  out << "nsd " << 4 << "\n";
  // ndf probably not needed
  // NOLINTNEXTLINE(*-magic-numbers)
  // out << "ndf " << 4 << "\n";

  out << "mxyz " << mxyz_file.string() << "\n";
  out << "mien " << mien_file.string() << "\n";
  out << "mrng " << mrng_file.string() << "\n";
}

std::filesystem::path writeMxyz(std::filesystem::path file, const std::vector<Vector4F> &vertices) {
  file.replace_extension(".mxyz");
  std::ofstream out(file);
  for (const auto &vertex : vertices) {
    for (const auto &coord : vertex) {
      auto coordinate_bytes = std::bit_cast<std::array<char, sizeof(FLOAT_T)>>(coord);
      if constexpr (std::endian::native != std::endian::big)
        std::reverse(coordinate_bytes.begin(), coordinate_bytes.end());
      out.write(coordinate_bytes.data(), sizeof(FLOAT_T));
    }
  }
  return file;
}

std::filesystem::path writeIntMixd(std::filesystem::path file, std::string_view extension,
                                   const std::vector<std::array<size_t, 5>> &full_cell_vertex_ids) {
  file.replace_extension(extension);
  std::ofstream out(file);
  for (const auto &vertex_ids : full_cell_vertex_ids) {
    for (const auto &id : vertex_ids) {
      auto id_bytes = std::bit_cast<std::array<char, sizeof(int)>>(static_cast<int>(id));
      if constexpr (std::endian::native != std::endian::big)
        std::reverse(id_bytes.begin(), id_bytes.end());
      out.write(id_bytes.data(), sizeof(int));
    }
  }
  return file;
}

bool positivePentatopeElementDet(const std::array<size_t, 5> &vertex_ids, const std::vector<Vector4F> &vertices) {
  Vector4F xr1 = vertices[vertex_ids[0]] - vertices[vertex_ids[4]];
  Vector4F xr2 = vertices[vertex_ids[1]] - vertices[vertex_ids[4]];
  Vector4F xr3 = vertices[vertex_ids[2]] - vertices[vertex_ids[4]];
  Vector4F xr4 = vertices[vertex_ids[3]] - vertices[vertex_ids[4]];

  const double cf11 = +xr4[3] * (xr2[1] * xr3[2] - xr3[1] * xr2[2]) - xr4[2] * (xr2[1] * xr3[3] - xr3[1] * xr2[3]) +
                      xr4[1] * (xr2[2] * xr3[3] - xr3[2] * xr2[3]);
  const double cf12 = -xr4[3] * (xr1[1] * xr3[2] - xr3[1] * xr1[2]) + xr4[2] * (xr1[1] * xr3[3] - xr3[1] * xr1[3]) -
                      xr4[1] * (xr1[2] * xr3[3] - xr3[2] * xr1[3]);
  const double cf13 = +xr4[3] * (xr1[1] * xr2[2] - xr2[1] * xr1[2]) - xr4[2] * (xr1[1] * xr2[3] - xr2[1] * xr1[3]) +
                      xr4[1] * (xr1[2] * xr2[3] - xr2[2] * xr1[3]);
  const double cf14 = -xr3[3] * (xr1[1] * xr2[2] - xr2[1] * xr1[2]) + xr3[2] * (xr1[1] * xr2[3] - xr2[1] * xr1[3]) -
                      xr3[1] * (xr1[2] * xr2[3] - xr2[2] * xr1[3]);

  const double dettmp = xr1[0] * cf11 + xr2[0] * cf12 + xr3[0] * cf13 + xr4[0] * cf14;

  return dettmp > FLOAT_T();
}
} // namespace stmesh::detail