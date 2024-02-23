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
} // namespace stmesh::detail
