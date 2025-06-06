#ifndef STMESH_MIXD_HPP
#define STMESH_MIXD_HPP

#include <array>
#include <cstddef>
#include <filesystem>
#include <string_view>
#include <vector>

#include "utility.hpp"

namespace stmesh::mixd {
// kToOmit.at(i) tells us which vertex is not part of the face at index i in the mrng file line
constexpr std::array<size_t, 5> kToOmit{4, 3, 2, 0, 1}; ///< The indices of the vertices to omit from the face
// kOmitted.at(j) tells us which mrng entry contains the boundary region index of the face which does not include the
// vertex of index j
constexpr std::array<size_t, 5> kOmitted{3, 4, 2, 1, 0}; ///< The index of the face omitting the vertex
// thus: kToOmit.at(kOmitted.at(j)) == j and kOmitted.at(kToOmit.at(i)) == i

/// Writes the minf file
/**
 * Writes the minf file. This function writes the minf file, which contains the paths to the mxyz, mien, and mrng files
 * as well as the number of elements and nodes in the mesh.
 *
 * @param minf_file The path to the minf file to write
 * @param mxyz_file The path to the mxyz file
 * @param mien_file The path to the mien file
 * @param mrng_file The path to the mrng file
 * @param number_elements The number of elements in the mesh
 * @param number_nodes The number of nodes in the mesh
 */
void writeMinf(const std::filesystem::path &minf_file, const std::filesystem::path &mxyz_file,
               const std::filesystem::path &mien_file, const std::filesystem::path &mrng_file, size_t number_elements,
               size_t number_nodes);

/// Writes the neim file
/**
 * Writes the neim file. This function writes the neim file, which contains the node element index mapping. The mapping
 * is padded with zeros to ensure all rows have the same length.
 *
 * @param neim_file The path to the neim file to write
 * @param neim The node element index mapping
 */
void writeNeim(const std::filesystem::path &neim_file, const std::vector<std::vector<int>> &neim);

/// Writes the mxyz file
/**
 * Writes the mxyz file. This function writes the mxyz file, which contains the vertices of the mesh. The vertices are
 * written in the format "x y z t" for each vertex, in big-endian double format.
 *
 * @param file The path to the mxyz file to write
 * @param vertices The vertices of the mesh
 * @param scale The scale of the mesh
 * @param min_time The minimum time of the mesh
 * @return The path to the written mxyz file
 */
std::filesystem::path writeMxyz(std::filesystem::path file, const std::vector<Vector4F> &vertices,
                                stmesh::FLOAT_T scale, stmesh::FLOAT_T min_time);

/// Writes an integer based mixd file
/**
 * Writes an integer based mixd file. This function writes an integer based mixd file, in big-endian int format.
 * This may be used for writing mien or mrng files.
 *
 * @param file The path to the mixd file to write
 * @param extension The extension of the mixd file to replace in the path
 * @param full_cell_vertex_ids The full cell vertex ids of the mesh
 * @return The path to the written mixd file
 */
std::filesystem::path writeIntMixd(std::filesystem::path file, std::string_view extension,
                                   const std::vector<std::array<int, 5>> &full_cell_vertex_ids);

/// Checks if a pentatope element is positive
/**
 * Checks if a pentatope element is positive. This function checks if a pentatope element is positive, i.e. if the
 * vertices of the pentatope are in positive orientation. This is required for XNS.
 *
 * @param vertex_ids The vertex ids of the pentatope
 * @param vertices The vertices of the mesh
 * @return Whether the pentatope element is positive
 */
[[nodiscard]] bool positivePentatopeElementDet(const std::array<int, 5> &vertex_ids,
                                               const std::vector<Vector4F> &vertices);

/// Reads the mxyz file
/**
 * Reads the mxyz file. This function reads the mxyz file, which contains the vertices of the mesh. The vertices are
 * read in the format "x y z t" for each vertex, in big-endian double format.
 *
 * @param mxyz_file The path to the mxyz file to read
 * @return The vertices of the mesh
 */
std::vector<Vector4F> readMxyz(const std::filesystem::path &mxyz_file);

/// Reads an integer based mixd file
/**
 * Reads an integer based mixd file. This function reads an integer based mixd file, in big-endian int format.
 * This may be used for reading mien or mrng files.
 *
 * @param file The path to the mixd file to read
 * @return The full cell vertex ids of the mesh
 */
std::vector<std::array<int, 5>> readIntMixd(const std::filesystem::path &file);

/// Reads an arbitrary mixd data file
/**
 * Reads an arbitrary mixd data file. This function reads an arbitrary mixd data file, in big-endian double format.
 * This may be used for reading output files from XNS.
 *
 * @param data_file The path to the data file to read
 * @param n The number of values per node
 * @return The data read from the file
 */
std::vector<std::vector<FLOAT_T>> readData(const std::filesystem::path &data_file, size_t n);

/// Output struct for readMinf
/**
 * Contains the paths to the mxyz, mien, and mrng files.
 */
struct MinfData {
  std::filesystem::path mxyz_file; ///< The path to the mxyz file
  std::filesystem::path mien_file; ///< The path to the mien file
  std::filesystem::path mrng_file; ///< The path to the mrng file
};

/// Reads the minf file
/**
 * Reads the minf file. This function reads the minf file, which contains the paths to the mxyz, mien, and mrng files.
 * It ignores all other information in the minf file.
 *
 * @param minf_file The path to the minf file to read
 * @return The paths to the mxyz, mien, and mrng files
 */
MinfData readMinf(const std::filesystem::path &minf_file);
} // namespace stmesh::mixd

#endif
