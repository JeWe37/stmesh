#ifndef STMESH_TRIANGULATION_HPP
#define STMESH_TRIANGULATION_HPP

#include <array>
#include <cstddef>
#include <filesystem>
#include <iterator>
#include <optional>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Dimension.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation_data_structure.h>
#include <CGAL/Triangulation_full_cell.h>
#include <CGAL/Triangulation_vertex.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <boost/geometry.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/parameters.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include "boundary_region_manager.hpp"
#include "edt.hpp"
#include "geometric_simplex.hpp"
#include "geometry_helpers.hpp"
#include "problem_types.hpp"
#include "surface_adapters.hpp"
#include "utility.hpp"
#include "writable_triangulation.hpp"

namespace stmesh {
template <typename MixdCell> class TriangulationFromMixdGeneric;
class TriangulationFromMixdWithData;

namespace detail {
class MixdDataCell;

template <typename Cell = std::monostate> class MixdCell {
  using CellType = std::conditional_t<std::is_same_v<Cell, std::monostate>, MixdCell, Cell>;

protected:
  const TriangulationFromMixdGeneric<CellType> *triangulation_;
  size_t index_;

  friend TriangulationFromMixdGeneric<CellType>;

public:
  MixdCell(const TriangulationFromMixdGeneric<CellType> *triangulation, size_t index);

  [[nodiscard]] GeometricSimplex<4> geometricSimplex() const noexcept;

  [[nodiscard]] bool isSurfaceSide(size_t i) const noexcept;
};

class MixdDataCell : public MixdCell<MixdDataCell> {
  friend class TriangulationFromMixdWithData;

public:
  using MixdCell::MixdCell;

  [[nodiscard]] Eigen::Matrix<FLOAT_T, Eigen::Dynamic, 5> data() const noexcept;
};
} // namespace detail

/// Triangulation created by loading a set of mixd files
/**
 * Triangulation created by loading a set of mixd files. This class represents a triangulation created by loading a
 * the files specified in the minf file. It can act as a WritableTriangulation and BoundaryRegionManager.
 */
template <typename MixdCell = detail::MixdCell<>> class TriangulationFromMixdGeneric {
protected:
  std::vector<Vector4F> mxyz_;
  std::vector<std::array<int, 5>> mien_;
  std::vector<std::array<int, 5>> mrng_;

  friend MixdCell;
  friend class detail::MixdCell<detail::MixdDataCell>;

public:
  /// An iterator for the elements in the file.
  /**
   * An iterator for the elements in the file. This class is an iterator for the elements in the file, and is used to
   * iterate over this class, returning WritableCell objects. Usable as a random access iterator.
   */
  // NOLINTNEXTLINE(readability-identifier-naming)
  class iterator : public boost::iterator_adaptor<iterator, std::vector<std::array<int, 5>>::const_iterator, MixdCell,
                                                  std::random_access_iterator_tag, MixdCell> {
    const TriangulationFromMixdGeneric *triangulation_;

    friend class boost::iterator_core_access;

    [[nodiscard]] MixdCell dereference() const noexcept;

  public:
    /// A constructor from a triangulation and an iterator
    /**
     * A constructor from a triangulation and an iterator. This constructor creates an iterator from a triangulation
     * pointer and an iterator.
     *
     * @param triangulation The triangulation to create the iterator for
     * @param it The iterator to create the iterator for
     */
    iterator(const TriangulationFromMixdGeneric *triangulation,
             const std::vector<std::array<int, 5>>::const_iterator &it);

    /// A default constructor
    /**
     * A default constructor. This constructor creates an iterator with a null triangulation pointer and a default
     * iterator, serving as the end iterator.
     */
    iterator();
  };

  /// A constructor from a minf file path
  /**
   * A constructor from a minf file path. This constructor creates a triangulation from a minf file path. The minf file
   * is read and the vertices, elements, and boundary regions are loaded into the triangulation.
   *
   * @param minf_file The path to the minf file
   */
  explicit TriangulationFromMixdGeneric(const std::filesystem::path &minf_file);

  /// Gets the bounding box of the triangulation
  /**
   * Gets the bounding box of the triangulation. This function returns the bounding box of the triangulation, which is
   * the bounding box of all the vertices in the triangulation.
   *
   * @return The bounding box
   */
  [[nodiscard]] Eigen::AlignedBox<double, 4> boundingBox() const noexcept;

  /// Gets the begin iterator for the elements in the file
  /**
   * Gets the begin iterator for the elements in the file. This function returns an iterator to the beginning of the
   * elements in the file, as WritableCell objects.
   *
   * @return The begin iterator for the elements in the file
   */
  [[nodiscard]] iterator begin() const noexcept;

  /// Gets the end iterator for the elements in the file
  /**
   * Gets the end iterator for the elements in the file. This function returns an iterator to the end of the elements
   * in the file, as WritableCell objects.
   *
   * @return The end iterator for the elements in the file
   */
  [[nodiscard]] iterator end() const noexcept;

  /// Gets the boundary region of a face
  /**
   * Gets the boundary region of a face. This function gets the boundary region of a face, i.e. the index of the
   * boundary region that the face belongs to. These are obtained from the mrng file.
   *
   * @param cell The cell a face of which to get the boundary region for
   * @param j The index of the face(i.e. the index of its covertex)
   * @return The boundary region
   */
  [[nodiscard]] size_t findBoundaryRegion(const MixdCell &cell, size_t j) const noexcept;
};

using TriangulationFromMixd = TriangulationFromMixdGeneric<>; ///< Standard TriangulationFromMixd

static_assert(WritableTriangulation<TriangulationFromMixd>);
static_assert(BoundaryRegionManager<TriangulationFromMixd>);

/// Triangulation created by loading a set of mixd files with additional data
/**
 * Triangulation created by loading a set of mixd files with additional data. This class represents a triangulation
 * created by loading a the files specified in the minf file, along with additional data from a data file, which can be
 * accessed through the data() method on its cells. It can act as a WritableTriangulation and BoundaryRegionManager.
 * Usable as a random access iterator.
 */
class TriangulationFromMixdWithData : public TriangulationFromMixdGeneric<detail::MixdDataCell> {
  std::vector<std::vector<FLOAT_T>> data_;

  friend class detail::MixdDataCell;

public:
  ProblemType problem_type_; ///< The problem type of the triangulation

  /// A constructor from a minf file path and a data file path
  /**
   * A constructor from a minf file path and a data file path. This constructor creates a triangulation from a minf
   * file path and a data file path. The minf file is read and the vertices, elements, and boundary regions are loaded
   * into the triangulation. The data file directly is read and the data is loaded into the triangulation.
   *
   * @param minf_file The path to the minf file
   * @param data_file The path to the data file
   * @param problem_type The problem type data to read
   */
  TriangulationFromMixdWithData(const std::filesystem::path &minf_file, const std::filesystem::path &data_file,
                                const std::optional<ProblemType> &problem_type);
};

static_assert(WritableTriangulation<TriangulationFromMixdWithData>);
static_assert(BoundaryRegionManager<TriangulationFromMixdWithData>);

namespace detail {
template <typename ExtraData, SurfaceAdapter4 Surface> class WritableMainTriangulation;
}

/// A Delaunay triangulation in 4D
/**
 * A Delaunay triangulation in 4D. This class is a wrapper around CGAL's Delaunay_triangulation_4 class, and provides
 * additional functionality for working with the triangulation.
 * In particular, it provides functionality for determining which full cells are added or removed when a vertex is
 * inserted or removed. Additionally, it provides functionality for determining if newly inserted vertices are good
 * points.
 * The triangulation is templated on the type of extra data stored in the triangulation. This extra data is
 * stored in the full cells of the triangulation, and can be used to store additional information about the full
 * cells.
 *
 * @tparam ExtraData The type of extra data stored in the triangulation
 */
template <typename ExtraData = std::monostate> class Triangulation {
  Eigen::AlignedBox<FLOAT_T, 4> bounding_box_;

  struct FullCellData {
    bool committed;
    ExtraData extra_data;
  };

  struct VertexData {
    bool nonfree_vertex;
  };

public:
  using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<4>>; ///< The CGAL kernel used by the triangulation
  using Point = Kernel::Point_d;                        ///< A CGAL point in 4D

  using TriangulationFullCell =
      CGAL::Triangulation_full_cell<Kernel, FullCellData>;                    ///< A full cell in the triangulation
  using TriangulationVertex = CGAL::Triangulation_vertex<Kernel, VertexData>; ///< A vertex in the triangulation
  using TriangulationDataStructure =
      CGAL::Triangulation_data_structure<Kernel::Dimension, TriangulationVertex,
                                         TriangulationFullCell>; ///< The data structure used by the triangulation
  using DelaunayTriangulation =
      CGAL::Delaunay_triangulation<Kernel, TriangulationDataStructure>; ///< The CGAL Delaunay triangulation

  using FullCell = DelaunayTriangulation::Full_cell;              ///< A full cell in the triangulation
  using VertexHandle = DelaunayTriangulation::Vertex_handle;      ///< A handle to a vertex in the triangulation
  using FullCellHandle = DelaunayTriangulation::Full_cell_handle; ///< A handle to a full cell in the triangulation
  using FullCellConstHandle =
      DelaunayTriangulation::Full_cell_const_handle; ///< A handle to a const full cell in the triangulation
  using Facet =
      DelaunayTriangulation::Facet; ///< A facet of the triangulation, consisting of a full cell and an opposite index
  using Face = DelaunayTriangulation::Face;              ///< A face of the triangulation
  using LocateType = DelaunayTriangulation::Locate_type; ///< The result of a locate operation

  using iterator = DelaunayTriangulation::Finite_full_cell_iterator; ///< An iterator for the finite full cells
  using const_iterator =
      DelaunayTriangulation::Finite_full_cell_const_iterator; ///< A const iterator for the finite full cells

private:
  DelaunayTriangulation triangulation_;

  using Tree = bg::index::rtree<BGPoint,
                                // NOLINTNEXTLINE(*-magic-numbers)
                                bg::index::rstar<16>>;

  Tree tree_;

  std::unordered_map<Vector4F, VertexHandle, Vector4FHash> vertex_handle_map_;

public:
  /// Converts a CGAL point to an Eigen vector
  /**
   * Converts a CGAL point to an Eigen vector. Useful for converting points from the triangulation to Eigen vectors.
   *
   * @param pt The point to convert
   * @return The vector
   */
  [[nodiscard]] static Vector4F pointToVec(const Point &pt);

  /// A constructor from a bounding box
  /**
   * A constructor from a bounding box. This constructor creates a Delaunay triangulation from a bounding box. The
   * corners of the bounding box are inserted as initial vertices into the triangulation.
   *
   * @param bounding_box The bounding box to create the triangulation from
   */
  explicit Triangulation(const Eigen::AlignedBox<FLOAT_T, 4> &bounding_box);

  /// Inserts a vertex into the triangulation
  /**
   * Inserts a vertex into the triangulation. This function inserts a vertex into the triangulation, and returns a
   * handle to the inserted vertex. It also updates the rtree with the new vertex.
   * The vertex is inserted at the location specified by the hint.
   * A vertex can be marked as nonfree, which means it lies on the box or surface of the triangulation.
   * During this process, the conflicting full cells are removed and full cells surrounding the vertex are added.
   *
   * @param point The point to insert
   * @param hint A hint for the location to insert the vertex
   * @param nonfree_vertex Whether the vertex is nonfree
   * @return The handle to the inserted vertex
   */
  VertexHandle insert(const Vector4F &point, FullCellHandle hint = {}, bool nonfree_vertex = false);

  /// Removes a vertex from the triangulation
  /**
   * Removes a vertex from the triangulation. This function removes a vertex from the triangulation, and returns a
   * handle to the new full cell containing the removed vertex. It also updates the rtree with the removed vertex.
   * During this process, the full cells surrounding the vertex are removed and surrounding the returned
   * FullCellHandle, uncommitted full cells are added.
   *
   * @param vertex The vertex to remove
   * @return The handle to the removed vertex
   */
  FullCellHandle remove(VertexHandle vertex);

  /// Finds all full cells in conflict with a point
  /**
   * Finds all full cells in conflict with a point. This function finds all full cells in conflict with a point, i.e.
   * all full cells whose circumsphere contains the point. It returns a vector of handles to the conflicting full
   * cells.
   *
   * @param point The point to find the conflicting full cells for
   * @param hint A hint for the location to search for the conflicting full cells
   * @return A vector of handles to the conflicting full cells
   */
  [[nodiscard]] std::vector<FullCellHandle> conflictZone(const Vector4F &point, FullCellHandle hint = {}) const;

  /// Find surrounding uncommitted full cells and commit them
  /**
   * Find surrounding uncommitted full cells and commit them. This function finds all uncommitted full cells
   * surrounding a full cell, and commits them. It returns a vector of handles to the committed full cells.
   *
   * @param start The full cell to start the search from
   * @return A vector of handles to the committed full cells
   */
  std::vector<FullCellHandle> commitUncommitted(FullCellHandle start);

  /// Find surrounding full cells
  /**
   * Find surrounding full cells. This function finds all full cells surrounding a vertex, and returns a vector of
   * handles to the surrounding full cells.
   *
   * @param vertex The vertex to find the surrounding full cells for
   * @param commit Whether to commit the surrounding full cells
   * @return A vector of handles to the surrounding full cells
   */
  std::vector<FullCellHandle> surroundingFullCells(VertexHandle vertex, bool commit = true);

  /// Find the mirror index of a full cell
  /**
   * Find the mirror index of a full cell. This function finds the mirror index of a full cell, i.e. the index of the
   * opposite vertex in the cell on the other side of the face opposite to the vertex with the given index.
   *
   * @param full_cell The full cell to find the mirror index for
   * @param index The index of the vertex in the full cell
   * @return The mirror index
   */
  [[nodiscard]] int mirrorIndex(const FullCellHandle full_cell, int index) const;

  /// Gets all vertices of a facet
  /**
   * Gets all vertices of a facet. This function gets all vertices of a facet, and returns them as a matrix, where
   * each column is a vertex.
   *
   * @param facet The facet to get the vertices for
   * @return The vertices of the facet
   */
  [[nodiscard]] Eigen::Matrix<FLOAT_T, 4, 4> facetVertices(Facet facet) const noexcept;

  /// Gets the facet from a matrix of vertices
  /**
   * Gets the facet from a matrix of vertices. This function gets the facet from a matrix of vertices, where each
   * column is a vertex. It returns the facet by locating it in the triangulation.
   *
   * @param vertices The matrix of vertices
   * @return The facet
   */
  [[nodiscard]] Facet facetFromVertices(const Eigen::Matrix<FLOAT_T, 4, 4> &vertices) const;

  /// Gets information about the cells on both sides of a facet
  /**
   * Gets information about the cells on both sides of a facet. This function gets information about the cells on both
   * sides of a facet, and returns it as a tuple. The first element of the tuple is a tuple containing the vertex
   * handle of the vertex opposite to the facet, the index of the vertex in the full cell, and the handle to the full
   * cell. The second element of the tuple is an optional tuple containing the same information for the cell on the
   * other side of the facet, if it is a finite cell.
   *
   * @param facet The facet to get the information for
   * @return A tuple containing information about the cells on both sides of the facet
   */
  [[nodiscard]] std::tuple<std::tuple<VertexHandle, int, FullCellHandle>,
                           std::optional<std::tuple<VertexHandle, int, FullCellHandle>>>
  facetMirrorVertices(Facet facet) const;

  /// Gets the vertex count of the triangulation
  /**
   * Gets the vertex count of the triangulation. This function returns the total number of vertices in the
   *triangulation.
   *
   * @return The vertex count
   */
  [[nodiscard]] size_t vertexCount() const noexcept;

  /// Gets the vertices within a radius of a point
  /**
   * Gets the vertices within a radius of a point. This function gets the vertices within a radius of a point, and
   * returns them as a vector of vertex handles. Internally it uses the rtree to query the vertices within the radius
   * efficiently.
   *
   * @param point The point to get the vertices within a radius for
   * @param radius The radius to search for vertices within
   * @return The vertices within the radius
   */
  [[nodiscard]] std::vector<VertexHandle> verticesInRadius(const Vector4F &point, FLOAT_T radius) const;

  /// Gets a GeometricSimplex<4> from a full cell
  /**
   * Gets a GeometricSimplex<4> from a full cell. This function gets a GeometricSimplex<4> from a full cell, and
   * returns it. The GeometricSimplex<4> is constructed from the vertices of the full cell.
   *
   * @param full_cell The full cell to get the GeometricSimplex<4> for
   * @return The GeometricSimplex<4>
   */
  [[nodiscard]] static GeometricSimplex<4> fullCellSimplex(FullCellConstHandle full_cell) noexcept;

  /// Checks if a point is a good point
  /**
   * Checks if a point is a good point. This function checks if a point is a good point, i.e. if the insertion of the
   * point did not cause the creation of any small sliver simplices. Equivalently, this checks that none of the full
   * cells surrounding the point are small sliver simplices.
   *
   * @param vertex The vertex to check
   * @param rho_bar The rho_bar parameter
   * @param tau_bar The tau_bar parameter
   * @param max_radius The maximum radius
   * @return Whether the point is a good point
   */
  [[nodiscard]] bool isGoodPoint(VertexHandle vertex, FLOAT_T rho_bar, FLOAT_T tau_bar, FLOAT_T max_radius) const;

  /// Checks if a full cell is infinite
  /**
   * Checks if a full cell is infinite. This function checks if a full cell is infinite, i.e. if it is a finite full
   * cell or an infinite full cell.
   *
   * @param full_cell The full cell to check
   * @return Whether the full cell is infinite
   */
  [[nodiscard]] bool isInfinite(FullCellConstHandle full_cell) const;

  /// Gets the bounding box of the triangulation
  /**
   * Gets the bounding box of the triangulation. This function returns the bounding box of the triangulation,
   * initially set in the constructor.
   *
   * @return The bounding box
   */
  [[nodiscard]] const Eigen::AlignedBox<FLOAT_T, 4> &boundingBox() const noexcept;

  /// Gets the const begin iterator for the finite full cells
  /**
   * Gets the const begin iterator for the finite full cells. This function returns a const iterator to the beginning
   * of the finite full cells of the triangulation.
   *
   * @return The const begin iterator for the finite full cells
   */
  [[nodiscard]] const_iterator cbegin() const;

  /// Gets the const end iterator for the finite full cells
  /**
   * Gets the const end iterator for the finite full cells. This function returns a const iterator to the end of the
   * finite full cells of the triangulation.
   *
   * @return The const end iterator for the finite full cells
   */
  [[nodiscard]] const_iterator cend() const;

  /// Gets the begin iterator for the finite full cells
  /**
   * Gets the begin iterator for the finite full cells. This function returns an iterator to the beginning of the
   * finite full cells of the triangulation.
   *
   * @return The begin iterator for the finite full cells
   */
  [[nodiscard]] iterator begin();

  /// Gets the end iterator for the finite full cells
  /**
   * Gets the end iterator for the finite full cells. This function returns an iterator to the end of the finite full
   * cells of the triangulation.
   *
   * @return The end iterator for the finite full cells
   */
  [[nodiscard]] iterator end();

  /// Gets the const begin iterator for the finite full cells
  /**
   * Gets the const begin iterator for the finite full cells. This function returns a const iterator to the beginning
   * of the finite full cells of the triangulation.
   *
   * @return The const begin iterator for the finite full cells
   */
  [[nodiscard]] const_iterator begin() const;

  /// Gets the const end iterator for the finite full cells
  /**
   * Gets the const end iterator for the finite full cells. This function returns a const iterator to the end of the
   * finite full cells of the triangulation.
   *
   * @return The const end iterator for the finite full cells
   */
  [[nodiscard]] const_iterator end() const;

  template <typename ExtraData2, SurfaceAdapter4 Surface> friend class detail::WritableMainTriangulation;

  /// Gets a writable triangulation from the triangulation
  /**
   * Gets a writable triangulation from the triangulation. This function returns a writable triangulation from the
   * triangulation, which can be used to write out the solution contained in it. As the solution only contains cells
   * whose circumcenters are inside the surface, the surface is needed to determine which cells are inside the surface.
   *
   * @param surface The surface to use for determining which cells are inside the surface
   * @return The writable triangulation
   */
  template <SurfaceAdapter4 Surface>
  [[nodiscard]] detail::WritableMainTriangulation<ExtraData, Surface>
  writableTriangulation(Surface *surface) const noexcept {
    return {this, surface};
  }
};

namespace detail {
template <typename ExtraData, SurfaceAdapter4 Surface> class WritableMainTriangulation {
  using TriangulationED = Triangulation<ExtraData>;
  const TriangulationED *triangulation_;
  const Surface *surface_;

  class WritableMainCell {
    const TriangulationED *triangulation_;
    const Surface *surface_;
    TriangulationED::FullCellConstHandle handle_;

  public:
    WritableMainCell(const TriangulationED *triangulation, const Surface *surface,
                     TriangulationED::FullCellConstHandle handle) noexcept
        : triangulation_(triangulation), surface_(surface), handle_(handle) {}

    [[nodiscard]] GeometricSimplex<4> geometricSimplex() const noexcept {
      return triangulation_->fullCellSimplex(handle_);
    }

    [[nodiscard]] bool isSurfaceSide(size_t i) const noexcept {
      auto neighbor = handle_->neighbor(i);
      return triangulation_->isInfinite(neighbor) ||
             !surface_->inside(triangulation_->fullCellSimplex(neighbor).circumsphere().center());
    }
  };

public:
  // NOLINTNEXTLINE(readability-identifier-naming)
  class iterator : public boost::iterator_adaptor<iterator, typename TriangulationED::FullCellConstHandle,
                                                  WritableMainCell, std::forward_iterator_tag, WritableMainCell> {
    const TriangulationED *triangulation_;
    const Surface *surface_;

    friend class boost::iterator_core_access;
    using iterator::iterator_adaptor_::base_reference;

    [[nodiscard]] WritableMainCell dereference() const noexcept { return {triangulation_, surface_, base_reference()}; }

    void increment() noexcept {
      ++base_reference();
      while (base_reference() != triangulation_->triangulation_.full_cells_end() &&
             (triangulation_->isInfinite(base_reference()) ||
              !surface_->inside(triangulation_->fullCellSimplex(base_reference()).circumsphere().center())))
        ++base_reference();
      if (base_reference() == triangulation_->triangulation_.full_cells_end()) {
        triangulation_ = nullptr;
        surface_ = nullptr;
      }
    }

  public:
    iterator() noexcept : triangulation_(nullptr), surface_(nullptr) {}

    iterator(TriangulationED::FullCellConstHandle it, const TriangulationED *triangulation,
             const Surface *surface) noexcept
        : iterator::iterator_adaptor_(it), triangulation_(triangulation), surface_(surface) {}
  };

  WritableMainTriangulation(const TriangulationED *triangulation, const Surface *surface) noexcept
      : triangulation_(triangulation), surface_(surface) {}

  [[nodiscard]] Eigen::AlignedBox<FLOAT_T, 4> boundingBox() const noexcept { return triangulation_->boundingBox(); }

  [[nodiscard]] iterator begin() const noexcept {
    auto it = triangulation_->triangulation_.full_cells_begin();
    while (it != triangulation_->triangulation_.full_cells_end() &&
           (triangulation_->isInfinite(it) ||
            !surface_->inside(triangulation_->fullCellSimplex(it).circumsphere().center())))
      ++it;
    return {it, triangulation_, surface_};
  }

  [[nodiscard]] iterator end() const noexcept {
    return {triangulation_->triangulation_.full_cells_end(), triangulation_, surface_};
  }
};
} // namespace detail
} // namespace stmesh

template <typename ExtraData, typename Surface>
// NOLINTNEXTLINE(readability-identifier-naming)
constexpr bool std::ranges::enable_borrowed_range<stmesh::detail::WritableMainTriangulation<ExtraData, Surface>> = true;

namespace stmesh {
static_assert(
    WritableTriangulation<detail::WritableMainTriangulation<std::monostate, EDTSurfaceAdapter<EDTReader<4>>>>);
} // namespace stmesh

#endif
