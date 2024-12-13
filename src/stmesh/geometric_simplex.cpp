#include "stmesh/geometric_simplex.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <CGAL/Dimension.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/convex_hull_3.h>
#include <Eigen/Core>

#include "stmesh/sdf.hpp"
#include "stmesh/utility.hpp"

namespace stmesh {
template <unsigned D, unsigned N>
GeometricSimplex<D, N>::GeometricSimplex(
    const Eigen::Matrix<FLOAT_T, static_cast<int>(D), static_cast<int>(N)> &vertices)
    : vertices_(vertices) {}

template <unsigned D, unsigned N> GeometricSimplex<D, N>::GeometricSimplex() = default;

template <unsigned D, unsigned N>
[[nodiscard]] const Eigen::Matrix<FLOAT_T, static_cast<int>(D), static_cast<int>(N)> &
GeometricSimplex<D, N>::vertices() const noexcept {
  return vertices_;
}

template <unsigned D, unsigned N>
[[nodiscard]] Eigen::Matrix<FLOAT_T, static_cast<int>(N + 1), static_cast<int>(N + 1)>
GeometricSimplex<D, N>::caleyMengerMatrix() const noexcept {
  Eigen::Matrix<FLOAT_T, static_cast<int>(N + 1), static_cast<int>(N + 1)> result =
      Eigen::Matrix<FLOAT_T, static_cast<int>(N + 1), static_cast<int>(N + 1)>::Ones();
  result.diagonal() = Eigen::Matrix<FLOAT_T, static_cast<int>(N + 1), 1>::Zero();
  for (unsigned i = 1; i < N + 1; ++i) {
    for (unsigned j = i + 1; j < N + 1; ++j)
      result(j, i) = result(i, j) = (vertices_.col(i - 1) - vertices_.col(j - 1)).squaredNorm();
  }
  return result;
}

// https://mathworld.wolfram.com/Cayley-MengerDeterminant.html
template <unsigned D, unsigned N>
[[nodiscard]] FLOAT_T GeometricSimplex<D, N>::caleyMengerDeterminant() const noexcept {
  return caleyMengerMatrix().determinant();
}

template <unsigned D, unsigned N> [[nodiscard]] FLOAT_T GeometricSimplex<D, N>::shortestEdgeLength() const noexcept {
  auto edge_lengths = allSquaredEdgeLengths();
  return std::sqrt(*std::ranges::min_element(edge_lengths));
}

template <unsigned D, unsigned N> [[nodiscard]] FLOAT_T GeometricSimplex<D, N>::content() const noexcept {
  return std::sqrt((N % 2 == 0 ? 1 : -1) * caleyMengerDeterminant() /
                   (std::pow(factorial(N - 1), FLOAT_T{2}) * std::pow(FLOAT_T{2}, FLOAT_T{N - 1})));
}

template <unsigned D, unsigned N> [[nodiscard]] HyperSphere<D> GeometricSimplex<D, N>::circumsphere() const noexcept {
  using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<static_cast<int>(D)>>;
  using Point = Kernel::Point_d;
  std::vector<Point> points;
  points.reserve(N + 1);
  std::ranges::transform(vertices_.colwise(), std::back_inserter(points),
                         [](const auto &vertex) { return Point(vertex.begin(), vertex.end()); });
  const Point point = Kernel().construct_circumcenter_d_object()(points.begin(), points.end());
  const Vector4F center = {CGAL::to_double(point[0]), CGAL::to_double(point[1]), CGAL::to_double(point[2]),
                           CGAL::to_double(point[3])};
  const FLOAT_T radius = (center - vertices_.col(0)).norm();
  return {radius, center};
}

template <unsigned D, unsigned N> [[nodiscard]] FLOAT_T GeometricSimplex<D, N>::radiusEdgeRatio() const noexcept {
  return circumsphere().radius() / shortestEdgeLength();
}

template <unsigned D, unsigned N> [[nodiscard]] FLOAT_T GeometricSimplex<D, N>::quality() const noexcept {
  return content() / std::pow(shortestEdgeLength(), FLOAT_T{N - 1});
}

template <unsigned D, unsigned N>
[[nodiscard]] Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(D)> GeometricSimplex<D, N>::normalRay() const noexcept
requires(D == N)
{
  const Eigen::Matrix<FLOAT_T, static_cast<int>(D), 1> kern =
      kernel<D, D - 1>(vertices_(Eigen::all, Eigen::seq(Eigen::fix<1>, Eigen::last)).colwise() - vertices_.col(0));
  return {circumsphere().center(), kern};
}

template <unsigned D, unsigned N>
[[nodiscard]] std::pair<Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(D)>, FLOAT_T>
GeometricSimplex<D, N>::toParameterizedLine() const noexcept
requires(N == 2)
{
  return {Eigen::ParametrizedLine<FLOAT_T, static_cast<int>(D)>::Through(vertices_.col(0), vertices_.col(1)),
          (vertices_.col(0) - vertices_.col(1)).norm()};
}

// TODO: this is cheating. also hard assuming planes with normal in direction of 4th axis
template <unsigned D, unsigned N>
[[nodiscard]] CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel>
GeometricSimplex<D, N>::planeCut(const Eigen::Hyperplane<FLOAT_T, static_cast<int>(D)> &plane) const
requires(D == 4 && N >= 4)
{
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  std::vector<Kernel::Point_3> points;
  auto edges = subSimplices<2>();
  for (const auto &edge : edges) {
    auto [line, factor] = edge.toParameterizedLine();
    const FLOAT_T lambda = line.intersectionParameter(plane);
    if (lambda <= factor && lambda >= FLOAT_T{0}) {
      VectorF<D> point = line.pointAt(lambda);
      points.emplace_back(point[0], point[1], point[2]);
    }
  }
  CGAL::Polyhedron_3<Kernel> result{};
  if constexpr (N == 5)
    CGAL::convex_hull_3(points.begin(), points.end(), result);
  else {
    if (points.size() >= 4) {
      for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = 0; j < points.size(); ++j) {
          constexpr double kEps = 1e-9;
          if (i != j && Kernel().compute_squared_distance_3_object()(points[i], points[j]) <
                            kEps * std::pow(shortestEdgeLength(), 2)) {
            points.erase(points.begin() + static_cast<std::vector<Kernel::Point_3>::difference_type>(j));
            --j;
            if (i > j)
              --i;
          }
        }
      }
    }
    if (points.size() == 4) {
      Eigen::Matrix<FLOAT_T, 3, 4> matrix;
      for (unsigned i = 0; i < 4; ++i)
        std::copy(points[i].cartesian_begin(), points[i].cartesian_end(), matrix.col(i).begin());
      const Eigen::JacobiSVD<Eigen::Matrix<FLOAT_T, 3, 3>> svd(
          matrix(Eigen::all, Eigen::seq(Eigen::fix<1>, Eigen::last)).colwise() - matrix.col(0), Eigen::ComputeFullU);
      const Eigen::Matrix<FLOAT_T, 3, 2> projection = svd.matrixU().leftCols(2);
      const Eigen::Vector<FLOAT_T, 2> center = (projection.transpose() * matrix).rowwise().mean();
      std::array<Eigen::Vector<FLOAT_T, 3>, 4> sorted;
      std::ranges::copy(matrix.colwise(), sorted.begin());
      std::ranges::sort(sorted, {}, [&](const Eigen::Vector<FLOAT_T, 3> &vec) {
        const Eigen::Vector<FLOAT_T, 2> centered = (projection.transpose() * vec).colwise() - center;
        return std::atan2(centered[1], centered[0]);
      });
      CGAL::Polyhedron_incremental_builder_3<CGAL::Polyhedron_3<Kernel>::HalfedgeDS> builder(result.hds());
      builder.begin_surface(4, 1);
      for (const auto &point : sorted)
        builder.add_vertex(Kernel::Point_3(point[0], point[1], point[2]));
      builder.begin_facet();
      for (size_t i = 0; i < 4; ++i)
        builder.add_vertex_to_facet(i);
      builder.end_facet();
      builder.end_surface();
    } else if (points.size() == 3)
      result.make_triangle(points[0], points[1], points[2]);
  }
  return result;
}

template <unsigned D, unsigned N>
[[nodiscard]] Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> GeometricSimplex<D, N>::boundingBox() const noexcept {
  Eigen::AlignedBox<FLOAT_T, static_cast<int>(D)> box;
  for (const auto &vert : vertices_.colwise())
    box.extend(vert);
  return box;
}

template <unsigned D, unsigned N>
[[nodiscard]] std::array<FLOAT_T, nChoosek(N, 2)> GeometricSimplex<D, N>::allSquaredEdgeLengths() const noexcept {
  std::array<FLOAT_T, nChoosek(N, 2)> result{};
  size_t i = 0;
  for (unsigned j = 0; j < N; ++j) {
    for (unsigned k = j + 1; k < N; ++k)
      result.at(i++) = (vertices_.col(j) - vertices_.col(k)).squaredNorm();
  }
  return result;
}

template <unsigned D, unsigned N>
[[nodiscard]] FLOAT_T GeometricSimplex<D, N>::omega() const noexcept
requires(D == 4 && N == 5)
{
  auto l = allSquaredEdgeLengths();
  auto square = [](auto x) { return x * x; };
  // NOLINTBEGIN(*-magic-numbers)
  return FLOAT_T(600.0) * square(l[0] - l[1]) + FLOAT_T(900.0) * square(l[4]) +
         FLOAT_T(100.0) * square(FLOAT_T(-2.0) * (l[0] + l[1]) + l[4]) +
         FLOAT_T(75.0) * square(l[0] - l[1] - FLOAT_T(3.0) * (l[5] + l[7])) +
         FLOAT_T(25.0) * square(l[0] + l[1] + l[4] - FLOAT_T(3.0) * (l[2] + l[5] + l[7])) +
         FLOAT_T(25.0) *
             square(l[0] + l[1] - FLOAT_T(6.0) * l[2] - FLOAT_T(2.0) * l[4] + FLOAT_T(3.0) * (l[5] + l[7])) +
         FLOAT_T(45.0) * square(l[0] - l[1] + l[5] - FLOAT_T(4.0) * l[6] - l[7] + FLOAT_T(4.0) * l[8]) +
         FLOAT_T(15.0) * square(l[0] + l[1] + FLOAT_T(2.0) * l[2] - FLOAT_T(8.0) * l[3] - FLOAT_T(2.0) * l[4] - l[5] +
                                FLOAT_T(4.0) * l[6] - l[7] + FLOAT_T(4.0) * l[8]) +
         FLOAT_T(30.0) * square(-l[0] - l[1] + l[2] + FLOAT_T(2.0) * l[3] - l[4] + l[5] + FLOAT_T(2.0) * l[6] + l[7] +
                                FLOAT_T(2.0) * l[8] - FLOAT_T(6.0) * l[9]) +
         FLOAT_T(9.0) * square(l[0] + l[1] + l[2] - FLOAT_T(4.0) * l[3] + l[4] + l[5] - FLOAT_T(4.0) * l[6] + l[7] -
                               FLOAT_T(4.0) * (l[8] + l[9]));
  // NOLINTEND(*-magic-numbers)
}

template <unsigned D, unsigned N>
[[nodiscard]] FLOAT_T GeometricSimplex<D, N>::metric1() const noexcept
requires(D == 4 && N == 5)
{
  auto l = allSquaredEdgeLengths();
  // NOLINTNEXTLINE(*-magic-numbers)
  return std::pow(FLOAT_T(5.0), FLOAT_T(3.0 / 4.0)) * std::sqrt(FLOAT_T(384.0) * content()) /
         std::accumulate(l.begin(), l.end(), FLOAT_T(0.0));
}

template <unsigned D, unsigned N>
[[nodiscard]] FLOAT_T GeometricSimplex<D, N>::metric2() const noexcept
requires(D == 4 && N == 5)
{
  auto l = allSquaredEdgeLengths();
  // NOLINTNEXTLINE(*-magic-numbers)
  return FLOAT_T(6.0) * std::accumulate(l.begin(), l.end(), FLOAT_T(0.0)) / std::sqrt(omega());
}

template <unsigned D, unsigned N>
[[nodiscard]] FLOAT_T GeometricSimplex<D, N>::metric3() const noexcept
requires(D == 4 && N == 5)
{
  return metric1() * metric2();
}

template <unsigned D, unsigned N>
[[nodiscard]] bool GeometricSimplex<D, N>::wellShaped(FLOAT_T rho_bar, FLOAT_T tau_bar) const noexcept {
  return radiusEdgeRatio() < rho_bar && quality() >= tau_bar;
}

template <unsigned D, unsigned N>
[[nodiscard]] bool GeometricSimplex<D, N>::sliver(FLOAT_T rho_bar, FLOAT_T tau_bar) const noexcept {
  if constexpr (N < 4)
    return false;
  else {
    if (radiusEdgeRatio() >= rho_bar || quality() >= tau_bar)
      return false;
    if constexpr (N > 4)
      // TODO: actually just need to check one dimension down for N > 5, but irrelevant for us anyway
      return std::apply(
          [&](auto &&...args) {
            return (std::accumulate(args.begin(), args.end(), true,
                                    [&](bool value, const auto &sub_simplex) {
                                      return value && sub_simplex.wellShaped(rho_bar, tau_bar);
                                    }) &&
                    ...);
          },
          allSubSimplicies<N - 1, 4>());
    return true;
  }
}

// already clear: edges/vertices are never bad, always good
// faces, facets and cells might be tho
// faces also always good for reasons explained in XY Li paper
// so for 4d: need to check if cell is bad and all 3d facets are good OR if any 3d facet is not good

// in general: for any k-facet in d-simplex(3<=k<=d) which is bad, check that all p-facets(3<=p<=k-1) contained in it
// are good

// weird: different definition of sliver in XY Li paper. there they check the goodness of the facets of the bad face,
// not of the entire cell at dimension

// actually equivalent for 4d case, it s fine
template <unsigned D, unsigned N>
[[nodiscard]] unsigned GeometricSimplex<D, N>::sliverSimplex(FLOAT_T rho_bar, FLOAT_T tau_bar) const noexcept {
  if constexpr (N < 4)
    return 0U;
  else {
    return std::apply(
        [&](auto &&...args) {
          unsigned i = 2;
          // increments i until a sliver is found, at which point it short circuits and returns i. if no sliver is
          // found, sets i to 0 at the end and thus returns 0
          ((i++, std::accumulate(args.begin(), args.end(), false,
                                 [&](bool value, const auto &sub_simplex) {
                                   return value || sub_simplex.sliver(rho_bar, tau_bar);
                                 })) ||
           ... || (i = 0, true));
          return i;
        },
        allSubSimplicies<N, 4>());
  }
  /*auto sub_simplices = allSubSimplicies();
  auto bad_simplices_test = [&](auto simplices) {
      bool good = true;
      bool bad = false;
      for (const auto &simplex : simplices) {
      FLOAT_T rho = simplex.radius_edge_ratio();
      FLOAT_T tau = simplex.quality();
      if (rho < rho_bar && tau < tau_bar)
          bad = true;
      if (rho > rho_bar || tau < tau_bar)
          good = false;
      }
      return std::make_pair(good, bad);
  };
  auto conditions =
      std::apply([&](auto &&...args) { return std::array{bad_simplices_test(args)...}; }, sub_simplices);
  int first_bad = 0;
  for (; first_bad < conditions.size(); ++first_bad) {
      if (conditions[first_bad].second)
      break;
  }
  int last_good = conditions.size() - 1;
  for (; last_good >= 0; --last_good) {
      if (conditions[last_good].first)
      break;
  }
  return first_bad > last_good;*/
}

template <unsigned D, unsigned N>
[[nodiscard]] bool GeometricSimplex<D, N>::smallSliverSimplex(FLOAT_T rho_bar, FLOAT_T tau_bar,
                                                              FLOAT_T max_radius) const noexcept {
  return circumsphere().radius() < max_radius && sliverSimplex(rho_bar, tau_bar);
}

template <unsigned D, unsigned N>
[[nodiscard]] bool GeometricSimplex<D, N>::operator==(const GeometricSimplex<D, N> &other) const noexcept = default;

template <unsigned D, unsigned N>
[[nodiscard]] bool GeometricSimplex<D, N>::operator!=(const GeometricSimplex<D, N> &other) const noexcept = default;

template class GeometricSimplex<4>;
template class GeometricSimplex<4, 4>;
template class GeometricSimplex<4, 3>;
template class GeometricSimplex<4, 2>;
} // namespace stmesh
