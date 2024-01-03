#ifndef SAMPLE_LIBRARY_HPP
#define SAMPLE_LIBRARY_HPP

#include "itkImage.h"
#include <CGAL/Dimension.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/point_generators_d.h>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <stmesh/stmesh_export.hpp>
#include <vector>

[[nodiscard]] STMESH_EXPORT int factorial(int) noexcept;

[[nodiscard]] constexpr int factorial_constexpr(int input) noexcept {
  if (input == 0) {
    return 1;
  }

  return input * factorial_constexpr(input - 1);
}

using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using Triangulation = CGAL::Triangulation<K>;

inline void test() {
  const int D = 5;   // we work in Euclidean 5-space
  const int N = 100; // we will insert 100 points
  // - - - - - - - - - - - - - - - - - - - - - - - - STEP 1
  const CGAL::Random_points_in_cube_d<Triangulation::Point> rand_it(D, 1.0);
  std::vector<Triangulation::Point> points;
  std::copy_n(rand_it, N, std::back_inserter(points));
  Triangulation t(D); // create triangulation
  assert(t.empty());
  t.insert(points.begin(), points.end()); // compute triangulation
  assert(t.is_valid());
  // - - - - - - - - - - - - - - - - - - - - - - - - STEP 2
  using Face = Triangulation::Face;
  using Faces = std::vector<Face>;
  Faces edges;
  const std::back_insert_iterator<Faces> out(edges);
  t.tds().incident_faces(t.infinite_vertex(), 1, out);
  // collect faces of dimension 1 (edges) incident to the infinite vertex
  std::cout << "There are " << edges.size() << " vertices on the convex hull.\n";
}

inline void test_itk() {
  using ImageType = itk::Image<unsigned short, 3>;

  ImageType::Pointer image = ImageType::New();

  std::cout << "Hello ITK World !" << std::endl;
}

void test_vtk();

#endif
