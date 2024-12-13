include(GenerateExportHeader)


add_library(stmesh 
  geometric_simplex.cpp
  sdf.cpp
  triangulation.cpp
  utility.cpp
  vtk_writer.cpp
  edt.cpp
  boundary_region_manager.cpp
  mixd.cpp
  marching_hypercubes.cpp
  voxel_complex.cpp
  rle_bitset.cpp
  bitset.cpp
  lfs_schemes.cpp
  radius_schemes.cpp
  ${CMAKE_BINARY_DIR}/src/generator/table.cpp
)

add_dependencies(stmesh RunGenerator)

set_source_files_properties(
  ${CMAKE_BINARY_DIR}/src/generator/table.cpp
  PROPERTIES SKIP_LINTING ON
)

add_library(stmesh::stmesh ALIAS stmesh)

target_link_libraries(stmesh PRIVATE 
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
  mph::mph
  stmesh_options 
  stmesh_warnings)

target_link_libraries(stmesh PUBLIC 
  fmt::fmt
  CGAL::CGAL
  Eigen
  spdlog::spdlog
  OpenMP::OpenMP_CXX
  ${Boost_INCLUDE_DIRS})

target_include_directories(stmesh PUBLIC "${CMAKE_BINARY_DIR}/configured_files/include")

vtk_module_autoinit(
  TARGETS stmesh
  MODULES ${VTK_LIBRARIES}
)

target_include_directories(stmesh ${WARNING_GUARD} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                                                                  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>)

target_compile_features(stmesh PUBLIC cxx_std_20)

set_target_properties(
  stmesh
  PROPERTIES VERSION ${PROJECT_VERSION}
             CXX_VISIBILITY_PRESET hidden
             VISIBILITY_INLINES_HIDDEN YES)

generate_export_header(stmesh EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/include/stmesh/stmesh_export.hpp)

if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(stmesh PUBLIC STMESH_STATIC_DEFINE)
endif()

if(BUILD_TESTING)
  target_compile_definitions(stmesh PRIVATE STMESH_TESTING)
endif()

install(
  TARGETS stmesh
  EXPORT stmesh-targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include)
