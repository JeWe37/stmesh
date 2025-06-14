cmake_minimum_required(VERSION 3.15...3.23)

project(CmakeConfigPackageTests LANGUAGES CXX)

# ---- Test as standalone project the exported config package ----

if(PROJECT_IS_TOP_LEVEL OR TEST_INSTALLED_VERSION)
  enable_testing()

  find_package(stmesh CONFIG REQUIRED) # for stmesher, project_options, ...

  if(NOT TARGET stmesh_options)
    message(FATAL_ERROR "Requiered config package not found!")
    return() # be strictly paranoid for Template Janitor github action! CK
  endif()
endif()

# ---- Dependencies ----

include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)

# Provide a simple smoke test to make sure that the CLI works and can display a --help message
add_test(NAME cli.has_help COMMAND stmesher --help)

# Provide a test to verify that the version being reported from the application
# matches the version given to CMake. This will be important once you package
# your program. Real world shows that this is the kind of simple mistake that is easy
# to make, but also easy to test for.
add_test(NAME cli.version_matches COMMAND stmesher --version)
set_tests_properties(cli.version_matches PROPERTIES PASS_REGULAR_EXPRESSION "${PROJECT_VERSION}")

add_executable(normal_tests 
  edt_tests.cpp
  geometric_simplex_tests.cpp
  marching_hypercubes_tests.cpp
  meshing_algorithm_tests.cpp
  sdf_tests.cpp
  surface_adapter_tests.cpp
  triangulation_tests.cpp
  lfs_schemes_tests.cpp
  utility_tests.cpp
  bitset_tests.cpp
  rle_bitset_tests.cpp
  radius_schemes_tests.cpp
  voxel_complex_tests.cpp
  mixd_tests.cpp
  mesh_projector_tests.cpp
)
target_link_libraries(
  normal_tests
  PRIVATE stmesh::stmesh_warnings
          stmesh::stmesh_options
          stmesh::stmesh
          Catch2::Catch2WithMain)

if(WIN32 AND BUILD_SHARED_LIBS)
  add_custom_command(
    TARGET normal_tests
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:normal_tests> $<TARGET_FILE_DIR:normal_tests>
    COMMAND_EXPAND_LISTS)
endif()

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# automatically discover tests that are defined in catch based test files you can modify the unittests. Set TEST_PREFIX
# to whatever you want, or use different for different binaries
catch_discover_tests(
  normal_tests
  TEST_PREFIX
  "unittests."
  PROPERTIES
  TIMEOUT
  3600
  REPORTER
  XML
  OUTPUT_DIR
  .
  OUTPUT_PREFIX
  "unittests."
  OUTPUT_SUFFIX
  .xml)

# Add a file containing a set of constexpr tests
add_executable(constexpr_tests constexpr_tests.cpp)
target_link_libraries(
  constexpr_tests
  PRIVATE stmesh::stmesh_warnings
          stmesh::stmesh_options
          stmesh::stmesh
          Catch2::Catch2WithMain)

catch_discover_tests(
  constexpr_tests
  TEST_PREFIX
  "constexpr."
  REPORTER
  XML
  OUTPUT_DIR
  .
  OUTPUT_PREFIX
  "constexpr."
  OUTPUT_SUFFIX
  .xml)

# Disable the constexpr portion of the test, and build again this allows us to have an executable that we can debug when
# things go wrong with the constexpr testing
add_executable(relaxed_constexpr_tests constexpr_tests.cpp)
target_link_libraries(
  relaxed_constexpr_tests
  PRIVATE stmesh::stmesh_warnings
          stmesh::stmesh_options
          stmesh::stmesh
          Catch2::Catch2WithMain)
target_compile_definitions(relaxed_constexpr_tests PRIVATE -DCATCH_CONFIG_RUNTIME_STATIC_REQUIRE)

catch_discover_tests(
  relaxed_constexpr_tests
  TEST_PREFIX
  "relaxed_constexpr."
  REPORTER
  XML
  OUTPUT_DIR
  .
  OUTPUT_PREFIX
  "relaxed_constexpr."
  OUTPUT_SUFFIX
  .xml)
