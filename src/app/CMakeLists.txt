add_executable(stmesher stmesher.cpp)

target_link_libraries(
  stmesher
  PRIVATE stmesh::stmesh_options
          stmesh::stmesh_warnings)

target_link_system_libraries(
  stmesher
  PRIVATE
          CLI11::CLI11
          fmt::fmt
          exprtk::exprtk
          spdlog::spdlog
          stmesh::stmesh)

target_include_directories(stmesher PRIVATE "${CMAKE_BINARY_DIR}/configured_files/include")

add_executable(mesh_analyzer mesh_analyzer.cpp)

target_link_libraries(
  mesh_analyzer
  PRIVATE stmesh::stmesh_options
          stmesh::stmesh_warnings)

target_link_system_libraries(
  mesh_analyzer
  PRIVATE
          CLI11::CLI11
          fmt::fmt
          spdlog::spdlog
          stmesh::stmesh)

target_include_directories(mesh_analyzer PRIVATE "${CMAKE_BINARY_DIR}/configured_files/include")

add_executable(data_visualizer data_visualizer.cpp)

target_link_libraries(
  data_visualizer
  PRIVATE stmesh::stmesh_options
          stmesh::stmesh_warnings)

target_link_system_libraries(
  data_visualizer
  PRIVATE
          CLI11::CLI11
          fmt::fmt
          spdlog::spdlog
          stmesh::stmesh)

target_include_directories(data_visualizer PRIVATE "${CMAKE_BINARY_DIR}/configured_files/include")

