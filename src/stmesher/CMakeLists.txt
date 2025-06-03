add_executable(stmesher main.cpp)

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

install(TARGETS stmesher DESTINATION bin)
