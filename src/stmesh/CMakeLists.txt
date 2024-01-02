include(GenerateExportHeader)


add_library(stmesh stmesh.cpp)



add_library(stmesh::stmesh ALIAS stmesh)

target_link_libraries(stmesh PRIVATE stmesh_options stmesh_warnings)

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
