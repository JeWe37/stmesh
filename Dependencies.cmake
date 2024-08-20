include(cmake/CPM.cmake)

# Done as a function so that updates to variables like
# CMAKE_CXX_FLAGS don't propagate out to other
# targets
function(stmesh_setup_dependencies)

  # For each dependency, see if it's
  # already been provided to us by a parent project

  set(CPM_USE_LOCAL_PACKAGES ON)

  if(NOT TARGET fmtlib::fmtlib)
    cpmaddpackage("gh:fmtlib/fmt#9.1.0")
  endif()

  if(NOT TARGET spdlog::spdlog)
    cpmaddpackage(
      NAME
      spdlog
      VERSION
      1.14.1
      GITHUB_REPOSITORY
      "gabime/spdlog"
      OPTIONS
      "SPDLOG_FMT_EXTERNAL ON")
  endif()

  if(NOT TARGET Catch2::Catch2WithMain)
    cpmaddpackage("gh:catchorg/Catch2@3.3.2")
  endif()

  if(NOT TARGET CLI11::CLI11)
    cpmaddpackage("gh:CLIUtils/CLI11@2.3.2")
  endif()

  if (NOT TARGET CGAL::CGAL)
    cpmaddpackage("gh:CGAL/CGAL@5.4")
  endif()

  set(ITK_USE_SYSTEM_EIGEN ON)
  find_package(ITK REQUIRED)
  include(${ITK_USE_FILE})
  set(ITK_LIBRARIES ${ITK_LIBRARIES} PARENT_SCOPE)

  find_package(VTK COMPONENTS 
    CommonColor
    CommonCore
    CommonDataModel
    FiltersGeneral
    IOXML
    InteractionStyle
    RenderingContextOpenGL2
    RenderingCore
    RenderingFreeType
    RenderingGL2PSOpenGL2
    RenderingOpenGL2
  )

  if(NOT VTK_FOUND)
      list(APPEND VTK_OPTIONS
          "BUILD_SHARED_LIBS OFF"
          "BUILD_TESTING OFF"
          "VTK_BUILD_EXAMPLES OFF"
          "VTK_BUILD_TESTING OFF"
          "VTK_ENABLE_WRAPPING OFF"
          "VTK_Group_Rendering OFF"
          "VTK_Group_StandAlone OFF"
          "VTK_USE_64BIT_IDS ON"
      )
      foreach(comp ${VTK_COMPONENTS})
          list(APPEND VTK_OPTIONS "Module_${comp} ON")
      endforeach()

      CPMAddPackage(
          NAME VTK
          GITHUB_REPOSITORY kitware/vtk
          VERSION 9.1.0
          OPTIONS ${VTK_OPTIONS}
          EXCLUDE_FROM_ALL YES
          GIT_SUBMODULES "" # Disable submodules
      )
      include(${VTK_BINARY_DIR}/VTKConfig.cmake)
  endif()

  set(VTK_LIBRARIES ${VTK_LIBRARIES} PARENT_SCOPE)

  CPMAddPackage(
    NAME Eigen
    VERSION 3.4.0
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    # Eigen's CMakelists are not intended for library use
    DOWNLOAD_ONLY YES 
  )

  if(Eigen_ADDED)
    add_library(Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
  endif()

  find_package(Boost 1.74.0 REQUIRED)

  find_package(OpenMP 4.5 REQUIRED)

  CPMAddPackage(
    Name mph
    GITHUB_REPOSITORY boost-ext/mph
    GIT_TAG v5.0.1
  )
  add_library(mph INTERFACE)
  target_include_directories(mph SYSTEM INTERFACE ${mph_SOURCE_DIR})
  if (CMAKE_CXX_COMPILER_ID MATCHES ".*GNU")
    target_compile_options(mph INTERFACE -fconstexpr-ops-limit=10000000000 -mbmi2)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    target_compile_options(mph INTERFACE -fconstexpr-steps=100000000 -mbmi2)
  endif()
  add_library(mph::mph ALIAS mph)

  CPMAddPackage(
    Name exprtk
    GITHUB_REPOSITORY ArashPartow/exprtk
    GIT_TAG f46bffcd6966d38a09023fb37ba9335214c9b959
    DOWNLOAD_ONLY YES
  )
  add_library(exprtk INTERFACE)
  target_include_directories(exprtk SYSTEM INTERFACE ${exprtk_SOURCE_DIR})
  add_library(exprtk::exprtk ALIAS exprtk)
endfunction()
