include(cmake/SystemLink.cmake)
include(cmake/LibFuzzer.cmake)
include(CMakeDependentOption)
include(CheckCXXCompilerFlag)


macro(stmesh_supports_sanitizers)
  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND NOT WIN32)
    set(SUPPORTS_UBSAN ON)
  else()
    message(WARNING "Undefined behavior sanitizer is not supported on this platform")
    set(SUPPORTS_UBSAN OFF)
  endif()

  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND WIN32)
    message(WARNING "Address sanitizer is not supported on this platform")
    set(SUPPORTS_ASAN OFF)
  else()
    set(SUPPORTS_ASAN ON)
  endif()
endmacro()

macro(stmesh_setup_options)
  option(stmesh_ENABLE_HARDENING "Enable hardening" ON)
  option(stmesh_ENABLE_COVERAGE "Enable coverage reporting" OFF)
  cmake_dependent_option(
    stmesh_ENABLE_GLOBAL_HARDENING
    "Attempt to push hardening options to built dependencies"
    ON
    stmesh_ENABLE_HARDENING
    OFF)

  stmesh_supports_sanitizers()

  if(NOT PROJECT_IS_TOP_LEVEL OR stmesh_PACKAGING_MAINTAINER_MODE)
    option(stmesh_ENABLE_IPO "Enable IPO/LTO" OFF)
    option(stmesh_WARNINGS_AS_ERRORS "Treat Warnings As Errors" OFF)
    option(stmesh_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(stmesh_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" OFF)
    option(stmesh_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(stmesh_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" OFF)
    option(stmesh_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(stmesh_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(stmesh_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    option(stmesh_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
    option(stmesh_ENABLE_CPPCHECK "Enable cpp-check analysis" OFF)
    option(stmesh_ENABLE_PCH "Enable precompiled headers" OFF)
    option(stmesh_ENABLE_CACHE "Enable ccache" OFF)
  else()
    option(stmesh_ENABLE_IPO "Enable IPO/LTO" ON)
    option(stmesh_WARNINGS_AS_ERRORS "Treat Warnings As Errors" ON)
    option(stmesh_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(stmesh_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" ${SUPPORTS_ASAN})
    option(stmesh_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(stmesh_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" ${SUPPORTS_UBSAN})
    option(stmesh_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(stmesh_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(stmesh_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    option(stmesh_ENABLE_CLANG_TIDY "Enable clang-tidy" ON)
    option(stmesh_ENABLE_CPPCHECK "Enable cpp-check analysis" ON)
    option(stmesh_ENABLE_PCH "Enable precompiled headers" OFF)
    option(stmesh_ENABLE_CACHE "Enable ccache" ON)
  endif()

  if(NOT PROJECT_IS_TOP_LEVEL)
    mark_as_advanced(
      stmesh_ENABLE_IPO
      stmesh_WARNINGS_AS_ERRORS
      stmesh_ENABLE_USER_LINKER
      stmesh_ENABLE_SANITIZER_ADDRESS
      stmesh_ENABLE_SANITIZER_LEAK
      stmesh_ENABLE_SANITIZER_UNDEFINED
      stmesh_ENABLE_SANITIZER_THREAD
      stmesh_ENABLE_SANITIZER_MEMORY
      stmesh_ENABLE_UNITY_BUILD
      stmesh_ENABLE_CLANG_TIDY
      stmesh_ENABLE_CPPCHECK
      stmesh_ENABLE_COVERAGE
      stmesh_ENABLE_PCH
      stmesh_ENABLE_CACHE)
  endif()

  stmesh_check_libfuzzer_support(LIBFUZZER_SUPPORTED)
  if(LIBFUZZER_SUPPORTED AND (stmesh_ENABLE_SANITIZER_ADDRESS OR stmesh_ENABLE_SANITIZER_THREAD OR stmesh_ENABLE_SANITIZER_UNDEFINED))
    set(DEFAULT_FUZZER ON)
  else()
    set(DEFAULT_FUZZER OFF)
  endif()

  option(stmesh_BUILD_FUZZ_TESTS "Enable fuzz testing executable" ${DEFAULT_FUZZER})

endmacro()

macro(stmesh_global_options)
  if(stmesh_ENABLE_IPO)
    include(cmake/InterproceduralOptimization.cmake)
    stmesh_enable_ipo()
  endif()

  stmesh_supports_sanitizers()

  if(stmesh_ENABLE_HARDENING AND stmesh_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN 
       OR stmesh_ENABLE_SANITIZER_UNDEFINED
       OR stmesh_ENABLE_SANITIZER_ADDRESS
       OR stmesh_ENABLE_SANITIZER_THREAD
       OR stmesh_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    message("${stmesh_ENABLE_HARDENING} ${ENABLE_UBSAN_MINIMAL_RUNTIME} ${stmesh_ENABLE_SANITIZER_UNDEFINED}")
    stmesh_enable_hardening(stmesh_options ON ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()
endmacro()

macro(stmesh_local_options)
  if(PROJECT_IS_TOP_LEVEL)
    include(cmake/StandardProjectSettings.cmake)
  endif()

  add_library(stmesh_warnings INTERFACE)
  add_library(stmesh_options INTERFACE)

  include(cmake/CompilerWarnings.cmake)
  stmesh_set_project_warnings(
    stmesh_warnings
    ${stmesh_WARNINGS_AS_ERRORS}
    ""
    ""
    ""
    "")

  if(stmesh_ENABLE_USER_LINKER)
    include(cmake/Linker.cmake)
    configure_linker(stmesh_options)
  endif()

  include(cmake/Sanitizers.cmake)
  stmesh_enable_sanitizers(
    stmesh_options
    ${stmesh_ENABLE_SANITIZER_ADDRESS}
    ${stmesh_ENABLE_SANITIZER_LEAK}
    ${stmesh_ENABLE_SANITIZER_UNDEFINED}
    ${stmesh_ENABLE_SANITIZER_THREAD}
    ${stmesh_ENABLE_SANITIZER_MEMORY})

  set_target_properties(stmesh_options PROPERTIES UNITY_BUILD ${stmesh_ENABLE_UNITY_BUILD})

  if(stmesh_ENABLE_PCH)
    target_precompile_headers(
      stmesh_options
      INTERFACE
      <vector>
      <string>
      <utility>)
  endif()

  if(stmesh_ENABLE_CACHE)
    include(cmake/Cache.cmake)
    stmesh_enable_cache()
  endif()

  include(cmake/StaticAnalyzers.cmake)
  if(stmesh_ENABLE_CLANG_TIDY)
    stmesh_enable_clang_tidy(stmesh_options ${stmesh_WARNINGS_AS_ERRORS})
  endif()

  if(stmesh_ENABLE_CPPCHECK)
    stmesh_enable_cppcheck(${stmesh_WARNINGS_AS_ERRORS} "" # override cppcheck options
    )
  endif()

  if(stmesh_ENABLE_COVERAGE)
    include(cmake/Tests.cmake)
    stmesh_enable_coverage(stmesh_options)
  endif()

  if(stmesh_WARNINGS_AS_ERRORS)
    check_cxx_compiler_flag("-Wl,--fatal-warnings" LINKER_FATAL_WARNINGS)
    if(LINKER_FATAL_WARNINGS)
      # This is not working consistently, so disabling for now
      # target_link_options(stmesh_options INTERFACE -Wl,--fatal-warnings)
    endif()
  endif()

  if(stmesh_ENABLE_HARDENING AND NOT stmesh_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN 
       OR stmesh_ENABLE_SANITIZER_UNDEFINED
       OR stmesh_ENABLE_SANITIZER_ADDRESS
       OR stmesh_ENABLE_SANITIZER_THREAD
       OR stmesh_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    stmesh_enable_hardening(stmesh_options OFF ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()

endmacro()
