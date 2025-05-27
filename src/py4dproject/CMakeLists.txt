nanobind_add_module(py4dproject NB_STATIC NB_SUPPRESS_WARNINGS STABLE_ABI py4dproject.cpp)
set_target_properties(nanobind-static-abi3 PROPERTIES CXX_CLANG_TIDY "")

set(OLD_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a") 
find_library(
    LIBC_LIB
    NAMES c
)
find_library(
    LIBGMP_LIB
    NAMES gmp
)
message(STATUS "libc: ${LIBC_LIB}")
message(STATUS "libgmp: ${LIBGMP_LIB}")
set(CMAKE_FIND_LIBRARY_SUFFIXES ${OLD_CMAKE_FIND_LIBRARY_SUFFIXES})

target_link_libraries(py4dproject PRIVATE 
    -nolibc -static-libgcc -static-libstdc++ ${LIBC_LIB}
    $<COMPILE_ONLY:stmesh::stmesh>
    $<TARGET_FILE:stmesh::stmesh>
    fmt::fmt
    ${LIBGMP_LIB}
)

install(TARGETS py4dproject DESTINATION lib/python3.12/site-packages)
