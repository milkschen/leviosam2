cmake_minimum_required(VERSION 3.10)

project(leviosam CXX C)

set(PROJECT_URL "https://github.com/alshai/leviosam")
set(PROJECT_DESCRIPTION "leviosam: library for moving alignments between aligned genomes")
set(CMAKE_CXX_STANDARD 11)

# Note: when I installed htslib and sdsl-lite on my Mac using brew,
# it created pkgconfig directories but no .cmake files or equivalent.
# I'm not a CMake pro, but I think that means we have to do this in
# two steps:
# 1. Set CMAKE_PREFIX_PATH to include the pkgconfig subdirectories
#    for htslib and sdsl-lite.  Mine is:
#    CMAKE_PREFIX_PATH=/usr/local/Cellar/sdsl-lite/2.1.1_2/lib/pkgconfig
# 2. Use the PkgConfig tool in CMake, as below

# htslib
find_package(PkgConfig REQUIRED)
pkg_search_module(PKG_HTSLIB REQUIRED htslib)

# sdsl
pkg_search_module(PKG_LIBDIVSUFSORT REQUIRED libdivsufsort)

add_executable(leviosam leviosam.cpp)
target_link_libraries(leviosam ${PKG_HTSLIB_LDFLAGS})
target_link_libraries(leviosam sdsl)
target_link_libraries(leviosam pthread)

INSTALL(TARGETS leviosam DESTINATION bin)
