cmake_minimum_required(VERSION 3.10)

project(liftover CXX C)

set(PROJECT_URL "https://github.com/alshai/liftover")
set(PROJECT_DESCRIPTION "liftover: library for moving alignments between aligned genomes")
set(CMAKE_CXX_STANDARD 14)

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



# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

enable_testing()

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

add_executable(liftover liftover.cpp)
target_link_libraries(liftover ${PKG_HTSLIB_LDFLAGS})
target_link_libraries(liftover sdsl)

add_executable(liftover_test liftover_test.cpp)
target_link_libraries(liftover_test ${PKG_HTSLIB_LDFLAGS})
target_link_libraries(liftover_test sdsl)
target_link_libraries(liftover_test gtest gtest_main)

add_test(liftover_test liftover_test)
