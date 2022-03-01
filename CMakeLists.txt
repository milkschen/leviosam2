cmake_minimum_required(VERSION 3.10)

project(leviosam CXX C)

set(PROJECT_URL "https://github.com/alshai/leviosam")
set(PROJECT_DESCRIPTION "leviosam: lifting over alignments between genomes")
set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_CXX_STANDARD 11)

find_library(HTS_LIB hts)
if(NOT HTS_LIB)
    message(FATAL_ERROR "htslib library not found. Please specify -D CMAKE_LIBRARY_PATH and -D CMAKE_INCLUDE_PATH in the CMake command")
endif()

find_library(SDSL_LIB sdsl)
if(NOT SDSL_LIB)
    message(FATAL_ERROR "libsdsl not found. Please specify -D CMAKE_LIBRARY_PATH and -D CMAKE_INCLUDE_PATH in the CMake command")
endif()


# From https://www.mattkeeter.com/blog/2018-01-06-versioning/
execute_process(COMMAND git log --pretty=format:'%h' -n 1
                OUTPUT_VARIABLE GIT_REV
                ERROR_QUIET)

# Check whether we got any revision (which isn't
# always the case, e.g. when someone downloaded a zip
# file from Github instead of a checkout)
if ("${GIT_REV}" STREQUAL "")
    set(GIT_REV "")
    set(GIT_DIFF "")
    set(GIT_TAG "")
    set(GIT_BRANCH "")
else()
    execute_process(
        COMMAND bash -c "git diff --quiet --exit-code || echo +"
        OUTPUT_VARIABLE GIT_DIFF)
    #execute_process(
    #    COMMAND git describe --exact-match --tags
    #    OUTPUT_VARIABLE GIT_TAG ERROR_QUIET)
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        OUTPUT_VARIABLE GIT_BRANCH)

    string(STRIP "${GIT_REV}" GIT_REV)
    string(SUBSTRING "${GIT_REV}" 1 7 GIT_REV)
    string(STRIP "-${GIT_REV}" GIT_REV)
    string(STRIP "${GIT_DIFF}" GIT_DIFF)
    #string(STRIP "-${GIT_TAG}" GIT_TAG)
    string(STRIP "-${GIT_BRANCH}" GIT_BRANCH)
endif()

set(VERSION
"#ifndef VERSION
#define VERSION \"0.5.2${GIT_BRANCH}${GIT_REV}${GIT_DIFF}\"
#endif
"
)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/version.hpp)
    file(READ ${CMAKE_CURRENT_SOURCE_DIR}/src/version.hpp VERSION_)
else()
    set(VERSION_ "")
endif()

if (NOT "${VERSION}" STREQUAL "${VERSION_}")
    file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/src/version.hpp "${VERSION}")
endif()


include_directories(${CMAKE_INCLUDE_PATH})

file(GLOB LIB_SOURCES   src/reconcile.cpp src/lift_bed.cpp src/collate.cpp src/chain.cpp
                        src/leviosam_utils.cpp src/bed.cpp src/lift_bed.cpp src/aln.cpp
                        src/ksw2_extd2_sse.c src/ksw2_extz.c
                        src/bam_md.c src/bam_aux.c src/gzstream.C)
file(GLOB LIB_HEADERS   src/leviosam.hpp
                        src/reconcile.hpp src/lift_bed.hpp src/collate.hpp src/chain.hpp
                        src/leviosam_utils.hpp src/bed.hpp src/lift_bed.hpp src/aln.hpp
                        src/ksw.h
                        src/bam.h src/IITree.h src/gzstream.h)
add_library(lvsam ${LIB_SOURCES} ${LIB_HEADERS})

add_executable(leviosam src/leviosam.cpp)
target_link_libraries(leviosam lvsam)
target_link_libraries(leviosam ${HTS_LIB})
target_link_libraries(leviosam ${SDSL_LIB})
target_link_libraries(leviosam pthread)
target_link_libraries(leviosam z)


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


add_executable(leviosam_test src/leviosam_test.cpp)
target_link_libraries(leviosam_test lvsam)
target_link_libraries(leviosam_test ${HTS_LIB})
target_link_libraries(leviosam_test ${SDSL_LIB})
target_link_libraries(leviosam_test gtest gtest_main)

add_test(NAME leviosam_test COMMAND leviosam_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testdata
)


add_executable(chain_test src/chain_test.cpp)
target_link_libraries(chain_test lvsam)
target_link_libraries(chain_test ${HTS_LIB})
target_link_libraries(chain_test ${SDSL_LIB})
target_link_libraries(chain_test gtest gtest_main)

add_test(NAME chain_test COMMAND chain_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testdata
)

INSTALL(TARGETS leviosam DESTINATION bin)
INSTALL(FILES src/leviosam.hpp DESTINATION include)

