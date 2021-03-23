cmake_minimum_required(VERSION 3.10)

project(leviosam CXX C)

set(PROJECT_URL "https://github.com/alshai/leviosam")
set(PROJECT_DESCRIPTION "leviosam: library for moving alignments between aligned genomes")
set(CMAKE_CXX_STANDARD 11)


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

add_executable(leviosam leviosam.cpp bam_md.c bam_aux.c)
target_link_libraries(leviosam ${PKG_HTSLIB_LDFLAGS})
target_link_libraries(leviosam sdsl)
target_link_libraries(leviosam pthread)
target_link_libraries(leviosam z)


add_executable(leviosam_test leviosam_test.cpp)
target_link_libraries(leviosam_test ${PKG_HTSLIB_LDFLAGS})
target_link_libraries(leviosam_test sdsl)
target_link_libraries(leviosam_test gtest gtest_main)

add_test(NAME leviosam_test COMMAND leviosam_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testdata
)

INSTALL(TARGETS leviosam DESTINATION bin)


