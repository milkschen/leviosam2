cmake_minimum_required(VERSION 3.10)

project(leviosam CXX C)

set(PROJECT_URL "https://github.com/alshai/leviosam")
set(PROJECT_DESCRIPTION "leviosam: library for moving alignments between aligned genomes")
set(CMAKE_CXX_STANDARD 11)

find_library(HTS_LIB hts)
if(NOT HTS_LIB)
    message(FATAL_ERROR "htslib library not found. Please specify -D CMAKE_LIBRARY_PATH and -D CMAKE_INCLUDE_PATH in the CMake command")
    # message(FATAL_ERROR "htslib library not found. Please specify -D HTS_LIB_DIR and -D HTS_INC_DIR in the CMake command")
endif()
message("${HTS_LIB}")

find_library(SDSL_LIB sdsl)
if(NOT SDSL_LIB)
    message(FATAL_ERROR "libsdsl not found. Please specify -D CMAKE_LIBRARY_PATH and -D CMAKE_INCLUDE_PATH in the CMake command")
    # message(FATAL_ERROR "sdsl-lite library not found (lib_path= ${SDSL_LIB_DIR}, inc_path= ${SDSL_INC_DIR}). Please specify -D SDSL_LIB_DIR and -D SDSL_INC_DIR in the CMake command")
endif()


include_directories(${CMAKE_INCLUDE_PATH})

add_executable(leviosam src/leviosam.cpp src/chain.cpp src/leviosam_utils.cpp src/bam_md.c src/bam_aux.c)
target_link_libraries(leviosam ${HTS_LIB})
target_link_libraries(leviosam ${SDSL_LIB})
target_link_libraries(leviosam pthread)
target_link_libraries(leviosam z)

add_executable(extract_unpaired src/extract_unpaired.cpp)
target_link_libraries(extract_unpaired ${HTS_LIB})
target_link_libraries(extract_unpaired ${SDSL_LIB})
target_link_libraries(extract_unpaired pthread)
target_link_libraries(extract_unpaired z)

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

add_executable(leviosam_test src/leviosam_test.cpp src/chain.cpp src/leviosam_utils.cpp)
target_link_libraries(leviosam_test ${HTS_LIB})
target_link_libraries(leviosam_test ${SDSL_LIB})
target_link_libraries(leviosam_test gtest gtest_main)

add_test(NAME leviosam_test COMMAND leviosam_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testdata
)

INSTALL(TARGETS leviosam DESTINATION bin)
INSTALL(FILES src/leviosam.hpp DESTINATION include)

