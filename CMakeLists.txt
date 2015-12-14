cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set (CMAKE_DISABLE_SOURCE_CHANGES  ON)

project (epa C CXX)

# set (epa_VERSION_MAJOR 0)
# set (epa_VERSION_MINOR 1)

# set (CMAKE_BUILD_TYPE DEBUG)
set (CMAKE_BUILD_TYPE RELEASE)
# set (CMAKE_VERBOSE_MAKEFILE ON)

message (STATUS "Building ${CMAKE_BUILD_TYPE}")

set (WARN_FLAGS              "-Wall -Wextra -pedantic -pedantic-errors")
set (CMAKE_CXX_FLAGS         "-std=c++11 ${WARN_FLAGS}")

set (CMAKE_CXX_FLAGS_DEBUG   "-g")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")

#find_library(libpll)

# include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

enable_testing()
add_subdirectory(${PROJECT_SOURCE_DIR}/test/src)
