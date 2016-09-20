cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set (CMAKE_DISABLE_SOURCE_CHANGES  ON)


# set( ENV{OMPI_CXX}            "clang++" PARENT_SCOPE )
# set( CMAKE_CXX_COMPILER       "mpic++")

project (epa C CXX)

#find_package(MPI)
#include_directories(${MPI_INCLUDE_PATH})

# set (epa_VERSION_MAJOR 0)
# set (epa_VERSION_MINOR 1)

set (CMAKE_BUILD_TYPE DEBUG)
# set (CMAKE_BUILD_TYPE RELEASE)
# set (CMAKE_VERBOSE_MAKEFILE ON)

message (STATUS "Building ${CMAKE_BUILD_TYPE}")

set (WARN_FLAGS               "-Wall -Wextra")
set (CMAKE_CXX_FLAGS          "-std=c++11 ${WARN_FLAGS}")

set (CMAKE_CXX_FLAGS_DEBUG    "-g")
set (CMAKE_CXX_FLAGS_RELEASE  "-O3")
#set (MPI_COMPILE_FLAGS        "-D__MPI")

include_directories(${PROJECT_SOURCE_DIR}/../cereal/include)


set(pll_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/libs/libpll/src)
set(pll_binary_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/src/binary)
set(pll_optimize_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/src/optimize)
set(pll_msa_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/src/msa)

add_custom_target(
  pll
  COMMAND make
  WORKING_DIRECTORY ${pll_dir}
)

add_custom_target(
  pll_binary
  COMMAND make
  WORKING_DIRECTORY ${pll_binary_dir}
)

add_custom_target(
  pll_optimize
  COMMAND make
  WORKING_DIRECTORY ${pll_optimize_dir}
)

add_custom_target(
  pll_msa
  COMMAND make
  WORKING_DIRECTORY ${pll_msa_dir}
)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

enable_testing()
add_subdirectory(${PROJECT_SOURCE_DIR}/test/src)
