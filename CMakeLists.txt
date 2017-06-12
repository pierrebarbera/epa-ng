cmake_minimum_required (VERSION 2.8.10 FATAL_ERROR)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)

set (ENABLE_MPI       OFF)
set (ENABLE_OMP       OFF)
set (ENABLE_PREFETCH  OFF)

project (epa CXX)

#check for minimum compiler versions
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(EPA_COMPILER_TARGET_VERSION "4.9")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS EPA_COMPILER_TARGET_VERSION)
    message (FATAL_ERROR "${CMAKE_CXX_COMPILER_ID} compiler too old! Minimum required: ${EPA_COMPILER_TARGET_VERSION}")
  endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(EPA_COMPILER_TARGET_VERSION "3.8")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS EPA_COMPILER_TARGET_VERSION)
    message (FATAL_ERROR "${CMAKE_CXX_COMPILER_ID} compiler too old! Minimum required: ${EPA_COMPILER_TARGET_VERSION}")
  endif()
endif()

# set (epa_VERSION_MAJOR 0)
# set (epa_VERSION_MINOR 1)

# set (CMAKE_BUILD_TYPE DEBUG)
set (CMAKE_BUILD_TYPE RELEASE)
# set (CMAKE_VERBOSE_MAKEFILE ON)

message (STATUS "Building ${CMAKE_BUILD_TYPE}")

set (WARN_FLAGS               "-Wall -Wextra")
set (CMAKE_CXX_FLAGS          "-std=c++14 ${WARN_FLAGS}")

set (CMAKE_CXX_FLAGS_DEBUG    "-g")
set (CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG")

if(ENABLE_PREFETCH)
find_package (Threads)
message(STATUS "Enabling Prefetching")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__PREFETCH")
endif()

if(ENABLE_OMP)
message(STATUS "Checking for OpenMP")
find_package( OpenMP QUIET )
if(OPENMP_FOUND)
message(STATUS "Checking for OpenMP -- found")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -D__OMP")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
message(STATUS "Checking for OpenMP -- NOT FOUND")
endif()
endif()

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-mavx HAS_AVX)
CHECK_CXX_COMPILER_FLAG(-msse3 HAS_SSE3)
if(HAS_AVX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -D__AVX")
elseif(HAS_SSE3)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3 -D__SSE3")
endif()

if(ENABLE_MPI)
find_package(MPI)
if(MPI_CXX_FOUND)
set( ENV{OMPI_CXX} ${CMAKE_CXX_COMPILER} PARENT_SCOPE )
#set( CMAKE_CXX_COMPILER       "mpicxx")
include_directories(${MPI_CXX_INCLUDE_PATH})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILER_FLAGS} -D__MPI")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
endif()
endif()

message(STATUS "Using flags: ${CMAKE_CXX_FLAGS}")

include_directories(${PROJECT_SOURCE_DIR}/libs/cereal/include)
include_directories(${PROJECT_SOURCE_DIR}/libs/cxxopts/src)
include_directories(${PROJECT_SOURCE_DIR}/libs/include/libpll)


set(pll_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/libs/libpll/src)
set(pll_binary_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/src/binary)
set(pll_optimize_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/src/optimize)
set(pll_msa_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/src/msa)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

enable_testing()
add_subdirectory(${PROJECT_SOURCE_DIR}/test/src)
