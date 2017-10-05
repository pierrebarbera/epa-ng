cmake_minimum_required (VERSION 2.8.7 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/libs/seqan/util/cmake")

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)

# main compile alternators and their defaults
set (ENABLE_MPI       OFF)
set (ENABLE_OMP       ON)
set (ENABLE_PREFETCH  ON)

if ( $ENV{EPA_HYBRID} )
	set (ENABLE_MPI       ON)
	set (ENABLE_OMP       ON)
	set (ENABLE_PREFETCH  ON)
endif ()

if ( $ENV{EPA_SERIAL} )
    set (ENABLE_MPI       OFF)
    set (ENABLE_OMP       OFF)
    set (ENABLE_PREFETCH  OFF)
endif ()

project ( epa CXX C )

set (epa_VERSION_MAJOR 0)
set (epa_VERSION_MINOR 1)
set (epa_VERSION_PATCH 0)

set (epa_PROJECT_VERSION "${epa_VERSION_MAJOR}")
set (epa_PROJECT_VERSION "${epa_PROJECT_VERSION}.${epa_VERSION_MINOR}")
set (epa_PROJECT_VERSION "${epa_PROJECT_VERSION}.${epa_VERSION_PATCH}")

#check for minimum compiler versions
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(EPA_COMPILER_TARGET_VERSION "4.9.2")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS EPA_COMPILER_TARGET_VERSION)
    message (FATAL_ERROR "${CMAKE_CXX_COMPILER_ID} compiler too old! Minimum required: ${EPA_COMPILER_TARGET_VERSION}")
  endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(EPA_COMPILER_TARGET_VERSION "3.8")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS EPA_COMPILER_TARGET_VERSION)
    message (FATAL_ERROR "${CMAKE_CXX_COMPILER_ID} compiler too old! Minimum required: ${EPA_COMPILER_TARGET_VERSION}")
  endif()
endif()

message(STATUS "EPA-ng version: ${epa_PROJECT_VERSION}")

set (CMAKE_BUILD_TYPE RELEASE)

if ( $ENV{EPA_DEBUG} )
  set (CMAKE_BUILD_TYPE DEBUG)
endif ()

message (STATUS "Building ${CMAKE_BUILD_TYPE}")

set (WARN_FLAGS               "-Wall -Wextra")
set (CMAKE_CXX_FLAGS          "-std=c++14 ${WARN_FLAGS}")

set (CMAKE_CXX_FLAGS_DEBUG    "-O0 -g -ggdb3 -DDEBUG -D_GLIBCXX_DEBUG")
set (CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG")

if(ENABLE_PREFETCH)
  set (CMAKE_THREAD_PREFER_PTHREAD ON)
  set (THREADS_PREFER_PTHREAD_FLAG ON)
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
  find_package(MPI REQUIRED)
  if(MPI_CXX_FOUND)
    set( ENV{OMPI_CXX} ${CMAKE_CXX_COMPILER} PARENT_SCOPE )
    #set( CMAKE_CXX_COMPILER       "mpicxx")
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILER_FLAGS} -D__MPI")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
  endif()
endif()

# pass special info to genesis
if(ENABLE_PREFETCH)
  set(GENESIS_USE_THREADS ON CACHE BOOL "Passing GENESIS_USE_THREADS=ON to genesis")
else()
  set(GENESIS_USE_THREADS OFF CACHE BOOL "Passing GENESIS_USE_THREADS=OFF to genesis")
endif()

set(GENESIS_USE_OPENMP OFF CACHE BOOL "Passing GENESIS_USE_OPENMP=OFF to genesis")

# Add Genesis as dependency
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/genesis)

# Use everything that Genesis exports, just to be sure that we use the same setup.
add_definitions( ${GENESIS_DEFINITIONS} )
include_directories( ${GENESIS_INCLUDE_DIR} )
set( CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${GENESIS_CXX_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GENESIS_EXE_LINKER_FLAGS}" )

message(STATUS "Using flags: ${CMAKE_CXX_FLAGS}")

# Add include directories.
include_directories(${PROJECT_SOURCE_DIR}/libs/cereal/include)
include_directories(${PROJECT_SOURCE_DIR}/libs/cxxopts/include)
include_directories(${PROJECT_SOURCE_DIR}/libs/include/libpll)

set(pll_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/libs/libpll/src)
set(pll_binary_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/src/binary)
set(pll_optimize_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/src/optimize)
set(pll_msa_dir ${PROJECT_SOURCE_DIR}/libs/pll-modules/src/msa)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

enable_testing()
add_subdirectory(${PROJECT_SOURCE_DIR}/test/src)
