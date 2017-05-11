find_package (GTest)

if(NOT GTEST_FOUND)
    message (STATUS "GTest not found")
    message (WARNING "Skipping building tests.")
    return()
endif()

message (STATUS "Building tests")

file (GLOB_RECURSE epa_test_sources ${PROJECT_SOURCE_DIR}/test/src/*.cpp ${PROJECT_SOURCE_DIR}/src/*.cpp)

# sources list now has 2 Main.cpp, old has to be removed
list(REMOVE_ITEM epa_test_sources "${PROJECT_SOURCE_DIR}/src/main.cpp")

include_directories (${PROJECT_SOURCE_DIR})

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/test/bin)

add_executable        (epa_test_module ${epa_test_sources})

target_link_libraries (epa_test_module ${PROJECT_SOURCE_DIR}/libs/lib/libpll_binary.a)
target_link_libraries (epa_test_module ${PROJECT_SOURCE_DIR}/libs/lib/libpll_optimize.a)
target_link_libraries (epa_test_module ${PROJECT_SOURCE_DIR}/libs/lib/libpll_msa.a)
target_link_libraries (epa_test_module ${PROJECT_SOURCE_DIR}/libs/lib/libpll_util.a)
target_link_libraries (epa_test_module ${PROJECT_SOURCE_DIR}/libs/lib/libpll.a)
target_link_libraries (epa_test_module m)

target_link_libraries (epa_test_module ${GTEST_BOTH_LIBRARIES} -pthread)

target_link_libraries (epa_test_module ${MPI_CXX_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(epa_test_module PROPERTIES
  COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(epa_test_module PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

set_target_properties (epa_test_module PROPERTIES OUTPUT_NAME epa_test)
set_target_properties (epa_test_module PROPERTIES PREFIX "")


add_test (epa_test ${PROJECT_SOURCE_DIR}/test/bin/epa_test)
