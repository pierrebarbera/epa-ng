file (GLOB_RECURSE epa_sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

include_directories (${pll_dir} ${pll_binary_dir} ${pll_optimize_dir} ${pll_msa_dir})

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_executable        (epa_module ${epa_sources})

add_dependencies(epa_module pll pll_binary pll_optimize pll_msa)

target_link_libraries (epa_module m)
target_link_libraries (epa_module ${pll_dir}/libpll.so)
target_link_libraries (epa_module ${pll_binary_dir}/libpll_binary.so)
target_link_libraries (epa_module ${pll_optimize_dir}/libpll_optimize.so)
target_link_libraries (epa_module ${pll_msa_dir}/libpll_msa.so)

target_link_libraries (epa_module ${MPI_CXX_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(epa_module PROPERTIES
  COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(epa_module PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

set_target_properties (epa_module PROPERTIES OUTPUT_NAME epa)
set_target_properties (epa_module PROPERTIES PREFIX "")
