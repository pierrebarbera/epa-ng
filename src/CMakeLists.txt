file (GLOB_RECURSE epa_sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_executable        (epa_module ${epa_sources})

target_link_libraries (epa_module pll)

set_target_properties (epa_module PROPERTIES OUTPUT_NAME epa)
set_target_properties (epa_module PROPERTIES PREFIX "")
