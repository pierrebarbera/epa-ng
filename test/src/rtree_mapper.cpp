#include "Epatest.hpp"

#include "core/pll/rtree_mapper.hpp"

#include "core/raxml/Model.hpp"
#include "io/file_io.hpp"
#include "tree/Tree_Numbers.hpp"

#include <string>
#include <vector>

using placepair = std::pair< unsigned int, double >;

static void test_placement_mapping( std::string const& tree_file,
                                    std::vector< placepair > const& utree_placements,
                                    std::vector< placepair > const& rtree_placements )
{
  // buildup
  pll_utree_t* tree;
  raxml::Model model;

  rtree_mapper mapper;
  Tree_Numbers nums;
  tree = build_tree_from_file( tree_file, nums, mapper );

  assert( utree_placements.size() == rtree_placements.size() );

  for( size_t i = 0; i < utree_placements.size(); ++i ) {

    auto const& u_p = utree_placements[ i ];
    auto const& r_p = rtree_placements[ i ];

    // translate
    unsigned int edge = 0;
    double distal     = -1.0;

    std::tie( edge, distal ) = mapper.in_rtree( u_p.first, u_p.second );

    EXPECT_EQ( r_p.first, edge );
    EXPECT_NEAR( r_p.second, distal, 1e-10 );
  }

  // printf("utree_root_edge_: %u\n", mapper.utree_root_edge_);
  // printf("rtree_proximal_edge_: %u\n", mapper.rtree_proximal_edge_);
  // printf("rtree_distal_edge_: %u\n", mapper.rtree_distal_edge_);

  // printf("map:\n");
  // auto const& map = mapper.map();
  // for (unsigned int i = 0; i < map.size(); ++i) {
  //   printf("%u %u\n", map[i], i);
  // }

  // teardown
  pll_utree_destroy( tree, nullptr );
}

/*  basic idea of this test: init the mapper by different rooted trees
    for which we know the resulting distal lengths of placements
*/
TEST( rtree_mapper, test_mapping_1 )
{
  test_placement_mapping(
      env->tree_file_rooted,
      {
          { 8, 1.0 },
          { 8, 1.5 },
          { 6, 0.5 },
          { 7, 0.001 },
      },
      {
          { 9, 1.0 },
          { 6, 0.63 },
          { 7, 0.5 },
          { 8, 0.001 },
      } );
}
TEST( rtree_mapper, test_mapping_2 )
{
  test_placement_mapping(
      env->tree_file_rooted_2,
      {
          { 0, 1.34 },
          { 0, 1.345 },
          { 8, 0.5 },
          { 2, 0.001 },
      },
      {
          { 0, 1.34 },
          { 9, 0.005 },
          { 8, 0.5 },
          { 2, 0.001 },
      } );
}
TEST( rtree_mapper, test_mapping_3 )
{
  test_placement_mapping(
      env->tree_file_rooted_3,
      {
          { 8, 0.5 },
          { 8, 0.005 },
          { 0, 0.5 },
          { 2, 0.001 },
      },
      {
          { 8, 1.41 },
          { 9, 0.005 },
          { 0, 0.5 },
          { 2, 0.001 },
      } );
}

//todo general rooting mapping test