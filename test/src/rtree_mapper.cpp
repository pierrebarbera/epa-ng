#include "Epatest.hpp"

#include "core/pll/rtree_mapper.hpp"

#include "io/file_io.hpp"
#include "tree/Tree_Numbers.hpp"
#include "core/raxml/Model.hpp"

#include <vector>
#include <string>

using placepair = std::pair<unsigned int, double>;

static void test_placement_mapping( std::string const& tree_file,
                                    std::vector<placepair> const& utree_placements,
                                    std::vector<placepair> const& rtree_placements )
{
  // buildup
  pll_utree_t * tree;
  raxml::Model model;

  rtree_mapper mapper;
  Tree_Numbers nums;
  tree = build_tree_from_file( tree_file, nums, mapper );

  assert( utree_placements.size() == rtree_placements.size() );

  for (size_t i = 0; i < utree_placements.size(); ++i) {

    auto const& u_p = utree_placements[i];
    auto const& r_p = rtree_placements[i];

    // translate
    unsigned int edge = 0;
    double distal = -1.0;

    std::tie(edge, distal) = mapper.in_rtree( u_p.first, u_p.second );

    EXPECT_EQ(    r_p.first,  edge );
    EXPECT_NEAR(  r_p.second, distal, 1e-10 );
  }

  // teardown
  pll_utree_destroy(tree, nullptr);
}

TEST(rtree_mapper, translate_distal)
{
  /*  basic idea of this test: init the mapper by different rooted trees
      for which we know the resulting distal lengths of placements
  */

  test_placement_mapping(
    env->tree_file_rooted,
    { {0, 1.0}, {7, 1.5}, {8, 0.5},   {8, 0.001}, },
    { {0, 1.0}, {8, 1.5}, {6, 0.52},  {9, 0.001}, }
  );

  test_placement_mapping(
    env->tree_file_rooted_2,
    { {0, 1.34}, {0, 1.345}, {8, 0.5},  {2, 0.001}, },
    { {0, 1.34}, {9, 0.005}, {8, 0.5},  {2, 0.001}, }
  );
}
