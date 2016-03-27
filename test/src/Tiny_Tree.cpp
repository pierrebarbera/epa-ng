#include "Epatest.hpp"

#include "src/pllhead.hpp"
#include "src/file_io.hpp"
#include "src/Tree_Numbers.hpp"
#include "src/Model.hpp"
#include "src/Tiny_Tree.hpp"
#include "src/MSA.hpp"
#include "src/Range.hpp"
#include "src/pll_util.hpp"
#include "src/epa_pll_util.hpp"

#include <tuple>

using namespace std;

TEST(Tiny_Tree, place_heuristic)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  MSA queries = build_MSA_from_file(env->query_file);
  Tree_Numbers nums = Tree_Numbers();
  Range range;
  range.span = msa.num_sites();
  pll_partition_t * part;
  pll_utree_t * tree;

  tree = build_tree_from_file(env->tree_file, nums);
  part = build_partition_from_file( env->model, nums, msa.num_sites());

  auto valid_map = vector<Range>(nums.tip_nodes);
  link_tree_msa(tree, part, msa, nums.tip_nodes, valid_map);
  precompute_clvs(tree, part, nums);

  // tests
  Tiny_Tree tt(tree, 0, part, env->model, false);

  for (auto const &x : queries)
  {
    auto place = tt.place(x);
    EXPECT_NE(place.likelihood(), 0.0);
    EXPECT_NE(place.distal_length(), 0.0);
    EXPECT_NE(place.pendant_length(), 0.0);
  }

  // teardown
  pll_partition_destroy(part);
  pll_utree_destroy(tree);
}

TEST(Tiny_Tree, place_BLO)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  MSA queries = build_MSA_from_file(env->query_file);
  Tree_Numbers nums = Tree_Numbers();
  pll_partition_t * part;
  pll_utree_t * tree;
  Range range;
  range.span = msa.num_sites();

  tree = build_tree_from_file(env->tree_file, nums);
  part = build_partition_from_file( env->model, nums, msa.num_sites());

  auto valid_map = vector<Range>(nums.tip_nodes);
  link_tree_msa(tree, part, msa, nums.tip_nodes, valid_map);
  precompute_clvs(tree, part, nums);

  // tests
  Tiny_Tree tt(tree, 0, part, env->model, true);

  for (auto const &x : queries)
  {
    auto place = tt.place(x);
    EXPECT_NE(place.likelihood(), 0.0);
    EXPECT_NE(place.distal_length(), 0.0);
    EXPECT_NE(place.pendant_length(), 0.0);
  }

  // teardown
  pll_partition_destroy(part);
  pll_utree_destroy(tree);
}
