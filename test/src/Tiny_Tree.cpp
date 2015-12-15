#include "Epatest.hpp"

#include "src/pllhead.hpp"
#include "src/file_io.hpp"
#include "src/Tree_Numbers.hpp"
#include "src/Model.hpp"
#include "src/Tiny_Tree.hpp"
#include "src/MSA.hpp"

#include <tuple>

using namespace std;

TEST(Tiny_Tree, place_heuristic)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  MSA queries = build_MSA_from_file(env->query_file);
  Tree_Numbers nums = Tree_Numbers();
  pll_partition_t * part;
  pll_utree_t * tree;

  tie(part, tree) = build_partition_from_file(env->tree_file, env->model, nums, msa.num_sites());

  // tests
  Tiny_Tree tt(tree, part, env->model, true);

  double logl, pendant, distal;
  for (auto const &x : queries)
  {
    tie(logl, distal, pendant) = tt.place(x);
    EXPECT_NE(logl, 0.0);
    EXPECT_NE(distal, 0.0);
    EXPECT_NE(pendant, 0.0);
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

  tie(part, tree) = build_partition_from_file(env->tree_file, env->model, nums, msa.num_sites());

  // tests
  Tiny_Tree tt(tree, part, env->model, true);

  double logl, pendant, distal;
  for (auto const &x : queries)
  {
    tie(logl, distal, pendant) = tt.place(x);
    EXPECT_NE(logl, 0.0);
    EXPECT_NE(distal, 0.0);
    EXPECT_NE(pendant, 0.0);
  }


  // teardown
  pll_partition_destroy(part);
  pll_utree_destroy(tree);
}
