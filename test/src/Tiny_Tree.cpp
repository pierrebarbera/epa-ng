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

  Options options;
  auto ref_tree = Tree(env->tree_file, msa, env->model, options);

  auto tree = ref_tree.tree();

  // tests
  Tiny_Tree tt(tree, 0, ref_tree, false);

  for (auto const &x : queries)
  {
    auto place = tt.place(x);
    EXPECT_NE(place.likelihood(), 0.0);
    EXPECT_NE(place.distal_length(), 0.0);
    EXPECT_NE(place.pendant_length(), 0.0);
  }

  // teardown
}

TEST(Tiny_Tree, place_BLO)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  MSA queries = build_MSA_from_file(env->query_file);
  Tree_Numbers nums = Tree_Numbers();

  Range range;
  range.span = msa.num_sites();

  Options options;
  auto ref_tree = Tree(env->tree_file, msa, env->model, options);

  auto tree = ref_tree.tree();

  // tests
  Tiny_Tree tt(tree, 0, ref_tree, true);

  for (auto const &x : queries)
  {
    auto place = tt.place(x);
    EXPECT_NE(place.likelihood(), 0.0);
    EXPECT_NE(place.distal_length(), 0.0);
    EXPECT_NE(place.pendant_length(), 0.0);
  }

  // teardown
}
