#include "Epatest.hpp"

#include "src/pllhead.hpp"
#include "src/pll_util.hpp"
#include "src/file_io.hpp"
#include "src/Tree_Numbers.hpp"
#include "src/Model.hpp"
#include "src/MSA.hpp"

#include <string>
#include <tuple>

using namespace std;

TEST(pll_util, utree_query_branches)
{

  MSA msa = build_MSA_from_file(env->reference_file);
  Tree_Numbers nums = Tree_Numbers();
  Model model = Model({0.25, 0.25, 0.25, 0.25}, {1,1,1,1,1,1}, 1.0);
  pll_partition_t * part;
  pll_utree_t * tree;

  tie(part, tree) = build_partition_from_file(env->tree_file, model, nums, msa.num_sites());

  vector<pll_utree_t *> node_list(nums.branches);
  auto num_traversed = utree_query_branches(tree, &node_list[0]);

  EXPECT_EQ(nums.branches, num_traversed);

  for (auto x : node_list) {
    EXPECT_NE(x, nullptr);
  }

  // check for duplicates
  for (auto x : node_list) {
    int count = 0;
    for (auto y : node_list) {
      if (x == y || x == y->back)
        count++;
    }
    EXPECT_EQ(count, 1);
  }

  pll_partition_destroy(part);
  pll_utree_destroy(tree);
}
