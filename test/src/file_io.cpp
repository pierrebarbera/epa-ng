#include "Epatest.hpp"

#include "src/pllhead.hpp"
#include "src/file_io.hpp"
#include "src/Tree_Numbers.hpp"
#include "src/Model.hpp"
#include "src/MSA.hpp"

#include <string>
#include <tuple>

using namespace std;

TEST(file_io, build_partition_from_file)
{
  MSA msa = build_MSA_from_file(env->reference_file);
  Tree_Numbers nums = Tree_Numbers();
  Model model = Model({0.25, 0.25, 0.25, 0.25}, {1,1,1,1,1,1}, 1.0);
  pll_partition_t * part;
  pll_utree_t * tree;

  tie(part, tree) = build_partition_from_file(env->tree_file, model, nums, msa.num_sites());

  EXPECT_EQ(nums.tip_nodes, 8);
  EXPECT_EQ(nums.nodes, 14);
  EXPECT_EQ(nums.inner_nodes, 6);
  EXPECT_EQ(nums.branches, 13);

  pll_partition_destroy(part);
}
