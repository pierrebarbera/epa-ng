#include "Epatest.hpp"

#include "src/pllhead.hpp"
#include "src/epa_pll_util.hpp"
#include "src/pll_util.hpp"
#include "src/file_io.hpp"
#include "src/Tree_Numbers.hpp"
#include "src/Model.hpp"
#include "src/Options.hpp"
#include "src/Tree.hpp"
#include "src/MSA.hpp"
#include "src/Range.hpp"

#include <string>

using namespace std;

TEST(epa_pll_util, link_tree_msa)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  Tree_Numbers nums = Tree_Numbers();
  pll_partition_t * part;
  pll_utree_t * tree;
  Model model;

  tree = build_tree_from_file(env->tree_file, nums);
  part = build_partition_from_file( env->model, nums, msa.num_sites());
  auto valid_map = vector<Range>(nums.tip_nodes);
  link_tree_msa(tree, part, model, msa, nums.tip_nodes, valid_map);

  // tests
  vector<pll_unode_t *> tip_nodes(nums.tip_nodes);
  tip_nodes.assign(tree->nodes, tree->nodes + nums.tip_nodes);

  for (auto n : tip_nodes) {
    ASSERT_NE(n, nullptr);
    if (part->attributes & PLL_ATTRIB_PATTERN_TIP)
      EXPECT_NE(part->tipchars[n->clv_index][0], 0);
    else
    {
      EXPECT_NE(part->clv[n->clv_index][0], 0);
    }
  }

  // teardown
  pll_partition_destroy(part);
  pll_utree_destroy(tree, nullptr);

}

TEST(epa_pll_util, precompute_clvs)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  Tree_Numbers nums = Tree_Numbers();
  pll_partition_t * part;
  pll_utree_t * tree;
  Model model;

  tree = build_tree_from_file(env->tree_file, nums);
  part = build_partition_from_file( env->model, nums, msa.num_sites());
  auto valid_map = vector<Range>(nums.tip_nodes);
  link_tree_msa(tree, part, model, msa, nums.tip_nodes, valid_map);
  precompute_clvs(tree, part, nums);

  // tests
  vector<pll_unode_t *> node_list(nums.branches);
  utree_query_branches(tree, &node_list[0]);

  // all edge logl should be the same
  auto first = true;
  double log_old = 0.0;
  double log_new;
  unsigned int param_indices[4] = {0};
  for (auto x : node_list)
  {
    log_new = pll_compute_edge_loglikelihood(part,
                                         x->clv_index,
                                         x->scaler_index,
                                         x->back->clv_index,
                                         x->back->scaler_index,
                                         x->pmatrix_index,
                                         param_indices, nullptr);
    if (!first)
    {
      EXPECT_DOUBLE_EQ(log_old, log_new);
    }

    log_old = log_new;
    first = false;
  }

  EXPECT_NE(log_new, 0.0);

  // teardown
  utree_free_node_data(get_root(tree));
  pll_partition_destroy(part);
  pll_utree_destroy(tree, nullptr);
}

TEST(epa_pll_util, split_combined_msa)
{
  // buildup
  auto combined_msa = build_MSA_from_file(env->combined_file);
  Model model;
  Options options;

  Tree tree(env->tree_file, combined_msa, model, options);

  // tests
  MSA query_msa;
  split_combined_msa(combined_msa, query_msa, tree);

  EXPECT_EQ(combined_msa.num_sites(), query_msa.num_sites());
  EXPECT_EQ(combined_msa.size(), 8);
  EXPECT_EQ(query_msa.size(), 2);

  for (auto x : combined_msa) {
    for (auto y : query_msa) {
      EXPECT_STRNE(x.header().c_str(), y.header().c_str());
    }
  }

  // teardown
}
