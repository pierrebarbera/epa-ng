#include "Epatest.hpp"

#include "src/pllhead.hpp"
#include "src/pll_util.hpp"
#include "src/epa_pll_util.hpp"
#include "src/file_io.hpp"
#include "src/Tree_Numbers.hpp"
#include "src/Model.hpp"
#include "src/MSA.hpp"
#include "src/Range.hpp"

#include <string>
#include <tuple>

using namespace std;

TEST(pll_util, utree_query_branches)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  Tree_Numbers nums = Tree_Numbers();
  pll_partition_t * part;
  pll_utree_t * tree;

  tie(part, tree) = build_partition_from_file(env->tree_file, env->model, nums, msa.num_sites());

  // tests
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

  // teardown
  pll_partition_destroy(part);
  pll_utree_destroy(tree);
}

TEST(pll_util, set_unique_clv_indices)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  Tree_Numbers nums = Tree_Numbers();
  pll_partition_t * part;
  pll_utree_t * tree;

  tie(part, tree) = build_partition_from_file(env->tree_file, env->model, nums, msa.num_sites());
  set_unique_clv_indices(tree, nums.tip_nodes);

  // tests
  vector<pll_utree_t *> inner_nodes(nums.inner_nodes);
  pll_utree_query_innernodes(tree, &inner_nodes[0]);
  vector<pll_utree_t *> tip_nodes(nums.tip_nodes);
  pll_utree_query_tipnodes(tree, &tip_nodes[0]);

  // check for duplicate clv indices
  vector<unsigned int> index_list;
  for (auto x : tip_nodes)
    index_list.push_back(x->clv_index);

  for (auto x : inner_nodes) {
    index_list.push_back(x->clv_index);
    index_list.push_back(x->next->clv_index);
    index_list.push_back(x->next->next->clv_index);
  }

  for (auto x : index_list) {
    int count = 0;
    for (auto y : index_list) {
      if (x == y)
        count++;
    }
    EXPECT_EQ(count, 1);
  }

  // teardown
  pll_partition_destroy(part);
  pll_utree_destroy(tree);
}

static int cb_set_branchlengths_one(pll_utree_t * node)
{
  node->length = 1;
  if(node->next)
  {
    node->next->length = 1;
    node->next->next->length = 1;
  }
  return 1;
};

TEST(pll_util, get_numbered_newick_string)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  Tree_Numbers nums = Tree_Numbers();
  pll_partition_t * part;
  pll_utree_t * tree;

  tie(part, tree) = build_partition_from_file(env->tree_file, env->model, nums, msa.num_sites());
  auto valid_map = vector<Range>(nums.tip_nodes);
  link_tree_msa(tree, part, msa, nums.tip_nodes, valid_map);

  // tests
  // valid output as returned by RAxML, with reset branch lengths, as we only want to test format
  string valid(
  "(Seal:1{0},(Whale:1{1},(Mouse:1{2},(Human:1{3},(Chicken:1{4},(Frog:1{5},Loach:1{6}):1{7}):1{8}):1{9}):1{10}):1{11},Cow:1{12});");

  vector<pll_utree_t*> travbuffer(nums.nodes);
  unsigned int traversal_size;
  pll_utree_traverse(tree, cb_set_branchlengths_one, &travbuffer[0], &traversal_size);

  auto ret_string =  get_numbered_newick_string(tree);

  EXPECT_STREQ(valid.c_str(), ret_string.c_str());

  // teardown

  pll_partition_destroy(part);
  pll_utree_destroy(tree);
}
