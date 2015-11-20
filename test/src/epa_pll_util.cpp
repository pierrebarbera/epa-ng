#include "Epatest.hpp"

#include "src/pllhead.hpp"
#include "src/epa_pll_util.hpp"
#include "src/pll_util.hpp"
#include "src/file_io.hpp"
#include "src/Tree_Numbers.hpp"
#include "src/Model.hpp"
#include "src/MSA.hpp"

#include <string>

using namespace std;

TEST(epa_pll_util, link_tree_msa)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  Tree_Numbers nums = Tree_Numbers();
  Model model = Model({0.25, 0.25, 0.25, 0.25}, {1,1,1,1,1,1}, 1.0);
  pll_partition_t * part;
  pll_utree_t * tree;

  tie(part, tree) = build_partition_from_file(env->tree_file, model, nums, msa.num_sites());
  link_tree_msa(tree, part, msa, nums.tip_nodes);

  // tests
  vector<pll_utree_t *> tip_nodes(nums.tip_nodes);
  auto size = pll_utree_query_tipnodes(tree, &tip_nodes[0]);

  EXPECT_EQ(size, 8);

  for (auto n : tip_nodes) {
    ASSERT_NE(n, nullptr);
    EXPECT_NE(part->clv[n->clv_index][0], 0.0);
  }

  // teardown
  pll_partition_destroy(part);
  pll_utree_destroy(tree);

}

TEST(epa_pll_util, precompute_clvs)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  Tree_Numbers nums = Tree_Numbers();
  Model model = Model({0.25, 0.25, 0.25, 0.25}, {1,1,1,1,1,1}, 1.0);
  pll_partition_t * part;
  pll_utree_t * tree;

  tie(part, tree) = build_partition_from_file(env->tree_file, model, nums, msa.num_sites());
  link_tree_msa(tree, part, msa, nums.tip_nodes);
  precompute_clvs(tree, part, nums);

  // tests
  vector<pll_utree_t *> node_list(nums.branches);
  utree_query_branches(tree, &node_list[0]);

  // all edge logl should be the same
  double log_old = 0.0;
  double log_new;
  for (auto x : node_list)
  {
    log_new = pll_compute_edge_loglikelihood(part,
                                         x->clv_index,
                                         x->scaler_index,
                                         x->back->clv_index,
                                         x->back->scaler_index,
                                         x->pmatrix_index,
                                         0);
    if (log_old != 0.0)
      EXPECT_DOUBLE_EQ(log_old, log_new);

    log_old = log_new;

  }

  EXPECT_NE(log_new, 0.0);

  // teardown
  utree_free_node_data(tree);
  pll_partition_destroy(part);
  pll_utree_destroy(tree);
}

void check_clvs_computed(pll_partition_t * p, pll_utree_t * t)
{
  int tip_nodes_count = 8;
  int inner_nodes_count = 6;
  int i;
  pll_utree_t ** tipnodes = (pll_utree_t  **)calloc(tip_nodes_count,
                                                    sizeof(pll_utree_t *));
  pll_utree_query_tipnodes(t, tipnodes);

  pll_utree_t ** innernodes = (pll_utree_t  **)calloc(inner_nodes_count,
                                                    sizeof(pll_utree_t *));

  pll_utree_query_innernodes(t, innernodes);

  double logl = 0.0;
  for (i = 0; i <tip_nodes_count + inner_nodes_count; ++i)
  {
    if (i < tip_nodes_count)
    {
      t = tipnodes[i]->back;
    }
    else
    {
      t = innernodes[i - tip_nodes_count];
    }

    logl = pll_compute_edge_loglikelihood(p,
                                         t->clv_index,
                                         t->scaler_index,
                                         t->back->clv_index,
                                         t->back->scaler_index,
                                         t->pmatrix_index,
                                         0);
    std::cout << logl << std::endl;
  }
};
