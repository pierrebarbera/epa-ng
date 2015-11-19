#ifndef EPA_EPA_PLL_UTIL_H_
#define EPA_EPA_PLL_UTIL_H_

#include <unordered_map>
#include <string>

#include "Tree_Numbers.hpp"

using namespace std;

static void link_tree_msa(pll_utree_t * tree, pll_partition_t * partition,
              const MSA& msa, const int num_tip_nodes)
{
  // obtain pointers to all tip nodes
  auto tip_nodes = new pll_utree_t*[num_tip_nodes];
  pll_utree_query_tipnodes(tree, tip_nodes);

  // and associate the sequences from the MSA file with the correct tips
  /* create a hash table of size num_tip_nodes */
  unordered_map<string, int> map; // mapping labels to tip clv indices

  /* populate the hash table with tree tip labels */
  for (int i = 0; i < num_tip_nodes; ++i)
    map[tip_nodes[i]->label] = i;

  /* find sequences in hash table and link them with the corresponding taxa */
  for (auto s : msa)
  {
    auto clv_index = map.find(s.header());

    if (clv_index == map.end())
      throw runtime_error{string("Sequence with header does not appear in the tree: ") + s.header()};

    pll_set_tip_states(partition, clv_index->second, pll_map_nt, s.sequence().c_str());
  }

  /* free the tip node list */
  delete [] tip_nodes;
};

static void precompute_clvs(pll_utree_t * tree, pll_partition_t * partition, const Tree_Numbers& nums)
{
  int num_matrices, num_ops;

  /* buffer for creating a postorder traversal structure */
  auto travbuffer = new pll_utree_t*[nums.nodes];
  auto branch_lengths = new double[nums.branches];
  auto matrix_indices = new int[nums.branches];
  auto operations = new pll_operation_t[nums.nodes];
  // TODO num_nodes too many? whats the upper bound considering tip traversal?

  // get a list of all tip nodes
  auto tip_nodes = new pll_utree_t*[nums.tip_nodes];
  pll_utree_query_tipnodes(tree, tip_nodes);

  /* adjust clv indices such that every direction has its own */
  set_unique_clv_indices(tree, nums.tip_nodes);

  for (int i = 0; i < nums.tip_nodes; ++i)
  {
    /* perform a partial postorder traversal of the unrooted tree */

    int traversal_size = pll_utree_traverse(tip_nodes[i]->back,
                                            cb_partial_traversal,
                                            travbuffer);
    if (traversal_size == -1)
      throw runtime_error{"Function pll_utree_traverse() requires inner nodes as parameters"};

    /* given the computed traversal descriptor, generate the operations
       structure, and the corresponding probability matrix indices that
       may need recomputing */
    pll_utree_create_operations(travbuffer,
                                traversal_size,
                                branch_lengths,
                                matrix_indices,
                                operations,
                                &num_matrices,
                                &num_ops);

    pll_update_prob_matrices(partition,
                             0,             // use model 0
                             matrix_indices,// matrices to update
                             branch_lengths,
                             num_matrices); // how many should be updated

    /* use the operations array to compute all num_ops inner CLVs. Operations
       will be carried out sequentially starting from operation 0 towrds num_ops-1 */
    pll_update_partials(partition, operations, num_ops);
  }

  delete [] travbuffer;
  delete [] tip_nodes;
  delete [] branch_lengths;
  delete [] matrix_indices;
  delete [] operations;

};

#endif
