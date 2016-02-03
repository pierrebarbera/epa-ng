#include "epa_pll_util.hpp"

#include <unordered_map>
#include <algorithm>

#include "pll_util.hpp"
#include "set_manipulators.hpp"

using namespace std;

void link_tree_msa(pll_utree_t * tree, pll_partition_t * partition,
              const MSA& msa, const unsigned int num_tip_nodes,
              vector<Range> &valid_map)
{
  // obtain pointers to all tip nodes
  vector<pll_utree_t*> tip_nodes(num_tip_nodes);
  pll_utree_query_tipnodes(tree, &tip_nodes[0]);

  // and associate the sequences from the MSA file with the correct tips
  /* create a hash table of size num_tip_nodes */
  unordered_map<string, unsigned int> map; // mapping labels to tip clv indices

  /* populate the hash table with tree tip labels */
  for (unsigned int i = 0; i < num_tip_nodes; ++i)
    map[tip_nodes[i]->label] = i;

  /* find sequences in hash table and link them with the corresponding taxa */
  for (auto const &s : msa)
  {
    auto map_value = map.find(s.header());
    auto clv_index = map_value->second;

    if (map_value == map.end())
      throw runtime_error{string("Sequence with header does not appear in the tree: ") + s.header()};

    // associates the sequence with the tip by calculating the tips clv buffers
    pll_set_tip_states(partition, clv_index, pll_map_nt, s.sequence().c_str());

    // TODO improvement?
    // remember the valid-range of the sequence, indexed by tip clv index
    valid_map[clv_index] = get_valid_range(s.sequence());
  }
}

void precompute_clvs(pll_utree_t * tree, pll_partition_t * partition, const Tree_Numbers& nums)
{
  unsigned int num_matrices, num_ops;

  /* various buffers for creating a postorder traversal and operations structures */
  vector<pll_utree_t*> travbuffer(nums.nodes);
  vector<double> branch_lengths(nums.branches);
  vector<unsigned int> matrix_indices(nums.branches);
  vector<pll_operation_t> operations(nums.nodes);

  // get a list of all tip nodes
  vector<pll_utree_t*> tip_nodes(nums.tip_nodes);
  pll_utree_query_tipnodes(tree, &tip_nodes[0]);

  /* adjust clv indices such that every direction has its own */
  set_unique_clv_indices(tree, nums.tip_nodes);

  for (auto node : tip_nodes)
  {
    /* perform a partial postorder traversal of the unrooted tree  starting at the current tip
      and returning every node whose clv in the direction of the tip hasn't been calculated yet*/
    unsigned int traversal_size;
    if(pll_utree_traverse(node->back, cb_partial_traversal, &travbuffer[0], &traversal_size)
                != PLL_SUCCESS)
      throw runtime_error{"Function pll_utree_traverse() requires inner nodes as parameters"};

    /* given the computed traversal descriptor, generate the operations
       structure, and the corresponding probability matrix indices that
       may need recomputing */
    pll_utree_create_operations(&travbuffer[0],
                                traversal_size,
                                &branch_lengths[0],
                                &matrix_indices[0],
                                &operations[0],
                                &num_matrices,
                                &num_ops);

    pll_update_prob_matrices(partition,
                             0,             // use model 0
                             &matrix_indices[0],// matrices to update
                             &branch_lengths[0],
                             num_matrices); // how many should be updated

    /* use the operations array to compute all num_ops inner CLVs. Operations
       will be carried out sequentially starting from operation 0 towrds num_ops-1 */
    pll_update_partials(partition, &operations[0], num_ops);
  }
}

void split_combined_msa(MSA& source, MSA& target, pll_utree_t * tree, unsigned int num_tip_nodes)
{
  vector<pll_utree_t*> tip_nodes(num_tip_nodes);
  pll_utree_query_tipnodes(tree, &tip_nodes[0]);

  auto falsegroup_begin = partition(source.begin(), source.end(),
    [&tip_nodes](const Sequence& em)
    {
      return find(tip_nodes.begin(), tip_nodes.end(), em) != tip_nodes.end();
    });
  target.num_sites(source.num_sites());
  target.move_sequences(falsegroup_begin, source.end());
  source.erase(falsegroup_begin, source.end());
}

bool operator==(const pll_utree_t * node, const Sequence& s)
{
  return s.header().compare(node->label) == 0;
}

bool operator==(const Sequence& s, const pll_utree_t * node)
{
  return operator==(node, s);
}
