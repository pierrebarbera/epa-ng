#ifndef EPA_EPA_PLL_UTIL_H_
#define EPA_EPA_PLL_UTIL_H_

#include <unordered_map>
#include <string>

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


#endif
