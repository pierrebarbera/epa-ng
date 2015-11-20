#ifndef EPA_EPA_PLL_UTIL_H_
#define EPA_EPA_PLL_UTIL_H_

#include "pllhead.hpp"
#include "Tree_Numbers.hpp"
#include "MSA.hpp"

void link_tree_msa(pll_utree_t * tree, pll_partition_t * partition,
              const MSA& msa, const int num_tip_nodes);
void precompute_clvs(pll_utree_t * tree, pll_partition_t * partition, const Tree_Numbers& nums);

#endif
