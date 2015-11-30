#ifndef EPA_EPA_PLL_UTIL_H_
#define EPA_EPA_PLL_UTIL_H_

#include "pllhead.hpp"
#include "Tree_Numbers.hpp"
#include "MSA.hpp"
#include "Sequence.hpp"

void link_tree_msa(pll_utree_t * tree, pll_partition_t * partition,
              const MSA& msa, const unsigned int num_tip_nodes);
void precompute_clvs(pll_utree_t * tree, pll_partition_t * partition, const Tree_Numbers& nums);
void bisect(MSA& source, MSA& target, pll_utree_t * tree, unsigned int num_tip_nodes);

// operator overloads
bool operator==(const pll_utree_t * node, const Sequence& s);
bool operator==(const Sequence& s, const pll_utree_t * node);

#endif
