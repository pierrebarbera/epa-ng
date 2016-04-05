#pragma once

#include <string>

#include "pllhead.hpp"
#include "Tree_Numbers.hpp"
#include "MSA.hpp"
#include "Sequence.hpp"
#include "Range.hpp"
#include "Model.hpp"

void link_tree_msa(pll_utree_t * tree, pll_partition_t * partition,
              const MSA& msa, const unsigned int num_tip_nodes,
              std::vector<Range> &valid_map);
void precompute_clvs(pll_utree_t * tree, pll_partition_t * partition, const Tree_Numbers& nums);
void split_combined_msa(MSA& source, MSA& target, pll_utree_t * tree, unsigned int num_tip_nodes);
Model get_model(pll_partition_t* partition);

// operator overloads
bool operator==(const pll_utree_t * node, const Sequence& s);
bool operator==(const Sequence& s, const pll_utree_t * node);
