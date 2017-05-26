#pragma once

#include <string>
#include <vector>

#include "pllhead.hpp"
#include "Tree_Numbers.hpp"
#include "MSA.hpp"
#include "MSA_Stream.hpp"
#include "Sequence.hpp"
#include "Range.hpp"
#include "Model.hpp"
#include "Tree.hpp"

void link_tree_msa( pll_utree_t * tree, 
                    pll_partition_t * partition, 
                    Model& model, 
                    const MSA& msa, 
                    const unsigned int num_tip_nodes);
void precompute_clvs( pll_utree_t const * const tree, 
                      pll_partition_t * partition, 
                      const Tree_Numbers& nums);
void split_combined_msa(MSA& source, 
                        MSA& target, 
                        Tree& tree);
Model get_model(pll_partition_t* partition);

// operator overloads
bool operator==(const pll_unode_t * node, const Sequence& s);
bool operator==(const Sequence& s, const pll_unode_t * node);
