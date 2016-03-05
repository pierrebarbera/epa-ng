#pragma once

#include "pllhead.hpp"

#include <string>
#include <sstream>

typedef struct
{
  int clv_valid;
} node_info_t;

/* some code adapted from libpll examples */

// interface
int cb_partial_traversal(pll_utree_t * node);
int cb_full_traversal(pll_utree_t * node);
int utree_free_node_data(pll_utree_t * node);
unsigned int utree_query_branches(pll_utree_t * node, pll_utree_t ** node_list);
void set_unique_clv_indices(pll_utree_t * tree, const int num_tip_nodes);
void set_missing_branch_length(pll_utree_t * tree, double length);
void set_branch_length(pll_utree_t * tree, double length);
std::string get_numbered_newick_string(pll_utree_t * root);
pll_utree_t * make_tiny_tree_structure(const pll_utree_t * old_proximal, const pll_utree_t * old_distal,
  const bool tip_tip_case);
pll_partition_t * make_tiny_partition(const pll_partition_t * old_partition, const pll_utree_t * tree,
  const pll_utree_t * old_proximal, const pll_utree_t * old_distal, const bool tip_tip_case);
pll_utree_t * get_tip_node(pll_utree_t * node);
