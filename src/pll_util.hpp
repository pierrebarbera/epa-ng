#ifndef EPA_PLL_UTIL_H_
#define EPA_PLL_UTIL_H_

#include "pllhead.hpp"

#include <string>
#include <sstream>

typedef struct
{
  int clv_valid;
} node_info_t;

/* some code adapted from libpll examples */

// helper functions
void set_missing_branch_length_recursive(pll_utree_t * tree, double length);
void set_unique_clv_indices_recursive(pll_utree_t * tree, const int num_tip_nodes);
void utree_query_branches_recursive(pll_utree_t * node, pll_utree_t ** node_list, int * index);
void free_node_data(pll_utree_t * node);
void get_numbered_newick_string_recursive(pll_utree_t * node, std::ostringstream &ss);

// interface
int cb_partial_traversal(pll_utree_t * node);
int cb_full_traversal(pll_utree_t * node);
int utree_free_node_data(pll_utree_t * node);
unsigned int utree_query_branches(pll_utree_t * node, pll_utree_t ** node_list);
void set_unique_clv_indices(pll_utree_t * tree, const int num_tip_nodes);
void set_missing_branch_length(pll_utree_t * tree, double length);
std::string get_numbered_newick_string(pll_utree_t * root);
pll_utree_t * make_tiny_tree_structure(const pll_utree_t * old_left, const pll_utree_t * old_right);

#endif
