#ifndef EPA_PLL_UTIL_H_
#define EPA_PLL_UTIL_H_

#include "pllhead.hpp"

#define DEFAULT_BRANCH_LENGTH 0.000001

typedef struct
{
  int clv_valid;
} node_info_t;

/* some code adapted from libpll examples */

void set_missing_branch_length_recursive(pll_utree_t * tree, double length);
/* set missing branch lengths to length */
void set_missing_branch_length(pll_utree_t * tree, double length);
void set_unique_clv_indices_recursive(pll_utree_t * tree, const int num_tip_nodes);
void set_unique_clv_indices(pll_utree_t * tree, const int num_tip_nodes);
/* a callback function for performing a partial traversal */
int cb_partial_traversal(pll_utree_t * node);
void free_node_data(pll_utree_t * node);
int utree_free_node_data(pll_utree_t * node);
void utree_query_branches_recursive(pll_utree_t * node, pll_utree_t ** node_list, int * index);
int utree_query_branches(pll_utree_t * node, pll_utree_t ** node_list);

#endif
