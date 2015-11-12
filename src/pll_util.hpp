#include "pllhead.hpp"

#define DEFAULT_BRANCH_LENGTH 0.000001

typedef struct
{
  int clv_valid;
} node_info_t;

// TODO attribution!
/* code adapted from libpll examples */

static void set_missing_branch_length_recursive(pll_utree_t * tree,
                                                double length)
{
  if (tree)
  {
    /* set branch length to length if not set */
    if (!tree->length)
      tree->length = length;

    if (tree->next)
    {
      if (!tree->next->length)
        tree->next->length = length;

      if (!tree->next->next->length)
        tree->next->next->length = length;

      set_missing_branch_length_recursive(tree->next->back, length);
      set_missing_branch_length_recursive(tree->next->next->back, length);
    }
  }
};

/* set missing branch lengths to length */
static void set_missing_branch_length(pll_utree_t * tree, double length)
{
  set_missing_branch_length_recursive(tree, length);
  set_missing_branch_length_recursive(tree->back, length);
};

static void set_unique_clv_indices_recursive(pll_utree_t * tree, const int num_tip_nodes)
{
  if (tree && tree->next)
  {
    int idx = tree->clv_index;
    /* new index is in principle old index * 3 + 0 for the first traversed, + 1 for the
      second etc., however we need to account for the first num_tip_nodes entries, as
      the tip nodes only have one clv  */
    idx = (idx - num_tip_nodes) * 3 + num_tip_nodes;
    tree->clv_index = idx;
    //tree->scaler_index = idx;
    tree->next->clv_index = ++idx;
    //tree->next->scaler_index = idx;
    tree->next->next->clv_index = ++idx;
    //tree->next->next->scaler_index = idx;

    // recurse
    set_unique_clv_indices_recursive(tree->next->back, num_tip_nodes);
    set_unique_clv_indices_recursive(tree->next->next->back, num_tip_nodes);
  }
};

static void set_unique_clv_indices(pll_utree_t * tree, const int num_tip_nodes)
{
  set_unique_clv_indices_recursive(tree, num_tip_nodes);
  set_unique_clv_indices_recursive(tree->back, num_tip_nodes);
};

/* a callback function for performing a partial traversal */
static int cb_partial_traversal(pll_utree_t * node)
{
  node_info_t * node_info;

  /* if we don't want tips in the traversal we must return 0 here. For now,
     allow tips */
  if (!node->next) return 1;

  /* get the data element from the node and check if the CLV vector is
     oriented in the direction that we want to traverse. If the data
     element is not yet allocated then we allocate it, set the direction
     and instruct the traversal routine to place the node in the traversal array
     by returning 1 */
  node_info = (node_info_t *)(node->data);
  if (!node_info)
  {
    /* allocate data element */
    node->data             = (node_info_t *)calloc(1,sizeof(node_info_t));
    node->next->data       = (node_info_t *)calloc(1,sizeof(node_info_t));
    node->next->next->data = (node_info_t *)calloc(1,sizeof(node_info_t));

    /* set orientation on selected direction and traverse the subtree */
    node_info = (node_info_t*) (node->data);
    node_info->clv_valid = 1;
    return 1;
  }

  /* if the data element was already there and the CLV on this direction is
     set, i.e. the CLV is valid, we instruct the traversal routine not to
     traverse the subtree rooted in this node/direction by returning 0 */
  if (node_info->clv_valid) return 0;

  /* otherwise, set orientation on selected direction */
  node_info->clv_valid = 1;

  return 1;
};

static void free_node_data(pll_utree_t * node)
{
  free(node->data);

  if (node->next) // we are at a inner node
  {
    // free all memory behind data of current node triple
    free(node->next->data);
    free(node->next->next->data);
    // recurse to sub trees
    free_node_data(node->next->back);
    free_node_data(node->next->next->back);
  }
};
