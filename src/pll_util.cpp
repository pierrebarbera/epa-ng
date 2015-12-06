#include "pll_util.hpp"

#include <iomanip>

#include "constants.hpp"

using namespace std;

void set_missing_branch_length_recursive(pll_utree_t * tree,
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
}

/* set missing branch lengths to length */
void set_missing_branch_length(pll_utree_t * tree, double length)
{
  set_missing_branch_length_recursive(tree, length);
  set_missing_branch_length_recursive(tree->back, length);
}

void set_unique_clv_indices_recursive(pll_utree_t * tree, const int num_tip_nodes)
{
  if (tree && tree->next)
  {
    unsigned int idx = tree->clv_index;
    /* new index is in principle old index * 3 + 0 for the first traversed, + 1 for the
      second etc., however we need to account for the first num_tip_nodes entries, as
      the tip nodes only have one clv  */
    idx = (idx - num_tip_nodes) * 3 + num_tip_nodes;
    tree->clv_index = idx;
    tree->scaler_index = idx;
    tree->next->clv_index = ++idx;
    tree->next->scaler_index = idx;
    tree->next->next->clv_index = ++idx;
    tree->next->next->scaler_index = idx;

    // recurse
    set_unique_clv_indices_recursive(tree->next->back, num_tip_nodes);
    set_unique_clv_indices_recursive(tree->next->next->back, num_tip_nodes);
  }
}

void set_unique_clv_indices(pll_utree_t * tree, const int num_tip_nodes)
{
  set_unique_clv_indices_recursive(tree, num_tip_nodes);
  set_unique_clv_indices_recursive(tree->back, num_tip_nodes);
}

/* a callback function for performing a partial traversal */
int cb_partial_traversal(pll_utree_t * node)
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
}

int cb_full_traversal(pll_utree_t * node)
{
  (void) node;
  return 1;
};

void free_node_data(pll_utree_t * node)
{

  // currently we don't allocate a data struct at the tips

  if (node->next) // we are at a inner node
  {
    // free all memory behind data of current node triplet
    free(node->data);
    free(node->next->data);
    free(node->next->next->data);
    // recurse to sub trees
    free_node_data(node->next->back);
    free_node_data(node->next->next->back);
  }
}

int utree_free_node_data(pll_utree_t * node)
{
  if (!node->next) return 0; // not a inner node!

  // we start at a trifurcation: call explicitly for this node and its adjacent node
  free_node_data(node);
  free_node_data(node->back);

  return 1;
}

void utree_query_branches_recursive(pll_utree_t * node, pll_utree_t ** node_list, unsigned int * index)
{
  // Postorder traversal

  if (node->next) // inner node
  {
    utree_query_branches_recursive(node->next->back, node_list, index);
    utree_query_branches_recursive(node->next->next->back, node_list, index);
  }
  node_list[*index] = node;
  *index = *index + 1;
}

unsigned int utree_query_branches(pll_utree_t * node, pll_utree_t ** node_list)
{
  unsigned int index = 0;

  // assure that we start at inner node
  if (!node->next) node = node->back;

  // utree-function: we start at a trifucation
  utree_query_branches_recursive(node->back, node_list, &index);
  utree_query_branches_recursive(node->next->back, node_list, &index);
  utree_query_branches_recursive(node->next->next->back, node_list, &index);

  return index;
}

void get_numbered_newick_string_recursive(pll_utree_t * node, ostringstream &ss, unsigned int * index)
{

  if (node->next) // inner node
  {
    ss << "(";
    get_numbered_newick_string_recursive(node->next->back, ss, index);
    ss << ",";
    get_numbered_newick_string_recursive(node->next->next->back, ss, index);
    ss << "):" << node->length << "{" << *index << "}";
  } else {
    ss << node->label << ":" << setprecision(20) << node->length << "{" << *index << "}";
  }
  *index = *index + 1;

}

string get_numbered_newick_string(pll_utree_t * root)
{
  ostringstream ss;
  unsigned int index = 0;
  //ss.precision(20);

  if (!root->next) root = root->back; // ensure that we start at inner node

  ss << "(";

  get_numbered_newick_string_recursive(root->back, ss, &index);
  ss << ",";
  get_numbered_newick_string_recursive(root->next->back, ss, &index);
  ss << ",";
  get_numbered_newick_string_recursive(root->next->next->back, ss, &index);

  ss << ")";
  ss << ";";

  return ss.str();
}

// doin it the hard way!
pll_utree_t * make_tiny_tree_structure(const pll_utree_t * old_left, const pll_utree_t * old_right)
{
  pll_utree_t * inner = (pll_utree_t *) calloc(1,sizeof(pll_utree_t));
  inner->next = (pll_utree_t *) calloc(1,sizeof(pll_utree_t));
  inner->next->next = (pll_utree_t *) calloc(1,sizeof(pll_utree_t));
  inner->next->next->next = inner;

  pll_utree_t * new_tip = (pll_utree_t *) calloc(1,sizeof(pll_utree_t));
  pll_utree_t * proximal = (pll_utree_t *) calloc(1,sizeof(pll_utree_t));
  pll_utree_t * distal = (pll_utree_t *) calloc(1,sizeof(pll_utree_t));

  // connect the nodes to each other
  inner->back = new_tip;
  new_tip->back = inner;
  inner->next->back = distal;
  distal->back = inner->next;
  inner->next->next->back = proximal;
  proximal->back = inner->next->next;

  // set up clv indices
  inner->clv_index = TINY_INNER_CLV_INDEX;
  inner->next->clv_index = TINY_INNER_CLV_INDEX;
  inner->next->next->clv_index = TINY_INNER_CLV_INDEX;
  proximal->clv_index = TINY_PROXIMAL_CLV_INDEX;
  distal->clv_index = TINY_DISTAL_CLV_INDEX;
  new_tip->clv_index = TINY_NEW_TIP_CLV_INDEX;

  // set up scaler indices
  new_tip->scaler_index = PLL_SCALE_BUFFER_NONE;
  inner->scaler_index = TINY_INNER_CLV_INDEX;
  inner->next->scaler_index = TINY_INNER_CLV_INDEX;
  inner->next->next->scaler_index = TINY_INNER_CLV_INDEX;
  proximal->scaler_index = (old_left->scaler_index == PLL_SCALE_BUFFER_NONE) ?
    PLL_SCALE_BUFFER_NONE : TINY_PROXIMAL_CLV_INDEX;
  distal->scaler_index = (old_right->scaler_index == PLL_SCALE_BUFFER_NONE) ?
    PLL_SCALE_BUFFER_NONE : TINY_DISTAL_CLV_INDEX;

  // set up branch lengths
  double half_old = old_right->length / 2;
  new_tip->length = DEFAULT_BRANCH_LENGTH;
  new_tip->back->length = DEFAULT_BRANCH_LENGTH;
  proximal->length = half_old;
  proximal->back->length = half_old;
  distal->length = half_old;
  distal->back->length = half_old;

  // set up pmatrix indices
  inner->pmatrix_index = 2;
  new_tip->pmatrix_index = 2;
  inner->next->pmatrix_index = 1;
  distal->pmatrix_index = 1;
  inner->next->next->pmatrix_index = 0;
  proximal->pmatrix_index = 0;

  return inner;
}
