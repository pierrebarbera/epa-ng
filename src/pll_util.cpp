#include "pll_util.hpp"

#include <iomanip>
#include <stdexcept>

#include "constants.hpp"

using namespace std;

// helper functions
void set_missing_branch_length_recursive(pll_utree_t * tree, double length);
void set_branch_length_recursive(pll_utree_t * tree, double length);
void set_unique_clv_indices_recursive(pll_utree_t * tree, const int num_tip_nodes);
void utree_query_branches_recursive(pll_utree_t * node, pll_utree_t ** node_list, int * index);
void free_node_data(pll_utree_t * node);
void get_numbered_newick_string_recursive(pll_utree_t * node, std::ostringstream &ss);

void set_missing_branch_length_recursive(pll_utree_t * tree, double length)
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

void set_missing_branch_length(pll_utree_t * tree, double length)
{
  set_missing_branch_length_recursive(tree, length);
  set_missing_branch_length_recursive(tree->back, length);
}

void set_branch_length_recursive(pll_utree_t * tree,
                                                double length)
{
  if (tree)
  {
    tree->length = length;

    if (tree->next)
    {
      tree->next->length = length;
      tree->next->next->length = length;

      set_branch_length_recursive(tree->next->back, length);
      set_branch_length_recursive(tree->next->next->back, length);
    }
  }
}

void set_branch_length(pll_utree_t * tree, double length)
{
  set_branch_length_recursive(tree, length);
  set_branch_length_recursive(tree->back, length);
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
}

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

pll_partition_t * make_tiny_partition(const pll_partition_t * old_partition, const pll_utree_t * tree,
  const pll_utree_t * old_proximal, const pll_utree_t * old_distal, const bool tip_tip_case)
{
  /**
    As we work with PLL_PATTERN_TIP functionality, special care has to be taken in regards to the tree and partition
    structure: PLL assumes that any node with clv index < number of tips is in fact a real tip, that is
    a tip that uses a character array instead of a real clv. Here we need to set up the tree/partition to fool pll:
    the tips that actually contain CLVs copied over from the reference tree have their index set to greater than
    number of tips. This results in a acceptable amount of wasted memory that is never used (num_sites * bytes
    * number of clv-tips)
  */
  unsigned int num_clv_tips = 2; // tip_inner case: both reference nodes are inner nodes
  if (tip_tip_case)
    num_clv_tips = 1; // one for the "proximal" clv tip

  auto distal = tree->next->back;
  auto proximal = tree->next->next->back;

  pll_partition_t * tiny = pll_partition_create(
    3, // tips
    1 + num_clv_tips, // extra clv's
    old_partition->states, old_partition->sites,
    0, // number of mixture models
    old_partition->rate_matrices,
    3, // number of prob. matrices (one per possible unique branch length)
    old_partition->rate_cats,
    3, // number of scale buffers (one per possible inner node)
    pll_map_nt,
    old_partition->attributes);

  assert(tiny);

  // TODO release the buffers allocated for these
  tiny->rates = old_partition->rates;
  tiny->subst_params = old_partition->subst_params;
  tiny->frequencies = old_partition->frequencies;
  tiny->eigenvecs = old_partition->eigenvecs;
  tiny->inv_eigenvecs = old_partition->inv_eigenvecs;
  tiny->eigenvals = old_partition->eigenvals;
  tiny->prop_invar = old_partition->prop_invar;
  tiny->eigen_decomp_valid = old_partition->eigen_decomp_valid;
  tiny->pattern_weights = old_partition->pattern_weights;

  // shalow/deep copy tip_tip_pattern specific things
  // shallow
  tiny->lh_statepair = old_partition->lh_statepair;
  tiny->charmap = old_partition->charmap;
  tiny->revmap = old_partition->revmap;
  // deep
  tiny->maxstates = old_partition->maxstates;
  tiny->log2_maxstates = old_partition->log2_maxstates;
  tiny->log2_rates = old_partition->log2_rates;
  tiny->log2_states = old_partition->log2_states;

  assert(old_partition->clv[old_proximal->clv_index] != NULL);
  assert(old_partition->clv[old_distal->clv_index] != NULL);

  unsigned int clv_size = sizeof(double) * old_partition->sites * old_partition->rate_cats * old_partition->states;

  // deep copy clv's
  memcpy(tiny->clv[proximal->clv_index],
    old_partition->clv[old_proximal->clv_index],
    clv_size);

  if(tip_tip_case)
    memcpy(tiny->tipchars[distal->clv_index],
      old_partition->tipchars[old_distal->clv_index],
      sizeof(char) * old_partition->sites );
  else
    memcpy(tiny->clv[distal->clv_index],
      old_partition->clv[old_distal->clv_index],
      clv_size);

  // deep copy scalers
  if (old_proximal->scaler_index != PLL_SCALE_BUFFER_NONE)
    memcpy(tiny->scale_buffer[proximal->scaler_index],
      old_partition->scale_buffer[old_proximal->scaler_index],
      sizeof(unsigned int) * old_partition->sites);
  if (old_distal->scaler_index != PLL_SCALE_BUFFER_NONE)
    memcpy(tiny->scale_buffer[distal->scaler_index],
      old_partition->scale_buffer[old_distal->scaler_index],
      sizeof(unsigned int) * old_partition->sites);

  return tiny;
}

pll_utree_t * make_tiny_tree_structure(const pll_utree_t * old_proximal, const pll_utree_t * old_distal,
  const bool tip_tip_case)
{
  const unsigned int inner_scaler_index = 1;
  const unsigned int proximal_scaler_index = 0;
  const unsigned int distal_scaler_index = 2;

  /**
    As we work with PLL_PATTERN_TIP functionality, special care has to be taken in regards to the tree and partition
    structure: PLL assumes that any node with clv index < number of tips is in fact a real tip, that is
    a tip that uses a character array instead of a real clv. Here we need to set up the tree/partition to fool pll:
    the tips that actually contain CLVs copied over from the reference tree have their index set to greater than
    number of tips. This results in a acceptable amount of wasted memory that is never used (num_sites * bytes
    * number of clv-tips)
  */
  // if tip-inner case
  unsigned int distal_clv_index = 5;
  const unsigned int proximal_clv_index = 4;
  const unsigned int inner_clv_index = 3;
  const unsigned int new_tip_clv_index = 1;

  if (tip_tip_case)
    distal_clv_index = 2;

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
  inner->clv_index = inner_clv_index;
  inner->next->clv_index = inner_clv_index;
  inner->next->next->clv_index = inner_clv_index;
  proximal->clv_index = proximal_clv_index;
  distal->clv_index = distal_clv_index;
  new_tip->clv_index = new_tip_clv_index;

  // set up scaler indices
  new_tip->scaler_index = PLL_SCALE_BUFFER_NONE;
  inner->scaler_index = inner_scaler_index;
  inner->next->scaler_index = inner_scaler_index;
  inner->next->next->scaler_index = inner_scaler_index;
  proximal->scaler_index = (old_proximal->scaler_index == PLL_SCALE_BUFFER_NONE) ?
    PLL_SCALE_BUFFER_NONE : proximal_scaler_index;
  distal->scaler_index = (old_distal->scaler_index == PLL_SCALE_BUFFER_NONE) ?
    PLL_SCALE_BUFFER_NONE : distal_scaler_index;

  // set up branch lengths
  double half_old = old_distal->length / 2.0;
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

/* Function to return the tip node if either <node> or <node->back> is one. Otherweise
  returns null. */
pll_utree_t * get_tip_node(pll_utree_t * node)
{
  pll_utree_t * tip_node = nullptr;
  // node is the tip
  if (!node->next)
    tip_node = node;
  else if (!node->back->next)
    tip_node = node->back;

  return tip_node;
}
