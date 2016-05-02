#include "pll_util.hpp"

#include <iomanip>
#include <stdexcept>

#include "constants.hpp"

using namespace std;

void utree_destroy(pll_utree_t * tree)
{
  if (tree)
  {
    utree_free_node_data(tree);
    pll_utree_destroy(tree);
  }
}
void fasta_close(pll_fasta_t* fptr) { if(fptr) pll_fasta_close(fptr); }

static void set_missing_branch_lengths_recursive(pll_utree_t * tree, double length)
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

      set_missing_branch_lengths_recursive(tree->next->back, length);
      set_missing_branch_lengths_recursive(tree->next->next->back, length);
    }
  }
}

void set_missing_branch_lengths(pll_utree_t * tree, double length)
{
  set_missing_branch_lengths_recursive(tree, length);
  set_missing_branch_lengths_recursive(tree->back, length);
}

static double sum_branch_lengths_recursive(const pll_utree_t * const tree)
{
  double length = 0.0;
  if (tree)
  {
    if (tree->next) // inner node
    {
      length = sum_branch_lengths_recursive(tree->next->back);
      length += sum_branch_lengths_recursive(tree->next->next->back);
    }
    length += tree->length;
  }
  return length;
}

double sum_branch_lengths(const pll_utree_t * const tree)
{
  double length = sum_branch_lengths_recursive(tree);
  length += sum_branch_lengths_recursive(tree->back);
  return length - tree->length;
}

static void set_branch_lengths_recursive(pll_utree_t * tree, double length)
{
  if (tree)
  {
    tree->length = length;

    if (tree->next)
    {
      tree->next->length = length;
      tree->next->next->length = length;

      set_branch_lengths_recursive(tree->next->back, length);
      set_branch_lengths_recursive(tree->next->next->back, length);
    }
  }
}

void set_branch_lengths(pll_utree_t * tree, double length)
{
  set_branch_lengths_recursive(tree, length);
  set_branch_lengths_recursive(tree->back, length);
}

static void set_unique_clv_indices_recursive(pll_utree_t * tree, const int num_tip_nodes)
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

static void free_node_data(pll_utree_t * node)
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

static void utree_query_branches_recursive(pll_utree_t * node, pll_utree_t ** node_list, unsigned int * index)
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

static void get_numbered_newick_string_recursive(pll_utree_t * node, ostringstream &ss, unsigned int * index)
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

void reset_triplet_lengths(pll_utree_t * toward_pendant, pll_partition_t * partition, const double old_length)
{
  double half_original = old_length / 2.0;

  if (toward_pendant)
  {
    // set up branch lengths
    toward_pendant->length = DEFAULT_BRANCH_LENGTH;
    toward_pendant->back->length = DEFAULT_BRANCH_LENGTH;
    toward_pendant->next->length = half_original;
    toward_pendant->next->back->length = half_original;
    toward_pendant->next->next->length = half_original;
    toward_pendant->next->next->back->length = half_original;

    // set up pmatrix indices
    toward_pendant->pmatrix_index = 2;
    toward_pendant->back->pmatrix_index = 2;
    toward_pendant->next->pmatrix_index = 1;
    toward_pendant->next->back->pmatrix_index = 1;
    toward_pendant->next->next->pmatrix_index = 0;
    toward_pendant->next->next->back->pmatrix_index = 0;
  }

  if (partition)
  {
    double branch_lengths[3] = {half_original, half_original, DEFAULT_BRANCH_LENGTH};
    unsigned int matrix_indices[3] = {0, 1, 2};
    unsigned int param_indices[RATE_CATS] = {0};
    pll_update_prob_matrices(partition, param_indices, matrix_indices, branch_lengths, 3);
  }
}

// TODO adjust when using pattern compression
void shift_partition_focus(pll_partition_t * partition, const int offset, const unsigned int span)
{
  const auto clv_size = partition->rate_cats * partition->states_padded;
  const auto num_tips = partition->tips;
  const auto max_index = num_tips + partition->clv_buffers;

  // shift the tip chars
  for (unsigned int i = 0; i < num_tips; i++)
    partition->tipchars[i] += offset;

  // shift the clvs
  for (unsigned int i = 0; i < max_index; i++)
    partition->clv[i] += offset * (int)clv_size;

  // shift the scalers
  for (unsigned int i = 0; i < partition->scale_buffers; i++)
    partition->scale_buffer[i] += offset;

  // shift the pattern weights
  partition->pattern_weights += offset;

  // adjust the number of sites
  partition->sites = span;

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
