#include "Tiny_Tree.hpp"

#include <vector>
#include <cassert>

#include "pll_util.hpp"
#include "optimize.hpp"
#include "Tree_Numbers.hpp"

using namespace std;

Tiny_Tree::Tiny_Tree(pll_utree_t * edge_node, pll_partition_t * old_partition, Model model,
   bool heuristic)
  : heuristic_(heuristic), model_(model)
{
  assert(edge_node != NULL);
  assert(old_partition != NULL);

  pll_utree_t * old_proximal = edge_node->back;
  pll_utree_t * old_distal = edge_node;

  tree_ = make_tiny_tree_structure(old_proximal, old_distal);

  partition_ =
  pll_partition_create(   3, // tips
                          1, // extra clv's
                          old_partition->states,
                          old_partition->sites,
                          0, // number of mixture models
                          old_partition->rate_matrices,
                          3, // number of prob. matrices (one per unique branch length)
                          old_partition->rate_cats,
                          4, // number of scale buffers (one per node)
                          old_partition->attributes);

  // shallow copy model params
  partition_->rates = old_partition->rates;
  partition_->subst_params = old_partition->subst_params;
  partition_->frequencies = old_partition->frequencies;

  // shallow copy 2 existing nodes clvs
  assert(old_partition->clv[old_proximal->clv_index] != NULL);
  assert(old_partition->clv[old_distal->clv_index] != NULL);

  //deep copy clv's
  memcpy(partition_->clv[TINY_PROXIMAL_CLV_INDEX], old_partition->clv[old_proximal->clv_index],
      sizeof(double) * old_partition->sites * old_partition->rate_cats * old_partition->states);
  memcpy(partition_->clv[TINY_DISTAL_CLV_INDEX], old_partition->clv[old_distal->clv_index],
      sizeof(double) * old_partition->sites * old_partition->rate_cats * old_partition->states);

  // deep copy scalers
  if (old_proximal->scaler_index != PLL_SCALE_BUFFER_NONE)
    memcpy(partition_->scale_buffer[TINY_PROXIMAL_CLV_INDEX],
        old_partition->scale_buffer[old_proximal->scaler_index],
        sizeof(unsigned int) * old_partition->sites);
  if (old_distal->scaler_index != PLL_SCALE_BUFFER_NONE)
    memcpy(partition_->scale_buffer[TINY_DISTAL_CLV_INDEX],
        old_partition->scale_buffer[old_distal->scaler_index],
        sizeof(unsigned int) * old_partition->sites);


  // precreate some of the operation fields that are static throughout the trees lifetime
  ops_.parent_clv_index    = TINY_INNER_CLV_INDEX;
  ops_.child1_clv_index    = TINY_PROXIMAL_CLV_INDEX;
  ops_.child1_scaler_index = tree_->next->next->back->scaler_index;
  ops_.child2_clv_index    = TINY_DISTAL_CLV_INDEX;
  ops_.child2_scaler_index = tree_->next->back->scaler_index;
  ops_.parent_scaler_index = TINY_INNER_CLV_INDEX;

  // if the heuristic is to be used, the branch lengths only need to be set once
  if (heuristic_)
  {
    /* heuristic insertion as described in EPA paper from 2011 (Berger et al.):
      original branch, now split by "inner", or base, node of the inserted sequence,
      defines the new branch lengths between inner and old left/right respectively
      as old branch length / 2.
      The new branch leading from inner to the new tip is initialized with length 0.9,
      which is the default branch length in RAxML.
    */
    double branch_lengths[2] = { old_proximal->length / 2, DEFAULT_BRANCH_LENGTH};
    unsigned int matrix_indices[2] = { 0, 1 };

    // use branch lengths to compute the probability matrices
    pll_update_prob_matrices(partition_, 0, matrix_indices, branch_lengths, 2);
  }

}

Tiny_Tree::~Tiny_Tree()
{
  if (partition_ != nullptr)
  {
    // unset shallow copied things
    partition_->rates = nullptr;
    partition_->subst_params = nullptr;
    partition_->frequencies = nullptr;

    pll_partition_destroy(partition_);
  }

  if (tree_ != nullptr)
    pll_utree_destroy(tree_);
}

std::tuple<double, double, double> Tiny_Tree::place(const Sequence &s)
{
  assert(partition_ != NULL);

  // init the new tip with s.sequence(), branch length
  pll_set_tip_states(partition_, TINY_NEW_TIP_CLV_INDEX, pll_map_nt, s.sequence().c_str());

  unsigned int inner_matrix_index = 1;

  // optimize branch lengths
  if (!heuristic_)
  {
    Tree_Numbers nums;
    nums.init(3);
    optimize(model_, tree_, partition_, nums, true);

    ops_.child1_matrix_index = tree_->next->next->pmatrix_index;
    ops_.child2_matrix_index = tree_->next->pmatrix_index;
    inner_matrix_index = tree_->pmatrix_index;
  }
  else
  {
    ops_.child1_matrix_index = 0;
    ops_.child2_matrix_index = 0;
  }

  // use update_partials to compute the clv pointing toward the new tip
  pll_update_partials(partition_, &ops_, TINY_NUM_OPS);

  // compute the loglikelihood using inner node and new tip
  auto logl = pll_compute_edge_loglikelihood(partition_,
                                        TINY_NEW_TIP_CLV_INDEX,
                                        PLL_SCALE_BUFFER_NONE,// scaler_index
                                        TINY_INNER_CLV_INDEX,
                                        TINY_INNER_CLV_INDEX,  // scaler_index
                                        inner_matrix_index,
                                        0);// freq index

  return make_tuple(logl, tree_->next->length, tree_->length);
}
