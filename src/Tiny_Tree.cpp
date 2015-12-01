#include "Tiny_Tree.hpp"

#include <vector>
#include <cassert>

#include "pll_util.hpp"

using namespace std;

Tiny_Tree::Tiny_Tree(pll_utree_t * edge_node, pll_partition_t * old_partition, bool heuristic)
  : heuristic_(heuristic)
{
  assert(edge_node != NULL);
  assert(old_partition != NULL);

  pll_utree_t * old_left = edge_node;
  pll_utree_t * old_right = edge_node->back;

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
  assert(old_partition->clv[old_left->clv_index] != NULL);
  assert(old_partition->clv[old_right->clv_index] != NULL);

  partition_->clv[OLD_LEFT_CLV_INDEX] =
                          old_partition->clv[old_left->clv_index];
  partition_->clv[OLD_RIGHT_CLV_INDEX] =
                          old_partition->clv[old_right->clv_index];

  // shallow copy scalers
  if (old_left->scaler_index != PLL_SCALE_BUFFER_NONE)
  {
    partition_->scale_buffer[OLD_LEFT_CLV_INDEX] =
                          old_partition->scale_buffer[old_left->scaler_index];
    ops_.child1_scaler_index = OLD_LEFT_CLV_INDEX;
  }
  else
    ops_.child1_scaler_index = PLL_SCALE_BUFFER_NONE;

  if (old_right->scaler_index != PLL_SCALE_BUFFER_NONE)
  {
    partition_->scale_buffer[OLD_RIGHT_CLV_INDEX] =
                          old_partition->scale_buffer[old_right->scaler_index];
    ops_.child2_scaler_index = OLD_RIGHT_CLV_INDEX;
  }
  else
    ops_.child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  // deep copy clv's
  // memcpy(partition_->clv[OLD_LEFT_CLV_INDEX], old_partition->clv[old_left->clv_index],
  //     sizeof(double) * old_partition->sites);
  // memcpy(partition_->clv[OLD_RIGHT_CLV_INDEX], old_partition->clv[old_right->clv_index],
  //     sizeof(double) * old_partition->sites);

  // deep copy scalers
  // if (old_left->scaler_index != PLL_SCALE_BUFFER_NONE)
  // {
  //   memcpy(partition_->scale_buffer[OLD_LEFT_CLV_INDEX],
  //       old_partition->scale_buffer[old_left->scaler_index],
  //       sizeof(unsigned int) * old_partition->sites);
  //   ops_.child1_scaler_index = OLD_LEFT_CLV_INDEX;
  // }
  // else
  //   ops_.child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  //
  // if (old_right->scaler_index != PLL_SCALE_BUFFER_NONE)
  // {
  //   memcpy(partition_->scale_buffer[OLD_RIGHT_CLV_INDEX],
  //       old_partition->scale_buffer[old_right->scaler_index],
  //       sizeof(unsigned int) * old_partition->sites);
  //   ops_.child2_scaler_index = OLD_RIGHT_CLV_INDEX;
  // }
  // else
  //   ops_.child2_scaler_index = PLL_SCALE_BUFFER_NONE;


  // precreate some of the operation fields that are static throughout the trees lifetime
  ops_.parent_clv_index    = INNER_CLV_INDEX;
  ops_.child1_clv_index    = OLD_LEFT_CLV_INDEX;
  ops_.child2_clv_index    = OLD_RIGHT_CLV_INDEX;
  ops_.parent_scaler_index = PLL_SCALE_BUFFER_NONE;//INNER_CLV_INDEX;

  // TODO might not be needed depending on how BLO is done
  old_branch_length_ = old_left->length;

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
    // TODO free
    auto branch_lengths = new double[2];// = { old_branch_length_ / 2, DEFAULT_BRANCH_LENGTH};
    branch_lengths[0] = old_branch_length_ / 2;
    branch_lengths[1] = DEFAULT_BRANCH_LENGTH;
    auto matrix_indices = new unsigned int[2];
    matrix_indices[0] = 0;
    matrix_indices[1] = 1;


  // use branch lengths to compute the probability matrices
  pll_update_prob_matrices(partition_, 0, matrix_indices, branch_lengths, 2);
  }

}

double Tiny_Tree::place(const Sequence &s)
{
  assert(partition_ != NULL);

  // init the new tip with s.sequence(), branch length
  pll_set_tip_states(partition_, NEW_TIP_CLV_INDEX, pll_map_nt, s.sequence().c_str());

  // optimize branch lengths
  if (!heuristic_)
  {
    assert(true);
    // TODO Newton-Raphson
    // TODO needs full tree? only done once?
    // something like: setup params structure, call opt. branch lengths local with radius 1
    // on the inner node
  } else {
    ops_.child1_matrix_index = 0;// TODO depends on NR vs heuristic
    ops_.child2_matrix_index = 0;// TODO depends on NR vs heuristic
  }


  // use update_partials to compute the clv pointing toward the new tip
  pll_update_partials(partition_, &ops_, NUM_OPS);

  // compute the loglikelihood using inner node and new tip
  return  pll_compute_edge_loglikelihood(partition_,
                                        NEW_TIP_CLV_INDEX,
                                        PLL_SCALE_BUFFER_NONE,// scaler_index
                                        INNER_CLV_INDEX,
                                        INNER_CLV_INDEX,  // scaler_index
                                        1,// matrix index of branch TODO depends on NR
                                        0);// freq index

}


Tiny_Tree::~Tiny_Tree()
{
  if (partition_ != NULL)
  {
    // unset model params
    partition_->rates = NULL;
    partition_->subst_params = NULL;
    partition_->frequencies = NULL;
    if (partition_->clv != NULL)
    {
      // unset existing nodes clvs
      partition_->clv[OLD_LEFT_CLV_INDEX] = NULL;
      partition_->clv[OLD_RIGHT_CLV_INDEX] = NULL;
    }
    if (partition_->scale_buffer != NULL)
    {
      // unset scalers
      partition_->scale_buffer[OLD_LEFT_CLV_INDEX] = NULL;
      partition_->scale_buffer[OLD_RIGHT_CLV_INDEX] = NULL;
    }
    pll_partition_destroy(partition_);
  }
}
