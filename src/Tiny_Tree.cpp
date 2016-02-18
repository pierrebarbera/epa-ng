#include "Tiny_Tree.hpp"

#include <vector>
#include <cassert>

#include "pll_util.hpp"
#include "optimize.hpp"
#include "Tree_Numbers.hpp"
#include "Range.hpp"
#include "set_manipulators.hpp"

using namespace std;

Tiny_Tree::Tiny_Tree(pll_utree_t *edge_node, pll_partition_t *old_partition,
                     Model model, bool opt_branches, Range reference_tip_range)
    : opt_branches_(opt_branches), original_branch_length_(edge_node->length),
      model_(model), reference_tip_range_(reference_tip_range)

{
  assert(edge_node != NULL);
  assert(old_partition != NULL);

  pll_utree_t *old_proximal = edge_node->back;
  pll_utree_t *old_distal = edge_node;

  // detect the tip-tip case. In the tip-tip case, the reference tip should
  // always be the DISTAL
  if (!old_distal->next) {
    tip_tip_case_ = true;
  } else if (!old_proximal->next) {
    tip_tip_case_ = true;
    // do the switcheroo
    old_distal = old_proximal;
    old_proximal = old_distal->back;
  }

  tree_ = make_tiny_tree_structure(old_proximal, old_distal);

  partition_ = pll_partition_create(
      3, // tips
      1, // extra clv's
      old_partition->states, old_partition->sites,
      0, // number of mixture models
      old_partition->rate_matrices,
      3, // number of prob. matrices (one per possible unique branch length)
      old_partition->rate_cats,
      4, // number of scale buffers (one per node)
      pll_map_nt,
      old_partition->attributes);

  // shallow copy model params
  partition_->rates = old_partition->rates;
  partition_->subst_params = old_partition->subst_params;
  partition_->frequencies = old_partition->frequencies;

  // shallow copy 2 existing nodes clvs
  assert(old_partition->clv[old_proximal->clv_index] != NULL);
  assert(old_partition->clv[old_distal->clv_index] != NULL);

  // deep copy clv's
  memcpy(partition_->clv[TINY_PROXIMAL_CLV_INDEX],
         old_partition->clv[old_proximal->clv_index],
         sizeof(double) * old_partition->sites * old_partition->rate_cats *
             old_partition->states);
  memcpy(partition_->clv[TINY_DISTAL_CLV_INDEX],
         old_partition->clv[old_distal->clv_index],
         sizeof(double) * old_partition->sites * old_partition->rate_cats *
             old_partition->states);

  // deep copy scalers
  if (old_proximal->scaler_index != PLL_SCALE_BUFFER_NONE)
    memcpy(partition_->scale_buffer[TINY_PROXIMAL_CLV_INDEX],
           old_partition->scale_buffer[old_proximal->scaler_index],
           sizeof(unsigned int) * old_partition->sites);
  if (old_distal->scaler_index != PLL_SCALE_BUFFER_NONE)
    memcpy(partition_->scale_buffer[TINY_DISTAL_CLV_INDEX],
           old_partition->scale_buffer[old_distal->scaler_index],
           sizeof(unsigned int) * old_partition->sites);

  // precreate some of the operation fields that are static throughout the trees
  // lifetime
  ops_.parent_clv_index = TINY_INNER_CLV_INDEX;
  ops_.child1_clv_index = TINY_PROXIMAL_CLV_INDEX;
  ops_.child1_scaler_index = tree_->next->next->back->scaler_index;
  ops_.child2_clv_index = TINY_DISTAL_CLV_INDEX;
  ops_.child2_scaler_index = tree_->next->back->scaler_index;
  ops_.parent_scaler_index = TINY_INNER_CLV_INDEX;

  /* heuristic insertion as described in EPA paper from 2011 (Berger et al.):
    original branch, now split by "inner", or base, node of the inserted sequence,
    defines the new branch lengths between inner and old left/right respectively
    as old branch length / 2.
    The new branch leading from inner to the new tip is initialized with length 0.9,
    which is the default branch length in RAxML.
  */ // TODO change this
  // wether heuristic is used or not, this is the initial branch length
  // configuration
  double branch_lengths[2] = {old_proximal->length / 2, DEFAULT_BRANCH_LENGTH};
  unsigned int matrix_indices[2] = {0, 1};

  // use branch lengths to compute the probability matrices
  // TODO replace with updating single matrix function
  // alternatively have modified version that recognizes the default branch length
  pll_update_prob_matrices(partition_, 0, matrix_indices, branch_lengths, 2);

  ops_.child1_matrix_index = 0;
  ops_.child2_matrix_index = 0;

  if (!opt_branches_)
    pll_update_partials(partition_, &ops_, TINY_NUM_OPS);
  // use update_partials to compute the clv pointing toward the new tip
}

Tiny_Tree::~Tiny_Tree() {
  if (partition_ != nullptr) {
    // unset shallow copied things
    partition_->rates = nullptr;
    partition_->subst_params = nullptr;
    partition_->frequencies = nullptr;

    pll_partition_destroy(partition_);
  }

  if (tree_ != nullptr)
    pll_utree_destroy(tree_);
}

std::tuple<double, double, double> Tiny_Tree::place(const Sequence &s) {
  assert(partition_ != NULL);

  // init the new tip with s.sequence(), branch length
  pll_set_tip_states(partition_, TINY_NEW_TIP_CLV_INDEX, pll_map_nt,
                     s.sequence().c_str());

  unsigned int inner_matrix_index = 1;
  auto distal_length = tree_->next->length;
  auto pendant_length = tree_->length;
  double logl = 0.0;

  if (opt_branches_) {
    Range range(0, partition_->sites);

    /* differentiate between the normal case and the tip tip case:
      in the normal case we want to compute the partial toward the newly palced sequence.
      In other words, we set the virtual root as the node whose back-neighbour is the new
      sequence, which is tree_. (*)
      */
    auto virtual_root = tree_;

    if (tip_tip_case_)
    {
      /* (cont. from (*))... however in the tip-tip case we want that virtual root to be toward the non-tip
        node of the reference tree. Thus virtual_root needs to reflect that.
        optimize_branch_triplet then takes care of the single computation operation for us
        using that node.*/
      virtual_root = tree_->next->next;

      /* setting the actual range of the ranged CLV computation to be the
        superset of the two tips ranges*/
      range = superset(get_valid_range(s.sequence()), reference_tip_range_);
      // range = get_valid_range(s.sequence());
    }

    logl = optimize_branch_triplet_ranged(partition_, virtual_root, range);

    assert(tree_->length >= 0);
    assert(tree_->next->length >= 0);
    assert(tree_->next->next->length >= 0);

    // rescale the distal length, as it has likely changed during optimization
    // done as in raxml
    double proximal_length = tree_->next->next->length;
    double new_total_branch_length = distal_length + proximal_length;
    distal_length = (original_branch_length_ / new_total_branch_length) * distal_length;
    pendant_length = tree_->length;

    // reset branch lenths to their default for the next insertion
    tree_->length = DEFAULT_BRANCH_LENGTH;
    tree_->back->length = DEFAULT_BRANCH_LENGTH;
    auto half_original = original_branch_length_ / 2;
    tree_->next->length = half_original;
    tree_->next->back->length = half_original;
    tree_->next->next->length = half_original;
    tree_->next->next->back->length = half_original;

    tree_->pmatrix_index = 1;
    tree_->back->pmatrix_index = 1;
    tree_->next->pmatrix_index = 0;
    tree_->next->back->pmatrix_index = 0;
    tree_->next->next->pmatrix_index = 0;
    tree_->next->next->back->pmatrix_index = 0;

    double branch_lengths[2] = {half_original, DEFAULT_BRANCH_LENGTH};
    unsigned int matrix_indices[2] = {0, 1};
    pll_update_prob_matrices(partition_, 0, matrix_indices, branch_lengths, 2);
  } else
    logl = pll_compute_edge_loglikelihood(partition_,
                                          TINY_NEW_TIP_CLV_INDEX,
                                          PLL_SCALE_BUFFER_NONE, // scaler_index
                                          TINY_INNER_CLV_INDEX,
                                          TINY_INNER_CLV_INDEX, // scaler_index
                                          inner_matrix_index,
                                          0); // freq index

  assert(distal_length <= original_branch_length_);

  return make_tuple(logl, distal_length, pendant_length);
}
