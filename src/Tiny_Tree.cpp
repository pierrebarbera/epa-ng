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

  const pll_utree_t *old_proximal = edge_node->back;
  const pll_utree_t *old_distal = edge_node;

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

  tree_ = make_tiny_tree_structure(old_proximal, old_distal, tip_tip_case_);

  partition_ = make_tiny_partition(old_partition, tree_, old_proximal, old_distal, tip_tip_case_);

  // operation for computing the clv toward the new tip (for initialization and logl in non-blo case)
  auto distal = tree_->next->back;
  auto proximal = tree_->next->next->back;
  pll_operation_t op;
  op.parent_clv_index = tree_->clv_index;
  op.child1_clv_index = distal->clv_index;
  op.child1_scaler_index = distal->scaler_index;
  op.child2_clv_index = proximal->clv_index;
  op.child2_scaler_index = proximal->scaler_index;
  op.parent_scaler_index = tree_->scaler_index;
  op.child1_matrix_index = distal->pmatrix_index;
  op.child2_matrix_index = proximal->pmatrix_index;

  // wether heuristic is used or not, this is the initial branch length configuration
  double branch_lengths[3] = {proximal->length, distal->length, tree_->length};
  unsigned int matrix_indices[3] = {proximal->pmatrix_index, distal->pmatrix_index, tree_->pmatrix_index};

  // use branch lengths to compute the probability matrices
  pll_update_prob_matrices(partition_, 0, matrix_indices, branch_lengths, 3);


  if (!opt_branches_)
    pll_update_partials(partition_, &op, 1);
}
// use update_partials to compute the clv pointing toward the new tip

Tiny_Tree::~Tiny_Tree() {
  if (partition_ != nullptr) {
    // unset shallow copied things
    partition_->rates = nullptr;
    partition_->subst_params = nullptr;
    partition_->frequencies = nullptr;
    partition_->eigenvecs = nullptr;
    partition_->inv_eigenvecs = nullptr;
    partition_->eigenvals = nullptr;
    partition_->prop_invar = nullptr;
    partition_->eigen_decomp_valid = nullptr;
    partition_->pattern_weights = nullptr;

    partition_->lh_statepair = nullptr;
    partition_->charmap = nullptr;
    partition_->revmap = nullptr;

    pll_partition_destroy(partition_);
  }

  if (tree_ != nullptr)
    pll_utree_destroy(tree_);
}

std::tuple<double, double, double> Tiny_Tree::place(const Sequence &s) {
  assert(partition_ != NULL);

  // init the new tip with s.sequence(), branch length
  pll_set_tip_states(partition_, tree_->back->clv_index, pll_map_nt,
                     s.sequence().c_str());

  auto distal_length = tree_->next->length;
  auto pendant_length = tree_->length;
  double logl = 0.0;

  if (opt_branches_) {
    Range range(0, partition_->sites);

    /* differentiate between the normal case and the tip tip case:
      in the normal case we want to compute the partial toward the newly placed sequence.
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
      // range = superset(get_valid_range(s.sequence()), reference_tip_range_);
      range = get_valid_range(s.sequence());
    }

    // optimize the branches using pnly the portion of the sites specified by range
    logl = call_focused(partition_, range, optimize_branch_triplet, virtual_root);


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
    double half_original = original_branch_length_ / 2.0;
    tree_->next->length = half_original;
    tree_->next->back->length = half_original;
    tree_->next->next->length = half_original;
    tree_->next->next->back->length = half_original;

    tree_->pmatrix_index = 2;
    tree_->back->pmatrix_index = 2;
    tree_->next->pmatrix_index = 1;
    tree_->next->back->pmatrix_index = 1;
    tree_->next->next->pmatrix_index = 0;
    tree_->next->next->back->pmatrix_index = 0;

    double branch_lengths[3] = {half_original, half_original, DEFAULT_BRANCH_LENGTH};
    unsigned int matrix_indices[3] = {0, 1, 2};
    pll_update_prob_matrices(partition_, 0, matrix_indices, branch_lengths, 3);
  } else
    logl = pll_compute_edge_loglikelihood(partition_,
                                          tree_->back->clv_index,
                                          PLL_SCALE_BUFFER_NONE, // scaler_index
                                          tree_->clv_index,
                                          tree_->scaler_index, // scaler_index
                                          tree_->pmatrix_index,
                                          0); // freq index

  assert(distal_length <= original_branch_length_);
  assert(distal_length >= 0.0);

  return make_tuple(logl, distal_length, pendant_length);
}
