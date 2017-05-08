#include "Tiny_Tree.hpp"

#ifdef __OMP
#include <omp.h>
#endif

#include <vector>
#include <numeric>

#include "tiny_util.hpp"
#include "pll_util.hpp"
#include "optimize.hpp"
#include "Tree_Numbers.hpp"
#include "Range.hpp"
#include "set_manipulators.hpp"
#include "logging.hpp"

static void precompute_sites_static(char nt,
                                    std::vector<double>& result,
                                    pll_partition_t* partition,
                                    pll_utree_t* tree,
                                    Model& model)
{
  const size_t sites = partition->sites;
  result.clear();
  result.resize(sites);
  std::string seq(sites, nt);

  std::vector<unsigned int> param_indices(model.rate_cats(), 0);

  auto err_check = pll_set_tip_states(partition, tree->back->clv_index, model.char_map(),
                     seq.c_str());

  if (err_check == PLL_FAILURE) {
    throw std::runtime_error{
      std::string("Set tip states during sites precompution failed! pll_errmsg: ")
      + pll_errmsg
    };
  }

  pll_compute_edge_loglikelihood( partition,
                                  tree->back->clv_index,
                                  PLL_SCALE_BUFFER_NONE, 
                                  tree->clv_index,
                                  tree->scaler_index,
                                  tree->pmatrix_index,
                                  &param_indices[0], &result[0]);
}


Tiny_Tree::Tiny_Tree( pll_utree_t * edge_node, 
                      const unsigned int branch_id, 
                      Tree& reference_tree, 
                      const bool opt_branches, 
                      const Options& options, 
                      std::shared_ptr<Lookup_Store>& lookup_store)
  : partition_(nullptr, tiny_partition_destroy)
  , tree_(nullptr, utree_destroy)
  , opt_branches_(opt_branches)
  , model_(reference_tree.model())
  , ranged_computation_(options.ranged)
  , sliding_blo_(options.sliding_blo)
  , branch_id_(branch_id)
  , lookup_(lookup_store)
{
  original_branch_length_ = (edge_node->length < 2*PLLMOD_OPT_MIN_BRANCH_LEN) ?
    2*PLLMOD_OPT_MIN_BRANCH_LEN : edge_node->length;
  assert(edge_node);

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

  tree_ = std::unique_ptr<pll_utree_t, utree_deleter>(
      	                    make_tiny_tree_structure(old_proximal, old_distal, tip_tip_case_),
                            utree_destroy);

  partition_ = std::unique_ptr<pll_partition_t, partition_deleter>(
                                make_tiny_partition(reference_tree, tree_.get(), old_proximal, old_distal, tip_tip_case_),
                                tiny_partition_destroy);

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
  std::vector<unsigned int> param_indices(model_.rate_cats(), 0);
  pll_update_prob_matrices(partition_.get(), &param_indices[0], matrix_indices, branch_lengths, 3);

  // use update_partials to compute the clv pointing toward the new tip
  pll_update_partials(partition_.get(), &op, 1);

  if (!opt_branches) {
    // LOG_DBG << "precomputation for branch " << branch_id << std::endl;
    const std::lock_guard<std::mutex> lock(lookup_store->get_mutex(branch_id));

    if( not lookup_store->has_branch(branch_id) ) {
      const auto size = lookup_store->char_map_size();

      // precompute all possible site likelihoods
      std::vector<std::vector<double>> precomputed_sites(size);
      for (size_t i = 0; i < size; ++i) {
        precompute_sites_static(lookup_store->char_map(i),
                                precomputed_sites[i],
                                partition_.get(),
                                tree_.get(),
                                model_);
      }
      lookup_store->init_branch(branch_id, precomputed_sites);
    }
  }
  
}

Placement Tiny_Tree::place(const Sequence &s) 
{
  assert(partition_);
  assert(tree_);

  auto distal_length = tree_->next->length;
  auto pendant_length = tree_->length;
  double logl = 0.0;
  std::vector<unsigned int> param_indices(model_.rate_cats(), 0);

  Range range(0, partition_->sites);
  if (ranged_computation_) {
    range = get_valid_range(s.sequence());
    // range = superset(get_valid_range(s.sequence()), reference_tip_range_);
  }

  if (opt_branches_) {

    /* differentiate between the normal case and the tip tip case:
      in the normal case we want to compute the partial toward the newly placed sequence.
      In other words, we set the virtual root as the node whose back-neighbour is the new
      sequence, which is tree_. (*)
      */
    auto virtual_root = tree_.get();

    if (tip_tip_case_) // TODO cant be used with pplacer blo
    {
      /* (cont. from (*))... however in the tip-tip case we want that virtual root to be toward the non-tip
        node of the reference tree. Thus virtual_root needs to reflect that.
        optimize_branch_triplet then takes care of the single computation operation for us
        using that node.*/
      // virtual_root = tree_->next->next;
    }

    // init the new tip with s.sequence(), branch length
    auto err_check = pll_set_tip_states(partition_.get(), tree_->back->clv_index, model_.char_map(),
                       s.sequence().c_str());

    if (err_check == PLL_FAILURE) {
      throw std::runtime_error{"Set tip states during placement failed!"};
    }

    // optimize the branches using pnly the portion of the sites specified by range
    logl = call_focused(partition_.get(), range, optimize_branch_triplet, virtual_root, sliding_blo_);
    // logl = optimize_branch_triplet(partition_.get(), virtual_root, sliding_blo_);

    assert(tree_->length >= 0);
    assert(tree_->next->length >= 0);
    assert(tree_->next->next->length >= 0);

    // rescale the distal length, as it has likely changed during optimization
    // done as in raxml
    double proximal_length = tree_->next->next->length;
    double new_total_branch_length = distal_length + proximal_length;
    distal_length = (original_branch_length_ / new_total_branch_length) * distal_length;
    pendant_length = tree_->length;

    reset_triplet_lengths(tree_.get(), partition_.get(), original_branch_length_);
    
    // re-update the partial
    auto child1 = virtual_root->next->back;
    auto child2 = virtual_root->next->next->back;

    pll_operation_t op;
    op.parent_clv_index = virtual_root->clv_index;
    op.parent_scaler_index = virtual_root->scaler_index;
    op.child1_clv_index = child1->clv_index;
    op.child1_scaler_index = child1->scaler_index;
    op.child1_matrix_index = child1->pmatrix_index;
    op.child2_clv_index = child2->clv_index;
    op.child2_scaler_index = child2->scaler_index;
    op.child2_matrix_index = child2->pmatrix_index;

    pll_update_partials(partition_.get(), &op, 1);

  } else {
    logl = lookup_->sum_precomputed_sitelk(branch_id_, s.sequence());
  }

  assert(distal_length <= original_branch_length_);
  assert(distal_length >= 0.0);

  return Placement(branch_id_, logl, pendant_length, distal_length);
}
