#include "Tiny_Tree.hpp"

#include <vector>
#include <algorithm>

#include "tiny_util.hpp"
#include "pll_util.hpp"
#include "optimize.hpp"
#include "Tree_Numbers.hpp"
#include "Range.hpp"
#include "set_manipulators.hpp"

using namespace std;

static constexpr size_t NT_A = 0;
static constexpr size_t NT_C = 1;
static constexpr size_t NT_G = 2;
static constexpr size_t NT_T = 3;
static constexpr size_t NT_GAP = 4;
static constexpr char NT_MAP[5] = {'A', 'C', 'G', 'T', '-'};

static void precompute_sites_static(char nt, vector<double>& result, 
  pll_partition_t* partition, pll_utree_t* tree, Model& model)
{
  const size_t sites = partition->sites;
  result.clear();
  result.resize(sites);
  string seq(sites, nt);

  vector<unsigned int> param_indices(model.rate_cats(), 0);

  auto err_check = pll_set_tip_states(partition, tree->back->clv_index, model.char_map(),
                     seq.c_str());

  if (err_check == PLL_FAILURE)
    throw runtime_error{"Set tip states during sites precompution failed!"};

  pll_compute_edge_loglikelihood(partition,
                                  tree->back->clv_index,
                                  PLL_SCALE_BUFFER_NONE, 
                                  tree->clv_index,
                                  tree->scaler_index,
                                  tree->pmatrix_index,
                                  &param_indices[0], &result[0]);
}

static double sum_precomputed_sitelk(vector<vector<double>>& lookup, const Sequence& s)
{
  string seq = s.sequence();
  assert(seq.length() == lookup[NT_A].size());
  assert(lookup[NT_G].size() == lookup[NT_A].size());
  assert(lookup[NT_C].size() == lookup[NT_A].size());
  assert(lookup[NT_T].size() == lookup[NT_A].size());
  assert(lookup[NT_GAP].size() == lookup[NT_A].size());

  // transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
  double sum = 0;
  for (size_t i = 0; i < seq.length(); ++i)
  {
    size_t c;
    switch (seq[i])
    {
      case 'A':
      case 'a':
        c = NT_A;
        break;
      case 'C':
      case 'c':
        c = NT_C;
        break;
      case 'G':
      case 'g':
        c = NT_G;
        break;
      case 'T':
      case 't':
        c = NT_T;
        break;
      case '-':
        c = NT_GAP;
        break;
      default:
        throw runtime_error{"derp"};
    }
    sum += lookup[c][i];
  }
  return sum;
}

Tiny_Tree::Tiny_Tree(pll_utree_t *edge_node, unsigned int branch_id, Tree& reference_tree,
                     bool opt_branches, Range reference_tip_range, bool ranged)
    : partition_(nullptr, tiny_partition_destroy), tree_(nullptr, utree_destroy), opt_branches_(opt_branches)
    , model_(reference_tree.model())
    , reference_tip_range_(reference_tip_range), ranged_computation_(ranged), branch_id_(branch_id)
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

  tree_ = unique_ptr<pll_utree_t, utree_deleter>(
      	                    make_tiny_tree_structure(old_proximal, old_distal, tip_tip_case_),
                            utree_destroy);

  partition_ = unique_ptr<pll_partition_t, partition_deleter>(
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
  vector<unsigned int> param_indices(model_.rate_cats(), 0);
  pll_update_prob_matrices(partition_.get(), &param_indices[0], matrix_indices, branch_lengths, 3);


  if (!opt_branches)
  {
    // use update_partials to compute the clv pointing toward the new tip
    pll_update_partials(partition_.get(), &op, 1);
   
    lookup_.clear();
    lookup_.resize(5);

    // precompute all possible site likelihoods
    size_t i = 0;
    for (char nt : NT_MAP)
      precompute_sites_static(nt, lookup_[i++], partition_.get(), tree_.get(), model_); 
  }
  
}

Placement Tiny_Tree::place(const Sequence &s) 
{
  assert(partition_);
  assert(tree_);

  auto distal_length = tree_->next->length;
  auto pendant_length = tree_->length;
  double logl = 0.0;
  vector<unsigned int> param_indices(model_.rate_cats(), 0);

  Range range(0, partition_->sites);
  if (ranged_computation_)
  {
    range = get_valid_range(s.sequence());
    // range = superset(get_valid_range(s.sequence()), reference_tip_range_);
  }

  if (opt_branches_)
  {

    /* differentiate between the normal case and the tip tip case:
      in the normal case we want to compute the partial toward the newly placed sequence.
      In other words, we set the virtual root as the node whose back-neighbour is the new
      sequence, which is tree_. (*)
      */
    auto virtual_root = tree_.get();

    if (tip_tip_case_)
    {
      /* (cont. from (*))... however in the tip-tip case we want that virtual root to be toward the non-tip
        node of the reference tree. Thus virtual_root needs to reflect that.
        optimize_branch_triplet then takes care of the single computation operation for us
        using that node.*/
      virtual_root = tree_->next->next;
    }

    // init the new tip with s.sequence(), branch length
    auto err_check = pll_set_tip_states(partition_.get(), tree_->back->clv_index, model_.char_map(),
                       s.sequence().c_str());

    if (err_check == PLL_FAILURE)
      throw runtime_error{"Set tip states during placement failed!"};

    // optimize the branches using pnly the portion of the sites specified by range
    logl = call_focused(partition_.get(), range, optimize_branch_triplet, virtual_root);

    logl = pll_compute_edge_loglikelihood(partition_,
                                  virtual_root->back->clv_index,
                                  PLL_SCALE_BUFFER_NONE, 
                                  virtual_root->clv_index,
                                  virtual_root->scaler_index,
                                  virtual_root->pmatrix_index,
                                  &param_indices[0], nullptr);

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
  }
  else
  {
    logl = sum_precomputed_sitelk(lookup_, s);
  }
  // logl = call_focused(partition_.get(), range, pll_compute_edge_loglikelihood,
  //                                       tree_->back->clv_index,
  //                                       PLL_SCALE_BUFFER_NONE, // scaler_index
  //                                       tree_->clv_index,
  //                                       tree_->scaler_index, // scaler_index
  //                                       tree_->pmatrix_index,
  //                                       &param_indices[0], nullptr); // freq index

  assert(distal_length <= original_branch_length_);
  assert(distal_length >= 0.0);

  return Placement(branch_id_, logl, pendant_length, distal_length);
}
