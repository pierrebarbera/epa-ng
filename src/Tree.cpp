#include "Tree.hpp"

#include <stdexcept>
#include <iostream>
#include <cstdio>

#include "epa_pll_util.hpp"
#include "file_io.hpp"
#include "Sequence.hpp"
#include "optimize.hpp"
#include "set_manipulators.hpp"
#include "Log.hpp"
#include "stringify.hpp"

using namespace std;

Tree::Tree(const string &tree_file, const MSA &msa, Model &model, Options& options)
  : partition_(nullptr, pll_partition_destroy), tree_(nullptr, utree_destroy)
  , ref_msa_(msa), model_(model), options_(options)
{
  tree_ = unique_ptr<pll_utree_t, utree_deleter>(
                      build_tree_from_file(tree_file, nums_),
                      utree_destroy);
  partition_ = unique_ptr<pll_partition_t, partition_deleter>(
                            build_partition_from_file(model_, nums_, ref_msa_.num_sites()),
                            pll_partition_destroy);

  locks_ = Mutex_List(partition_->tips + partition_->clv_buffers);

  valid_map_ = vector<Range>(nums_.tip_nodes);
  link_tree_msa(tree_.get(), partition_.get(), model_, ref_msa_, nums_.tip_nodes, valid_map_);

  // find_collapse_equal_sequences(query_msa_);

  // perform branch length and model optimization on the reference tree
  optimize(model_, tree_.get(), partition_.get(), nums_, options_.opt_branches, options_.opt_model);

  lgr.dbg() << to_string(model_);

  lgr.dbg() << "Tree length: " << sum_branch_lengths(tree_.get()) << endl;

  precompute_clvs(tree_.get(), partition_.get(), nums_);

  lgr.dbg() << "\nPost-optimization reference tree log-likelihood: ";
  lgr.dbg() << to_string(this->ref_tree_logl()) << endl;
}

/**
  Constructs the structures from binary file.
*/
Tree::Tree(const string& bin_file, Model &model, Options& options)
  : partition_(nullptr, pll_partition_destroy), tree_(nullptr, utree_destroy), model_(model)
  , options_(options), binary_(bin_file)
{
  tree_ = unique_ptr<pll_utree_t, utree_deleter>(
                      binary_.load_utree(),
                      utree_destroy);
  partition_ = unique_ptr<pll_partition_t, partition_deleter>(
                            binary_.load_partition(),
                            pll_partition_destroy);

  locks_ = Mutex_List(partition_->tips + partition_->clv_buffers);

  nums_.init(partition_->tips);
}

/**
  Returns a pointer either to the CLV or tipchar buffer, depending on the index.
  If they are not currently in memory, fetches them from file.
  Ensures that associated scalers are allocated and ready on return.
*/
void * Tree::get_clv(const pll_utree_t* node)
{
  auto i = node->clv_index;

  // prevent race condition from concurrent access to this function
  Scoped_Mutex lock_by_clv_id(locks_[i]);
  // printf("Got lock %d!\n", i);

  auto scaler = node->scaler_index;
  bool use_tipchars = partition_->attributes & PLL_ATTRIB_PATTERN_TIP;

  if(i >= partition_->tips + partition_->clv_buffers)
    throw runtime_error{"Node index out of bounds"};

  void* clv_ptr;
  if (use_tipchars && i < partition_->tips)
  {
    clv_ptr = partition_->tipchars[i];
    // dynamically load from disk if not in memory
    if(!clv_ptr)
    {
      binary_.load_tipchars(partition_.get(), i);
      clv_ptr = partition_->tipchars[i];
    }
  }
  else
  {
    clv_ptr = partition_->clv[i];
    // dynamically load from disk if not in memory
    if(!clv_ptr)
    {
      binary_.load_clv(partition_.get(), i);
      clv_ptr = partition_->clv[i];
    }
  }

  // dynamically load the scaler if needed
  if(scaler != PLL_SCALE_BUFFER_NONE && !(partition_->scale_buffer[scaler]))
    binary_.load_scaler(partition_.get(), scaler);

  assert(clv_ptr);
  
  // printf("Release lock %d!\n", i);

  return clv_ptr;
}

double Tree::ref_tree_logl()
{
  vector<unsigned int> param_indices(model_.rate_cats(), 0);
  // ensure clvs are there
  this->get_clv(tree_.get());
  this->get_clv(tree_.get()->back);

  return pll_compute_edge_loglikelihood(
      partition_.get(), tree_->clv_index, tree_->scaler_index, tree_->back->clv_index,
      tree_->back->scaler_index, tree_->pmatrix_index, &param_indices[0], nullptr);
}
