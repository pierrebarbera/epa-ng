#include "Tree.hpp"

#include <stdexcept>
#include <iostream>
#include <cstdio>

#include "pll_util.hpp"
#include "epa_pll_util.hpp"
#include "file_io.hpp"
#include "Sequence.hpp"
#include "optimize.hpp"
#include "set_manipulators.hpp"
#include "logging.hpp"
#include "stringify.hpp"

using namespace std;

Tree::Tree(const string &tree_file, const MSA &msa, Model &model, Options& options)
    : ref_msa_(msa), model_(model), options_(options)
{
  tree_ = std::unique_ptr<pll_utree_t>(build_tree_from_file(tree_file, nums_));
  partition_ = std::unique_ptr<pll_partition_t>(build_partition_from_file(model_, nums_, ref_msa_.num_sites()));

  // split msa if no separate query msa was supplied
  // if (query.num_sites() == 0)
  //   split_combined_msa(ref_msa_, query_msa_, tree_.get(), nums_.tip_nodes);

  valid_map_ = vector<Range>(nums_.tip_nodes);
  link_tree_msa(tree_.get(), partition_.get(), ref_msa_, nums_.tip_nodes, valid_map_);

  // find_collapse_equal_sequences(query_msa_);

  // compute_and_set_empirical_frequencies(partition_.get(), model_);

  // perform branch length and model optimization on the reference tree
  optimize(model_, tree_.get(), partition_.get(), nums_, options_.opt_branches, options_.opt_model);

  lgr << to_string(model_);

  lgr << "Tree length: " << sum_branch_lengths(tree_.get()) << endl;

  precompute_clvs(tree_.get(), partition_.get(), nums_);

  lgr << "\nPost-optimization reference tree log-likelihood: ";
  lgr << to_string(this->ref_tree_logl()) << endl;
}

/**
  Constructs the structures from binary file.
*/
Tree::Tree(const string& bin_file, const string&, Options& options)
  : options_(options), binary_(bin_file)
{
  // tree_ = build_tree_from_file(tree_file, nums_);

  tree_ = std::unique_ptr<pll_utree_t>(binary_.load_utree());
  partition_ = std::unique_ptr<pll_partition_t>(binary_.load_partition());

  // mirror the model from the partition to the model_ object
  model_ = get_model(partition_.get());
}

Tree::~Tree()
{
  // free data segment of tree nodes
  if (tree_)
  {
    utree_free_node_data(tree_.get());
    pll_utree_destroy(tree_.get());
  }
  if (partition_)
    pll_partition_destroy(partition_.get());
}

void * Tree::get_clv(unsigned int i)
{
  assert(i < partition_.get()->tips + partition_.get()->clv_buffers);
  void* clv_ptr;
  if (i < partition_.get()->tips)
  {
    clv_ptr = partition_.get()->tipchars[i];
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
  assert(clv_ptr);
  return clv_ptr;
}

double Tree::ref_tree_logl()
{
  // TODO useless if out-of-core
  return pll_compute_edge_loglikelihood(
      partition_.get(), tree_->clv_index, tree_->scaler_index, tree_->back->clv_index,
      tree_->back->scaler_index, tree_->pmatrix_index, 0);
}
