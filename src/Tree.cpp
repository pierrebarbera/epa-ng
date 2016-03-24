#include "Tree.hpp"

#include <stdexcept>
#include <iostream>
#include <omp.h>

#include "pll_util.hpp"
#include "epa_pll_util.hpp"
#include "file_io.hpp"
#include "Sequence.hpp"
#include "optimize.hpp"
#include "set_manipulators.hpp"
#include "logging.hpp"
#include "stringify.hpp"

using namespace std;

Tree::Tree(const string &tree_file, const MSA &msa, Model &model,
           Options options, MSA &query)
    : ref_msa_(msa), query_msa_(query), model_(model), options_(options), out_of_core_(false)
{
  tree_ = build_tree_from_file(tree_file, nums_);
  partition_ = build_partition_from_file(model_, nums_, ref_msa_.num_sites());

  // split msa if no separate query msa was supplied
  if (query.num_sites() == 0)
    split_combined_msa(ref_msa_, query_msa_, tree_, nums_.tip_nodes);

  valid_map_ = vector<Range>(nums_.tip_nodes);
  link_tree_msa(tree_, partition_, ref_msa_, nums_.tip_nodes, valid_map_);

  find_collapse_equal_sequences(query_msa_);

  // compute_and_set_empirical_frequencies(partition_, model_);

  // perform branch length and model optimization on the reference tree
  optimize(model_, tree_, partition_, nums_, options_.opt_branches, options_.opt_model);

  lgr << to_string(model_);

  lgr << "Tree length: " << sum_branch_lengths(tree_) << endl;

  precompute_clvs(tree_, partition_, nums_);

  lgr << "\nPost-optimization reference tree log-likelihood: ";
  lgr << to_string(this->ref_tree_logl()) << endl;
}

/**
  Constructs the structures from binary files and sets the relevant out-of-core structures.
  Mainly this involves using MSA_Stream and having the partition only be as big as needed
*/
Tree::Tree(const string& bin_file, const string& tree_file, Options& options) : options_(options), out_of_core_(true)
{
  tree_ = build_tree_from_file(tree_file, nums_);
  partition_ = build_partition_from_binary(bin_file, out_of_core_);

}

double Tree::ref_tree_logl()
{
  // TODO useless if out-of-core
  return pll_compute_edge_loglikelihood(
      partition_, tree_->clv_index, tree_->scaler_index, tree_->back->clv_index,
      tree_->back->scaler_index, tree_->pmatrix_index, 0);
}

Tree::~Tree()
{
  // free data segment of tree nodes
  utree_free_node_data(tree_);
  pll_partition_destroy(partition_);
  pll_utree_destroy(tree_);
}

Sample Tree::place()
{
  const auto num_branches = nums_.branches;
  const auto num_queries = query_msa_.size();
  // get all edges
  vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(tree_, &branches[0]);
  assert(num_traversed_branches == num_branches);

  lgr << "\nPlacing "<< to_string(num_queries) << " reads on " <<
    to_string(num_branches) << " branches." << endl;

  // build all tiny trees with corresponding edges
  vector<Tiny_Tree> insertion_trees;
  for (unsigned int branch_id = 0; branch_id < num_branches; ++branch_id)
    insertion_trees.emplace_back(branches[branch_id], branch_id, partition_, model_, !options_.prescoring);
    /* clarification: last arg here is a flag specifying whether to optimize the branches.
      we don't want that if the mode is prescoring */

  // output class
  Sample sample(get_numbered_newick_string(tree_));
  for (unsigned int sequence_id = 0; sequence_id < num_queries; sequence_id++)
    sample.emplace_back(sequence_id, num_branches);


  // place all s on every edge
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int branch_id = 0; branch_id < num_branches; ++branch_id)
  {
    for (unsigned int sequence_id = 0; sequence_id < num_queries; ++sequence_id)
    {
      sample[sequence_id][branch_id] = insertion_trees[branch_id].place(query_msa_[sequence_id]);
    }
  }
  // now that everything has been placed, we can compute the likelihood weight ratio
  compute_and_set_lwr(sample);

  /* prescoring was chosen: perform a second round, but only on candidate edges identified
    during the first run */
  if (options_.prescoring)
  {
    lgr << "Entering second phase of placement. \n";
    if (options_.prescoring_by_percentage)
      discard_bottom_x_percent(sample, (1.0 - options_.prescoring_threshold));
    else
      discard_by_accumulated_threshold(sample, options_.prescoring_threshold);

    // build a list of placements per edge that need to be recomputed
    vector<vector<tuple<Placement *, const unsigned int>>> recompute_list(num_branches);
    for (auto & pq : sample)
      for (auto & placement : pq)
        recompute_list[placement.branch_id()].push_back(make_tuple(&placement, pq.sequence_id()));

    #pragma omp parallel for schedule(dynamic)
    for (unsigned int branch_id = 0; branch_id < num_branches; branch_id++)
    {
      Placement * placement;
      // Sequence * sequence;
      auto& branch = insertion_trees[branch_id];
      branch.opt_branches(true); // TODO only needs to be done once
      for (auto recomp_tuple : recompute_list[branch_id])
      {
        placement = get<0>(recomp_tuple);
        *placement = branch.place(query_msa_[get<1>(recomp_tuple)]);

      }
    }
    compute_and_set_lwr(sample);
  }

  // finally, trim the output
  if (options_.acc_threshold)
    discard_by_accumulated_threshold(sample, options_.support_threshold);
  else
    discard_by_support_threshold(sample, options_.support_threshold);

  return sample;
}
