#include "Tree.hpp"

#include <stdexcept>
#include <iostream>
#include <omp.h>

#include "pll_util.hpp"
#include "epa_pll_util.hpp"
#include "file_io.hpp"
#include "Sequence.hpp"
#include "optimize.hpp"
#include "calculation.hpp"

using namespace std;

Tree::Tree(const string& tree_file, const MSA& msa, Model& model, Options options,
  const MSA& query)
  : ref_msa_(msa), query_msa_(query), model_(model), options_(options)
{
  //parse, build tree
  nums_ = Tree_Numbers();
  tie(partition_, tree_) = build_partition_from_file(tree_file, model_, nums_, ref_msa_.num_sites());

  // split msa if no separate query msa was supplied
  if (query.num_sites() == 0)
    split_combined_msa(ref_msa_, query_msa_, tree_, nums_.tip_nodes);

  valid_map_ = vector<tuple<unsigned int, unsigned int>>(nums_.tip_nodes);
  link_tree_msa(tree_, partition_, ref_msa_, nums_.tip_nodes, valid_map_);

  compute_and_set_empirical_frequencies(partition_);

  // perform branch length and model optimization on the reference tree
  optimize(model_, tree_, partition_, nums_, options_.opt_branches, options_.opt_model);

  precompute_clvs(tree_, partition_, nums_);
}

double Tree::ref_tree_logl()
{
  return pll_compute_edge_loglikelihood (partition_, tree_->clv_index,
                                                tree_->scaler_index,
                                                tree_->back->clv_index,
                                                tree_->back->scaler_index,
                                                tree_->pmatrix_index, 0);
}

Tree::~Tree()
{
  // free data segment of tree nodes
  utree_free_node_data(tree_);
	pll_partition_destroy(partition_);
  pll_utree_destroy(tree_);

}

PQuery_Set Tree::place() const
{
  const auto num_branches = nums_.branches;
  // get all edges
  vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(tree_, &branches[0]);
  assert(num_traversed_branches == num_branches);

  // build all tiny trees with corresponding edges
  vector<Tiny_Tree> insertion_trees;
  for (auto node : branches)
    insertion_trees.emplace_back(node, partition_, model_, !options_.prescoring);
    /* clarification: last arg here is a flag specifying whether to optimize the branches.
      we don't want that if the mode is prescoring */

  // output class
  auto num_queries = query_msa_.size();
  PQuery_Set pquerys(get_numbered_newick_string(tree_));
  for (const auto & s : query_msa_)
    pquerys.emplace_back(s, num_branches);

  // place all s on every edge
  double logl, distal, pendant;

  // omp_set_num_threads(4);
  #pragma omp parallel for
  for (unsigned int branch_id = 0; branch_id < num_branches; ++branch_id)
  {
    for (unsigned int sequence_id = 0; sequence_id < num_queries; ++sequence_id)
    {
      tie(logl, distal, pendant) = insertion_trees[branch_id].place(query_msa_.get(sequence_id));
      pquerys[sequence_id][branch_id] = Placement(
        branch_id,
        logl,
        pendant,
        distal);
    }
  }
  // now that everything has been placed, we can compute the likelihood weight ratio
  compute_and_set_lwr(pquerys);

  /* prescoring was chosen: perform a second round, but only on candidate edges identified
    during the first run */
  // TODO for MPI: think about how to split up the work, probably build list of sequences
  // per branch, then pass
  if (options_.prescoring)
  {
    discard_by_accumulated_threshold(pquerys, 0.95);// TODO outside input? constant?
    for (auto &pq : pquerys)
      for (auto &placement : pq)
      {
        unsigned int id = placement.branch_id();
        insertion_trees[id].opt_branches(true); // TODO only needs to be done once
        tie(logl, distal, pendant) = insertion_trees[id].place(pq.sequence());
        placement.likelihood(logl);
        placement.pendant_length(pendant);
        placement.distal_length(distal);
      }
    compute_and_set_lwr(pquerys);
  }

  // finally, trim the output
  if (options_.acc_threshold)
    discard_by_accumulated_threshold(pquerys, options_.support_threshold);
  else
    discard_by_support_threshold(pquerys, options_.support_threshold);

  return pquerys;
}
