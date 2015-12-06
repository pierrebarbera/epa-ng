#include "Tree.hpp"

#include <stdexcept>
#include <tuple>
#include <vector>

#include "pll_util.hpp"
#include "epa_pll_util.hpp"
#include "file_io.hpp"
#include "Sequence.hpp"
#include "optimize.hpp"
#include "calculation.hpp"

using namespace std;

Tree::Tree(const string& tree_file, const MSA& msa, const Model& model, const bool heuristic,
  const MSA& query)
  : ref_msa_(msa), query_msa_(query), model_(model), heuristic_(heuristic)
{
  //parse, build tree
  nums_ = Tree_Numbers();
  tie(partition_, tree_) = build_partition_from_file(tree_file, model_, nums_, ref_msa_.num_sites());

  // split msa if no separate query msa was supplied
  if (query.num_sites() == 0)
    split_combined_msa(ref_msa_, query_msa_, tree_, nums_.tip_nodes);

  link_tree_msa(tree_, partition_, ref_msa_, nums_.tip_nodes);

  // perform branch length optimization on the reference tree
  optimize(tree_, partition_, nums_);

  precompute_clvs(tree_, partition_, nums_);
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
  // get all edges
  vector<pll_utree_t *> branches(nums_.branches);
  auto num_traversed = utree_query_branches(tree_, &branches[0]);

  // build all tiny trees with corresponding edges
  vector<Tiny_Tree> insertion_trees;
  for (auto node : branches)
    insertion_trees.emplace_back(node, partition_, model_, heuristic_);

  // output class
  PQuery_Set pquerys(get_numbered_newick_string(tree_));

  // place all s on every edge
  double logl, distal, pendant;
  for (auto const &s : query_msa_)// make sure a reference, not a copy, is returned
  {
    pquerys.emplace_back(s);
    for (unsigned int i = 0; i < num_traversed; ++i)
    {
      tie(logl, distal, pendant) = insertion_trees[i].place(s);
      pquerys.back().emplace_back(
        i, // branch_id
        logl, // likelihood
        pendant, // pendant length
        distal // distal length
        );
    }
  }
  // now that everything has been placed, we can compute the likelihood weight ratio
  compute_and_set_lwr(pquerys);
  return pquerys;
}
