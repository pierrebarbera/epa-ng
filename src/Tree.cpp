#include "Tree.hpp"

#include <stdexcept>
#include <tuple>
#include <vector>

#include "pll_util.hpp"
#include "epa_pll_util.hpp"
#include "file_io.hpp"
#include "Sequence.hpp"

using namespace std;

Tree::Tree(const string& tree_file, const MSA& msa, const Model& model, const MSA& query)
  : ref_msa_(msa), query_msa_(query), model_(model)
{
  //parse, build tree
  nums_ = Tree_Numbers();
  tie(partition_, tree_) = build_partition_from_file(tree_file, model_, nums_, ref_msa_.num_sites());

  // split msa if no separate query msa was supplied
  if (query.num_sites() == 0)
    split_combined_msa(ref_msa_, query_msa_, tree_, nums_.tip_nodes);

  link_tree_msa(tree_, partition_, ref_msa_, nums_.tip_nodes);

  // perform branch length optimization on the reference tree
  optimize(tree_, partition_, nums_, model_);

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
    insertion_trees.emplace_back(node, partition_);

  // output class
  PQuery_Set pquerys(get_numbered_newick_string(tree_));

  // place all s on every edge
  for (auto const &s : query_msa_)// make sure a reference, not a copy, is returned
  {
    pquerys.emplace_back(nums_.branches, s);
    for (unsigned int i = 0; i < num_traversed; ++i)
    {
      pquerys.back().set(i, insertion_trees[i].place(s));
    }
  }

  // insertion_trees[0]->place(query_msa_.get(0));
  //
  // free(insertion_trees);
  return pquerys;
}
