#ifndef EPA_TREE_H_
#define EPA_TREE_H_

#include <string>

#include "pllhead.hpp"
#include "MSA.hpp"
#include "Model.hpp"

/* Encapsulates the pll data structures for ML computation */
class Tree
{
public:
  Tree(const std::string& tree_file, const std::string& msa_file, Model& model);
  ~Tree();
  // TODO should return placement object
  // TODO doesnt follow parallelization scheme: overload?
  void place(const Sequence& s) const;

private:
  pll_partition_t * partition;
  pll_utree_t * tree;
  MSA msa;
  Model model;

  double place_on_edge(Sequence& s, pll_utree_t * node);
  void build_partition_from_file(const std::string& tree_file);
  void link_tree_msa(const int num_tip_nodes, pll_utree_t* tree);
  void precompute_clvs(int num_tip_nodes, pll_utree_t* tree);

};

#endif
