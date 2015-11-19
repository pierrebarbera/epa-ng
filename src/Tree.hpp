#ifndef EPA_TREE_H_
#define EPA_TREE_H_

#include <string>
#include <functional>

#include "pllhead.hpp"
#include "MSA.hpp"
#include "Model.hpp"

/* Encapsulates the pll data structures for ML computation */
class Tree
{
public:
  Tree(const std::string& tree_file, const std::string& msa_file, const Model& model);
  ~Tree();
  // TODO should return placement object
  // TODO doesnt follow parallelization scheme: overload?
  void place(const MSA& msa) const;
  // TODO bad: breaks encapsulation
  void visit(std::function<void(pll_partition_t *, pll_utree_t *)> f);

private:
  // pll structures
  pll_partition_t * partition_;
  pll_utree_t * tree_; // must be top level node as parsed in newick! (for jplace)

  // tree related numbers
  int num_branches_;
  int num_tip_nodes_;

  // epa related classes
  MSA ref_msa_;
  Model model_;

  double place_on_edge(Sequence& s, pll_utree_t * node) const;
  void build_partition_from_file(const std::string& tree_file);
  void precompute_clvs();

};

#endif
