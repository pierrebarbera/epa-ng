#ifndef EPA_TREE_H_
#define EPA_TREE_H_

#include <string>

#include "pllhead.hpp"
#include "MSA.hpp"
#include "Model.hpp"
#include "Tree_Numbers.hpp"
#include "Placement.hpp"


/* Encapsulates the pll data structures for ML computation */
class Tree
{
public:
  Tree(const std::string& tree_file, const std::string& msa_file, const Model& model);
  ~Tree();
  // TODO should return placement object
  // TODO doesnt follow parallelization scheme: overload?
  std::vector<Placement> place(const MSA& msa) const;

  // member access
  inline Tree_Numbers nums() const {return nums_;};

private:
  // pll structures
  pll_partition_t * partition_;
  pll_utree_t * tree_; // must be top level node as parsed in newick! (for jplace)

  // tree related numbers
  Tree_Numbers nums_;

  // epa related classes
  MSA ref_msa_;
  Model model_;

  double place_on_edge(const Sequence& s, pll_utree_t * node) const;

};

#endif
