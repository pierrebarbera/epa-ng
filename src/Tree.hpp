#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "pllhead.hpp"
#include "MSA.hpp"
#include "Model.hpp"
#include "Tree_Numbers.hpp"
#include "Sample.hpp"
#include "Tiny_Tree.hpp"
#include "Options.hpp"
#include "Range.hpp"

/* Encapsulates the pll data structures for ML computation */
class Tree
{
public:
  Tree(const std::string& tree_file, const MSA& msa, Model& model, Options options,
      const MSA& query = MSA());
  ~Tree();
  Sample place();

  // member access
  inline Tree_Numbers nums() const {return nums_;};
  inline Model model() {return model_;};

  double ref_tree_logl();

private:
  // pll structures
  pll_partition_t * partition_;
  pll_utree_t * tree_; // must be top level node as parsed in newick! (for jplace)

  // tree related numbers
  Tree_Numbers nums_;

  // epa related classes
  MSA ref_msa_;
  MSA query_msa_;
  Model model_;
  Options options_;

  // useful internal strucutres
  std::vector<Range> valid_map_;

  double place_on_edge(const Sequence& s, pll_utree_t * node, bool optimize=false) const;

};
