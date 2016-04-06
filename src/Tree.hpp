#pragma once

#include <string>
#include <vector>

#include "pllhead.hpp"
#include "MSA.hpp"
#include "Model.hpp"
#include "Tree_Numbers.hpp"
#include "Sample.hpp"
#include "Tiny_Tree.hpp"
#include "Options.hpp"
#include "Range.hpp"
#include "Binary.hpp"

/* Encapsulates the pll data structures for ML computation */
class Tree
{
public:
  Tree(const std::string& tree_file, const MSA& msa, Model& model, Options options, MSA& query);
  Tree(const std::string& bin_file, const std::string& tree_file, Options& options);
  Tree() : partition_(nullptr), tree_(nullptr) { }
  ~Tree();
  Sample place();

  // member access
  Tree_Numbers& nums() { return nums_; }
  Model& model() { return model_; }
  MSA& query_msa() { return query_msa_; }
  pll_partition_t * partition() { return partition_; }
  pll_utree_t * tree() { return tree_; }

  void * get_clv(unsigned int i);

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
  Binary binary_;

  // useful internal strucutres
  std::vector<Range> valid_map_;

};
