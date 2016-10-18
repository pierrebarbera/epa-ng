#pragma once

#include <string>
#include <vector>
#include <memory>

#include "pllhead.hpp"
#include "MSA.hpp"
#include "Model.hpp"
#include "Tree_Numbers.hpp"
#include "Options.hpp"
#include "Range.hpp"
#include "Binary.hpp"
#include "Named_Lock.hpp"
#include "pll_util.hpp"

/* Encapsulates the pll data structures for ML computation */
class Tree
{
public:
  Tree(const std::string& tree_file, const MSA& msa, Model& model, Options& options);
  Tree(const std::string& bin_file, Options& options);
  Tree() : partition_(nullptr, pll_partition_destroy), tree_(nullptr, pll_utree_destroy) { }
  ~Tree() = default;

  Tree(Tree const& other) = delete;
  Tree(Tree&& other) = default;

  Tree& operator= (Tree const& other) = delete;
  Tree& operator= (Tree && other) = default;

  // member access
  Tree_Numbers& nums() { return nums_; }
  Model& model() { return model_; }
  Options& options() { return options_; }
  pll_partition_t * partition() { return partition_.get(); }
  pll_utree_t * tree() { return tree_.get(); }

  void * get_clv(const pll_utree_t*);

  double ref_tree_logl();

private:
  // pll structures
  std::unique_ptr<pll_partition_t, partition_deleter> partition_;
  std::unique_ptr<pll_utree_t, utree_deleter> tree_; // must be top level node as parsed in newick! (for jplace)

  // tree related numbers
  Tree_Numbers nums_;

  // epa related classes
  MSA ref_msa_;
  Model model_;
  Options options_;
  Binary binary_;

  // useful internal strucutres
  std::vector<Range> valid_map_;

  // thread safety
  Named_Lock locks_;

};
