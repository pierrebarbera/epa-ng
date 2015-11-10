#include "stdinc.h"

#ifndef __PLL__
#define __PLL__
extern "C" {
#include "pll.h"
}
#endif

#include "util.h"
#include "MSA.h"
#include "Model.h"


/* Encapsulates the pll data structures for ML computation */
class Tree
{
public:
  Tree(const std::string& tree_file, const std::string& msa_file, Model& model);
  ~Tree();

private:
  pll_partition_t * partition;
  MSA* msa;
  Model model;

  void build_partition_from_file(const std::string& tree_file);
  void link_tree_msa(const int num_tip_nodes, pll_utree_t* tree);
  void precompute_clvs();
};
