#ifndef EPA_TREE_H_
#define EPA_TREE_H_

#include <string>

#ifndef PLL_H_
#define PLL_H_
extern "C" {
#include "pll.h"
}
#endif

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
  void precompute_clvs(int num_tip_nodes, pll_utree_t* tree);
};

#endif
