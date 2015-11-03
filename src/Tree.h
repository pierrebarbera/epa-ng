#include <string>
#include <cassert>

#ifndef __PLL__
#define __PLL__
#include "pll.h"
#endif

#include "util.h"

class Tree
{
public:
  Tree(std::string* tree_file, std::string* msa_file);
  ~Tree();

private:
  int num_tip_nodes;
  pll_partition_t * partition;

  void build_tree_from_file(std::string* tree_file, pll_utree_t * (*tree_parse_f) (const char*, int*) );
  void get_msa_from_file(std::string* msa_file);
};
