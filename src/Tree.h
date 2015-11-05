#include <string>
#include <assert.h>
#include <stdexcept>

#ifndef __PLL__
#define __PLL__
#include "pll.h"
#endif

#include "util.h"
#include "MSA.h"


/* Encapsulates the pll data structures for ML computation */
class Tree
{
public:
  Tree(const std::string& tree_file, const std::string& msa_file);
  ~Tree();

private:
  int num_tip_nodes;
  pll_partition_t * partition;
  MSA* msa;

  void build_partition_from_file(const std::string& tree_file, pll_utree_t * (*tree_parse_f) (const char*, int*) );
  MSA* read_msa(const std::string& msa_file); 
};
