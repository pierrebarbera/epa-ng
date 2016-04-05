#pragma once

#include <string>
#include <tuple>

#include "Tree_Numbers.hpp"
#include "Model.hpp"
#include "pllhead.hpp"

// forward declarations
class MSA;

MSA build_MSA_from_file(const std::string& msa_file);
pll_utree_t * build_tree_from_file(const std::string& tree_file, Tree_Numbers& nums);
pll_partition_t *  build_partition_from_file(const Model& model,
  Tree_Numbers& nums,  const int num_sites);
void file_check(const std::string& file_path);
// pll_partition_t *  build_partition_from_binary(const std::string& bin_file, const bool out_of_core);
