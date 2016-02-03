#pragma once

#include <string>
#include <tuple>

#include "Tree_Numbers.hpp"
#include "Model.hpp"
#include "pllhead.hpp"

// forward declarations
class MSA;

/* reads in sequences from a file
  msa_file: string specifying the file path
*/
MSA build_MSA_from_file(const std::string& msa_file);
std::tuple<pll_partition_t *, pll_utree_t *> build_partition_from_file(
  const std::string& tree_file, const Model& model, Tree_Numbers& nums, const int num_sites);
void file_check(const std::string& file_path);
