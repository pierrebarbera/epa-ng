#pragma once

#include <string>
#include <vector>

#include "Tree.hpp"
#include "MSA.hpp"
#include "Model.hpp"
#include "Options.hpp"

void epa(const std::string& tree_file, const std::string& reference_msa_file, const std::string& query_msa_file,
  std::string& outdir, Model model, const Options& options, const std::string& invocation);

// void epa(const std::string& binary_file, const string& tree_file, const std::string &query_msa_file,
//   const std::string& outdir, Options& options);
