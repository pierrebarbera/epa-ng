#ifndef EPA_EPA_H_
#define EPA_EPA_H_

#include <string>
#include <vector>

#include "Tree.hpp"
#include "MSA.hpp"
#include "Model.hpp"
#include "Options.hpp"

void epa(std::string& tree_file, std::string& reference_msa_file, std::string& query_msa_file,
                Model model, Options options, std::string invocation="");

#endif
