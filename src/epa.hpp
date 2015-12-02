#ifndef EPA_EPA_H_
#define EPA_EPA_H_

#include <string>
#include <vector>

#include "Tree.hpp"
#include "MSA.hpp"

void epa(std::string& tree_file, std::string& reference_msa_file, std::string& query_msa_file,
                std::vector<double> base_frequencies, std::vector<double> substitution_rates,
                double alpha, bool heuristic, std::string invocation="");

#endif
