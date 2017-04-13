#pragma once

#include <fstream>
#include <vector>

#include "stringify.hpp"
#include "PQuery.hpp"
#include "Sample.hpp"
#include "Placement.hpp"
#include "MSA.hpp"

std::string placement_to_jplace_string(const Placement& p);
std::string pquery_to_jplace_string(const PQuery& p, const MSA& msa);
std::string full_jplace_string( const Sample& sample, 
                                const std::string& invocation, 
                                const MSA& msa);
std::string init_jplace_string(const std::string& numbered_newick);
std::string finalize_jplace_string(const std::string& invocation);
std::string sample_to_jplace_string(const Sample& sample, const MSA& msa);
void merge_into(std::ofstream& dest, const std::vector<std::string>& sources);
