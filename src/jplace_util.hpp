#pragma once

#include "stringify.hpp"
#include "PQuery.hpp"
#include "Sample.hpp"
#include "Placement.hpp"

std::string placement_to_jplace_string(const Placement& p);
std::string pquery_to_jplace_string(const PQuery& p);
std::string sample_to_jplace_string(const Sample& ps, const std::string& invocation);
std::string init_jplace_string(const std::string& numbered_newick);
std::string finalize_jplace_string(const std::string& invocation);
