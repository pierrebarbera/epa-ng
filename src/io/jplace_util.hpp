#pragma once

#include <fstream>
#include <vector>

#include "util/stringify.hpp"
#include "sample/PQuery.hpp"
#include "sample/Sample.hpp"
#include "sample/Placement.hpp"
#include "seq/MSA.hpp"

void placement_to_jplace_string(const Placement& p, std::ostream& os);
void pquery_to_jplace_string(const PQuery<Placement>& p, std::ostream& os);
std::string full_jplace_string( const Sample<Placement>& sample,
                                const std::string& invocation);
void init_jplace_string(const std::string& numbered_newick, std::ostream& os);
void finalize_jplace_string(const std::string& invocation, std::ostream& os);
void sample_to_jplace_string(const Sample<Placement>& sample, std::ostream& os);
void merge_into(std::ofstream& dest, const std::vector<std::string>& sources);
