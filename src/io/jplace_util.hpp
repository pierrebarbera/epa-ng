#pragma once

#include <fstream>
#include <vector>

#include "util/stringify.hpp"
#include "sample/PQuery.hpp"
#include "sample/Sample.hpp"
#include "sample/Placement.hpp"
#include "seq/MSA.hpp"
#include "core/pll/rtree_mapper.hpp"

void sample_to_jplace_string( Sample<Placement> const& sample, std::ostream& os, rtree_mapper const& mapper );
void pquery_to_jplace_string( PQuery<Placement> const& p, std::ostream& os, rtree_mapper const& mapper );
void placement_to_jplace_string( Placement const& p, std::ostream& os, rtree_mapper const& mapper );
std::string full_jplace_string( Sample<Placement> const& sample,
                                std::string const& invocation,
                                rtree_mapper const& mapper );
void init_jplace_string( std::string const& numbered_newick, std::ostream& os );
void finalize_jplace_string( std::string const& invocation, std::ostream& os );
void merge_into( std::ofstream& dest, std::vector<std::string> const& sources );
