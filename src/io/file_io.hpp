#pragma once

#include <string>
#include <tuple>
#include <limits>

#include "core/raxml/Model.hpp"
#include "core/pll/pllhead.hpp"
#include "core/pll/rtree_mapper.hpp"
#include "tree/Tree_Numbers.hpp"
#include "tree/Tree.hpp"
#include "seq/MSA.hpp"
#include "seq/MSA_Info.hpp"
#include "util/Options.hpp"

// forward declarations
class MSA;

MSA build_MSA_from_file( 	std::string const& msa_file,
													MSA_Info const& info,
													bool const premasking = false);
pll_utree_s * build_tree_from_file( std::string const& tree_file,
                                    Tree_Numbers& nums,
                                    rtree_mapper& mapper,
                                    bool const preserve_rooting = true );
pll_partition_t * make_partition( raxml::Model const& model,
                                  Tree_Numbers& nums,
                                  int const num_sites,
                                  Options const& options );
void file_check(const std::string& file_path);
std::vector<size_t> get_offsets( std::string const& file, MSA& msa );
int pll_fasta_fseek( 	pll_fasta_t* fd,
											long int const offset,
											int const whence );
