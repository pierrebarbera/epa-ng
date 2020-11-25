#pragma once

#include "core/raxml/Model.hpp"
#include "seq/MSA_Info.hpp"
#include "util/Options.hpp"

/*
	Try to estimate the memory footprint of the program in bytes.
*/
size_t estimate_footprint( MSA_Info const& ref_info,
                           MSA_Info const& qry_info,
                           raxml::Model const& model,
                           Options const& options );

/*
	convert number of bytes into a formatte string of the highest sensible magnitude
 */
std::string format_byte_num( size_t size );

size_t get_max_memory();
