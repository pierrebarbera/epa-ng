#pragma once

#include "core/pll/pllhead.hpp"
#include "core/raxml/Model.hpp"
#include "io/file_io.hpp"
#include "seq/MSA_Info.hpp"
#include "util/Options.hpp"
#include "util/logging.hpp"

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