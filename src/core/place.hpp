#pragma once

#include "core/raxml/Model.hpp"
#include "seq/MSA_Info.hpp"
#include "seq/MSA_Stream.hpp"
#include "tree/Tree.hpp"
#include "util/Options.hpp"

#include <string>

void simple_mpi( Tree& tree,
                 std::string const& query_file,
                 MSA_Info const& msa_info,
                 std::string const& outdir,
                 Options const& options,
                 std::string const& invocation );
