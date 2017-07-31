#pragma once

#include "MSA_Stream.hpp"
#include "Options.hpp"
#include "Tree.hpp"
#include "Model.hpp"

#include <string>

void pipeline_place(Tree& tree,
                    const std::string& query_file,
                    const std::string& outdir,
                    const Options& options,
                    const std::string& invocation);

void simple_mpi(Tree& tree,
                const std::string& query_file,
                const std::string& outdir,
                const Options& options,
                const std::string& invocation);

