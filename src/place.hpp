#pragma once

#include "MSA_Stream.hpp"
#include "Options.hpp"
#include "Tree.hpp"
#include "Model.hpp"

#include <string>

void process( Tree& tree, 
              MSA_Stream& msa_stream, 
              const std::string& outdir, 
              const Options& options,
              const std::string& invocation);

void tmp_pipeline_test( Tree& tree, 
                        MSA_Stream& msa_stream, 
                        const std::string& outdir, 
                        const Options& options,
                        const std::string& invocation);

void simple_mpi(Tree& tree, 
                const std::string& query_file, 
                const std::string& outdir, 
                const Options& options,
                const std::string& invocation);

