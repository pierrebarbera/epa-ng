#pragma once

#include "MSA_Stream.hpp"
#include "Options.hpp"
#include "Tree.hpp"
#include "Model.hpp"
#include "Sample.hpp"
#include "MSA.hpp"

#include <string>

void process(Tree& tree, MSA_Stream& msa_stream, const std::string& outdir, const std::string& invocation);

Sample place(Tree& epa_tree, MSA& query_msa_);
