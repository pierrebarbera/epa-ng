#pragma once

#include <memory>

#include "seq/MSA_Stream.hpp"
#include "io/Binary_Fasta.hpp"
#include "io/file_io.hpp"
#include "util/stringify.hpp"
#include "util/logging.hpp"
#include "util/Options.hpp"
#include "io/msa_reader_interface.hpp"

auto make_msa_reader( const std::string& file_name,
                      const Options& options,
                      const size_t sites = 0)
{
  std::unique_ptr<msa_reader> result(nullptr);

  try {
    result = std::make_unique<Binary_Fasta_Reader>(file_name);
  } catch(const std::exception&) {
    LOG_DBG << "Failed to parse input as binary fasta (bfast), trying `fasta` instead.";
    result = std::make_unique<MSA_Stream>(file_name, options.chunk_size);
  }

  return std::move(result);
}
